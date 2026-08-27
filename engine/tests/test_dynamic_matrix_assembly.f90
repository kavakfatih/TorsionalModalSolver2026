program test_dynamic_matrix_assembly
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
    ieee_positive_inf, ieee_negative_inf
  use tms_kinds, only : dp
  use tms_dof_types, only : TORSIONAL_ROTATION
  use tms_local_matrix, only : local_matrix_2x2
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t, &
    validate_torsional_element, get_local_stiffness, &
    get_local_loss_stiffness, get_local_damping
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element, get_torsional_element, &
    validate_torsional_system
  use tms_dof_map, only : dof_map_t, initialize_full_dof_map
  use tms_constraint_types, only : constraint_t, FIXED_CONSTRAINT, &
    PRESCRIBED_VALUE_CONSTRAINT
  use tms_constraint_manager, only : constraint_manager_t, &
    initialize_constraint_manager, add_constraint
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    get_stiffness_matrix_values
  use tms_loss_stiffness_matrix, only : loss_stiffness_matrix_t, &
    get_loss_stiffness_matrix_values
  use tms_damping_matrix, only : damping_matrix_t, &
    get_damping_matrix_values
  use tms_mass_matrix, only : mass_matrix_t, get_mass_matrix_values
  use tms_matrix_assembly, only : assemble_full_stiffness, &
    assemble_full_loss_stiffness, assemble_full_damping, &
    assemble_full_inertia
  use tms_reduced_dynamic_system, only : &
    reduced_dynamic_torsional_system_t, &
    build_reduced_dynamic_torsional_system, &
    get_reduced_dynamic_stiffness, &
    get_reduced_dynamic_loss_stiffness, &
    get_reduced_dynamic_damping, get_reduced_dynamic_mass, &
    get_reduced_dynamic_active_dof_count, recover_harmonic_response
  use tms_torsional_system, only : two_inertia_tvd_system_t, &
    build_generalized_two_inertia_system
  implicit none

  real(dp), parameter :: tolerance = 1.0e-12_dp
  character(len=64) :: validation_case

  ! Geçersiz K'' ve c girdileri aynı executable dosyasının CTest
  ! WILL_FAIL selector vakalarıyla üretim validator'ında sınanır.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_local_dynamic_contributions()
  call test_full_dynamic_matrix_assembly()
  call test_reduced_dynamic_matrices_and_recovery()
  call test_prescribed_offset_not_recovered()
  call test_two_inertia_loss_stiffness_bridge()
  call test_input_immutability()

  print *, "K', K'', C ve M dynamic matrix katmanı doğrulandı."

contains

  !> Tek elemanın K', K'' ve C lokal katkılarını exact referansla sınar.
  !! Fiziksel model: node_i->node_j yönünde K'=100 N*m/rad,
  !! K''=10 N*m/rad ve c=2 N*m*s/rad pasif lineer bağlantıdır.
  !! Matematiksel model: Her kanal a*[[1,-1],[-1,1]] üretir; simetri,
  !! sıfır satır toplamı ve rigid [1,1] null katkısı korunur.
  !! Girdi ve çıktı test içinde sabittir; production getter'ları kullanılır.
  subroutine test_local_dynamic_contributions()
    type(torsional_element_t) :: default_element
    type(torsional_element_t) :: element
    type(local_matrix_2x2) :: local_stiffness
    type(local_matrix_2x2) :: local_loss_stiffness
    type(local_matrix_2x2) :: local_damping
    real(dp) :: expected_stiffness(2, 2)
    real(dp) :: expected_loss_stiffness(2, 2)
    real(dp) :: expected_damping(2, 2)

    default_element = torsional_element_t( &
      id=2, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=0.0_dp)
    if (abs(default_element%loss_stiffness_nm_per_rad) > tolerance) then
      error stop "Torsional elemanın default K'' değeri sıfır değil."
    end if

    element = torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=2.0_dp, &
      loss_stiffness_nm_per_rad=10.0_dp)
    local_stiffness = get_local_stiffness(element)
    local_loss_stiffness = get_local_loss_stiffness(element)
    local_damping = get_local_damping(element)

    expected_stiffness = connection_matrix(100.0_dp)
    expected_loss_stiffness = connection_matrix(10.0_dp)
    expected_damping = connection_matrix(2.0_dp)
    call assert_matrix_close( &
      local_stiffness%value, expected_stiffness, &
      "Lokal K' katkısı exact referansla uyuşmuyor.")
    call assert_matrix_close( &
      local_loss_stiffness%value, expected_loss_stiffness, &
      "Lokal K'' katkısı exact referansla uyuşmuyor.")
    call assert_matrix_close( &
      local_damping%value, expected_damping, &
      "Lokal C katkısı exact referansla uyuşmuyor.")

    call assert_connection_invariants( &
      local_loss_stiffness%value, "Lokal K''")
    call assert_connection_invariants(local_damping%value, "Lokal C")
  end subroutine test_local_dynamic_contributions

  !> Üç düğümlü zincirin full K', K'', C ve M assembly sonucunu sınar.
  !! Fiziksel model: 10--(100,10,1)--20--(200,40,3)--30; parantezler
  !! sırasıyla K' [N*m/rad], K'' [N*m/rad], c [N*m*s/rad] değerleridir.
  !! Düğüm ataletleri [0.1,0.2,0.3] kg*m^2'dir.
  !! Matematiksel model: Her lokal katkı full equation uzayına scatter-add
  !! edilir. Test exact katsayıları, simetriyi ve K''/C rigid null modunu
  !! production assembly API'leri üzerinden doğrular.
  subroutine test_full_dynamic_matrix_assembly()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: stiffness
    type(loss_stiffness_matrix_t) :: loss_stiffness
    type(damping_matrix_t) :: damping
    type(mass_matrix_t) :: mass
    real(dp), allocatable :: actual_stiffness(:, :)
    real(dp), allocatable :: actual_loss_stiffness(:, :)
    real(dp), allocatable :: actual_damping(:, :)
    real(dp), allocatable :: actual_mass(:, :)
    real(dp) :: expected_stiffness(3, 3)
    real(dp) :: expected_loss_stiffness(3, 3)
    real(dp) :: expected_damping(3, 3)
    real(dp) :: expected_mass(3, 3)

    call build_three_node_dynamic_chain(system)
    call initialize_full_dof_map(mapping, system)
    stiffness = assemble_full_stiffness(system, mapping)
    loss_stiffness = assemble_full_loss_stiffness(system, mapping)
    damping = assemble_full_damping(system, mapping)
    mass = assemble_full_inertia(system, mapping)

    actual_stiffness = get_stiffness_matrix_values(stiffness)
    actual_loss_stiffness = &
      get_loss_stiffness_matrix_values(loss_stiffness)
    actual_damping = get_damping_matrix_values(damping)
    actual_mass = get_mass_matrix_values(mass)

    expected_stiffness = reshape([ &
      100.0_dp, -100.0_dp, 0.0_dp, &
      -100.0_dp, 300.0_dp, -200.0_dp, &
      0.0_dp, -200.0_dp, 200.0_dp], [3, 3])
    expected_loss_stiffness = reshape([ &
      10.0_dp, -10.0_dp, 0.0_dp, &
      -10.0_dp, 50.0_dp, -40.0_dp, &
      0.0_dp, -40.0_dp, 40.0_dp], [3, 3])
    expected_damping = reshape([ &
      1.0_dp, -1.0_dp, 0.0_dp, &
      -1.0_dp, 4.0_dp, -3.0_dp, &
      0.0_dp, -3.0_dp, 3.0_dp], [3, 3])
    expected_mass = reshape([ &
      0.1_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.2_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.3_dp], [3, 3])

    call assert_matrix_close(actual_stiffness, expected_stiffness, &
      "Full K' assembly exact referansla uyuşmuyor.")
    call assert_matrix_close( &
      actual_loss_stiffness, expected_loss_stiffness, &
      "Full K'' assembly exact referansla uyuşmuyor.")
    call assert_matrix_close(actual_damping, expected_damping, &
      "Full C assembly exact referansla uyuşmuyor.")
    call assert_matrix_close(actual_mass, expected_mass, &
      "Full M assembly exact referansla uyuşmuyor.")

    call assert_connection_invariants(actual_loss_stiffness, "Full K''")
    call assert_connection_invariants(actual_damping, "Full C")
  end subroutine test_full_dynamic_matrix_assembly

  !> İlk düğüm fixed iken dört matrisin aynı principal reduction'ını
  !! ve kompleks harmonik recovery'yi doğrular.
  !! Matematiksel model: A_r=P^T*A_full*P, A={K',K'',C,M}; retained sıra
  !! [node20,node30]'dur. theta_hat=P*theta_hat_r ile node10 sıfır kalır.
  !! Girdi birimleri matris sözleşmeleriyle aynı, response [rad] birimindedir.
  subroutine test_reduced_dynamic_matrices_and_recovery()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_dynamic_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: stiffness
    type(loss_stiffness_matrix_t) :: loss_stiffness
    type(damping_matrix_t) :: damping
    type(mass_matrix_t) :: mass
    complex(dp) :: reduced_response(2)
    complex(dp), allocatable :: physical_response(:)

    call build_three_node_dynamic_chain(system)
    call initialize_constraint_manager(manager)
    call add_constraint( &
      manager, make_constraint(1, 10, 0.0_dp, FIXED_CONSTRAINT), system)
    reduced_system = &
      build_reduced_dynamic_torsional_system(system, manager)

    if (get_reduced_dynamic_active_dof_count(reduced_system) /= 2) then
      error stop "Reduced dynamic active DOF sayısı iki değil."
    end if

    stiffness = get_reduced_dynamic_stiffness(reduced_system)
    loss_stiffness = &
      get_reduced_dynamic_loss_stiffness(reduced_system)
    damping = get_reduced_dynamic_damping(reduced_system)
    mass = get_reduced_dynamic_mass(reduced_system)

    call assert_matrix_close( &
      get_stiffness_matrix_values(stiffness), reshape([ &
        300.0_dp, -200.0_dp, -200.0_dp, 200.0_dp], [2, 2]), &
      "Reduced K' principal submatrix yanlış.")
    call assert_matrix_close( &
      get_loss_stiffness_matrix_values(loss_stiffness), reshape([ &
        50.0_dp, -40.0_dp, -40.0_dp, 40.0_dp], [2, 2]), &
      "Reduced K'' principal submatrix yanlış.")
    call assert_matrix_close( &
      get_damping_matrix_values(damping), reshape([ &
        4.0_dp, -3.0_dp, -3.0_dp, 3.0_dp], [2, 2]), &
      "Reduced C principal submatrix yanlış.")
    call assert_matrix_close( &
      get_mass_matrix_values(mass), reshape([ &
        0.2_dp, 0.0_dp, 0.0_dp, 0.3_dp], [2, 2]), &
      "Reduced M principal submatrix yanlış.")

    reduced_response = [ &
      cmplx(1.0_dp, 0.5_dp, kind=dp), &
      cmplx(-0.25_dp, 0.75_dp, kind=dp)]
    physical_response = &
      recover_harmonic_response(reduced_system, reduced_response)
    call assert_complex_vector_close(physical_response, [ &
      cmplx(0.0_dp, 0.0_dp, kind=dp), reduced_response(1), &
      reduced_response(2)], "Kompleks harmonik recovery yanlış.")
  end subroutine test_reduced_dynamic_matrices_and_recovery

  !> Saklanan statik prescribed offset'in harmonik phasor'a eklenmediğini
  !! doğrular. Fiziksel modelde node10 için q_p=0.125 rad saklanır, ancak
  !! V0.6 perturbation recovery theta_hat=P*theta_hat_r olduğundan constrained
  !! complex response tam olarak sıfır kalmalıdır.
  subroutine test_prescribed_offset_not_recovered()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_dynamic_torsional_system_t) :: reduced_system
    complex(dp) :: reduced_response(2)
    complex(dp), allocatable :: physical_response(:)

    call build_three_node_dynamic_chain(system)
    call initialize_constraint_manager(manager)
    call add_constraint(manager, make_constraint( &
      1, 10, 0.125_dp, PRESCRIBED_VALUE_CONSTRAINT), system)
    reduced_system = &
      build_reduced_dynamic_torsional_system(system, manager)

    reduced_response = [ &
      cmplx(0.2_dp, -0.1_dp, kind=dp), &
      cmplx(0.4_dp, 0.3_dp, kind=dp)]
    physical_response = &
      recover_harmonic_response(reduced_system, reduced_response)
    if (abs(physical_response(1)) > tolerance) then
      error stop "Statik prescribed offset harmonik phasor'a eklendi."
    end if
    call assert_complex_vector_close( &
      physical_response(2:3), reduced_response, &
      "Active harmonik cevap recovery sırasında değişti.")
  end subroutine test_prescribed_offset_not_recovered

  !> İki-atalet TVD köprüsünün K'' değerini ayrı element kanalına
  !! aktardığını ve c=0 bıraktığını doğrular.
  !! Girdiler: J_h/J_r [kg*m^2], K'/K'' [N*m/rad]. Çıktı generalized
  !! eleman alanlarıdır. Test K''=c veya K''/omega dönüşümünü engeller.
  subroutine test_two_inertia_loss_stiffness_bridge()
    type(two_inertia_tvd_system_t) :: source_system
    type(torsional_system_t) :: generalized_system
    type(torsional_element_t) :: element

    source_system = two_inertia_tvd_system_t( &
      hub_polar_inertia_kg_m2=0.1_dp, &
      ring_polar_inertia_kg_m2=0.2_dp, &
      storage_stiffness_nm_per_rad=1000.0_dp, &
      loss_stiffness_nm_per_rad=100.0_dp, &
      loss_factor=0.1_dp, &
      material_reference_frequency_hz=100.0_dp, &
      material_temperature_k=293.15_dp)
    generalized_system = &
      build_generalized_two_inertia_system(source_system)
    element = get_torsional_element(generalized_system, 1)

    if (abs(element%stiffness_nm_per_rad - 1000.0_dp) > tolerance .or. &
        abs(element%loss_stiffness_nm_per_rad - 100.0_dp) > tolerance) then
      error stop "TVD köprüsü K' veya K'' değerini aynen aktarmadı."
    end if
    if (abs(element%damping_nms_per_rad) > tolerance) then
      error stop "TVD köprüsü K'' değerini yanlışlıkla viskoz c'ye aktardı."
    end if
  end subroutine test_two_inertia_loss_stiffness_bridge

  !> Assembly, getter ve recovery yordamlarının authoritative girdileri
  !! değiştirmediğini sınar. Matris katsayıları kendi SI birimlerindedir;
  !! kompleks response [rad] birimindedir. Getter'dan dönen dense kopyanın
  !! değiştirilmesi private matrix storage'a geri yazmamalıdır.
  subroutine test_input_immutability()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(loss_stiffness_matrix_t) :: loss_stiffness
    type(torsional_element_t) :: element_before
    type(torsional_element_t) :: element_after
    real(dp), allocatable :: values_copy(:, :)
    real(dp), allocatable :: values_after(:, :)

    call build_three_node_dynamic_chain(system)
    element_before = get_torsional_element(system, 1)
    call initialize_full_dof_map(mapping, system)
    loss_stiffness = assemble_full_loss_stiffness(system, mapping)
    values_copy = get_loss_stiffness_matrix_values(loss_stiffness)
    values_copy(1, 1) = -999.0_dp
    values_after = get_loss_stiffness_matrix_values(loss_stiffness)
    element_after = get_torsional_element(system, 1)

    if (abs(values_after(1, 1) - 10.0_dp) > tolerance) then
      error stop "K'' getter kopyası private matrix storage'ı değiştirdi."
    end if
    if (element_after%id /= element_before%id .or. &
        abs(element_after%stiffness_nm_per_rad - &
        element_before%stiffness_nm_per_rad) > tolerance .or. &
        abs(element_after%loss_stiffness_nm_per_rad - &
        element_before%loss_stiffness_nm_per_rad) > tolerance .or. &
        abs(element_after%damping_nms_per_rad - &
        element_before%damping_nms_per_rad) > tolerance) then
      error stop "Dynamic matrix assembly torsional eleman girdisini değiştirdi."
    end if
  end subroutine test_input_immutability

  !> Üç düğümlü dynamic zinciri ortak test fixture olarak kurar.
  !! Düğüm ataletleri [kg*m^2], eleman K'/K'' değerleri [N*m/rad] ve
  !! c değerleri [N*m*s/rad] birimindedir. Çıktı geçerli torsional sistemdir.
  subroutine build_three_node_dynamic_chain(system)
    type(torsional_system_t), intent(out) :: system

    call add_torsional_node(system, torsional_node_t( &
      id=10, polar_inertia_kg_m2=0.1_dp, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=20, polar_inertia_kg_m2=0.2_dp, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=30, polar_inertia_kg_m2=0.3_dp, constrained=.false.))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=10, node_j_id=20, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=1.0_dp, &
      loss_stiffness_nm_per_rad=10.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      id=2, node_i_id=20, node_j_id=30, &
      stiffness_nm_per_rad=200.0_dp, damping_nms_per_rad=3.0_dp, &
      loss_stiffness_nm_per_rad=40.0_dp))
    call validate_torsional_system(system)
  end subroutine build_three_node_dynamic_chain

  !> Tek bir torsional constraint kaydını üretir. Kimlikler boyutsuz,
  !! value [rad] birimindedir. Çıktı FIXED veya PRESCRIBED_VALUE türündeki
  !! constraint_t'dir; fiziksel DOF TORSIONAL_ROTATION olarak sabittir.
  pure function make_constraint( &
      constraint_id, node_id, value, constraint_type) result(constraint)
    integer, intent(in) :: constraint_id
    integer, intent(in) :: node_id
    real(dp), intent(in) :: value
    integer, intent(in) :: constraint_type
    type(constraint_t) :: constraint

    constraint = constraint_t( &
      constraint_id=constraint_id, node_id=node_id, &
      dof_type=TORSIONAL_ROTATION, value=value, &
      constraint_type=constraint_type)
  end function make_constraint

  !> İki düğümlü bağlantı matrisinin analitik test referansını kurar.
  !! Matematiksel model A=a[[1,-1],[-1,1]]. Girdi a, çağıran kanalın
  !! SI birimindedir; çıktı aynı birimli 2x2 matristir. Bu helper yalnız
  !! exact assembly referansı için test kodunda kullanılır.
  pure function connection_matrix(coefficient) result(matrix)
    real(dp), intent(in) :: coefficient
    real(dp) :: matrix(2, 2)

    matrix = reshape([ &
      coefficient, -coefficient, &
      -coefficient, coefficient], [2, 2])
  end function connection_matrix

  !> K'' veya C bağlantı matrisinin simetri, sıfır satır toplamı ve
  !! rigid common motion null-space invariantlarını doğrular.
  !! Girdi square matrix kendi SI birimindedir; label tanı bilgisidir.
  subroutine assert_connection_invariants(matrix, label)
    real(dp), intent(in) :: matrix(:, :)
    character(len=*), intent(in) :: label
    real(dp), allocatable :: rigid_motion(:)

    if (size(matrix, 1) /= size(matrix, 2)) then
      error stop label // " matrisi kare değil."
    end if
    if (maxval(abs(matrix - transpose(matrix))) > tolerance) then
      error stop label // " matrisi simetrik değil."
    end if
    if (maxval(abs(sum(matrix, dim=2))) > tolerance) then
      error stop label // " matrisinin satır toplamı sıfır değil."
    end if

    allocate(rigid_motion(size(matrix, 1)))
    rigid_motion = 1.0_dp
    if (maxval(abs(matmul(matrix, rigid_motion))) > tolerance) then
      error stop label // " matrisi rigid common motion'u korumuyor."
    end if
  end subroutine assert_connection_invariants

  !> Gerçek matrisleri scale-aware mutlak/bağıl toleransla karşılaştırır.
  !! Girdiler aynı fiziksel birimde olmalıdır. Boyut veya katsayı farkı
  !! tolerance*(1+max|expected|) sınırını aşarsa error stop üretir.
  subroutine assert_matrix_close(actual, expected, message)
    real(dp), intent(in) :: actual(:, :)
    real(dp), intent(in) :: expected(:, :)
    character(len=*), intent(in) :: message
    real(dp) :: scale

    if (any(shape(actual) /= shape(expected))) error stop message
    scale = 1.0_dp + maxval(abs(expected))
    if (maxval(abs(actual - expected)) > tolerance * scale) then
      error stop message
    end if
  end subroutine assert_matrix_close

  !> Kompleks response vektörlerini scale-aware mutlak/bağıl toleransla
  !! karşılaştırır. Girdiler [rad] gibi aynı fiziksel birimde olmalıdır.
  subroutine assert_complex_vector_close(actual, expected, message)
    complex(dp), intent(in) :: actual(:)
    complex(dp), intent(in) :: expected(:)
    character(len=*), intent(in) :: message
    real(dp) :: scale

    if (size(actual) /= size(expected)) error stop message
    scale = 1.0_dp + maxval(abs(expected))
    if (maxval(abs(actual - expected)) > tolerance * scale) then
      error stop message
    end if
  end subroutine assert_complex_vector_close

  !> Geçersiz K'' veya c değerini production element validator'ına gönderir.
  !! K'' [N*m/rad], c [N*m*s/rad] birimindedir. Negatif, NaN ve sonsuz
  !! pasif olmayan/nonfinite değerler error stop ile reddedilmelidir.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name
    type(torsional_element_t) :: invalid_element

    invalid_element = torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=1.0_dp, &
      loss_stiffness_nm_per_rad=10.0_dp)

    select case (case_name)
      case ("negative_loss_stiffness")
        invalid_element%loss_stiffness_nm_per_rad = -1.0_dp
      case ("nan_loss_stiffness")
        invalid_element%loss_stiffness_nm_per_rad = &
          ieee_value(0.0_dp, ieee_quiet_nan)
      case ("positive_infinity_loss_stiffness")
        invalid_element%loss_stiffness_nm_per_rad = &
          ieee_value(0.0_dp, ieee_positive_inf)
      case ("negative_infinity_loss_stiffness")
        invalid_element%loss_stiffness_nm_per_rad = &
          ieee_value(0.0_dp, ieee_negative_inf)
      case ("negative_damping")
        invalid_element%damping_nms_per_rad = -1.0_dp
      case ("nan_damping")
        invalid_element%damping_nms_per_rad = &
          ieee_value(0.0_dp, ieee_quiet_nan)
      case ("positive_infinity_damping")
        invalid_element%damping_nms_per_rad = &
          ieee_value(0.0_dp, ieee_positive_inf)
      case ("negative_infinity_damping")
        invalid_element%damping_nms_per_rad = &
          ieee_value(0.0_dp, ieee_negative_inf)
      case default
        error stop "Bilinmeyen dynamic matrix validation selector istendi."
    end select

    call validate_torsional_element(invalid_element)
    print *, "Geçersiz dynamic eleman beklenmedik biçimde kabul edildi: ", &
      case_name
  end subroutine exercise_invalid_case

end program test_dynamic_matrix_assembly
