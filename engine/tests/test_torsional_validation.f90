module tms_validation_test_helpers
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  public :: assert_absolute_close
  public :: assert_relative_close
  public :: assert_matrix_close
  public :: assert_vector_close
  public :: assert_symmetric_matrix
  public :: assert_nonnegative_strain_energy

contains

  !> İki skaler sonucu aynı fiziksel birimde mutlak toleransla karşılaştırır.
  !! Matematiksel koşul: |actual-expected| <= tolerance. Girdiler sonlu ve
  !! aynı SI biriminde, tolerans ise aynı birimde ve negatif olmamalıdır.
  pure subroutine assert_absolute_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (tolerance < 0.0_dp .or. .not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected)) then
      error stop "Mutlak karşılaştırma girdileri sonlu ve tolerans geçerli olmalıdır."
    end if

    if (abs(actual - expected) > tolerance) error stop message
  end subroutine assert_absolute_close

  !> İki pozitif skaler sonucu boyutsuz bağıl toleransla karşılaştırır.
  !! Matematiksel koşul: |actual-expected|/|expected| <= tolerance. Girdiler
  !! aynı SI biriminde ve beklenen değer sıfırdan farklı olmalıdır.
  pure subroutine assert_relative_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (tolerance < 0.0_dp .or. abs(expected) <= tiny(1.0_dp) .or. &
        .not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) then
      error stop "Bağıl karşılaştırma girdileri sonlu ve referans sıfırdan farklı olmalıdır."
    end if

    if (abs(actual - expected) / abs(expected) > tolerance) then
      error stop message
    end if
  end subroutine assert_relative_close

  !> Aynı boyut ve fiziksel birimdeki iki matrisin tüm katsayılarını doğrular.
  !! Matematiksel koşul maksimum katsayı hatasının mutlak toleransı aşmamasıdır.
  pure subroutine assert_matrix_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:, :)
    real(dp), intent(in) :: expected(:, :)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (any(shape(actual) /= shape(expected))) error stop message
    if (tolerance < 0.0_dp .or. .not. all(ieee_is_finite(actual)) .or. &
        .not. all(ieee_is_finite(expected))) then
      error stop "Matris karşılaştırma girdileri sonlu ve tolerans geçerli olmalıdır."
    end if

    if (maxval(abs(actual - expected)) > tolerance) error stop message
  end subroutine assert_matrix_close

  !> Aynı boyut ve fiziksel birimdeki iki vektörü mutlak toleransla doğrular.
  !! Matematiksel ölçüt maksimum bileşen hatasıdır.
  pure subroutine assert_vector_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:)
    real(dp), intent(in) :: expected(:)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (size(actual) /= size(expected)) error stop message
    if (tolerance < 0.0_dp .or. .not. all(ieee_is_finite(actual)) .or. &
        .not. all(ieee_is_finite(expected))) then
      error stop "Vektör karşılaştırma girdileri sonlu ve tolerans geçerli olmalıdır."
    end if

    if (maxval(abs(actual - expected)) > tolerance) error stop message
  end subroutine assert_vector_close

  !> Kare matrisin simetri kalıntısını Frobenius normuyla sınar.
  !! Matematiksel koşul: ||A-A^T||_F <= tolerance. Matrisin fiziksel birimi
  !! ile tolerans aynıdır; yordam matrisi veya üretim durumunu değiştirmez.
  pure subroutine assert_symmetric_matrix(matrix, tolerance, message)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    real(dp) :: symmetry_norm

    if (size(matrix, 1) /= size(matrix, 2)) then
      error stop "Simetri kontrolü kare matris gerektirir."
    end if
    if (tolerance < 0.0_dp .or. .not. all(ieee_is_finite(matrix))) then
      error stop "Simetri kontrolü sonlu matris ve geçerli tolerans gerektirir."
    end if

    symmetry_norm = sqrt(sum((matrix - transpose(matrix))**2))
    if (.not. ieee_is_finite(symmetry_norm) .or. &
        symmetry_norm > tolerance) then
      error stop message
    end if
  end subroutine assert_symmetric_matrix

  !> Torsional rijitlik matrisinin negatif elastik enerji üretmediğini sınar.
  !! Fiziksel model: U=1/2 theta^T K theta. K [N*m/rad], theta [rad] ve U
  !! [N*m] birimindedir. Küçük negatif yuvarlama hatası mutlak enerji
  !! toleransına kadar kabul edilir; çıktı üretmez.
  pure subroutine assert_nonnegative_strain_energy( &
      stiffness, theta, energy_tolerance_nm, message)
    real(dp), intent(in) :: stiffness(:, :)
    real(dp), intent(in) :: theta(:)
    real(dp), intent(in) :: energy_tolerance_nm
    character(len=*), intent(in) :: message

    real(dp) :: strain_energy_nm

    if (size(stiffness, 1) /= size(stiffness, 2) .or. &
        size(stiffness, 1) /= size(theta)) then
      error stop "Enerji kontrolünde matris ve dönme vektörü boyutları uyumsuz."
    end if
    if (energy_tolerance_nm < 0.0_dp .or. &
        .not. all(ieee_is_finite(stiffness)) .or. &
        .not. all(ieee_is_finite(theta))) then
      error stop "Enerji kontrolü sonlu girdiler ve geçerli tolerans gerektirir."
    end if

    strain_energy_nm = 0.5_dp * &
      dot_product(theta, matmul(stiffness, theta))
    if (.not. ieee_is_finite(strain_energy_nm) .or. &
        strain_energy_nm < -energy_tolerance_nm) then
      error stop message
    end if
  end subroutine assert_nonnegative_strain_energy

end module tms_validation_test_helpers

program test_torsional_validation
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_local_matrix, only : local_matrix_2x2
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t, get_local_stiffness
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element
  use tms_dof_map, only : dof_map_t, initialize_dof_map, &
    lookup_equation_id, get_active_dof_count, get_dof_map_entry_count
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    get_stiffness_matrix_values
  use tms_mass_matrix, only : mass_matrix_t, get_mass_matrix_values
  use tms_matrix_assembly, only : assemble_stiffness, assemble_inertia
  use tms_torsional_system, only : two_inertia_tvd_system_t, &
    two_inertia_modal_result_t, build_generalized_two_inertia_system, &
    solve_free_free_two_inertia_modes
  use tms_validation_test_helpers, only : assert_absolute_close, &
    assert_relative_close, assert_matrix_close, assert_vector_close, &
    assert_symmetric_matrix, assert_nonnegative_strain_energy
  implicit none

  real(dp), parameter :: absolute_tolerance = 1.0e-10_dp
  real(dp), parameter :: relative_tolerance = 1.0e-10_dp

  call test_single_torsional_element()
  call test_rigid_body_mode()
  call test_global_assembly_regression()
  call test_two_inertia_analytic_benchmark()
  call test_dof_mapping_regression()
  call test_matrix_quality()

  print *, "V0.3.1 torsional analitik doğrulama kontrolleri geçti."

contains

  !> Tek lineer torsional elemanın lokal rijitlik işaret konvansiyonunu sınar.
  !! Model: k=100 N*m/rad ve K_e=k[[1,-1],[-1,1]]. Girdi üretim elemanıdır;
  !! çıktı yoktur. Küçük dönme, ortak eksen ve kütlesiz bağlantı varsayılır.
  subroutine test_single_torsional_element()
    type(torsional_element_t) :: element
    type(local_matrix_2x2) :: local_stiffness

    element = torsional_element_t( &
      id=1, node_i_id=10, node_j_id=20, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=0.0_dp)
    local_stiffness = get_local_stiffness(element)

    call assert_absolute_close( &
      local_stiffness%value(1, 1), 100.0_dp, absolute_tolerance, &
      "Tek eleman K11 katsayısı +k değil.")
    call assert_absolute_close( &
      local_stiffness%value(1, 2), -100.0_dp, absolute_tolerance, &
      "Tek eleman K12 katsayısı -k değil.")
    call assert_absolute_close( &
      local_stiffness%value(2, 1), -100.0_dp, absolute_tolerance, &
      "Tek eleman K21 katsayısı -k değil.")
    call assert_absolute_close( &
      local_stiffness%value(2, 2), 100.0_dp, absolute_tolerance, &
      "Tek eleman K22 katsayısı +k değil.")
    call assert_symmetric_matrix( &
      local_stiffness%value, absolute_tolerance, &
      "Tek eleman lokal rijitlik matrisi simetrik değil.")
  end subroutine test_single_torsional_element

  !> Serbest iki düğümlü sistemin ortak dönme rijit-cisim modunu doğrular.
  !! Model: theta=[1,1]^T ve K*theta=0. K [N*m/rad], theta [rad], kalıntı
  !! [N*m] birimindedir. Her iki düğüm serbest ve bağlantı lineerdir.
  subroutine test_rigid_body_mode()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: global_stiffness
    real(dp), allocatable :: stiffness_values(:, :)
    real(dp) :: rigid_rotation(2)
    real(dp) :: residual_nm(2)

    call build_two_node_system(system, 100.0_dp)
    call initialize_dof_map(mapping, system)
    global_stiffness = assemble_stiffness(system, mapping)
    stiffness_values = get_stiffness_matrix_values(global_stiffness)

    rigid_rotation = [1.0_dp, 1.0_dp]
    residual_nm = matmul(stiffness_values, rigid_rotation)
    call assert_vector_close( &
      residual_nm, [0.0_dp, 0.0_dp], absolute_tolerance, &
      "Serbest iki düğümlü sistem K*[1,1] rijit-cisim modunu korumuyor.")
  end subroutine test_rigid_body_mode

  !> İki lokal eleman katkısının üç düğümlü global K matrisinde toplandığını
  !! doğrular. Model: 30--k1=100--10--k2=200--70; fiziksel kimlikler matris
  !! indisi değildir. K [N*m/rad] olup K=sum(A_e^T*K_e*A_e) beklenir.
  subroutine test_global_assembly_regression()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: global_stiffness
    real(dp), allocatable :: actual(:, :)
    real(dp) :: expected(3, 3)

    call build_three_node_chain(system, first_node_constrained=.false.)
    call initialize_dof_map(mapping, system)
    global_stiffness = assemble_stiffness(system, mapping)
    actual = get_stiffness_matrix_values(global_stiffness)
    expected = reshape([ &
      100.0_dp, -100.0_dp, 0.0_dp, &
      -100.0_dp, 300.0_dp, -200.0_dp, &
      0.0_dp, -200.0_dp, 200.0_dp], [3, 3])

    call assert_matrix_close( &
      actual, expected, absolute_tolerance, &
      "Üç düğümlü global K assembly analitik referansla uyuşmuyor.")
    call assert_symmetric_matrix( &
      actual, absolute_tolerance, &
      "Üç düğümlü global K matrisi simetrik değil.")

    if (lookup_equation_id(mapping, 30) /= 1 .or. &
        lookup_equation_id(mapping, 10) /= 2 .or. &
        lookup_equation_id(mapping, 70) /= 3) then
      error stop "Global assembly fiziksel node ID ile equation ID ayrımını bozdu."
    end if
  end subroutine test_global_assembly_regression

  !> İki ataletli serbest-serbest TVD'nin analitik elastik frekansını mevcut
  !! solver ve assembled M/K matrisleriyle bağımsız olarak çapraz doğrular.
  !! Model: J1=0.1, J2=0.2 kg*m^2, K=1000 N*m/rad; sönümsüz, lineer,
  !! kütlesiz bağlantı ve bilinen phi_e=[1,-J1/J2]^T kabul edilir.
  subroutine test_two_inertia_analytic_benchmark()
    real(dp), parameter :: first_inertia_kg_m2 = 0.10_dp
    real(dp), parameter :: second_inertia_kg_m2 = 0.20_dp
    real(dp), parameter :: stiffness_nm_per_rad = 1000.0_dp
    real(dp), parameter :: expected_frequency_hz = &
      19.492420030841906_dp

    type(two_inertia_tvd_system_t) :: source_system
    type(two_inertia_modal_result_t) :: modal_result
    type(torsional_system_t) :: generalized_system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: global_stiffness
    type(mass_matrix_t) :: global_mass
    real(dp), allocatable :: stiffness_values(:, :)
    real(dp), allocatable :: mass_values(:, :)
    real(dp) :: analytical_frequency_hz
    real(dp) :: expected_omega_squared
    real(dp) :: elastic_mode(2)
    real(dp) :: rigid_mode(2)
    real(dp) :: stiffness_action(2)
    real(dp) :: inertia_action(2)
    real(dp) :: modal_residual(2)
    real(dp) :: residual_scale
    real(dp) :: rayleigh_omega_squared
    real(dp) :: rayleigh_frequency_hz
    real(dp) :: mass_orthogonality

    source_system = two_inertia_tvd_system_t( &
      hub_polar_inertia_kg_m2=first_inertia_kg_m2, &
      ring_polar_inertia_kg_m2=second_inertia_kg_m2, &
      storage_stiffness_nm_per_rad=stiffness_nm_per_rad)

    expected_omega_squared = stiffness_nm_per_rad * &
      (1.0_dp / first_inertia_kg_m2 + 1.0_dp / second_inertia_kg_m2)
    analytical_frequency_hz = &
      sqrt(expected_omega_squared) / (2.0_dp * pi)
    call assert_relative_close( &
      analytical_frequency_hz, expected_frequency_hz, relative_tolerance, &
      "İki ataletli analitik frekans sabit referansla uyuşmuyor.")

    modal_result = solve_free_free_two_inertia_modes(source_system)
    call assert_relative_close( &
      modal_result%elastic_frequency_hz, expected_frequency_hz, &
      relative_tolerance, &
      "Mevcut iki ataletli solver analitik elastik frekansı vermiyor.")

    generalized_system = build_generalized_two_inertia_system(source_system)
    call initialize_dof_map(mapping, generalized_system)
    global_stiffness = assemble_stiffness(generalized_system, mapping)
    global_mass = assemble_inertia(generalized_system, mapping)
    stiffness_values = get_stiffness_matrix_values(global_stiffness)
    mass_values = get_mass_matrix_values(global_mass)

    elastic_mode = [1.0_dp, &
      -first_inertia_kg_m2 / second_inertia_kg_m2]
    rigid_mode = [1.0_dp, 1.0_dp]
    stiffness_action = matmul(stiffness_values, elastic_mode)
    inertia_action = matmul(mass_values, elastic_mode)
    modal_residual = stiffness_action - &
      expected_omega_squared * inertia_action
    residual_scale = max( &
      1.0_dp, maxval(abs(stiffness_action)), &
      maxval(abs(expected_omega_squared * inertia_action)))
    if (maxval(abs(modal_residual)) / residual_scale > relative_tolerance) then
      error stop "Assembled M/K analitik iki-atalet elastik modunu sağlamıyor."
    end if

    rayleigh_omega_squared = &
      dot_product(elastic_mode, stiffness_action) / &
      dot_product(elastic_mode, inertia_action)
    rayleigh_frequency_hz = &
      sqrt(rayleigh_omega_squared) / (2.0_dp * pi)
    call assert_relative_close( &
      rayleigh_frequency_hz, expected_frequency_hz, relative_tolerance, &
      "Assembled M/K Rayleigh frekansı analitik referansla uyuşmuyor.")

    mass_orthogonality = dot_product( &
      rigid_mode, matmul(mass_values, elastic_mode))
    call assert_absolute_close( &
      mass_orthogonality, 0.0_dp, absolute_tolerance, &
      "İki ataletli rijit ve elastik modlar M-ortogonal değil.")
  end subroutine test_two_inertia_analytic_benchmark

  !> Fiziksel düğüm kimliklerinden aktif equation ID üretimini iki sınır
  !! koşulunda doğrular. Case A: [30,10,70]->[1,2,3]. Case B: ilk düğüm
  !! kısıtlıyken [30,10,70]->[0,1,2]. Kimlikler boyutsuzdur.
  subroutine test_dof_mapping_regression()
    type(torsional_system_t) :: all_free_system
    type(torsional_system_t) :: first_constrained_system
    type(dof_map_t) :: all_free_mapping
    type(dof_map_t) :: first_constrained_mapping

    call build_three_node_chain( &
      all_free_system, first_node_constrained=.false.)
    call initialize_dof_map(all_free_mapping, all_free_system)

    if (get_dof_map_entry_count(all_free_mapping) /= 3 .or. &
        get_active_dof_count(all_free_mapping) /= 3) then
      error stop "Case A DOF haritasının kayıt veya aktif denklem sayısı yanlış."
    end if
    if (lookup_equation_id(all_free_mapping, 30) /= 1 .or. &
        lookup_equation_id(all_free_mapping, 10) /= 2 .or. &
        lookup_equation_id(all_free_mapping, 70) /= 3) then
      error stop "Case A serbest düğüm equation ID sürekliliğini korumuyor."
    end if

    call build_three_node_chain( &
      first_constrained_system, first_node_constrained=.true.)
    call initialize_dof_map(first_constrained_mapping, &
      first_constrained_system)

    if (get_dof_map_entry_count(first_constrained_mapping) /= 3 .or. &
        get_active_dof_count(first_constrained_mapping) /= 2) then
      error stop "Case B DOF haritasının kayıt veya aktif denklem sayısı yanlış."
    end if
    if (lookup_equation_id(first_constrained_mapping, 30) /= 0 .or. &
        lookup_equation_id(first_constrained_mapping, 10) /= 1 .or. &
        lookup_equation_id(first_constrained_mapping, 70) /= 2) then
      error stop "Case B constraint sentinel veya equation ID sürekliliği yanlış."
    end if
  end subroutine test_dof_mapping_regression

  !> Serbest üç düğümlü zincirin global K kalite invariantlarını doğrular.
  !! Simetri ||K-K^T||_F ile, enerji ise U=1/2 theta^T K theta ile sınanır.
  !! K [N*m/rad], theta [rad], U [N*m]; bağlantılar lineer ve pasiftir.
  subroutine test_matrix_quality()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: global_stiffness
    real(dp), allocatable :: stiffness_values(:, :)
    real(dp) :: rotations(3, 4)
    integer :: case_index

    call build_three_node_chain(system, first_node_constrained=.false.)
    call initialize_dof_map(mapping, system)
    global_stiffness = assemble_stiffness(system, mapping)
    stiffness_values = get_stiffness_matrix_values(global_stiffness)

    call assert_symmetric_matrix( &
      stiffness_values, absolute_tolerance, &
      "Global K matrisinin Frobenius simetri normu toleransı aşıyor.")

    rotations(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
    rotations(:, 2) = [1.0_dp, 1.0_dp, 1.0_dp]
    rotations(:, 3) = [0.2_dp, -0.1_dp, 0.4_dp]
    rotations(:, 4) = [-0.3_dp, 0.7_dp, -0.2_dp]
    do case_index = 1, size(rotations, 2)
      call assert_nonnegative_strain_energy( &
        stiffness_values, rotations(:, case_index), absolute_tolerance, &
        "Global K matrisi seçili dönme halinde negatif enerji üretti.")
    end do
  end subroutine test_matrix_quality

  !> Serbest iki düğüm ve tek lineer torsional bağlantıdan test sistemi kurar.
  !! J=[0.1,0.2] kg*m^2, K girdi olarak [N*m/rad], c=0 N*m*s/rad'dır.
  subroutine build_two_node_system(system, stiffness_nm_per_rad)
    type(torsional_system_t), intent(out) :: system
    real(dp), intent(in) :: stiffness_nm_per_rad

    call add_torsional_node(system, torsional_node_t( &
      id=10, polar_inertia_kg_m2=0.10_dp, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=20, polar_inertia_kg_m2=0.20_dp, constrained=.false.))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=10, node_j_id=20, &
      stiffness_nm_per_rad=stiffness_nm_per_rad, &
      damping_nms_per_rad=0.0_dp))
  end subroutine build_two_node_system

  !> Fiziksel ID sırası [30,10,70] olan üç düğümlü zinciri kurar.
  !! J=[0.1,0.2,0.3] kg*m^2, K=[100,200] N*m/rad'dır. Opsiyonel fiziksel
  !! kısıt ilk düğümün equation ID değerini sıfır yapar; çıktı sistemdir.
  subroutine build_three_node_chain(system, first_node_constrained)
    type(torsional_system_t), intent(out) :: system
    logical, intent(in) :: first_node_constrained

    call add_torsional_node(system, torsional_node_t( &
      id=30, polar_inertia_kg_m2=0.10_dp, &
      constrained=first_node_constrained))
    call add_torsional_node(system, torsional_node_t( &
      id=10, polar_inertia_kg_m2=0.20_dp, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=70, polar_inertia_kg_m2=0.30_dp, constrained=.false.))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=30, node_j_id=10, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=0.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      id=2, node_i_id=10, node_j_id=70, &
      stiffness_nm_per_rad=200.0_dp, damping_nms_per_rad=0.0_dp))
  end subroutine build_three_node_chain

end program test_torsional_validation
