program test_matrix_assembly
  use tms_kinds, only : dp
  use tms_matrix_types, only : dense_matrix_t, initialize_dense_matrix, &
    add_dense_matrix_entry, get_dense_matrix_entry, &
    get_dense_matrix_row_count, get_dense_matrix_column_count
  use tms_local_matrix, only : local_matrix_2x2
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t, get_local_stiffness
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element
  use tms_dof_map, only : dof_map_t, dof_map_entry_t, &
    initialize_dof_map, lookup_equation_id, get_active_dof_count, &
    get_dof_map_entry_count, get_dof_map_entry
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    initialize_stiffness_matrix, add_local_stiffness, &
    get_stiffness_matrix_values
  use tms_mass_matrix, only : mass_matrix_t, initialize_mass_matrix, &
    add_nodal_inertia, get_mass_matrix_values
  use tms_matrix_assembly, only : assemble_stiffness, assemble_inertia
  implicit none

  real(dp), parameter :: tolerance = 1.0e-10_dp
  character(len=64) :: validation_case

  ! Boyut ve fizik önkoşulları aynı yürütülebilir dosyanın ayrı WILL_FAIL
  ! çağrılarıyla doğrudan üretim API'leri üzerinden sınanır.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_dense_matrix_container()
  call test_single_element_assembly()
  call test_three_node_chain_assembly()
  call test_constrained_dof_assembly()
  call test_fully_constrained_assembly()

  print *, "DOF eşleme ile global torsional M/K assembly doğrulandı."

contains

  !> Genel dense matris taşıyıcısının boyut ve katsayı erişimini doğrular.
  !! Girdiler/çıktılar test içinde sabittir; katsayı birimi genel taşıyıcıda
  !! üretici katmana bağlıdır. Fiziksel assembly hesabı yapmaz.
  subroutine test_dense_matrix_container()
    type(dense_matrix_t) :: matrix

    call initialize_dense_matrix(matrix, 2, 3)
    call add_dense_matrix_entry(matrix, 2, 3, 5.0_dp)

    if (get_dense_matrix_row_count(matrix) /= 2 .or. &
        get_dense_matrix_column_count(matrix) /= 3) then
      error stop "Dense matris boyutları doğru saklanmadı."
    end if

    if (abs(get_dense_matrix_entry(matrix, 2, 3) - 5.0_dp) > tolerance) then
      error stop "Dense matris katsayı katkısı doğru saklanmadı."
    end if
  end subroutine test_dense_matrix_container

  !> Tek elemanın lokal 2x2 rijitliğinin DOF haritası üzerinden aynı global
  !! matrise taşındığını doğrular.
  !! Model: k=100 N*m/rad, node ID'leri 10 ve 20, iki serbest torsional DOF.
  subroutine test_single_element_assembly()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(dof_map_entry_t) :: entry
    type(stiffness_matrix_t) :: global_stiffness
    real(dp), allocatable :: actual(:, :)
    real(dp) :: expected(2, 2)

    call build_two_node_system(system, 10, 20, .false., .false.)
    call initialize_dof_map(mapping, system)

    if (get_dof_map_entry_count(mapping) /= 2 .or. &
        get_active_dof_count(mapping) /= 2) then
      error stop "İki serbest düğüm için DOF haritası boyutu yanlış."
    end if

    entry = get_dof_map_entry(mapping, 1)
    if (entry%node_id /= 10 .or. entry%equation_id /= 1) then
      error stop "Fiziksel düğüm ile denklem kimliği ayrımı doğru kurulmadı."
    end if
    if (lookup_equation_id(mapping, 20) /= 2) then
      error stop "İkinci fiziksel düğüm doğru denklem kimliğine eşlenmedi."
    end if

    global_stiffness = assemble_stiffness(system, mapping)
    actual = get_stiffness_matrix_values(global_stiffness)
    expected = reshape([ &
      100.0_dp, -100.0_dp, &
      -100.0_dp, 100.0_dp], [2, 2])
    call assert_matrix_close( &
      actual, expected, "Tek eleman global K matrisi yanlış.")
  end subroutine test_single_element_assembly

  !> Üç düğümlü serbest zincirin global K ve yığılmış dönel M matrislerini
  !! analitik referanslarla doğrular.
  !! Model: 10--k1=100--20--k2=200--30; J=[0.1,0.2,0.3] kg*m^2.
  !! Serbest sistemde simetri, sıfır satır toplamı, rijit-cisim null modu ve
  !! negatif olmayan elastik enerji de sınanır.
  subroutine test_three_node_chain_assembly()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: global_stiffness
    type(mass_matrix_t) :: global_mass
    real(dp), allocatable :: stiffness_values(:, :)
    real(dp), allocatable :: mass_values(:, :)
    real(dp) :: expected_stiffness(3, 3)
    real(dp) :: expected_mass(3, 3)
    real(dp) :: rigid_rotation(3)
    real(dp) :: residual(3)
    real(dp) :: theta(3)
    real(dp) :: quadratic_form

    call build_three_node_chain(system)
    call initialize_dof_map(mapping, system)
    global_stiffness = assemble_stiffness(system, mapping)
    global_mass = assemble_inertia(system, mapping)
    stiffness_values = get_stiffness_matrix_values(global_stiffness)
    mass_values = get_mass_matrix_values(global_mass)

    expected_stiffness = reshape([ &
      100.0_dp, -100.0_dp, 0.0_dp, &
      -100.0_dp, 300.0_dp, -200.0_dp, &
      0.0_dp, -200.0_dp, 200.0_dp], [3, 3])
    expected_mass = reshape([ &
      0.10_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.20_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.30_dp], [3, 3])

    call assert_matrix_close( &
      stiffness_values, expected_stiffness, &
      "Üç düğümlü zincirin global K matrisi yanlış.")
    call assert_matrix_close( &
      mass_values, expected_mass, &
      "Global dönel atalet matrisi diagonal J değerlerini taşımıyor.")

    if (maxval(abs(stiffness_values - transpose(stiffness_values))) > &
        tolerance) then
      error stop "Global torsional K matrisi simetrik değil."
    end if

    if (maxval(abs(sum(stiffness_values, dim=2))) > tolerance) then
      error stop "Serbest global K matrisinin satır toplamları sıfır değil."
    end if

    rigid_rotation = 1.0_dp
    residual = matmul(stiffness_values, rigid_rotation)
    if (maxval(abs(residual)) > tolerance) then
      error stop "Serbest global K matrisi rijit-cisim modunu korumuyor."
    end if

    theta = [0.20_dp, -0.10_dp, 0.40_dp]
    quadratic_form = dot_product(theta, matmul(stiffness_values, theta))
    if (quadratic_form < -tolerance) then
      error stop "Global K matrisi negatif elastik enerji üretti."
    end if
  end subroutine test_three_node_chain_assembly

  !> Kısıtlı düğümün equation_id=0 ile homojen olarak elendiğini doğrular.
  !! Model: node 10 sabit, node 20 serbest, k=100 N*m/rad. İndirgenmiş
  !! K=[100 N*m/rad], M=[J20=0.2 kg*m^2] olmalıdır. Kısıtlı sistemde sıfır
  !! satır toplamı veya rijit-cisim modu beklenmez.
  subroutine test_constrained_dof_assembly()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: global_stiffness
    type(mass_matrix_t) :: global_mass
    real(dp), allocatable :: stiffness_values(:, :)
    real(dp), allocatable :: mass_values(:, :)

    call build_two_node_system(system, 10, 20, .true., .false.)
    call initialize_dof_map(mapping, system)

    if (lookup_equation_id(mapping, 10) /= 0 .or. &
        lookup_equation_id(mapping, 20) /= 1) then
      error stop "Kısıtlı ve serbest düğüm denklem kimlikleri yanlış."
    end if

    global_stiffness = assemble_stiffness(system, mapping)
    global_mass = assemble_inertia(system, mapping)
    stiffness_values = get_stiffness_matrix_values(global_stiffness)
    mass_values = get_mass_matrix_values(global_mass)

    if (any(shape(stiffness_values) /= [1, 1]) .or. &
        abs(stiffness_values(1, 1) - 100.0_dp) > tolerance) then
      error stop "Homojen kısıt sonrası indirgenmiş K matrisi yanlış."
    end if

    if (any(shape(mass_values) /= [1, 1]) .or. &
        abs(mass_values(1, 1) - 0.20_dp) > tolerance) then
      error stop "Kısıt sonrası indirgenmiş dönel M matrisi yanlış."
    end if
  end subroutine test_constrained_dof_assembly

  !> Tamamen kısıtlı geçerli topolojinin sıfır aktif DOF ve 0x0 global M/K
  !! matrisleri ürettiğini doğrular.
  subroutine test_fully_constrained_assembly()
    type(torsional_system_t) :: system
    type(dof_map_t) :: mapping
    type(stiffness_matrix_t) :: global_stiffness
    type(mass_matrix_t) :: global_mass
    real(dp), allocatable :: stiffness_values(:, :)
    real(dp), allocatable :: mass_values(:, :)

    call build_two_node_system(system, 10, 20, .true., .true.)
    call initialize_dof_map(mapping, system)
    if (get_active_dof_count(mapping) /= 0) then
      error stop "Tamamen kısıtlı sistemin aktif DOF sayısı sıfır değil."
    end if

    global_stiffness = assemble_stiffness(system, mapping)
    global_mass = assemble_inertia(system, mapping)
    stiffness_values = get_stiffness_matrix_values(global_stiffness)
    mass_values = get_mass_matrix_values(global_mass)

    if (any(shape(stiffness_values) /= [0, 0]) .or. &
        any(shape(mass_values) /= [0, 0])) then
      error stop "Tamamen kısıtlı sistem 0x0 global matris üretmedi."
    end if
  end subroutine test_fully_constrained_assembly

  !> İki düğüm ve k=100 N*m/rad değerli tek elemandan test sistemi kurar.
  !! Düğüm ataletleri sırasıyla 0.1 ve 0.2 kg*m^2'dir. Kısıt girdileri
  !! homojen sıfır dönme sınır koşulunu seçer; çıktı doğrulanabilir topolojidir.
  subroutine build_two_node_system( &
      system, first_id, second_id, first_constrained, second_constrained)
    type(torsional_system_t), intent(out) :: system
    integer, intent(in) :: first_id
    integer, intent(in) :: second_id
    logical, intent(in) :: first_constrained
    logical, intent(in) :: second_constrained

    call add_torsional_node(system, torsional_node_t( &
      id=first_id, polar_inertia_kg_m2=0.10_dp, &
      constrained=first_constrained))
    call add_torsional_node(system, torsional_node_t( &
      id=second_id, polar_inertia_kg_m2=0.20_dp, &
      constrained=second_constrained))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=first_id, node_j_id=second_id, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=0.0_dp))
  end subroutine build_two_node_system

  !> Üç serbest düğümlü analitik zinciri kurar.
  !! Node ID [10,20,30], J [0.1,0.2,0.3] kg*m^2 ve eleman rijitlikleri
  !! [100,200] N*m/rad değerlerindedir. Çıktı global assembly test sistemidir.
  subroutine build_three_node_chain(system)
    type(torsional_system_t), intent(out) :: system

    call add_torsional_node(system, torsional_node_t( &
      id=10, polar_inertia_kg_m2=0.10_dp, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=20, polar_inertia_kg_m2=0.20_dp, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=30, polar_inertia_kg_m2=0.30_dp, constrained=.false.))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=10, node_j_id=20, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=0.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      id=2, node_i_id=20, node_j_id=30, &
      stiffness_nm_per_rad=200.0_dp, damping_nms_per_rad=0.0_dp))
  end subroutine build_three_node_chain

  !> Gerçek matris ile aynı boyuttaki analitik referansın tüm katsayılarını
  !! mutlak 1e-10 toleransıyla karşılaştırır. Birimler çağıran fizik testinin
  !! K [N*m/rad] veya M [kg*m^2] sözleşmesinden gelir.
  subroutine assert_matrix_close(actual, expected, message)
    real(dp), intent(in) :: actual(:, :)
    real(dp), intent(in) :: expected(:, :)
    character(len=*), intent(in) :: message

    if (any(shape(actual) /= shape(expected))) then
      error stop message
    end if
    if (maxval(abs(actual - expected)) > tolerance) then
      error stop message
    end if
  end subroutine assert_matrix_close

  !> Dense matris, DOF arama ve global katkı API'lerinin geçersiz boyut veya
  !! fiziksel girdileri reddettiğini sınar. Her vaka error stop üretmelidir.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(dense_matrix_t) :: dense_matrix
    type(torsional_system_t) :: source_system
    type(torsional_system_t) :: other_system
    type(dof_map_t) :: mapping
    type(torsional_element_t) :: element
    type(local_matrix_2x2) :: local_stiffness
    type(stiffness_matrix_t) :: global_stiffness
    type(mass_matrix_t) :: global_mass
    integer :: equation_id

    select case (case_name)
      case ("negative_matrix_dimension")
        call initialize_dense_matrix(dense_matrix, -1, 2)
      case ("unknown_node_lookup")
        call build_two_node_system( &
          source_system, 10, 20, .false., .false.)
        call initialize_dof_map(mapping, source_system)
        equation_id = lookup_equation_id(mapping, 999)
        print *, "Bilinmeyen düğüm beklenmedik biçimde bulundu: ", equation_id
      case ("out_of_range_equation")
        element = torsional_element_t( &
          id=1, node_i_id=10, node_j_id=20, &
          stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=0.0_dp)
        local_stiffness = get_local_stiffness(element)
        call initialize_stiffness_matrix(global_stiffness, 1)
        call add_local_stiffness( &
          global_stiffness, [1, 2], local_stiffness)
      case ("nonpositive_inertia")
        call initialize_mass_matrix(global_mass, 1)
        call add_nodal_inertia(global_mass, 1, -0.1_dp)
      case ("incompatible_dof_map")
        call build_two_node_system( &
          source_system, 10, 20, .false., .false.)
        call build_two_node_system( &
          other_system, 30, 40, .false., .false.)
        call initialize_dof_map(mapping, source_system)
        global_stiffness = assemble_stiffness(other_system, mapping)
        print *, "Uyumsuz harita beklenmedik biçimde kabul edildi: ", &
          size(get_stiffness_matrix_values(global_stiffness), 1)
      case default
        error stop "Bilinmeyen geçersiz matrix assembly testi istendi."
    end select
  end subroutine exercise_invalid_case

end program test_matrix_assembly
