program test_constraint_foundation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan
  use tms_kinds, only : dp
  use tms_dof_types, only : TORSIONAL_ROTATION
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element
  use tms_dof_map, only : dof_map_t, initialize_full_dof_map
  use tms_constraint_types, only : constraint_t, FIXED_CONSTRAINT, &
    PRESCRIBED_VALUE_CONSTRAINT
  use tms_constraint_manager, only : constraint_manager_t, active_dof_map_t, &
    active_dof_map_entry_t, initialize_constraint_manager, add_constraint, &
    initialize_constraint_manager_from_system, &
    generate_active_dof_map, &
    get_active_equation_count, get_physical_dof_count, &
    get_active_dof_map_entry, lookup_active_equation_id, &
    lookup_full_equation_id_for_active_map
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    get_stiffness_matrix_values
  use tms_mass_matrix, only : mass_matrix_t, get_mass_matrix_values
  use tms_matrix_assembly, only : assemble_full_stiffness, &
    assemble_full_inertia
  use tms_reduced_system, only : reduced_torsional_system_t, &
    build_reduced_torsional_system, get_reduced_stiffness, &
    get_reduced_mass, get_reduced_active_dof_map, &
    get_reduced_active_dof_count, recover_physical_state, recover_mode_shape
  implicit none

  real(dp), parameter :: tolerance = 1.0e-12_dp
  character(len=64) :: validation_case

  ! Geçersiz girdiler üretim API'sinde error stop oluşturmalıdır. Her seçim
  ! CTest tarafından bağımsız WILL_FAIL vakası olarak çalıştırılır.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_two_node_fixed_reduction()
  call test_three_node_chain_reduction()
  call test_fully_constrained_reduction()
  call test_prescribed_value_recovery()
  call test_legacy_constraint_adapter()
  call test_node_permutation_equivalence()

  print *, "V0.4 constraint, reduction ve recovery doğrulamaları başarılı."

contains

  !> İki düğümlü tek elemanda full assembly'nin constraint-blind kaldığını ve
  !! node 10 fixed sonrası K_r=[k], M_r=[J_20] olduğunu doğrular.
  !! Model: k=100 N*m/rad, J=[0.1,0.2] kg*m^2. Direct elimination yalnız
  !! active node 20 satır/sütununu tutar; test üretim fonksiyonlarını kullanır.
  subroutine test_two_node_fixed_reduction()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(dof_map_t) :: full_mapping
    type(stiffness_matrix_t) :: full_stiffness
    type(mass_matrix_t) :: full_mass
    type(reduced_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: reduced_stiffness
    type(mass_matrix_t) :: reduced_mass
    type(active_dof_map_t) :: active_mapping
    real(dp), allocatable :: full_k(:, :)
    real(dp), allocatable :: full_m(:, :)
    real(dp), allocatable :: reduced_k(:, :)
    real(dp), allocatable :: reduced_m(:, :)
    real(dp), allocatable :: full_state(:)
    real(dp), allocatable :: full_mode(:)

    call build_two_node_system(system, 10, 20, .false., .false.)
    call initialize_full_dof_map(full_mapping, system)
    full_stiffness = assemble_full_stiffness(system, full_mapping)
    full_mass = assemble_full_inertia(system, full_mapping)
    full_k = get_stiffness_matrix_values(full_stiffness)
    full_m = get_mass_matrix_values(full_mass)

    call assert_matrix_close(full_k, reshape([ &
      100.0_dp, -100.0_dp, -100.0_dp, 100.0_dp], [2, 2]), &
      "Constraint öncesi full K matrisi korunmadı.")
    call assert_matrix_close(full_m, reshape([ &
      0.10_dp, 0.0_dp, 0.0_dp, 0.20_dp], [2, 2]), &
      "Constraint öncesi full M matrisi korunmadı.")

    call initialize_constraint_manager(manager)
    call add_constraint(manager, make_fixed_constraint(1, 10), system)
    reduced_system = build_reduced_torsional_system(system, manager)
    reduced_stiffness = get_reduced_stiffness(reduced_system)
    reduced_mass = get_reduced_mass(reduced_system)
    active_mapping = get_reduced_active_dof_map(reduced_system)
    reduced_k = get_stiffness_matrix_values(reduced_stiffness)
    reduced_m = get_mass_matrix_values(reduced_mass)

    if (get_reduced_active_dof_count(reduced_system) /= 1 .or. &
        get_active_equation_count(active_mapping) /= 1 .or. &
        get_physical_dof_count(active_mapping) /= 2) then
      error stop "İki düğümlü fixed modelin physical/active DOF sayısı yanlış."
    end if
    if (lookup_full_equation_id_for_active_map( &
        active_mapping, 10, TORSIONAL_ROTATION) /= 1 .or. &
        lookup_full_equation_id_for_active_map( &
        active_mapping, 20, TORSIONAL_ROTATION) /= 2 .or. &
        lookup_active_equation_id( &
        active_mapping, 10, TORSIONAL_ROTATION) /= 0 .or. &
        lookup_active_equation_id( &
        active_mapping, 20, TORSIONAL_ROTATION) /= 1) then
      error stop "Physical, full ve active equation eşlemesi yanlış."
    end if

    call assert_matrix_close( &
      reduced_k, reshape([100.0_dp], [1, 1]), "Reduced K=[k] değil.")
    call assert_matrix_close( &
      reduced_m, reshape([0.20_dp], [1, 1]), "Reduced M=[J] değil.")

    full_state = recover_physical_state(reduced_system, [0.25_dp])
    full_mode = recover_mode_shape(reduced_system, [0.25_dp])
    call assert_vector_close( &
      full_state, [0.0_dp, 0.25_dp], "Fixed state recovery yanlış.")
    call assert_vector_close( &
      full_mode, [0.0_dp, 0.25_dp], "Fixed modal recovery yanlış.")
  end subroutine test_two_node_fixed_reduction

  !> Üç düğümlü zincirde ilk DOF elendikten sonra active node 20 ve 30 için
  !! principal K/M alt matrislerini ve [1,0.5] -> [0,1,0.5] recovery'sini
  !! doğrular. K birimi [N*m/rad], M birimi [kg*m^2]'dir.
  subroutine test_three_node_chain_reduction()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: stiffness
    type(mass_matrix_t) :: mass
    real(dp), allocatable :: actual_k(:, :)
    real(dp), allocatable :: actual_m(:, :)
    real(dp), allocatable :: full_state(:)

    call build_three_node_chain(system, [10, 20, 30])
    call initialize_constraint_manager(manager)
    call add_constraint(manager, make_fixed_constraint(1, 10), system)
    reduced_system = build_reduced_torsional_system(system, manager)
    stiffness = get_reduced_stiffness(reduced_system)
    mass = get_reduced_mass(reduced_system)
    actual_k = get_stiffness_matrix_values(stiffness)
    actual_m = get_mass_matrix_values(mass)

    call assert_matrix_close(actual_k, reshape([ &
      300.0_dp, -200.0_dp, -200.0_dp, 200.0_dp], [2, 2]), &
      "Üç düğümlü zincirin reduced K matrisi yanlış.")
    call assert_matrix_close(actual_m, reshape([ &
      0.20_dp, 0.0_dp, 0.0_dp, 0.30_dp], [2, 2]), &
      "Üç düğümlü zincirin reduced M matrisi yanlış.")

    if (maxval(abs(actual_k - transpose(actual_k))) > tolerance) then
      error stop "Reduced K simetriyi korumuyor."
    end if
    full_state = recover_physical_state(reduced_system, [1.0_dp, 0.5_dp])
    call assert_vector_close(full_state, [0.0_dp, 1.0_dp, 0.5_dp], &
      "Üç düğümlü result recovery yanlış.")
  end subroutine test_three_node_chain_reduction

  !> Bütün fiziksel DOF'ların fixed olduğu geçerli sistemde active DOF
  !! sayısının sıfır, reduced K/M boyutunun 0x0 ve full recovery'nin iki sıfır
  !! dönme [rad] olmasını doğrular.
  subroutine test_fully_constrained_reduction()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: stiffness
    type(mass_matrix_t) :: mass
    real(dp), allocatable :: reduced_k(:, :)
    real(dp), allocatable :: reduced_m(:, :)
    real(dp), allocatable :: empty(:)
    real(dp), allocatable :: full_state(:)
    real(dp), allocatable :: full_mode(:)

    call build_two_node_system(system, 10, 20, .false., .false.)
    call initialize_constraint_manager(manager)
    call add_constraint(manager, make_fixed_constraint(1, 10), system)
    call add_constraint(manager, make_fixed_constraint(2, 20), system)
    reduced_system = build_reduced_torsional_system(system, manager)
    stiffness = get_reduced_stiffness(reduced_system)
    mass = get_reduced_mass(reduced_system)
    reduced_k = get_stiffness_matrix_values(stiffness)
    reduced_m = get_mass_matrix_values(mass)

    if (get_reduced_active_dof_count(reduced_system) /= 0 .or. &
        any(shape(reduced_k) /= [0, 0]) .or. &
        any(shape(reduced_m) /= [0, 0])) then
      error stop "Tam constraint sistemi 0x0 reduced K/M üretmedi."
    end if

    allocate(empty(0))
    full_state = recover_physical_state(reduced_system, empty)
    full_mode = recover_mode_shape(reduced_system, empty)
    call assert_vector_close( &
      full_state, [0.0_dp, 0.0_dp], "Tam constraint state recovery yanlış.")
    call assert_vector_close( &
      full_mode, [0.0_dp, 0.0_dp], "Tam constraint modal recovery yanlış.")
  end subroutine test_fully_constrained_reduction

  !> Prescribed sıfır olmayan dönmenin matris reduction'dan bağımsız metadata
  !! olarak saklandığını doğrular. q=P*q_r+q_p state recovery'de 0.125 rad
  !! eklenir; phi=P*phi_r modal recovery'de constrained bileşen sıfır kalır.
  subroutine test_prescribed_value_recovery()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(constraint_t) :: prescribed
    type(stiffness_matrix_t) :: stiffness
    type(mass_matrix_t) :: mass
    real(dp), allocatable :: stiffness_values(:, :)
    real(dp), allocatable :: mass_values(:, :)
    real(dp), allocatable :: state(:)
    real(dp), allocatable :: mode(:)

    call build_two_node_system(system, 10, 20, .false., .false.)
    call initialize_constraint_manager(manager)
    prescribed = constraint_t( &
      constraint_id=1, node_id=10, dof_type=TORSIONAL_ROTATION, &
      value=0.125_dp, constraint_type=PRESCRIBED_VALUE_CONSTRAINT)
    call add_constraint(manager, prescribed, system)
    reduced_system = build_reduced_torsional_system(system, manager)
    stiffness = get_reduced_stiffness(reduced_system)
    mass = get_reduced_mass(reduced_system)
    stiffness_values = get_stiffness_matrix_values(stiffness)
    mass_values = get_mass_matrix_values(mass)

    ! Prescribed değer q_p yalnız recovery metadata'sıdır; V0.4 yük/RHS
    ! çözmediği için principal Kr/Mr fixed-zero indirgemesiyle aynı kalır.
    call assert_matrix_close( &
      stiffness_values, reshape([100.0_dp], [1, 1]), &
      "Prescribed değer reduced K matrisine sızdı.")
    call assert_matrix_close( &
      mass_values, reshape([0.20_dp], [1, 1]), &
      "Prescribed değer reduced M matrisine sızdı.")

    state = recover_physical_state(reduced_system, [0.50_dp])
    mode = recover_mode_shape(reduced_system, [0.50_dp])
    call assert_vector_close( &
      state, [0.125_dp, 0.50_dp], "Prescribed state recovery yanlış.")
    call assert_vector_close( &
      mode, [0.0_dp, 0.50_dp], "Prescribed modal recovery offset içeriyor.")
  end subroutine test_prescribed_value_recovery

  !> V0.3 node%constrained işaretinin explicit fixed-zero manager kaydına
  !! dönüştürülebildiğini doğrular. Adaptör sonrası yeni full->reduced yolun
  !! K_r=[100 N*m/rad] üretmesi eski public davranışı korur.
  subroutine test_legacy_constraint_adapter()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: stiffness
    real(dp), allocatable :: values(:, :)

    call build_two_node_system(system, 10, 20, .true., .false.)
    call initialize_constraint_manager_from_system(manager, system)
    reduced_system = build_reduced_torsional_system(system, manager)
    stiffness = get_reduced_stiffness(reduced_system)
    values = get_stiffness_matrix_values(stiffness)
    call assert_matrix_close(values, reshape([100.0_dp], [1, 1]), &
      "Legacy constraint adaptörü reduced K değerini korumadı.")
  end subroutine test_legacy_constraint_adapter

  !> Aynı fiziksel üç düğümlü zincirin farklı node ekleme sıralarında reduced
  !! K/M katsayılarının fiziksel (node_id,dof_type) anahtarına göre aynı
  !! kaldığını doğrular. Ham matris sırası değil, Q^T*A*Q permütasyon eşdeğerliği
  !! ve aynı fiziksel dönme alanının elastik enerjisi karşılaştırılır.
  subroutine test_node_permutation_equivalence()
    type(torsional_system_t) :: system_a
    type(torsional_system_t) :: system_b
    type(constraint_manager_t) :: manager_a
    type(constraint_manager_t) :: manager_b
    type(reduced_torsional_system_t) :: reduced_a
    type(reduced_torsional_system_t) :: reduced_b
    type(active_dof_map_t) :: map_a
    type(active_dof_map_t) :: map_b
    type(stiffness_matrix_t) :: stiffness_a
    type(stiffness_matrix_t) :: stiffness_b
    type(mass_matrix_t) :: mass_a
    type(mass_matrix_t) :: mass_b
    real(dp), allocatable :: k_a(:, :)
    real(dp), allocatable :: k_b(:, :)
    real(dp), allocatable :: m_a(:, :)
    real(dp), allocatable :: m_b(:, :)
    real(dp) :: theta_a(2)
    real(dp) :: theta_b(2)
    real(dp) :: energy_a
    real(dp) :: energy_b
    real(dp), allocatable :: full_state_a(:)
    real(dp), allocatable :: full_state_b(:)
    integer :: active_nodes(2)
    integer :: physical_nodes(3)
    integer :: equation_a_i
    integer :: equation_a_j
    integer :: equation_b_i
    integer :: equation_b_j
    integer :: i
    integer :: j

    call build_three_node_chain(system_a, [10, 20, 30])
    call build_three_node_chain(system_b, [30, 10, 20])
    call initialize_constraint_manager(manager_a)
    call initialize_constraint_manager(manager_b)
    call add_constraint(manager_a, make_fixed_constraint(1, 10), system_a)
    call add_constraint(manager_b, make_fixed_constraint(1, 10), system_b)
    reduced_a = build_reduced_torsional_system(system_a, manager_a)
    reduced_b = build_reduced_torsional_system(system_b, manager_b)
    map_a = get_reduced_active_dof_map(reduced_a)
    map_b = get_reduced_active_dof_map(reduced_b)
    stiffness_a = get_reduced_stiffness(reduced_a)
    stiffness_b = get_reduced_stiffness(reduced_b)
    mass_a = get_reduced_mass(reduced_a)
    mass_b = get_reduced_mass(reduced_b)
    k_a = get_stiffness_matrix_values(stiffness_a)
    k_b = get_stiffness_matrix_values(stiffness_b)
    m_a = get_mass_matrix_values(mass_a)
    m_b = get_mass_matrix_values(mass_b)

    active_nodes = [20, 30]
    physical_nodes = [10, 20, 30]
    do i = 1, 2
      do j = 1, 2
        equation_a_i = lookup_active_equation_id( &
          map_a, active_nodes(i), TORSIONAL_ROTATION)
        equation_a_j = lookup_active_equation_id( &
          map_a, active_nodes(j), TORSIONAL_ROTATION)
        equation_b_i = lookup_active_equation_id( &
          map_b, active_nodes(i), TORSIONAL_ROTATION)
        equation_b_j = lookup_active_equation_id( &
          map_b, active_nodes(j), TORSIONAL_ROTATION)
        call assert_close(k_a(equation_a_i, equation_a_j), &
          k_b(equation_b_i, equation_b_j), &
          "Node permütasyonu fiziksel K katsayısını değiştirdi.")
        call assert_close(m_a(equation_a_i, equation_a_j), &
          m_b(equation_b_i, equation_b_j), &
          "Node permütasyonu fiziksel M katsayısını değiştirdi.")
      end do
    end do

    theta_a = 0.0_dp
    theta_b = 0.0_dp
    theta_a(lookup_active_equation_id( &
      map_a, 20, TORSIONAL_ROTATION)) = 1.0_dp
    theta_a(lookup_active_equation_id( &
      map_a, 30, TORSIONAL_ROTATION)) = 0.5_dp
    theta_b(lookup_active_equation_id( &
      map_b, 20, TORSIONAL_ROTATION)) = 1.0_dp
    theta_b(lookup_active_equation_id( &
      map_b, 30, TORSIONAL_ROTATION)) = 0.5_dp
    energy_a = 0.5_dp * dot_product(theta_a, matmul(k_a, theta_a))
    energy_b = 0.5_dp * dot_product(theta_b, matmul(k_b, theta_b))
    call assert_close( &
      energy_a, energy_b, "Node permütasyonu fiziksel enerjiyi değiştirdi.")

    full_state_a = recover_physical_state(reduced_a, theta_a)
    full_state_b = recover_physical_state(reduced_b, theta_b)
    do i = 1, 3
      call assert_close( &
        lookup_physical_state_value(full_state_a, map_a, &
          physical_nodes(i)), &
        lookup_physical_state_value(full_state_b, map_b, &
          physical_nodes(i)), &
        "Node permütasyonu physical result recovery değerini değiştirdi.")
    end do
  end subroutine test_node_permutation_equivalence

  !> Geçersiz constraint ve recovery girdilerinin üretim API'si tarafından
  !! reddedildiğini sınar. Kimlikler boyutsuz, value [rad] ve atalet/rjitlik
  !! birimleri test sistemi sözleşmesindeki SI değerleridir.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_system_t) :: system
    type(torsional_system_t) :: other_system
    type(constraint_manager_t) :: manager
    type(dof_map_t) :: full_mapping
    type(active_dof_map_t) :: active_mapping
    type(reduced_torsional_system_t) :: reduced_system
    type(constraint_t) :: constraint
    real(dp), allocatable :: unexpected_state(:)

    select case (case_name)
      case ("missing_node")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        call add_constraint(manager, make_fixed_constraint(1, 999), system)
      case ("duplicate_physical_dof")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        call add_constraint(manager, make_fixed_constraint(1, 10), system)
        call add_constraint(manager, make_fixed_constraint(2, 10), system)
      case ("duplicate_constraint_id")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        call add_constraint(manager, make_fixed_constraint(1, 10), system)
        call add_constraint(manager, make_fixed_constraint(1, 20), system)
      case ("invalid_dof_type")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        constraint = make_fixed_constraint(1, 10)
        constraint%dof_type = 999
        call add_constraint(manager, constraint, system)
      case ("invalid_constraint_type")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        constraint = make_fixed_constraint(1, 10)
        constraint%constraint_type = 999
        call add_constraint(manager, constraint, system)
      case ("nonfinite_value")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        constraint = constraint_t( &
          constraint_id=1, node_id=10, dof_type=TORSIONAL_ROTATION, &
          value=ieee_value(0.0_dp, ieee_quiet_nan), &
          constraint_type=PRESCRIBED_VALUE_CONSTRAINT)
        call add_constraint(manager, constraint, system)
      case ("fixed_nonzero_value")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        constraint = make_fixed_constraint(1, 10)
        constraint%value = 0.1_dp
        call add_constraint(manager, constraint, system)
      case ("legacy_constraint_mismatch")
        call build_two_node_system(system, 10, 20, .true., .false.)
        call initialize_constraint_manager(manager)
        reduced_system = build_reduced_torsional_system(system, manager)
      case ("foreign_full_dof_map")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call build_two_node_system( &
          other_system, 20, 10, .false., .false.)
        call initialize_full_dof_map(full_mapping, system)
        call initialize_constraint_manager(manager)
        call generate_active_dof_map( &
          active_mapping, manager, full_mapping, other_system)
      case ("invalid_recovery_size")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        call add_constraint(manager, make_fixed_constraint(1, 10), system)
        reduced_system = build_reduced_torsional_system(system, manager)
        unexpected_state = recover_physical_state( &
          reduced_system, [1.0_dp, 2.0_dp])
      case ("invalid_recovery_value")
        call build_two_node_system(system, 10, 20, .false., .false.)
        call initialize_constraint_manager(manager)
        call add_constraint(manager, make_fixed_constraint(1, 10), system)
        reduced_system = build_reduced_torsional_system(system, manager)
        unexpected_state = recover_physical_state(reduced_system, &
          [ieee_value(0.0_dp, ieee_quiet_nan)])
      case default
        print *, "Bilinmeyen validation seçimi: ", trim(case_name)
        return
    end select

    ! Üretim API'si geçersiz girdiyi kabul ederse normal çıkış yapılır;
    ! CTest WILL_FAIL özelliği böylece yanlış-pozitif başarı üretemez.
    print *, "Geçersiz constraint girdisi beklenmedik biçimde kabul edildi: ", &
      trim(case_name)
  end subroutine exercise_invalid_case

  !> Homojen fixed torsional constraint kaydı üretir. Girdiler boyutsuz
  !! constraint/node kimlikleri, çıktı theta=0 rad koşulunu taşıyan kayıttır.
  pure function make_fixed_constraint(constraint_id, node_id) &
      result(constraint)
    integer, intent(in) :: constraint_id
    integer, intent(in) :: node_id
    type(constraint_t) :: constraint

    constraint = constraint_t( &
      constraint_id=constraint_id, node_id=node_id, &
      dof_type=TORSIONAL_ROTATION, value=0.0_dp, &
      constraint_type=FIXED_CONSTRAINT)
  end function make_fixed_constraint

  !> k=100 N*m/rad tek eleman ve J=[0.1,0.2] kg*m^2 düğümlerinden
  !! doğrulanabilir iki-node torsional sistem kurar. Constraint bayrakları
  !! yalnız legacy uyumluluk/mismatch testlerinde kullanılır.
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

  !> Aynı fiziksel üç-node zinciri verilen ekleme sırasıyla kurar.
  !! Fiziksel kimliklere bağlı değerler J10/J20/J30=[0.1,0.2,0.3] kg*m^2;
  !! k10-20=100 ve k20-30=200 N*m/rad'dır. Bu ayrım permutation testinin
  !! topolojik eşdeğerliği matris indeksinden bağımsız sınamasını sağlar.
  subroutine build_three_node_chain(system, node_order)
    type(torsional_system_t), intent(out) :: system
    integer, intent(in) :: node_order(3)

    integer :: index
    real(dp) :: inertia

    do index = 1, 3
      select case (node_order(index))
        case (10)
          inertia = 0.10_dp
        case (20)
          inertia = 0.20_dp
        case (30)
          inertia = 0.30_dp
        case default
          error stop "Test zinciri yalnız 10, 20 ve 30 node ID kabul eder."
      end select
      call add_torsional_node(system, torsional_node_t( &
        id=node_order(index), polar_inertia_kg_m2=inertia, &
        constrained=.false.))
    end do

    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=10, node_j_id=20, &
      stiffness_nm_per_rad=100.0_dp, damping_nms_per_rad=0.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      id=2, node_i_id=20, node_j_id=30, &
      stiffness_nm_per_rad=200.0_dp, damping_nms_per_rad=0.0_dp))
  end subroutine build_three_node_chain

  !> Physical state vektöründeki değeri node/dof anahtarıyla döndürür.
  !! Girdi full state [rad], active mapping ve boyutsuz node kimliğidir; çıktı
  !! aynı fiziksel torsional dönme [rad] değeridir. Matris sırasına dayanmaz.
  function lookup_physical_state_value(full_state, mapping, node_id) &
      result(value)
    real(dp), intent(in) :: full_state(:)
    type(active_dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    real(dp) :: value

    type(active_dof_map_entry_t) :: entry
    integer :: entry_index

    if (size(full_state) /= get_physical_dof_count(mapping)) then
      error stop "Physical state boyutu active mapping ile uyumsuz."
    end if

    do entry_index = 1, get_physical_dof_count(mapping)
      entry = get_active_dof_map_entry(mapping, entry_index)
      if (entry%physical_dof%node_id == node_id .and. &
          entry%physical_dof%dof_type == TORSIONAL_ROTATION) then
        value = full_state(entry%physical_dof_id)
        return
      end if
    end do

    error stop "Physical state içinde istenen torsional DOF bulunamadı."
  end function lookup_physical_state_value

  !> İki aynı boyutlu fiziksel matrisin tüm SI katsayılarını ölçeğe duyarlı
  !! 1e-12 toleransıyla karşılaştırır. Matematiksel model:
  !! |a-b| <= tolerance*max(1,|b|).
  subroutine assert_matrix_close(actual, expected, message)
    real(dp), intent(in) :: actual(:, :)
    real(dp), intent(in) :: expected(:, :)
    character(len=*), intent(in) :: message

    if (any(shape(actual) /= shape(expected))) error stop message
    if (.not. all(ieee_is_finite(actual)) .or. &
        .not. all(ieee_is_finite(expected))) error stop message
    if (size(actual) == 0) return
    if (maxval(abs(actual - expected)) > &
        tolerance * max(1.0_dp, maxval(abs(expected)))) then
      error stop message
    end if
  end subroutine assert_matrix_close

  !> İki aynı boyutlu state veya mode vektörünü ölçeğe duyarlı 1e-12
  !! toleransıyla karşılaştırır. Girdi birimi çağıran teste göre [rad] veya
  !! modal genliktir; fizik formülünü yeniden hesaplamaz.
  subroutine assert_vector_close(actual, expected, message)
    real(dp), intent(in) :: actual(:)
    real(dp), intent(in) :: expected(:)
    character(len=*), intent(in) :: message

    if (size(actual) /= size(expected)) error stop message
    if (.not. all(ieee_is_finite(actual)) .or. &
        .not. all(ieee_is_finite(expected))) error stop message
    if (size(actual) == 0) return
    if (maxval(abs(actual - expected)) > &
        tolerance * max(1.0_dp, maxval(abs(expected)))) then
      error stop message
    end if
  end subroutine assert_vector_close

  !> İki skaler fiziksel değeri ölçeğe duyarlı 1e-12 toleransıyla sınar.
  !! Girdiler aynı SI birimindedir; sonuç üretmez, uyumsuzluk reddedilir.
  subroutine assert_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected)) error stop message
    if (abs(actual - expected) > &
        tolerance * max(1.0_dp, abs(expected))) error stop message
  end subroutine assert_close

end program test_constraint_foundation
