program test_modal_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_dof_types, only : TORSIONAL_ROTATION
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element
  use tms_constraint_types, only : constraint_t, FIXED_CONSTRAINT
  use tms_constraint_manager, only : constraint_manager_t, &
    initialize_constraint_manager, add_constraint
  use tms_reduced_system, only : reduced_torsional_system_t, &
    build_reduced_torsional_system, get_reduced_mass
  use tms_mass_matrix, only : mass_matrix_t, get_mass_matrix_values
  use tms_torsional_system, only : two_inertia_tvd_system_t, &
    two_inertia_modal_result_t, build_generalized_two_inertia_system, &
    solve_free_free_two_inertia_modes
  use tms_generalized_eigen_problem, only : generalized_eigen_problem_t, &
    create_generalized_eigen_problem
  use tms_eigen_solution, only : eigen_solution_t, get_eigenvalues
  use tms_generalized_eigen_solver, only : solve_generalized_eigen_problem
  use tms_modal_validation, only : AUTO_RIGID_TOLERANCE_MULTIPLIER, &
    calculate_auto_rigid_eigenvalue_tolerance, is_rigid_eigenvalue
  use tms_modal_result, only : modal_result_t, RIGID_MODE, ELASTIC_MODE, &
    get_modal_mode_count, get_modal_eigenvalues, &
    get_modal_frequencies_hz, get_modal_mode_classifications, &
    get_modal_relative_residuals, get_modal_reduced_mode_shapes, &
    get_modal_physical_mode_shapes, get_modal_mass_orthogonality_error, &
    get_modal_solver_backend_identity, is_modal_result_linear, &
    is_modal_result_undamped, is_modal_result_frozen_property
  use tms_modal_analysis, only : modal_analysis_options_t, &
    analyze_reduced_torsional_system
  implicit none

  real(dp), parameter :: relative_tolerance = 1.0e-10_dp
  real(dp), parameter :: absolute_tolerance = 1.0e-11_dp
  real(dp), parameter :: diagnostic_tolerance = 1.0e-10_dp
  character(len=64) :: validation_case

  ! Geçersiz modal-analysis vakaları CTest tarafından ayrı WILL_FAIL
  ! süreçlerinde çalıştırılır. Bilinmeyen selector normal çıkmalıdır.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_fixed_one_dof_modal_solution()
  call test_free_free_two_inertia_cross_validation()
  call test_three_node_free_free_chain()
  call test_constrained_three_node_chain()
  call test_multiple_rigid_body_modes()
  call test_scale_aware_rigid_classification()
  call test_small_negative_roundoff_policy()
  call test_partial_spectrum_tolerance_contract()
  call test_requested_mode_count()

  print *, "V0.5 reduced-system modal analiz doğrulamaları başarılı."

contains

  !> Fixed---k---J modelinin tek elastic eigenpair'ini uçtan uca doğrular.
  !! Beklenen lambda=k/J [1/s^2] ve f=sqrt(k/J)/(2*pi) [Hz]'dir. Test mass
  !! normalization, relative residual, backend-neutral metadata ve
  !! phi=P*phi_r recovery sonucunda fixed physical bileşenin sıfır kalmasını
  !! üretim API'leri üzerinden sınar.
  subroutine test_fixed_one_dof_modal_solution()
    real(dp), parameter :: stiffness = 1000.0_dp
    real(dp), parameter :: free_inertia = 0.20_dp
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_result_t) :: result
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: frequencies(:)
    real(dp), allocatable :: reduced_modes(:, :)
    real(dp), allocatable :: physical_modes(:, :)
    real(dp), allocatable :: residuals(:)
    integer, allocatable :: classifications(:)
    character(len=:), allocatable :: backend

    call build_two_node_system(system, 0.10_dp, free_inertia, stiffness)
    call initialize_constraint_manager(manager)
    call add_constraint(manager, make_fixed_constraint(1, 1), system)
    reduced_system = build_reduced_torsional_system(system, manager)
    result = analyze_reduced_torsional_system(reduced_system)

    eigenvalues = get_modal_eigenvalues(result)
    frequencies = get_modal_frequencies_hz(result)
    classifications = get_modal_mode_classifications(result)
    reduced_modes = get_modal_reduced_mode_shapes(result)
    physical_modes = get_modal_physical_mode_shapes(result)
    residuals = get_modal_relative_residuals(result)
    backend = get_modal_solver_backend_identity(result)

    if (get_modal_mode_count(result) /= 1 .or. &
        classifications(1) /= ELASTIC_MODE) then
      error stop "Fixed 1-DOF sistem tek elastic mode üretmedi."
    end if
    call assert_relative_close( &
      eigenvalues(1), stiffness / free_inertia, relative_tolerance, &
      "Fixed 1-DOF eigenvalue yanlış.")
    call assert_relative_close( &
      frequencies(1), sqrt(stiffness / free_inertia) / (2.0_dp * pi), &
      relative_tolerance, "Fixed 1-DOF doğal frekansı yanlış.")
    call assert_absolute_close( &
      free_inertia * reduced_modes(1, 1)**2, 1.0_dp, &
      diagnostic_tolerance, "Fixed 1-DOF mode M-normalized değil.")
    call assert_absolute_close(physical_modes(1, 1), 0.0_dp, &
      absolute_tolerance, "Fixed physical mode constrained DOF'ta sıfır değil.")
    if (abs(physical_modes(2, 1) - reduced_modes(1, 1)) > &
        absolute_tolerance) then
      error stop "Fixed mode physical recovery active genliği korumadı."
    end if
    call assert_diagnostics(result, residuals)
    if (index(backend, "DSYGV") == 0) then
      error stop "Modal sonuç reference DSYGV backend metadata'sını taşımıyor."
    end if
    if (.not. is_modal_result_linear(result) .or. &
        .not. is_modal_result_undamped(result) .or. &
        .not. is_modal_result_frozen_property(result)) then
      error stop "Modal sonuç linear/undamped/frozen-property etiketi taşımıyor."
    end if
  end subroutine test_fixed_one_dof_modal_solution

  !> Serbest-serbest iki ataletli sistemi mevcut analitik solver ile çapraz
  !! doğrular. Beklenen lambda=[0,k(1/Jh+1/Jr)] ve bir rigid/bir elastic
  !! mode'dur. DSYGV mass-normalized ve sign-arbitrary mode shape ürettiği için
  !! test exact bileşen değil M-ağırlıklı modal assurance criterion kullanır.
  subroutine test_free_free_two_inertia_cross_validation()
    real(dp), parameter :: hub_inertia = 0.10_dp
    real(dp), parameter :: ring_inertia = 0.20_dp
    real(dp), parameter :: stiffness = 1000.0_dp
    type(two_inertia_tvd_system_t) :: source_system
    type(two_inertia_modal_result_t) :: analytical_result
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_result_t) :: result
    real(dp) :: mass(2, 2)
    real(dp) :: expected_rigid(2)
    real(dp) :: expected_elastic(2)
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: frequencies(:)
    real(dp), allocatable :: modes(:, :)
    real(dp), allocatable :: physical_modes(:, :)
    real(dp), allocatable :: residuals(:)
    integer, allocatable :: classifications(:)

    source_system%hub_polar_inertia_kg_m2 = hub_inertia
    source_system%ring_polar_inertia_kg_m2 = ring_inertia
    source_system%storage_stiffness_nm_per_rad = stiffness
    source_system%loss_stiffness_nm_per_rad = 0.0_dp
    source_system%loss_factor = 0.0_dp
    source_system%material_reference_frequency_hz = 100.0_dp
    source_system%material_temperature_k = 293.15_dp

    analytical_result = solve_free_free_two_inertia_modes(source_system)
    system = build_generalized_two_inertia_system(source_system)
    call initialize_constraint_manager(manager)
    reduced_system = build_reduced_torsional_system(system, manager)
    result = analyze_reduced_torsional_system(reduced_system)

    eigenvalues = get_modal_eigenvalues(result)
    frequencies = get_modal_frequencies_hz(result)
    classifications = get_modal_mode_classifications(result)
    modes = get_modal_reduced_mode_shapes(result)
    physical_modes = get_modal_physical_mode_shapes(result)
    residuals = get_modal_relative_residuals(result)
    mass = diagonal_matrix([hub_inertia, ring_inertia])
    expected_rigid = [1.0_dp, 1.0_dp]
    expected_elastic = [1.0_dp, -hub_inertia / ring_inertia]

    if (get_modal_mode_count(result) /= 2 .or. &
        any(classifications /= [RIGID_MODE, ELASTIC_MODE])) then
      error stop "İki ataletli sistem rigid/elastic mode sırasını üretmedi."
    end if
    call assert_absolute_close(eigenvalues(1), 0.0_dp, &
      absolute_tolerance, "İki ataletli rigid eigenvalue sıfır değil.")
    call assert_relative_close(eigenvalues(2), stiffness * &
      (1.0_dp / hub_inertia + 1.0_dp / ring_inertia), &
      relative_tolerance, "İki ataletli elastic eigenvalue yanlış.")
    call assert_relative_close( &
      frequencies(2), analytical_result%elastic_frequency_hz, &
      relative_tolerance, "Generalized ve analitik iki-atalet frekansı farklı.")
    call assert_absolute_close( &
      mass_weighted_mac(modes(:, 1), expected_rigid, mass), &
      1.0_dp, diagnostic_tolerance, "Rigid mode eigenspace uyuşmuyor.")
    call assert_absolute_close( &
      mass_weighted_mac(modes(:, 2), expected_elastic, mass), &
      1.0_dp, diagnostic_tolerance, "Elastic mode eigenspace uyuşmuyor.")
    call assert_matrix_close(physical_modes, modes, absolute_tolerance, &
      "Constraint bulunmayan modelde physical ve reduced mode farklı.")
    call assert_diagnostics(result, residuals)
  end subroutine test_free_free_two_inertia_cross_validation

  !> Eş J ve k değerli üç-node free-free zincirde analitik
  !! lambda=[0,k/J,3k/J] sırasını, tek rigid mode'u, iki elastic mode'u,
  !! residual ve M-orthogonality sonuçlarını doğrular.
  subroutine test_three_node_free_free_chain()
    real(dp), parameter :: inertia = 2.0_dp
    real(dp), parameter :: stiffness = 30.0_dp
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_result_t) :: result
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: residuals(:)
    integer, allocatable :: classifications(:)

    call build_three_node_chain(system, inertia, stiffness)
    call initialize_constraint_manager(manager)
    reduced_system = build_reduced_torsional_system(system, manager)
    result = analyze_reduced_torsional_system(reduced_system)
    eigenvalues = get_modal_eigenvalues(result)
    classifications = get_modal_mode_classifications(result)
    residuals = get_modal_relative_residuals(result)

    call assert_absolute_close(eigenvalues(1), 0.0_dp, &
      absolute_tolerance, "Üç-node zincir rigid eigenvalue yanlış.")
    call assert_relative_close(eigenvalues(2), stiffness / inertia, &
      relative_tolerance, "Üç-node zincir ikinci eigenvalue yanlış.")
    call assert_relative_close(eigenvalues(3), 3.0_dp * stiffness / inertia, &
      relative_tolerance, "Üç-node zincir üçüncü eigenvalue yanlış.")
    if (count(classifications == RIGID_MODE) /= 1 .or. &
        count(classifications == ELASTIC_MODE) /= 2) then
      error stop "Üç-node zincir mode sınıflandırması yanlış."
    end if
    call assert_diagnostics(result, residuals)
  end subroutine test_three_node_free_free_chain

  !> İlk düğümü fixed eş üç-node zincirin reduced analitik eigenvalue'larını
  !! ve physical recovery'de constrained ilk bileşenin her mode için sıfır
  !! kaldığını doğrular. Beklenen lambda=(k/J)(3+-sqrt(5))/2 [1/s^2]'dir.
  subroutine test_constrained_three_node_chain()
    real(dp), parameter :: inertia = 2.0_dp
    real(dp), parameter :: stiffness = 30.0_dp
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_result_t) :: result
    real(dp) :: expected_eigenvalues(2)
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: frequencies(:)
    real(dp), allocatable :: physical_modes(:, :)
    real(dp), allocatable :: residuals(:)
    integer, allocatable :: classifications(:)

    call build_three_node_chain(system, inertia, stiffness)
    call initialize_constraint_manager(manager)
    call add_constraint(manager, make_fixed_constraint(1, 1), system)
    reduced_system = build_reduced_torsional_system(system, manager)
    result = analyze_reduced_torsional_system(reduced_system)
    eigenvalues = get_modal_eigenvalues(result)
    frequencies = get_modal_frequencies_hz(result)
    classifications = get_modal_mode_classifications(result)
    physical_modes = get_modal_physical_mode_shapes(result)
    residuals = get_modal_relative_residuals(result)

    expected_eigenvalues = stiffness / inertia * &
      [3.0_dp - sqrt(5.0_dp), 3.0_dp + sqrt(5.0_dp)] / 2.0_dp
    call assert_relative_close(eigenvalues(1), expected_eigenvalues(1), &
      relative_tolerance, "Constrained zincir ilk eigenvalue yanlış.")
    call assert_relative_close(eigenvalues(2), expected_eigenvalues(2), &
      relative_tolerance, "Constrained zincir ikinci eigenvalue yanlış.")
    call assert_relative_close(frequencies(1), &
      sqrt(expected_eigenvalues(1)) / (2.0_dp * pi), &
      relative_tolerance, "Constrained zincir ilk frekans yanlış.")
    call assert_relative_close(frequencies(2), &
      sqrt(expected_eigenvalues(2)) / (2.0_dp * pi), &
      relative_tolerance, "Constrained zincir ikinci frekans yanlış.")
    if (any(classifications /= ELASTIC_MODE)) then
      error stop "Constrained zincir bütün modları elastic sınıflandırmadı."
    end if
    if (maxval(abs(physical_modes(1, :))) > absolute_tolerance) then
      error stop "Constrained zincir physical mode recovery ilk DOF'u sıfırlamadı."
    end if
    call assert_diagnostics(result, residuals)
  end subroutine test_constrained_three_node_chain

  !> İki ayrık free-free iki-node alt sistemin iki bağımsız rigid-body mode
  !! ürettiğini doğrular. Individual rigid vectors unique olmadığı için exact
  !! mode shape kıyaslanmaz; sınıf sayısı, residual ve M-orthogonality sınanır.
  subroutine test_multiple_rigid_body_modes()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_result_t) :: result
    real(dp), allocatable :: residuals(:)
    integer, allocatable :: classifications(:)

    call build_disconnected_free_system(system)
    call initialize_constraint_manager(manager)
    reduced_system = build_reduced_torsional_system(system, manager)
    result = analyze_reduced_torsional_system(reduced_system)
    classifications = get_modal_mode_classifications(result)
    residuals = get_modal_relative_residuals(result)

    if (count(classifications == RIGID_MODE) /= 2 .or. &
        count(classifications == ELASTIC_MODE) /= 2) then
      error stop "Ayrık alt sistemler iki rigid ve iki elastic mode üretmedi."
    end if
    call assert_diagnostics(result, residuals)
  end subroutine test_multiple_rigid_body_modes

  !> Aynı iki-node free-free fiziğini 1e-12, 1 ve 1e12 rijitlik ölçeklerinde
  !! çözer. Sabit Hz eşiği kullanılmadığını; her ölçekte bir rigid ve bir
  !! elastic mode sınıfının korunduğunu doğrular.
  subroutine test_scale_aware_rigid_classification()
    real(dp), parameter :: scales(3) = [1.0e-12_dp, 1.0_dp, 1.0e12_dp]
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_result_t) :: result
    integer, allocatable :: classifications(:)
    integer :: scale_index

    do scale_index = 1, size(scales)
      call build_two_node_system(system, 1.0_dp, 1.0_dp, scales(scale_index))
      call initialize_constraint_manager(manager)
      reduced_system = build_reduced_torsional_system(system, manager)
      result = analyze_reduced_torsional_system(reduced_system)
      classifications = get_modal_mode_classifications(result)
      if (count(classifications == RIGID_MODE) /= 1 .or. &
          count(classifications == ELASTIC_MODE) /= 1) then
        error stop "AUTO rijit sınıflandırması model ölçeğine göre değişti."
      end if
    end do
  end subroutine test_scale_aware_rigid_classification

  !> Pozitif spectral scale yanında AUTO tolerans içindeki küçük negatif
  !! eigenvalue'nun DSYGV çözümünden sonra rigid kabul edildiğini doğrular.
  !! Bu sentetik regression sabit frekans eşiği veya production formül kopyası
  !! kullanmaz; üretim tolerance ve classification yordamlarını çağırır.
  subroutine test_small_negative_roundoff_policy()
    type(generalized_eigen_problem_t) :: problem
    type(eigen_solution_t) :: solution
    real(dp) :: stiffness(2, 2)
    real(dp) :: mass(2, 2)
    real(dp) :: tolerance
    real(dp), allocatable :: eigenvalues(:)

    stiffness = diagonal_matrix([ &
      -50.0_dp * epsilon(1.0_dp), 1.0_dp])
    mass = identity_matrix(2)
    problem = create_generalized_eigen_problem(stiffness, mass)
    solution = solve_generalized_eigen_problem(problem)
    eigenvalues = get_eigenvalues(solution)
    tolerance = calculate_auto_rigid_eigenvalue_tolerance( &
      stiffness, mass, eigenvalues, AUTO_RIGID_TOLERANCE_MULTIPLIER)

    if (.not. is_rigid_eigenvalue(eigenvalues(1), tolerance)) then
      error stop "AUTO tolerans içindeki küçük negatif eigenvalue rigid değil."
    end if
    if (is_rigid_eigenvalue(eigenvalues(2), tolerance)) then
      error stop "Pozitif elastic eigenvalue yanlışlıkla rigid sınıflandırıldı."
    end if
  end subroutine test_small_negative_roundoff_policy

  !> Gelecekteki partial-spectrum backend'in 1<=m<=n özdeğer döndürmesini
  !! ortak AUTO tolerans katmanının kabul ettiğini doğrular. V0.5 DSYGV tüm
  !! spectrum'u çözmeye devam eder; test yalnız backend-neutral post-processing
  !! sözleşmesini full-spectrum varsayımına kilitlememeyi amaçlar.
  subroutine test_partial_spectrum_tolerance_contract()
    real(dp) :: stiffness(3, 3)
    real(dp) :: mass(3, 3)
    real(dp) :: partial_eigenvalues(2)
    real(dp) :: tolerance

    stiffness = diagonal_matrix([0.0_dp, 10.0_dp, 30.0_dp])
    mass = identity_matrix(3)
    partial_eigenvalues = [0.0_dp, 10.0_dp]
    tolerance = calculate_auto_rigid_eigenvalue_tolerance( &
      stiffness, mass, partial_eigenvalues, &
      AUTO_RIGID_TOLERANCE_MULTIPLIER)

    if (tolerance <= 0.0_dp .or. &
        .not. is_rigid_eigenvalue(partial_eigenvalues(1), tolerance) .or. &
        is_rigid_eigenvalue(partial_eigenvalues(2), tolerance)) then
      error stop "Partial-spectrum AUTO rijit tolerans sözleşmesi başarısız."
    end if
  end subroutine test_partial_spectrum_tolerance_contract

  !> requested_mode_count seçeneğinin DSYGV bütün spectrum'u çözdükten sonra
  !! ortak post-processing katmanında artan sıradaki ilk iki mode'u birlikte
  !! seçtiğini doğrular. Public API belirli backend'in extraction davranışına
  !! kilitlenmemelidir.
  subroutine test_requested_mode_count()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_analysis_options_t) :: options
    type(modal_result_t) :: result
    real(dp), allocatable :: eigenvalues(:)

    call build_three_node_chain(system, 2.0_dp, 30.0_dp)
    call initialize_constraint_manager(manager)
    reduced_system = build_reduced_torsional_system(system, manager)
    options%requested_mode_count = 2
    result = analyze_reduced_torsional_system(reduced_system, options)
    eigenvalues = get_modal_eigenvalues(result)

    if (get_modal_mode_count(result) /= 2) then
      error stop "Requested mode count sonuç sözleşmesine uygulanmadı."
    end if
    call assert_absolute_close(eigenvalues(1), 0.0_dp, &
      absolute_tolerance, "Requested spectrum rigid mode'u korumadı.")
    call assert_relative_close(eigenvalues(2), 15.0_dp, &
      relative_tolerance, "Requested spectrum artan eigenpair sırasını bozdu.")
  end subroutine test_requested_mode_count

  !> Modal-analysis hata yollarını production ve synthetic girdilerle çalıştırır.
  !! Fully constrained 0x0 sistem LAPACK'e gitmemeli; significant negative
  !! lambda clean instability tanısı üretmeli; geçersiz options reddedilmelidir.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: reduced_system
    type(modal_analysis_options_t) :: options
    type(modal_result_t) :: result
    type(generalized_eigen_problem_t) :: problem
    type(eigen_solution_t) :: solution
    real(dp) :: stiffness(2, 2)
    real(dp) :: mass(2, 2)
    real(dp) :: tolerance
    real(dp), allocatable :: eigenvalues(:)
    logical :: unexpected_rigid

    select case (case_name)
      case ("fully_constrained")
        call build_two_node_system(system, 0.10_dp, 0.20_dp, 100.0_dp)
        call initialize_constraint_manager(manager)
        call add_constraint(manager, make_fixed_constraint(1, 1), system)
        call add_constraint(manager, make_fixed_constraint(2, 2), system)
        reduced_system = build_reduced_torsional_system(system, manager)
        result = analyze_reduced_torsional_system(reduced_system)
      case ("significant_negative_eigenvalue")
        stiffness = diagonal_matrix([-1.0_dp, 1.0_dp])
        mass = identity_matrix(2)
        problem = create_generalized_eigen_problem(stiffness, mass)
        solution = solve_generalized_eigen_problem(problem)
        eigenvalues = get_eigenvalues(solution)
        tolerance = calculate_auto_rigid_eigenvalue_tolerance( &
          stiffness, mass, eigenvalues, AUTO_RIGID_TOLERANCE_MULTIPLIER)
        unexpected_rigid = is_rigid_eigenvalue(eigenvalues(1), tolerance)
        if (unexpected_rigid) print *, "Negatif eigenvalue yanlış rigid kabul edildi."
      case ("negative_requested_mode_count")
        call build_two_node_system(system, 1.0_dp, 1.0_dp, 10.0_dp)
        call initialize_constraint_manager(manager)
        reduced_system = build_reduced_torsional_system(system, manager)
        options%requested_mode_count = -1
        result = analyze_reduced_torsional_system(reduced_system, options)
      case ("excessive_requested_mode_count")
        call build_two_node_system(system, 1.0_dp, 1.0_dp, 10.0_dp)
        call initialize_constraint_manager(manager)
        reduced_system = build_reduced_torsional_system(system, manager)
        options%requested_mode_count = 3
        result = analyze_reduced_torsional_system(reduced_system, options)
      case ("nonfinite_rigid_multiplier")
        call build_two_node_system(system, 1.0_dp, 1.0_dp, 10.0_dp)
        call initialize_constraint_manager(manager)
        reduced_system = build_reduced_torsional_system(system, manager)
        options%auto_rigid_tolerance_multiplier = &
          ieee_value(0.0_dp, ieee_quiet_nan)
        result = analyze_reduced_torsional_system(reduced_system, options)
      case default
        print *, "Bilinmeyen modal doğrulama selector'ı: ", case_name
        return
    end select

    print *, "Geçersiz modal-analysis girdisi beklenmedik biçimde kabul edildi."
  end subroutine exercise_invalid_case

  !> Homojen fixed torsional constraint kaydı üretir. Kimlikler boyutsuz,
  !! prescribed dönme değeri sıfır [rad]'dır.
  pure function make_fixed_constraint(constraint_id, node_id) result(constraint)
    integer, intent(in) :: constraint_id
    integer, intent(in) :: node_id
    type(constraint_t) :: constraint

    constraint = constraint_t( &
      constraint_id=constraint_id, node_id=node_id, &
      dof_type=TORSIONAL_ROTATION, value=0.0_dp, &
      constraint_type=FIXED_CONSTRAINT)
  end function make_fixed_constraint

  !> Verilen iki pozitif polar atalet [kg*m^2] ve bağlantı rijitliği
  !! [N*m/rad] ile iki-node free torsional sistem kurar.
  subroutine build_two_node_system(system, first_inertia, second_inertia, stiffness)
    type(torsional_system_t), intent(out) :: system
    real(dp), intent(in) :: first_inertia
    real(dp), intent(in) :: second_inertia
    real(dp), intent(in) :: stiffness

    call add_torsional_node(system, torsional_node_t( &
      id=1, polar_inertia_kg_m2=first_inertia, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=2, polar_inertia_kg_m2=second_inertia, constrained=.false.))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=stiffness, damping_nms_per_rad=0.0_dp))
  end subroutine build_two_node_system

  !> Eş polar atalet J [kg*m^2] ve eş torsional rijitlik k [N*m/rad] ile
  !! 1---2---3 free-free zincirini kurar.
  subroutine build_three_node_chain(system, inertia, stiffness)
    type(torsional_system_t), intent(out) :: system
    real(dp), intent(in) :: inertia
    real(dp), intent(in) :: stiffness
    integer :: node_id

    do node_id = 1, 3
      call add_torsional_node(system, torsional_node_t( &
        id=node_id, polar_inertia_kg_m2=inertia, constrained=.false.))
    end do
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=stiffness, damping_nms_per_rad=0.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      id=2, node_i_id=2, node_j_id=3, &
      stiffness_nm_per_rad=stiffness, damping_nms_per_rad=0.0_dp))
  end subroutine build_three_node_chain

  !> Dört eş ataletli düğümde yalnız 1---2 ve 3---4 bağlantılarını kurar.
  !! Ayrık iki free-free alt sistemin her biri bir rigid-body mode taşır.
  subroutine build_disconnected_free_system(system)
    type(torsional_system_t), intent(out) :: system
    integer :: node_id

    do node_id = 1, 4
      call add_torsional_node(system, torsional_node_t( &
        id=node_id, polar_inertia_kg_m2=1.0_dp, constrained=.false.))
    end do
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=10.0_dp, damping_nms_per_rad=0.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      id=2, node_i_id=3, node_j_id=4, &
      stiffness_nm_per_rad=20.0_dp, damping_nms_per_rad=0.0_dp))
  end subroutine build_disconnected_free_system

  !> Modal result residual ve M-orthogonality tanılarını finite üst sınırlarla
  !! doğrular. Residual ve ortogonallik boyutsuzdur.
  subroutine assert_diagnostics(result, residuals)
    type(modal_result_t), intent(in) :: result
    real(dp), intent(in) :: residuals(:)

    if (.not. all(ieee_is_finite(residuals)) .or. &
        any(residuals > diagnostic_tolerance)) then
      error stop "Modal relative residual mühendislik toleransını aşıyor."
    end if
    if (.not. ieee_is_finite(get_modal_mass_orthogonality_error(result)) .or. &
        get_modal_mass_orthogonality_error(result) > diagnostic_tolerance) then
      error stop "Modal M-orthogonality hatası toleransı aşıyor."
    end if
  end subroutine assert_diagnostics

  !> İki mode shape arasındaki sign ve scale invariant M-ağırlıklı MAC
  !! değerini hesaplar. Sonuç [0,1] aralığında boyutsuzdur.
  function mass_weighted_mac(actual, expected, mass) result(mac)
    real(dp), intent(in) :: actual(:)
    real(dp), intent(in) :: expected(:)
    real(dp), intent(in) :: mass(:, :)
    real(dp) :: mac
    real(dp) :: actual_norm
    real(dp) :: cross_inner_product
    real(dp) :: expected_norm

    actual_norm = dot_product(actual, matmul(mass, actual))
    expected_norm = dot_product(expected, matmul(mass, expected))
    cross_inner_product = dot_product(actual, matmul(mass, expected))
    if (actual_norm <= 0.0_dp .or. expected_norm <= 0.0_dp) then
      error stop "M-ağırlıklı MAC pozitif modal normlar gerektirir."
    end if
    mac = cross_inner_product**2 / (actual_norm * expected_norm)
  end function mass_weighted_mac

  !> Pozitif referanslı skaler sonuçları boyutsuz bağıl toleransla sınar.
  subroutine assert_relative_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp .or. &
        abs(expected) <= tiny(1.0_dp)) then
      error stop "Modal bağıl assertion girdileri geçersiz."
    end if
    if (abs(actual - expected) / abs(expected) > tolerance) error stop message
  end subroutine assert_relative_close

  !> Aynı fiziksel birimdeki skaler sonuçları mutlak toleransla sınar.
  subroutine assert_absolute_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) then
      error stop "Modal mutlak assertion girdileri geçersiz."
    end if
    if (abs(actual - expected) > tolerance) error stop message
  end subroutine assert_absolute_close

  !> Eş boyutlu sonlu matrislerin maksimum katsayı farkını sınar.
  subroutine assert_matrix_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:, :)
    real(dp), intent(in) :: expected(:, :)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (any(shape(actual) /= shape(expected))) error stop message
    if (.not. all(ieee_is_finite(actual)) .or. &
        .not. all(ieee_is_finite(expected)) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) then
      error stop "Modal matris assertion girdileri geçersiz."
    end if
    if (maxval(abs(actual - expected)) > tolerance) error stop message
  end subroutine assert_matrix_close

  !> Verilen köşegen değerlerinden n x n diagonal matris oluşturur.
  pure function diagonal_matrix(diagonal) result(matrix)
    real(dp), intent(in) :: diagonal(:)
    real(dp) :: matrix(size(diagonal), size(diagonal))
    integer :: index

    matrix = 0.0_dp
    do index = 1, size(diagonal)
      matrix(index, index) = diagonal(index)
    end do
  end function diagonal_matrix

  !> Verilen boyutta boyutsuz birim matris oluşturur.
  pure function identity_matrix(order) result(matrix)
    integer, intent(in) :: order
    real(dp) :: matrix(order, order)
    integer :: index

    matrix = 0.0_dp
    do index = 1, order
      matrix(index, index) = 1.0_dp
    end do
  end function identity_matrix

end program test_modal_analysis
