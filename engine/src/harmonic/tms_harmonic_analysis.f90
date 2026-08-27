module tms_harmonic_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
  use tms_kinds, only : dp
  use tms_generalized_torsional_system, only : torsional_system_t
  use tms_constraint_manager, only : constraint_manager_t, &
    active_dof_map_t, get_physical_dof_count
  use tms_stiffness_matrix, only : stiffness_matrix_t
  use tms_loss_stiffness_matrix, only : loss_stiffness_matrix_t
  use tms_damping_matrix, only : damping_matrix_t
  use tms_mass_matrix, only : mass_matrix_t
  use tms_reduced_dynamic_system, only : &
    reduced_dynamic_torsional_system_t, &
    build_reduced_dynamic_torsional_system, &
    get_reduced_dynamic_stiffness, &
    get_reduced_dynamic_loss_stiffness, &
    get_reduced_dynamic_damping, get_reduced_dynamic_mass, &
    get_reduced_dynamic_active_dof_map, &
    get_reduced_dynamic_active_dof_count, recover_harmonic_response
  use tms_harmonic_excitation, only : harmonic_excitation_t, &
    assemble_harmonic_load_vector
  use tms_dynamic_stiffness, only : dynamic_stiffness_matrix_t, &
    build_dynamic_stiffness, get_dynamic_stiffness_values
  use tms_complex_linear_problem, only : complex_linear_problem_t, &
    create_complex_linear_problem
  use tms_complex_linear_solution, only : complex_linear_solution_t, &
    get_complex_linear_solution_status, has_complex_linear_response, &
    get_complex_linear_response, &
    get_complex_linear_reciprocal_condition_number, &
    get_complex_linear_forward_error_bounds, &
    get_complex_linear_backward_errors, &
    get_complex_linear_relative_residuals, &
    get_complex_linear_backend_identity
  use tms_complex_linear_solver, only : solve_complex_linear_problem
  use tms_harmonic_response, only : harmonic_response_t, &
    create_harmonic_response, validate_harmonic_frequency_sweep
  implicit none
  private

  public :: analyze_harmonic_response

contains

  !> Genel torsional sistemin direct full-order frequency-domain response'unu
  !! explicit frequency grid üzerinde çözer.
  !!
  !! Fiziksel açıklama: Complex peak nodal torque uyarımı altında küçük
  !! genlikli steady-state torsional response hesaplanır. K', K'', C ve M
  !! sweep boyunca frozen/current değerlerinde tutulur.
  !! Matematiksel model: Her f>0 için omega=2*pi*f ve
  !! [K'_r-omega^2*M_r+i(K''_r+omega*C_r)]*theta_hat_r=T_hat_r çözülür;
  !! ardından theta_hat=P*theta_hat_r ile physical DOF uzayına açılır.
  !! Girdiler: torsional sistem, mevcut direct-elimination constraint manager,
  !! kesin artan frequency array [Hz] ve complex peak torque kayıtları [N*m].
  !! Çıktı: Status, RCOND, residual, FERR/BERR ve [rad] response taşıyan
  !! harmonic_response_t.
  !! Varsayımlar ve geçerlilik: Constrained harmonic DOF=0 kabul edilir;
  !! stored prescribed static offset eklenmez. Time-varying prescribed motion,
  !! material interpolation, mode superposition ve transient çözüm yoktur.
  !! Exact singular frequency bir analysis status'tür; sweep sürdürülür ve o
  !! noktada uydurma response oluşturulmaz. Geçersiz input error stop üretir.
  function analyze_harmonic_response( &
      system, manager, frequencies_hz, excitations) result(response)
    type(torsional_system_t), intent(in) :: system
    type(constraint_manager_t), intent(in) :: manager
    real(dp), intent(in) :: frequencies_hz(:)
    type(harmonic_excitation_t), intent(in) :: excitations(:)
    type(harmonic_response_t) :: response

    type(active_dof_map_t) :: active_mapping
    type(complex_linear_problem_t) :: linear_problem
    type(complex_linear_solution_t) :: linear_solution
    type(damping_matrix_t) :: reduced_damping
    type(dynamic_stiffness_matrix_t) :: dynamic_stiffness
    type(loss_stiffness_matrix_t) :: reduced_loss_stiffness
    type(mass_matrix_t) :: reduced_mass
    type(reduced_dynamic_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: reduced_storage_stiffness
    character(len=:), allocatable :: backend_identity
    complex(dp), allocatable :: coefficient_matrix(:, :)
    complex(dp), allocatable :: load_matrix(:, :)
    complex(dp), allocatable :: load_vector(:)
    complex(dp), allocatable :: physical_responses(:, :)
    complex(dp), allocatable :: reduced_responses(:, :)
    complex(dp), allocatable :: solution_values(:, :)
    integer, allocatable :: solution_statuses(:)
    real(dp), allocatable :: backward_errors(:)
    real(dp), allocatable :: forward_error_bounds(:)
    real(dp), allocatable :: reciprocal_condition_numbers(:)
    real(dp), allocatable :: relative_residuals(:)
    real(dp), allocatable :: solution_backward_errors(:)
    real(dp), allocatable :: solution_forward_error_bounds(:)
    real(dp), allocatable :: solution_relative_residuals(:)
    real(dp) :: quiet_nan
    integer :: active_dof_count
    integer :: frequency_count
    integer :: frequency_index
    integer :: physical_dof_count

    call validate_harmonic_frequency_sweep(frequencies_hz)
    reduced_system = build_reduced_dynamic_torsional_system(system, manager)
    active_dof_count = get_reduced_dynamic_active_dof_count(reduced_system)
    if (active_dof_count == 0) then
      error stop "Harmonic analysis için en az bir active DOF gereklidir."
    end if

    active_mapping = get_reduced_dynamic_active_dof_map(reduced_system)
    physical_dof_count = get_physical_dof_count(active_mapping)
    load_vector = assemble_harmonic_load_vector(active_mapping, excitations)
    allocate(load_matrix(active_dof_count, 1))
    load_matrix(:, 1) = load_vector

    reduced_storage_stiffness = &
      get_reduced_dynamic_stiffness(reduced_system)
    reduced_loss_stiffness = &
      get_reduced_dynamic_loss_stiffness(reduced_system)
    reduced_damping = get_reduced_dynamic_damping(reduced_system)
    reduced_mass = get_reduced_dynamic_mass(reduced_system)

    frequency_count = size(frequencies_hz)
    allocate(solution_statuses(frequency_count))
    allocate(reciprocal_condition_numbers(frequency_count))
    allocate(relative_residuals(frequency_count))
    allocate(forward_error_bounds(frequency_count))
    allocate(backward_errors(frequency_count))
    allocate(reduced_responses(active_dof_count, frequency_count))
    allocate(physical_responses(physical_dof_count, frequency_count))

    quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)
    solution_statuses = 0
    reciprocal_condition_numbers = 0.0_dp
    relative_residuals = quiet_nan
    forward_error_bounds = quiet_nan
    backward_errors = quiet_nan
    reduced_responses = cmplx(quiet_nan, quiet_nan, kind=dp)
    physical_responses = cmplx(quiet_nan, quiet_nan, kind=dp)

    do frequency_index = 1, frequency_count
      dynamic_stiffness = build_dynamic_stiffness( &
        reduced_storage_stiffness, reduced_loss_stiffness, reduced_damping, &
        reduced_mass, frequencies_hz(frequency_index))
      coefficient_matrix = get_dynamic_stiffness_values(dynamic_stiffness)
      linear_problem = create_complex_linear_problem( &
        coefficient_matrix, load_matrix)
      linear_solution = solve_complex_linear_problem(linear_problem)

      solution_statuses(frequency_index) = &
        get_complex_linear_solution_status(linear_solution)
      reciprocal_condition_numbers(frequency_index) = &
        get_complex_linear_reciprocal_condition_number(linear_solution)
      if (.not. allocated(backend_identity)) then
        backend_identity = get_complex_linear_backend_identity(linear_solution)
      end if

      if (.not. has_complex_linear_response(linear_solution)) cycle

      solution_values = get_complex_linear_response(linear_solution)
      reduced_responses(:, frequency_index) = solution_values(:, 1)
      physical_responses(:, frequency_index) = recover_harmonic_response( &
        reduced_system, reduced_responses(:, frequency_index))

      solution_forward_error_bounds = &
        get_complex_linear_forward_error_bounds(linear_solution)
      solution_backward_errors = &
        get_complex_linear_backward_errors(linear_solution)
      solution_relative_residuals = &
        get_complex_linear_relative_residuals(linear_solution)
      forward_error_bounds(frequency_index) = &
        solution_forward_error_bounds(1)
      backward_errors(frequency_index) = solution_backward_errors(1)
      relative_residuals(frequency_index) = solution_relative_residuals(1)
    end do

    response = create_harmonic_response( &
      frequencies_hz, solution_statuses, reciprocal_condition_numbers, &
      relative_residuals, forward_error_bounds, backward_errors, &
      reduced_responses, physical_responses, backend_identity)
  end function analyze_harmonic_response

end module tms_harmonic_analysis
