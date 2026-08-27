module tms_material_aware_harmonic_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
  use tms_kinds, only : dp
  use tms_dof_types, only : TORSIONAL_ROTATION
  use tms_generalized_torsional_system, only : torsional_system_t, &
    get_torsional_element_count, get_torsional_element
  use tms_torsional_element, only : torsional_element_t
  use tms_constraint_manager, only : constraint_manager_t, &
    active_dof_map_t, get_physical_dof_count, lookup_active_equation_id
  use tms_local_matrix, only : local_matrix_2x2
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    initialize_stiffness_matrix, add_local_stiffness
  use tms_loss_stiffness_matrix, only : loss_stiffness_matrix_t, &
    initialize_loss_stiffness_matrix, add_local_loss_stiffness
  use tms_damping_matrix, only : damping_matrix_t
  use tms_mass_matrix, only : mass_matrix_t
  use tms_reduced_dynamic_system, only : &
    reduced_dynamic_torsional_system_t, &
    build_reduced_dynamic_torsional_system, &
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
  use tms_dynamic_torsional_property_binding, only : &
    dynamic_torsional_property_binding_t, &
    dynamic_torsional_property_state_t, &
    validate_dynamic_torsional_property_binding, &
    evaluate_dynamic_torsional_property, get_dynamic_binding_element_id
  use tms_material_state_trace, only : material_binding_trace_t, &
    material_state_trace_t, create_material_binding_trace, &
    create_material_state_trace
  use tms_material_aware_harmonic_response, only : &
    material_aware_harmonic_response_t, &
    create_material_aware_harmonic_response
  implicit none
  private

  public :: analyze_material_aware_harmonic_response

contains

  !> Tabulated dynamic shear material binding'leriyle direct full-order
  !! frequency-domain torsional response çözer.
  !!
  !! Fiziksel model: exp(+i*omega*t) convention'ında passive G''>=0 ve
  !! K''>=0 kullanılır. Bound elemanlar için K'(f)=C_theta*G'(f),
  !! K''(f)=C_theta*G''(f); unbound elemanların K'/K'' ve bütün elemanların
  !! viscous c değerleri frequency-independent kalır.
  !! Matematiksel model: Her f için
  !! Z=K'_r(f)-omega^2*M_r+i*(K''_r(f)+omega*C_r) ve Z*theta=T çözülür.
  !! Bound elemanın stored nominal K'/K'' değeri yalnız frozen V0.6 yoluna
  !! aittir; bu yordamda base matristen çıkarılır ve dynamic değerle override
  !! edilir, ayrıca toplanmaz.
  !! Girdiler: SI birimli system, constraint manager, strictly increasing
  !! frequency grid [Hz], explicit operating temperature [K], complex peak
  !! torque excitations [N*m] ve en az bir binding. Çıktı: V0.6 harmonic
  !! response ile her binding/frequency için G*/K* ve dataset trace'i.
  !! Varsayımlar ve sınırlar: Provider domain/temperature/mode ve tüm sweep
  !! material state'leri ZSYSVX çağrısından önce prevalidate edilir. Topology,
  !! DOF map, constraints, M, C, base K ve C_theta yalnız bir kez hazırlanır.
  !! Modal, transient, extrapolation, TTS ve temperature interpolation yoktur.
  !! Ayrıntılar: docs/architecture/V0.7_dynamic_material_provider.md ve
  !! docs/mathematics/dynamic_modulus_interpolation.md.
  function analyze_material_aware_harmonic_response( &
      system, manager, frequencies_hz, operating_temperature_k, &
      excitations, bindings) result(response)
    type(torsional_system_t), intent(in) :: system
    type(constraint_manager_t), intent(in) :: manager
    real(dp), intent(in) :: frequencies_hz(:)
    real(dp), intent(in) :: operating_temperature_k
    type(harmonic_excitation_t), intent(in) :: excitations(:)
    type(dynamic_torsional_property_binding_t), intent(in) :: bindings(:)
    type(material_aware_harmonic_response_t) :: response

    type(active_dof_map_t) :: active_mapping
    type(complex_linear_problem_t) :: linear_problem
    type(complex_linear_solution_t) :: linear_solution
    type(damping_matrix_t) :: reduced_damping
    type(dynamic_stiffness_matrix_t) :: dynamic_stiffness
    type(harmonic_response_t) :: harmonic_response
    type(loss_stiffness_matrix_t) :: base_loss_stiffness
    type(loss_stiffness_matrix_t) :: frequency_loss_stiffness
    type(mass_matrix_t) :: reduced_mass
    type(reduced_dynamic_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: base_storage_stiffness
    type(stiffness_matrix_t) :: frequency_storage_stiffness
    type(torsional_element_t) :: current_element
    type(torsional_element_t), allocatable :: bound_elements(:)
    type(dynamic_torsional_property_state_t) :: property_state
    type(material_binding_trace_t), allocatable :: binding_traces(:)
    type(material_state_trace_t), allocatable :: state_traces(:, :)
    character(len=:), allocatable :: backend_identity
    complex(dp), allocatable :: coefficient_matrix(:, :)
    complex(dp), allocatable :: load_matrix(:, :)
    complex(dp), allocatable :: load_vector(:)
    complex(dp), allocatable :: physical_responses(:, :)
    complex(dp), allocatable :: reduced_responses(:, :)
    complex(dp), allocatable :: solution_values(:, :)
    integer, allocatable :: binding_equation_ids(:, :)
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
    integer :: binding_count
    integer :: binding_index
    integer :: frequency_count
    integer :: frequency_index
    integer :: element_equation_ids(2)
    integer :: element_index
    integer :: physical_dof_count

    call validate_harmonic_frequency_sweep(frequencies_hz)
    binding_count = size(bindings)
    if (binding_count == 0) then
      error stop "Material-aware harmonic analysis en az bir binding gerektirir."
    end if

    ! V0.6 reduction bir kez kullanılır; M/C ve recovery bağlamı frozen kalır.
    reduced_system = build_reduced_dynamic_torsional_system(system, manager)
    active_dof_count = get_reduced_dynamic_active_dof_count(reduced_system)
    if (active_dof_count == 0) then
      error stop "Material-aware harmonic analysis active DOF gerektirir."
    end if
    active_mapping = get_reduced_dynamic_active_dof_map(reduced_system)
    physical_dof_count = get_physical_dof_count(active_mapping)
    load_vector = assemble_harmonic_load_vector(active_mapping, excitations)
    allocate(load_matrix(active_dof_count, 1))
    load_matrix(:, 1) = load_vector

    call initialize_stiffness_matrix(base_storage_stiffness, active_dof_count)
    call initialize_loss_stiffness_matrix( &
      base_loss_stiffness, active_dof_count)
    reduced_damping = get_reduced_dynamic_damping(reduced_system)
    reduced_mass = get_reduced_dynamic_mass(reduced_system)

    allocate(bound_elements(binding_count))
    allocate(binding_equation_ids(2, binding_count))
    allocate(binding_traces(binding_count))
    do binding_index = 1, binding_count
      call validate_dynamic_torsional_property_binding( &
        bindings(binding_index))
      call require_unique_binding(bindings, binding_index)
      bound_elements(binding_index) = find_bound_element( &
        system, get_dynamic_binding_element_id(bindings(binding_index)))
      binding_equation_ids(:, binding_index) = [ &
        lookup_active_equation_id( &
          active_mapping, bound_elements(binding_index)%node_i_id, &
          TORSIONAL_ROTATION), &
        lookup_active_equation_id( &
          active_mapping, bound_elements(binding_index)%node_j_id, &
          TORSIONAL_ROTATION)]
      binding_traces(binding_index) = &
        create_material_binding_trace(bindings(binding_index))

    end do

    ! K'_base ve K''_base sıfırdan yalnız unbound elemanlarla assemble edilir.
    ! Bound stored nominal değerler base'e hiç girmez; frozen API'de ise mevcut
    ! assembly davranışı aynen sürer. Bütün elemanların viscous c katkısı yukarıda
    ! bir kez alınan reduced_damping içinde korunur.
    do element_index = 1, get_torsional_element_count(system)
      current_element = get_torsional_element(system, element_index)
      if (is_bound_element(bindings, current_element%id)) cycle
      element_equation_ids = [ &
        lookup_active_equation_id( &
          active_mapping, current_element%node_i_id, &
          TORSIONAL_ROTATION), &
        lookup_active_equation_id( &
          active_mapping, current_element%node_j_id, &
          TORSIONAL_ROTATION)]
      call add_local_stiffness( &
        base_storage_stiffness, element_equation_ids, &
        create_local_coefficient_matrix( &
          current_element%stiffness_nm_per_rad))
      call add_local_loss_stiffness( &
        base_loss_stiffness, element_equation_ids, &
        create_local_coefficient_matrix( &
          current_element%loss_stiffness_nm_per_rad))
    end do

    ! Gate: Bütün provider/binding x frequency matrisi herhangi bir complex
    ! solver çağrısından önce değerlendirilir. Domain hatası partial solve yapmaz.
    frequency_count = size(frequencies_hz)
    allocate(state_traces(binding_count, frequency_count))
    do binding_index = 1, binding_count
      do frequency_index = 1, frequency_count
        property_state = evaluate_dynamic_torsional_property( &
          bindings(binding_index), frequencies_hz(frequency_index), &
          operating_temperature_k)
        state_traces(binding_index, frequency_index) = &
          create_material_state_trace(property_state)
      end do
    end do

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
      frequency_storage_stiffness = base_storage_stiffness
      frequency_loss_stiffness = base_loss_stiffness
      do binding_index = 1, binding_count
        call add_local_stiffness( &
          frequency_storage_stiffness, &
          binding_equation_ids(:, binding_index), &
          create_local_coefficient_matrix( &
            state_traces(binding_index, frequency_index)% &
            storage_stiffness_nm_per_rad))
        call add_local_loss_stiffness( &
          frequency_loss_stiffness, &
          binding_equation_ids(:, binding_index), &
          create_local_coefficient_matrix( &
            state_traces(binding_index, frequency_index)% &
            loss_stiffness_nm_per_rad))
      end do

      dynamic_stiffness = build_dynamic_stiffness( &
        frequency_storage_stiffness, frequency_loss_stiffness, &
        reduced_damping, reduced_mass, frequencies_hz(frequency_index))
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

    harmonic_response = create_harmonic_response( &
      frequencies_hz, solution_statuses, reciprocal_condition_numbers, &
      relative_residuals, forward_error_bounds, backward_errors, &
      reduced_responses, physical_responses, backend_identity)
    response = create_material_aware_harmonic_response( &
      harmonic_response, binding_traces, state_traces)
  end function analyze_material_aware_harmonic_response

  !> k katsayısından [k,-k;-k,k] lokal torsional matrisini üretir.
  !! Girdi k, storage veya loss rijitliği için [N*m/rad]; çıktı aynı birimde
  !! simetrik, sıfır satır toplamlı local_matrix_2x2'dir.
  pure function create_local_coefficient_matrix(coefficient) result(matrix)
    real(dp), intent(in) :: coefficient
    type(local_matrix_2x2) :: matrix

    matrix%value(1, 1) = coefficient
    matrix%value(1, 2) = -coefficient
    matrix%value(2, 1) = -coefficient
    matrix%value(2, 2) = coefficient
  end function create_local_coefficient_matrix

  !> Binding element ID'sini sistem koleksiyonunda bulup bağımsız kopyasını
  !! döndürür. Bilinmeyen eleman ID clean error ile reddedilir.
  pure function find_bound_element(system, element_id) result(element)
    type(torsional_system_t), intent(in) :: system
    integer, intent(in) :: element_id
    type(torsional_element_t) :: element

    integer :: element_index

    do element_index = 1, get_torsional_element_count(system)
      element = get_torsional_element(system, element_index)
      if (element%id == element_id) return
    end do
    error stop "Dynamic material binding bilinmeyen torsional elemana ait."
  end function find_bound_element

  !> Her torsional elemanın en fazla bir authoritative dynamic binding
  !! taşıdığını doğrular. Fizik hesabı yapmaz; yinelenen element ID reddedilir.
  pure subroutine require_unique_binding(bindings, current_index)
    type(dynamic_torsional_property_binding_t), intent(in) :: bindings(:)
    integer, intent(in) :: current_index

    integer :: previous_index

    do previous_index = 1, current_index - 1
      if (get_dynamic_binding_element_id(bindings(previous_index)) == &
          get_dynamic_binding_element_id(bindings(current_index))) then
        error stop "Aynı torsional eleman için duplicate dynamic binding bulundu."
      end if
    end do
  end subroutine require_unique_binding

  !> Bir element ID'nin authoritative dynamic binding kümesinde bulunup
  !! bulunmadığını döndürür. Bu topoloji yardımcı yordamı fizik hesabı yapmaz.
  pure function is_bound_element(bindings, element_id) result(is_bound)
    type(dynamic_torsional_property_binding_t), intent(in) :: bindings(:)
    integer, intent(in) :: element_id
    logical :: is_bound

    integer :: binding_index

    is_bound = .false.
    do binding_index = 1, size(bindings)
      if (get_dynamic_binding_element_id(bindings(binding_index)) == &
          element_id) then
        is_bound = .true.
        return
      end if
    end do
  end function is_bound_element

end module tms_material_aware_harmonic_analysis
