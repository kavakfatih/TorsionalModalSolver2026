program test_thermorheological_harmonic_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_dof_types, only : TORSIONAL_ROTATION
  use tms_geometry, only : rubber_geometry_t
  use tms_material_frequency, only : material_frequency_point
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    DYNAMIC_DEFORMATION_MODE_SHEAR
  use tms_dynamic_modulus_provider, only : LINEAR_LOG_FREQUENCY
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t, &
    create_tabulated_dynamic_modulus_provider
  use tms_tabulated_temperature_shift, only : &
    tabulated_log10_shift_provider_t, &
    create_tabulated_temperature_shift_provider
  use tms_wlf_temperature_shift, only : wlf_temperature_shift_provider_t, &
    create_wlf_temperature_shift_provider
  use tms_arrhenius_temperature_shift, only : &
    arrhenius_temperature_shift_provider_t, &
    create_arrhenius_temperature_shift_provider
  use tms_temperature_shift_types, only : TABULATED_LOG10_SHIFT, &
    WLF_TEMPERATURE_SHIFT, ARRHENIUS_TEMPERATURE_SHIFT
  use tms_thermorheological_dynamic_modulus_provider, only : &
    thermorheological_dynamic_modulus_provider_t, &
    create_thermorheological_dynamic_modulus_provider
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element
  use tms_constraint_types, only : constraint_t, FIXED_CONSTRAINT
  use tms_constraint_manager, only : constraint_manager_t, &
    initialize_constraint_manager, add_constraint
  use tms_harmonic_excitation, only : harmonic_excitation_t
  use tms_harmonic_response, only : harmonic_response_t, &
    get_physical_complex_response, get_harmonic_solution_statuses, &
    get_harmonic_response_availability
  use tms_complex_linear_solution, only : COMPLEX_SOLVE_SINGULAR
  use tms_dynamic_torsional_property_binding, only : &
    dynamic_torsional_property_binding_t, &
    dynamic_torsional_property_state_t, &
    create_dynamic_torsional_property_binding, &
    evaluate_dynamic_torsional_property
  use tms_material_state_trace, only : material_binding_trace_t, &
    material_state_trace_t
  use tms_material_aware_harmonic_response, only : &
    material_aware_harmonic_response_t, &
    create_material_aware_harmonic_response, get_base_harmonic_response, &
    get_material_binding_count, get_material_binding_traces, &
    get_material_state_traces, get_material_state_trace
  use tms_material_aware_harmonic_analysis, only : &
    analyze_material_aware_harmonic_response
  implicit none

  real(dp), parameter :: tolerance = 5.0e-11_dp
  real(dp), parameter :: gas_constant_reference = 8.31446261815324_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_one_dof_thermorheological_chain()
  call test_multiple_shift_providers()
  call test_singular_trace_is_preserved()

  print *, "V0.8 thermorheological harmonic zinciri doğrulandı."

contains

  !> Existing analyze_material_aware_harmonic_response API'sini kullanarak
  !! T -> a_T -> f_r -> G*/K* -> Z -> theta tam zincirini fixed-hub 1-DOF
  !! sistemde bağımsız complex analytical sonuçla doğrular. Yeni harmonic
  !! orchestration gerekmediğini ve physical/reduced trace ayrımını kilitler.
  subroutine test_one_dof_thermorheological_chain()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    type(tabulated_log10_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: bindings(1)
    type(material_aware_harmonic_response_t) :: response
    type(harmonic_response_t) :: harmonic_response
    type(material_state_trace_t) :: trace
    type(harmonic_excitation_t) :: excitation(1)
    type(rubber_geometry_t) :: rubber
    complex(dp), allocatable :: physical(:, :)
    complex(dp) :: expected
    complex(dp) :: torque
    real(dp) :: c_theta
    real(dp) :: omega
    real(dp) :: frequencies(1)

    master_curve = make_master_curve( &
      "ONE-DOF-TTS", 293.15_dp, 2.0e6_dp, 4.0e6_dp, &
      0.2e6_dp, 0.4e6_dp)
    shift_provider = make_tabulated_shift()
    provider = create_thermorheological_dynamic_modulus_provider( &
      master_curve, shift_provider)
    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    bindings(1) = create_dynamic_torsional_property_binding( &
      10, provider, rubber)
    call build_fixed_system( &
      system, manager, 10, 9000.0_dp, 800.0_dp, 2.0_dp, 0.25_dp)
    torque = cmplx(3.0_dp, -1.0_dp, kind=dp)
    excitation(1) = make_excitation(2, torque)
    frequencies = [100.0_dp]

    response = analyze_material_aware_harmonic_response( &
      system, manager, frequencies, 313.15_dp, excitation, bindings)
    harmonic_response = get_base_harmonic_response(response)
    physical = get_physical_complex_response(harmonic_response)

    c_theta = annular_factor(rubber)
    omega = 2.0_dp*pi*frequencies(1)
    expected = torque/cmplx( &
      c_theta*3.0e6_dp-omega**2*0.25_dp, &
      c_theta*0.3e6_dp+omega*2.0_dp, kind=dp)
    call assert_complex_close(physical(2, 1), expected, tolerance, &
      "Thermorheological 1-DOF response analitik sonuçla uyuşmuyor.")

    trace = get_material_state_trace(response, 1, 1)
    call assert_close(trace%physical_frequency_hz, &
      100.0_dp, tolerance, "Harmonic trace physical f hatalı.")
    call assert_close(trace%lookup_frequency_hz, &
      10.0_dp, tolerance, "Harmonic trace reduced f_r hatalı.")
    call assert_close(trace%operating_temperature_k, &
      313.15_dp, tolerance, "Harmonic trace operating T hatalı.")
    call assert_close(trace%reference_temperature_k, &
      293.15_dp, tolerance, "Harmonic trace reference T hatalı.")
    call assert_close(trace%log10_a_t, &
      -1.0_dp, tolerance, "Harmonic trace log10(a_T) hatalı.")
    call assert_close(trace%a_t, 0.1_dp, tolerance, &
      "Harmonic trace a_T hatalı.")
    call assert_close(trace%storage_modulus_pa, &
      3.0e6_dp, tolerance, "Harmonic trace G' hatalı.")
    call assert_close(trace%loss_modulus_pa, &
      0.3e6_dp, tolerance, "Harmonic trace G'' hatalı.")
    call assert_close(trace%storage_stiffness_nm_per_rad, &
      c_theta*3.0e6_dp, tolerance, "Harmonic trace K' hatalı.")
    call assert_close(trace%loss_stiffness_nm_per_rad, &
      c_theta*0.3e6_dp, tolerance, "Harmonic trace K'' hatalı.")
    if (.not. trace%temperature_shift_applied .or. &
        trace%shift_model_kind /= TABULATED_LOG10_SHIFT .or. &
        .not. trace%shift_exact_temperature_point) then
      error stop "Thermorheological harmonic shift trace context'i eksik."
    end if
  end subroutine test_one_dof_thermorheological_chain

  !> Aynı torsional sistemde WLF ve Arrhenius provider'larının kendi T_ref,
  !! a_T ve reduced frequency değerleriyle bağımsız çalıştığını doğrular.
  !! Result existing 2-DOF complex matrix inverse ile karşılaştırılır.
  subroutine test_multiple_shift_providers()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(tabulated_dynamic_modulus_provider_t) :: master_wlf
    type(tabulated_dynamic_modulus_provider_t) :: master_arrhenius
    type(wlf_temperature_shift_provider_t) :: shift_wlf
    type(arrhenius_temperature_shift_provider_t) :: shift_arrhenius
    type(thermorheological_dynamic_modulus_provider_t) :: provider_wlf
    type(thermorheological_dynamic_modulus_provider_t) :: provider_arrhenius
    type(dynamic_torsional_property_binding_t) :: bindings(2)
    type(material_aware_harmonic_response_t) :: response
    type(harmonic_response_t) :: harmonic_response
    type(material_state_trace_t) :: trace_wlf
    type(material_state_trace_t) :: trace_arrhenius
    type(harmonic_excitation_t) :: excitation(1)
    type(rubber_geometry_t) :: rubber_wlf
    type(rubber_geometry_t) :: rubber_arrhenius
    complex(dp), allocatable :: physical(:, :)
    complex(dp) :: expected(2)
    complex(dp) :: torque
    real(dp) :: c_wlf
    real(dp) :: c_arrhenius
    real(dp) :: g_storage_wlf
    real(dp) :: g_loss_wlf
    real(dp) :: g_storage_arrhenius
    real(dp) :: g_loss_arrhenius
    real(dp) :: s_wlf
    real(dp) :: s_arrhenius
    real(dp) :: frequencies(1)

    master_wlf = make_master_curve( &
      "WLF-MASTER", 293.15_dp, 1.0e6_dp, 3.0e6_dp, &
      0.1e6_dp, 0.3e6_dp)
    master_arrhenius = make_master_curve( &
      "ARRHENIUS-MASTER", 298.15_dp, 0.5e6_dp, 1.5e6_dp, &
      0.05e6_dp, 0.15e6_dp)
    shift_wlf = create_wlf_temperature_shift_provider( &
      2.0_dp, 100.0_dp, 293.15_dp, 280.0_dp, 320.0_dp)
    shift_arrhenius = create_arrhenius_temperature_shift_provider( &
      12000.0_dp, 298.15_dp, 280.0_dp, 320.0_dp)
    provider_wlf = create_thermorheological_dynamic_modulus_provider( &
      master_wlf, shift_wlf)
    provider_arrhenius = create_thermorheological_dynamic_modulus_provider( &
      master_arrhenius, shift_arrhenius)

    rubber_wlf = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    rubber_arrhenius = rubber_geometry_t(0.015_dp, 0.04_dp, 0.012_dp)
    bindings(1) = create_dynamic_torsional_property_binding( &
      10, provider_wlf, rubber_wlf)
    bindings(2) = create_dynamic_torsional_property_binding( &
      20, provider_arrhenius, rubber_arrhenius)
    call build_three_node_system(system, manager)
    frequencies = [10.0_dp]
    torque = cmplx(2.0_dp, 0.4_dp, kind=dp)
    excitation(1) = make_excitation(3, torque)

    s_wlf = -2.0_dp*10.0_dp/(100.0_dp+10.0_dp)
    s_arrhenius = 12000.0_dp/(gas_constant_reference*log(10.0_dp)) * &
      (1.0_dp/303.15_dp-1.0_dp/298.15_dp)
    g_storage_wlf = log_interpolate( &
      1.0e6_dp, 3.0e6_dp, 1.0_dp+s_wlf)
    g_loss_wlf = log_interpolate( &
      0.1e6_dp, 0.3e6_dp, 1.0_dp+s_wlf)
    g_storage_arrhenius = log_interpolate( &
      0.5e6_dp, 1.5e6_dp, 1.0_dp+s_arrhenius)
    g_loss_arrhenius = log_interpolate( &
      0.05e6_dp, 0.15e6_dp, 1.0_dp+s_arrhenius)
    c_wlf = annular_factor(rubber_wlf)
    c_arrhenius = annular_factor(rubber_arrhenius)

    response = analyze_material_aware_harmonic_response( &
      system, manager, frequencies, 303.15_dp, excitation, bindings)
    harmonic_response = get_base_harmonic_response(response)
    physical = get_physical_complex_response(harmonic_response)
    expected = solve_two_dof_chain( &
      frequencies(1), torque, 0.2_dp, 0.3_dp, &
      c_wlf*g_storage_wlf, c_wlf*g_loss_wlf, 0.4_dp, &
      c_arrhenius*g_storage_arrhenius, &
      c_arrhenius*g_loss_arrhenius, 0.6_dp)
    call assert_complex_close(physical(2, 1), expected(1), tolerance, &
      "Multiple-provider node 2 response hatalı.")
    call assert_complex_close(physical(3, 1), expected(2), tolerance, &
      "Multiple-provider node 3 response hatalı.")
    if (get_material_binding_count(response) /= 2) then
      error stop "Multiple thermorheological provider trace sayısı hatalı."
    end if

    trace_wlf = get_material_state_trace(response, 1, 1)
    trace_arrhenius = get_material_state_trace(response, 2, 1)
    if (trace_wlf%shift_model_kind /= WLF_TEMPERATURE_SHIFT .or. &
        trace_arrhenius%shift_model_kind /= &
          ARRHENIUS_TEMPERATURE_SHIFT) then
      error stop "Multiple-provider shift model kimlikleri karıştı."
    end if
    call assert_close(trace_wlf%reference_temperature_k, &
      293.15_dp, tolerance, "WLF trace kendi T_ref değerini korumadı.")
    call assert_close(trace_arrhenius%reference_temperature_k, &
      298.15_dp, tolerance, "Arrhenius trace kendi T_ref değerini korumadı.")
    call assert_close(trace_wlf%lookup_frequency_hz, &
      10.0_dp**(1.0_dp+s_wlf), tolerance, "WLF reduced frequency hatalı.")
    call assert_close(trace_arrhenius%lookup_frequency_hz, &
      10.0_dp**(1.0_dp+s_arrhenius), tolerance, &
      "Arrhenius reduced frequency hatalı.")
  end subroutine test_multiple_shift_providers

  !> Material evaluation tamamlandıktan sonra Z exact singular olsa bile
  !! physical f/T, shift, reduced f_r, G* ve K* trace'inin korunduğunu sınar.
  subroutine test_singular_trace_is_preserved()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    type(tabulated_log10_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: bindings(1)
    type(dynamic_torsional_property_state_t) :: state
    type(material_aware_harmonic_response_t) :: response
    type(harmonic_response_t) :: harmonic_response
    type(material_state_trace_t) :: trace
    type(harmonic_excitation_t) :: excitation(1)
    type(rubber_geometry_t) :: rubber
    integer, allocatable :: statuses(:)
    logical, allocatable :: availability(:)
    real(dp) :: inertia
    real(dp) :: frequencies(1)

    master_curve = make_master_curve( &
      "SINGULAR-TTS", 293.15_dp, 1.0e6_dp, 1.2e6_dp, 0.0_dp, 0.0_dp)
    shift_provider = make_tabulated_shift()
    provider = create_thermorheological_dynamic_modulus_provider( &
      master_curve, shift_provider)
    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    bindings(1) = create_dynamic_torsional_property_binding( &
      10, provider, rubber)
    frequencies = [10.0_dp]
    state = evaluate_dynamic_torsional_property( &
      bindings(1), frequencies(1), 293.15_dp)
    inertia = find_exact_resonant_inertia( &
      state%storage_stiffness_nm_per_rad, &
      (2.0_dp*pi*frequencies(1))**2)
    call build_fixed_system( &
      system, manager, 10, 1234.0_dp, 0.0_dp, 0.0_dp, inertia)
    excitation(1) = make_excitation(2, cmplx(1.0_dp, 0.0_dp, kind=dp))

    response = analyze_material_aware_harmonic_response( &
      system, manager, frequencies, 293.15_dp, excitation, bindings)
    harmonic_response = get_base_harmonic_response(response)
    statuses = get_harmonic_solution_statuses(harmonic_response)
    availability = get_harmonic_response_availability(harmonic_response)
    if (statuses(1) /= COMPLEX_SOLVE_SINGULAR .or. availability(1)) then
      error stop "Thermorheological exact singular status hatalı."
    end if

    trace = get_material_state_trace(response, 1, 1)
    call assert_close(trace%physical_frequency_hz, &
      10.0_dp, 0.0_dp, "Singular trace physical f kayboldu.")
    call assert_close(trace%lookup_frequency_hz, &
      10.0_dp, 0.0_dp, "Singular trace reduced f_r kayboldu.")
    call assert_close(trace%a_t, &
      1.0_dp, 0.0_dp, "Singular trace a_T kayboldu.")
    call assert_close(trace%operating_temperature_k, &
      293.15_dp, 0.0_dp, "Singular trace operating T kayboldu.")
    call assert_close(trace%reference_temperature_k, &
      293.15_dp, 0.0_dp, "Singular trace reference T kayboldu.")
    call assert_close(trace%log10_a_t, &
      0.0_dp, 0.0_dp, "Singular trace log10(a_T) kayboldu.")
    call assert_close(trace%storage_modulus_pa, &
      1.1e6_dp, tolerance, "Singular trace G' kayboldu.")
    call assert_close(trace%loss_modulus_pa, &
      0.0_dp, 0.0_dp, "Singular trace G'' kayboldu.")
    call assert_close(trace%storage_stiffness_nm_per_rad, &
      state%storage_stiffness_nm_per_rad, 0.0_dp, &
      "Singular trace K' kayboldu.")
    call assert_close(trace%loss_stiffness_nm_per_rad, &
      0.0_dp, 0.0_dp, "Singular trace K'' kayboldu.")
    if (.not. trace%temperature_shift_applied .or. &
        trace%shift_model_kind /= TABULATED_LOG10_SHIFT .or. &
        .not. trace%shift_exact_temperature_point) then
      error stop "Singular trace shift model context'i kayboldu."
    end if
  end subroutine test_singular_trace_is_preserved

  !> Harmonic sweep prevalidation'ın shift-temperature ve reduced-frequency
  !! domain ihlallerini solver çağrısından önce reddettiğini tetikler.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(thermorheological_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: bindings(1)
    type(material_aware_harmonic_response_t) :: response
    type(harmonic_response_t) :: harmonic_response
    type(material_binding_trace_t), allocatable :: binding_traces(:)
    type(material_state_trace_t), allocatable :: state_traces(:, :)
    type(harmonic_excitation_t) :: excitation(1)
    type(rubber_geometry_t) :: rubber
    real(dp) :: frequencies(1)

    provider = make_thermorheological_provider()
    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    bindings(1) = create_dynamic_torsional_property_binding( &
      10, provider, rubber)
    call build_fixed_system( &
      system, manager, 10, 1000.0_dp, 100.0_dp, 0.0_dp, 0.25_dp)
    excitation(1) = make_excitation(2, cmplx(1.0_dp, 0.0_dp, kind=dp))

    select case (case_name)
    case ("uncovered_reduced_sweep")
      frequencies = [0.5_dp]
      response = analyze_material_aware_harmonic_response( &
        system, manager, frequencies, 313.15_dp, excitation, bindings)
    case ("temperature_outside_domain")
      frequencies = [10.0_dp]
      response = analyze_material_aware_harmonic_response( &
        system, manager, frequencies, 333.15_dp, excitation, bindings)
    case ("reference_trace_metadata_mismatch")
      frequencies = [10.0_dp]
      response = analyze_material_aware_harmonic_response( &
        system, manager, frequencies, 293.15_dp, excitation, bindings)
      harmonic_response = get_base_harmonic_response(response)
      binding_traces = get_material_binding_traces(response)
      state_traces = get_material_state_traces(response)
      state_traces(1, 1)%reference_temperature_k = 298.15_dp
      response = create_material_aware_harmonic_response( &
        harmonic_response, binding_traces, state_traces)
    case default
      error stop "Bilinmeyen thermorheological harmonic selector."
    end select
  end subroutine exercise_invalid_case

  subroutine build_fixed_system( &
      system, manager, element_id, stiffness, loss_stiffness, damping, inertia)
    type(torsional_system_t), intent(out) :: system
    type(constraint_manager_t), intent(out) :: manager
    integer, intent(in) :: element_id
    real(dp), intent(in) :: stiffness
    real(dp), intent(in) :: loss_stiffness
    real(dp), intent(in) :: damping
    real(dp), intent(in) :: inertia

    type(constraint_t) :: fixed_constraint_record

    call add_torsional_node( &
      system, torsional_node_t(1, 0.5_dp, 0.0_dp, .false.))
    call add_torsional_node( &
      system, torsional_node_t(2, inertia, 0.0_dp, .false.))
    call add_torsional_element(system, torsional_element_t( &
      element_id, 1, 2, stiffness, damping, loss_stiffness))
    call initialize_constraint_manager(manager)
    fixed_constraint_record = constraint_t( &
      1, 1, TORSIONAL_ROTATION, 0.0_dp, FIXED_CONSTRAINT)
    call add_constraint(manager, fixed_constraint_record, system)
  end subroutine build_fixed_system

  subroutine build_three_node_system(system, manager)
    type(torsional_system_t), intent(out) :: system
    type(constraint_manager_t), intent(out) :: manager

    type(constraint_t) :: fixed_constraint_record

    call add_torsional_node( &
      system, torsional_node_t(1, 0.5_dp, 0.0_dp, .false.))
    call add_torsional_node( &
      system, torsional_node_t(2, 0.2_dp, 0.0_dp, .false.))
    call add_torsional_node( &
      system, torsional_node_t(3, 0.3_dp, 0.0_dp, .false.))
    call add_torsional_element(system, torsional_element_t( &
      10, 1, 2, 7000.0_dp, 0.4_dp, 700.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      20, 2, 3, 300.0_dp, 0.6_dp, 30.0_dp))
    call initialize_constraint_manager(manager)
    fixed_constraint_record = constraint_t( &
      1, 1, TORSIONAL_ROTATION, 0.0_dp, FIXED_CONSTRAINT)
    call add_constraint(manager, fixed_constraint_record, system)
  end subroutine build_three_node_system

  function make_thermorheological_provider() result(provider)
    type(thermorheological_dynamic_modulus_provider_t) :: provider

    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    type(tabulated_log10_shift_provider_t) :: shift_provider

    master_curve = make_master_curve( &
      "INVALID-HARMONIC", 293.15_dp, 1.0e6_dp, 2.0e6_dp, &
      0.1e6_dp, 0.2e6_dp)
    shift_provider = make_tabulated_shift()
    provider = create_thermorheological_dynamic_modulus_provider( &
      master_curve, shift_provider)
  end function make_thermorheological_provider

  pure function make_tabulated_shift() result(provider)
    type(tabulated_log10_shift_provider_t) :: provider
    real(dp) :: temperature_k(3)
    real(dp) :: log10_a_t(3)

    temperature_k = [273.15_dp, 293.15_dp, 313.15_dp]
    log10_a_t = [1.0_dp, 0.0_dp, -1.0_dp]
    provider = create_tabulated_temperature_shift_provider( &
      temperature_k, log10_a_t, 293.15_dp)
  end function make_tabulated_shift

  pure function make_master_curve( &
      dataset_id, temperature_k, storage_1, storage_2, &
      loss_1, loss_2) result(provider)
    character(len=*), intent(in) :: dataset_id
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: storage_1
    real(dp), intent(in) :: storage_2
    real(dp), intent(in) :: loss_1
    real(dp), intent(in) :: loss_2
    type(tabulated_dynamic_modulus_provider_t) :: provider

    type(material_frequency_point) :: points(2)

    points(1) = make_point( &
      1.0_dp, storage_1, loss_1, temperature_k)
    points(2) = make_point( &
      100.0_dp, storage_2, loss_2, temperature_k)
    provider = create_tabulated_dynamic_modulus_provider( &
      points, make_metadata(dataset_id, temperature_k), &
      LINEAR_LOG_FREQUENCY)
  end function make_master_curve

  pure function make_point( &
      frequency_hz, storage_pa, loss_pa, temperature_k) result(point)
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: storage_pa
    real(dp), intent(in) :: loss_pa
    real(dp), intent(in) :: temperature_k
    type(material_frequency_point) :: point

    point%frequency = frequency_hz
    point%temperature = temperature_k
    point%storage_modulus = storage_pa
    point%loss_modulus = loss_pa
  end function make_point

  pure function make_metadata(dataset_id, temperature_k) result(metadata)
    character(len=*), intent(in) :: dataset_id
    real(dp), intent(in) :: temperature_k
    type(dynamic_material_metadata_t) :: metadata

    metadata%dataset_identifier = dataset_id
    metadata%material_identifier = "EPDM-TTS"
    metadata%dataset_temperature_k = temperature_k
    metadata%has_dynamic_shear_strain_amplitude = .true.
    metadata%dynamic_shear_strain_amplitude = 0.005_dp
    metadata%has_static_shear_prestrain = .true.
    metadata%static_shear_prestrain = 0.02_dp
    metadata%deformation_mode = DYNAMIC_DEFORMATION_MODE_SHEAR
  end function make_metadata

  pure function make_excitation(node_id, torque) result(excitation)
    integer, intent(in) :: node_id
    complex(dp), intent(in) :: torque
    type(harmonic_excitation_t) :: excitation

    excitation = harmonic_excitation_t( &
      node_id, TORSIONAL_ROTATION, torque)
  end function make_excitation

  pure function annular_factor(rubber) result(c_theta)
    type(rubber_geometry_t), intent(in) :: rubber
    real(dp) :: c_theta

    c_theta = 4.0_dp*pi*rubber%axial_length_m * &
      rubber%inner_radius_m**2*rubber%outer_radius_m**2 / &
      (rubber%outer_radius_m**2-rubber%inner_radius_m**2)
  end function annular_factor

  !> 1 ve 100 Hz endpoints arasında LINEAR_LOG_FREQUENCY interpolation
  !! sonucu verir. lookup log10 coordinate'i 1+s, alpha=(1+s)/2'dir.
  pure function log_interpolate(lower_value, upper_value, log_frequency) &
      result(value)
    real(dp), intent(in) :: lower_value
    real(dp), intent(in) :: upper_value
    real(dp), intent(in) :: log_frequency
    real(dp) :: value

    value = (1.0_dp-log_frequency/2.0_dp)*lower_value + &
      (log_frequency/2.0_dp)*upper_value
  end function log_interpolate

  !> İki active torsional DOF için complex symmetric dynamic stiffness
  !! matrisinin bağımsız kapalı-form inverse çözümüdür.
  pure function solve_two_dof_chain( &
      frequency_hz, torque, inertia_2, inertia_3, &
      storage_1, loss_1, damping_1, &
      storage_2, loss_2, damping_2) result(theta)
    real(dp), intent(in) :: frequency_hz
    complex(dp), intent(in) :: torque
    real(dp), intent(in) :: inertia_2
    real(dp), intent(in) :: inertia_3
    real(dp), intent(in) :: storage_1
    real(dp), intent(in) :: loss_1
    real(dp), intent(in) :: damping_1
    real(dp), intent(in) :: storage_2
    real(dp), intent(in) :: loss_2
    real(dp), intent(in) :: damping_2
    complex(dp) :: theta(2)

    complex(dp) :: z11
    complex(dp) :: z12
    complex(dp) :: z22
    complex(dp) :: determinant
    real(dp) :: omega

    omega = 2.0_dp*pi*frequency_hz
    z11 = cmplx(storage_1+storage_2-omega**2*inertia_2, &
      loss_1+loss_2+omega*(damping_1+damping_2), kind=dp)
    z12 = -cmplx(storage_2, loss_2+omega*damping_2, kind=dp)
    z22 = cmplx(storage_2-omega**2*inertia_3, &
      loss_2+omega*damping_2, kind=dp)
    determinant = z11*z22-z12*z12
    theta(1) = -z12*torque/determinant
    theta(2) = z11*torque/determinant
  end function solve_two_dof_chain

  pure function find_exact_resonant_inertia( &
      stiffness, omega_squared) result(inertia)
    real(dp), intent(in) :: stiffness
    real(dp), intent(in) :: omega_squared
    real(dp) :: inertia

    integer :: attempt
    real(dp) :: difference

    inertia = stiffness/omega_squared
    do attempt = 1, 256
      difference = stiffness-omega_squared*inertia
      if (abs(difference) <= 0.0_dp) return
      if (difference > 0.0_dp) then
        inertia = nearest(inertia, 1.0_dp)
      else
        inertia = nearest(inertia, -1.0_dp)
      end if
    end do
    error stop "Exact thermorheological resonance inertia bulunamadı."
  end function find_exact_resonant_inertia

  subroutine assert_complex_close(actual, expected, relative_tolerance, message)
    complex(dp), intent(in) :: actual
    complex(dp), intent(in) :: expected
    real(dp), intent(in) :: relative_tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(real(actual, kind=dp)) .or. &
        .not. ieee_is_finite(aimag(actual)) .or. &
        .not. ieee_is_finite(real(expected, kind=dp)) .or. &
        .not. ieee_is_finite(aimag(expected)) .or. &
        .not. ieee_is_finite(relative_tolerance) .or. &
        relative_tolerance < 0.0_dp .or. abs(actual-expected) > &
        relative_tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_complex_close

  subroutine assert_close(actual, expected, relative_tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: relative_tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(relative_tolerance) .or. &
        relative_tolerance < 0.0_dp .or. abs(actual-expected) > &
        relative_tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_close

end program test_thermorheological_harmonic_analysis
