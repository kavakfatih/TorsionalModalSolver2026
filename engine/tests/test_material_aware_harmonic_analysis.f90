program test_material_aware_harmonic_analysis
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_dof_types, only : TORSIONAL_ROTATION
  use tms_geometry, only : rubber_geometry_t
  use tms_material_frequency, only : material_frequency_point
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    DYNAMIC_DEFORMATION_MODE_SHEAR
  use tms_dynamic_modulus_provider, only : LINEAR_FREQUENCY
  use tms_temperature_shift_types, only : TEMPERATURE_SHIFT_NONE
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t, &
    create_tabulated_dynamic_modulus_provider
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element, get_torsional_element
  use tms_constraint_types, only : constraint_t, FIXED_CONSTRAINT
  use tms_constraint_manager, only : constraint_manager_t, &
    initialize_constraint_manager, add_constraint
  use tms_harmonic_excitation, only : harmonic_excitation_t
  use tms_harmonic_analysis, only : analyze_harmonic_response
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
    material_aware_harmonic_response_t, get_base_harmonic_response, &
    get_material_binding_count, get_material_binding_traces, &
    get_material_state_traces, get_material_state_trace
  use tms_material_aware_harmonic_analysis, only : &
    analyze_material_aware_harmonic_response
  implicit none

  real(dp), parameter :: tight_tolerance = 2.0e-11_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_frequency_dependent_one_dof_chain()
  call test_mixed_and_multiple_dynamic_bindings()
  call test_material_trace_at_singular_point()

  print *, "V0.7 material-aware harmonic response zinciri doğrulandı."

contains

  !> Interpolate G'/G'' -> map K'/K'' -> build Z -> solve theta tam zincirini
  !! fixed-hub 1-DOF TVD için bağımsız analitik ifadeyle doğrular.
  !! Aynı test dynamic override/no-double-counting, frequency-independent c,
  !! frozen V0.6 geriye uyumluluğu, trace ve input immutability'yi de kilitler.
  subroutine test_frequency_dependent_one_dof_chain()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: bindings(1)
    type(material_aware_harmonic_response_t) :: material_response
    type(harmonic_response_t) :: harmonic_response
    type(harmonic_response_t) :: frozen_response
    type(material_binding_trace_t), allocatable :: binding_traces(:)
    type(material_state_trace_t), allocatable :: traces(:, :)
    type(material_state_trace_t) :: trace
    type(harmonic_excitation_t) :: excitation(1)
    type(torsional_element_t) :: element_before
    type(torsional_element_t) :: element_after
    type(rubber_geometry_t) :: rubber
    complex(dp), allocatable :: physical(:, :)
    complex(dp) :: expected
    complex(dp) :: frozen_expected
    complex(dp) :: torque
    real(dp) :: angular_frequency
    real(dp) :: c_theta
    real(dp) :: expected_loss_modulus
    real(dp) :: expected_storage_modulus
    real(dp) :: frequencies(1)
    real(dp) :: dissipated_energy

    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    provider = make_provider( &
      "ONE-DOF", 1.0e6_dp, 2.0e6_dp, 0.1e6_dp, 0.2e6_dp)
    bindings(1) = create_dynamic_torsional_property_binding( &
      10, provider, rubber)
    call build_fixed_system( &
      system, manager, 10, 9000.0_dp, 800.0_dp, 2.0_dp, 0.25_dp)
    excitation(1) = make_excitation(2, cmplx(3.0_dp, -1.0_dp, kind=dp))
    frequencies = [15.0_dp]
    torque = excitation(1)%torque_amplitude_nm
    element_before = get_torsional_element(system, 1)

    material_response = analyze_material_aware_harmonic_response( &
      system, manager, frequencies, 293.15_dp, excitation, bindings)
    harmonic_response = get_base_harmonic_response(material_response)
    physical = get_physical_complex_response(harmonic_response)

    c_theta = annular_factor(rubber)
    expected_storage_modulus = 1.5e6_dp
    expected_loss_modulus = 0.15e6_dp
    angular_frequency = 2.0_dp*pi*frequencies(1)
    expected = torque / cmplx( &
      c_theta*expected_storage_modulus - &
        angular_frequency**2*0.25_dp, &
      c_theta*expected_loss_modulus + angular_frequency*2.0_dp, kind=dp)
    call assert_complex_close(physical(2, 1), expected, tight_tolerance, &
      "Material-aware 1-DOF tam zincir sonucu analitik değerle uyuşmuyor.")

    ! Nominal 9000/800 değerleri material-aware Z'ye eklenmiş olsaydı bu
    ! karşılaştırma açıkça başarısız olur; stored K' yine pozitif tutulur.
    if (element_before%stiffness_nm_per_rad <= 0.0_dp) then
      error stop "Dynamic-bound elemanın nominal K' değeri pozitif kalmadı."
    end if
    frozen_response = analyze_harmonic_response( &
      system, manager, frequencies, excitation)
    physical = get_physical_complex_response(frozen_response)
    frozen_expected = torque / cmplx( &
      9000.0_dp-angular_frequency**2*0.25_dp, &
      800.0_dp+angular_frequency*2.0_dp, kind=dp)
    call assert_complex_close(physical(2, 1), frozen_expected, &
      tight_tolerance, "V0.6 frozen nominal K'/K'' davranışı değişti.")

    if (get_material_binding_count(material_response) /= 1) then
      error stop "Material-aware result binding sayısını korumadı."
    end if
    binding_traces = get_material_binding_traces(material_response)
    traces = get_material_state_traces(material_response)
    trace = get_material_state_trace(material_response, 1, 1)
    if (size(traces, 1) /= 1 .or. size(traces, 2) /= 1 .or. &
        trim(binding_traces(1)%metadata%dataset_identifier) /= "ONE-DOF") then
      error stop "Material/dataset trace boyutu veya kimliği hatalı."
    end if
    call assert_close(trace%storage_modulus_pa, &
      expected_storage_modulus, tight_tolerance, "Trace G' hatalı.")
    call assert_close(trace%loss_modulus_pa, &
      expected_loss_modulus, tight_tolerance, "Trace G'' hatalı.")
    call assert_close(trace%storage_stiffness_nm_per_rad, &
      c_theta*expected_storage_modulus, tight_tolerance, "Trace K' hatalı.")
    call assert_close(trace%loss_stiffness_nm_per_rad, &
      c_theta*expected_loss_modulus, tight_tolerance, "Trace K'' hatalı.")
    if (trace%interpolation_policy /= LINEAR_FREQUENCY .or. &
        trace%exact_table_point .or. &
        abs(trace%interpolation_alpha-0.5_dp) > tight_tolerance) then
      error stop "Trace interpolation policy/bracket bilgisi hatalı."
    end if
    if (trace%temperature_shift_applied .or. &
        trace%shift_model_kind /= TEMPERATURE_SHIFT_NONE .or. &
        abs(trace%physical_frequency_hz-trace%frequency_hz) > 0.0_dp .or. &
        abs(trace%lookup_frequency_hz-trace%frequency_hz) > 0.0_dp .or. &
        abs(trace%operating_temperature_k-trace%temperature_k) > 0.0_dp .or. &
        abs(trace%log10_a_t) > 0.0_dp .or. &
        abs(trace%a_t-1.0_dp) > 0.0_dp .or. &
        trace%has_temperature_bracket) then
      error stop "V0.7 material trace unshifted semantiği değişti."
    end if
    dissipated_energy = pi*trace%loss_stiffness_nm_per_rad*abs(expected)**2
    if (dissipated_energy < 0.0_dp) then
      error stop "Passive dynamic material negatif dissipated energy üretti."
    end if

    element_after = get_torsional_element(system, 1)
    if (abs(element_after%stiffness_nm_per_rad- &
        element_before%stiffness_nm_per_rad) > 0.0_dp .or. &
        abs(element_after%loss_stiffness_nm_per_rad- &
        element_before%loss_stiffness_nm_per_rad) > 0.0_dp .or. &
        abs(element_after%damping_nms_per_rad- &
        element_before%damping_nms_per_rad) > 0.0_dp .or. &
        abs(frequencies(1)-15.0_dp) > 0.0_dp) then
      error stop "Material-aware analysis authoritative input'u değiştirdi."
    end if
  end subroutine test_frequency_dependent_one_dof_chain

  !> Constant ve dynamic elemanların aynı sistemde doğru ayrıldığını, ayrıca
  !! iki farklı provider'ın yalnız kendi element ID'sine uygulandığını 2-DOF
  !! complex matrix inverse ile doğrular. Stored c iki durumda da korunur.
  subroutine test_mixed_and_multiple_dynamic_bindings()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(tabulated_dynamic_modulus_provider_t) :: provider_a
    type(tabulated_dynamic_modulus_provider_t) :: provider_b
    type(dynamic_torsional_property_binding_t) :: mixed_bindings(1)
    type(dynamic_torsional_property_binding_t) :: all_bindings(2)
    type(material_aware_harmonic_response_t) :: response
    type(harmonic_response_t) :: harmonic_response
    type(harmonic_excitation_t) :: excitation(1)
    type(rubber_geometry_t) :: rubber_a
    type(rubber_geometry_t) :: rubber_b
    complex(dp), allocatable :: physical(:, :)
    complex(dp) :: expected(2)
    complex(dp) :: torque
    real(dp) :: c_theta_a
    real(dp) :: c_theta_b
    real(dp) :: frequencies(1)

    rubber_a = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    rubber_b = rubber_geometry_t(0.015_dp, 0.04_dp, 0.012_dp)
    provider_a = make_provider( &
      "MATERIAL-A", 1.0e6_dp, 2.0e6_dp, 0.1e6_dp, 0.2e6_dp)
    provider_b = make_provider( &
      "MATERIAL-B", 0.5e6_dp, 0.9e6_dp, 0.05e6_dp, 0.09e6_dp)
    call build_three_node_system(system, manager)
    mixed_bindings(1) = create_dynamic_torsional_property_binding( &
      10, provider_a, rubber_a)
    all_bindings(1) = mixed_bindings(1)
    all_bindings(2) = create_dynamic_torsional_property_binding( &
      20, provider_b, rubber_b)
    frequencies = [15.0_dp]
    torque = cmplx(2.0_dp, 0.4_dp, kind=dp)
    excitation(1) = make_excitation(3, torque)
    c_theta_a = annular_factor(rubber_a)
    c_theta_b = annular_factor(rubber_b)

    response = analyze_material_aware_harmonic_response( &
      system, manager, frequencies, 293.15_dp, excitation, mixed_bindings)
    harmonic_response = get_base_harmonic_response(response)
    physical = get_physical_complex_response(harmonic_response)
    expected = solve_two_dof_chain( &
      frequencies(1), torque, 0.2_dp, 0.3_dp, &
      c_theta_a*1.5e6_dp, c_theta_a*0.15e6_dp, 1.0_dp, &
      300.0_dp, 30.0_dp, 0.5_dp)
    call assert_complex_close(physical(2, 1), expected(1), &
      tight_tolerance, "Mixed system node 2 response hatalı.")
    call assert_complex_close(physical(3, 1), expected(2), &
      tight_tolerance, "Mixed system node 3 response hatalı.")

    response = analyze_material_aware_harmonic_response( &
      system, manager, frequencies, 293.15_dp, excitation, all_bindings)
    harmonic_response = get_base_harmonic_response(response)
    physical = get_physical_complex_response(harmonic_response)
    expected = solve_two_dof_chain( &
      frequencies(1), torque, 0.2_dp, 0.3_dp, &
      c_theta_a*1.5e6_dp, c_theta_a*0.15e6_dp, 1.0_dp, &
      c_theta_b*0.7e6_dp, c_theta_b*0.07e6_dp, 0.5_dp)
    call assert_complex_close(physical(2, 1), expected(1), &
      tight_tolerance, "Multiple-provider node 2 response hatalı.")
    call assert_complex_close(physical(3, 1), expected(2), &
      tight_tolerance, "Multiple-provider node 3 response hatalı.")
    if (get_material_binding_count(response) /= 2) then
      error stop "Multiple dynamic provider trace sayısı hatalı."
    end if
  end subroutine test_mixed_and_multiple_dynamic_bindings

  !> Material evaluation geçerli olduktan sonra Z exact singular olsa bile
  !! G'/G'' ve K'/K'' trace'inin her requested frequency için korunduğunu
  !! doğrular. Harmonic response unavailable, material trace available'dır.
  subroutine test_material_trace_at_singular_point()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: bindings(1)
    type(dynamic_torsional_property_state_t) :: state
    type(material_aware_harmonic_response_t) :: response
    type(harmonic_response_t) :: harmonic_response
    type(material_state_trace_t) :: trace
    type(harmonic_excitation_t) :: excitation(1)
    type(material_frequency_point) :: points(2)
    type(dynamic_material_metadata_t) :: metadata
    type(rubber_geometry_t) :: rubber
    integer, allocatable :: statuses(:)
    logical, allocatable :: availability(:)
    real(dp) :: frequencies(1)
    real(dp) :: inertia
    real(dp) :: omega_squared

    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    frequencies = [10.0_dp]
    omega_squared = (2.0_dp*pi*frequencies(1))**2
    points(1) = make_point(10.0_dp, 1.0e6_dp, 0.0_dp)
    points(2) = make_point(20.0_dp, 1.2e6_dp, 0.0_dp)
    metadata = make_metadata("SINGULAR-TRACE")
    provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    bindings(1) = create_dynamic_torsional_property_binding( &
      10, provider, rubber)
    state = evaluate_dynamic_torsional_property( &
      bindings(1), frequencies(1), 293.15_dp)
    inertia = find_exact_resonant_inertia( &
      state%storage_stiffness_nm_per_rad, omega_squared)
    call build_fixed_system( &
      system, manager, 10, 1234.0_dp, 0.0_dp, 0.0_dp, inertia)
    excitation(1) = make_excitation(2, cmplx(1.0_dp, 0.0_dp, kind=dp))

    response = analyze_material_aware_harmonic_response( &
      system, manager, frequencies, 293.15_dp, excitation, bindings)
    harmonic_response = get_base_harmonic_response(response)
    statuses = get_harmonic_solution_statuses(harmonic_response)
    availability = get_harmonic_response_availability(harmonic_response)
    if (statuses(1) /= COMPLEX_SOLVE_SINGULAR .or. availability(1)) then
      error stop "Exact singular material-aware nokta doğru status üretmedi."
    end if
    trace = get_material_state_trace(response, 1, 1)
    call assert_close(trace%storage_modulus_pa, &
      1.0e6_dp, 0.0_dp, "Singular noktada material G' trace kayboldu.")
    call assert_close(trace%storage_stiffness_nm_per_rad, &
      state%storage_stiffness_nm_per_rad, 0.0_dp, &
      "Singular noktada mapped K' trace kayboldu.")
  end subroutine test_material_trace_at_singular_point

  !> Duplicate/unknown binding, full-sweep domain, isotherm ve boş binding
  !! prevalidation hata yollarını ayrı CTest süreçlerinde tetikler.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t), allocatable :: bindings(:)
    type(material_aware_harmonic_response_t) :: response
    type(harmonic_excitation_t) :: excitation(1)
    type(rubber_geometry_t) :: rubber
    real(dp), allocatable :: frequencies(:)

    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    provider = make_provider( &
      "INVALID-ANALYSIS", 1.0e6_dp, 2.0e6_dp, 0.1e6_dp, 0.2e6_dp)
    call build_fixed_system( &
      system, manager, 10, 1000.0_dp, 100.0_dp, 0.0_dp, 0.25_dp)
    excitation(1) = make_excitation(2, cmplx(1.0_dp, 0.0_dp, kind=dp))
    frequencies = [10.0_dp]

    select case (case_name)
    case ("duplicate_binding")
      allocate(bindings(2))
      bindings(1) = create_dynamic_torsional_property_binding( &
        10, provider, rubber)
      bindings(2) = bindings(1)
    case ("unknown_element")
      allocate(bindings(1))
      bindings(1) = create_dynamic_torsional_property_binding( &
        99, provider, rubber)
    case ("uncovered_sweep")
      allocate(bindings(1))
      bindings(1) = create_dynamic_torsional_property_binding( &
        10, provider, rubber)
      frequencies = [10.0_dp, 30.0_dp]
    case ("temperature_mismatch")
      allocate(bindings(1))
      bindings(1) = create_dynamic_torsional_property_binding( &
        10, provider, rubber)
    case ("empty_bindings")
      allocate(bindings(0))
    case default
      error stop "Bilinmeyen material-aware validation selector."
    end select

    if (case_name == "temperature_mismatch") then
      response = analyze_material_aware_harmonic_response( &
        system, manager, frequencies, 303.15_dp, excitation, bindings)
    else
      response = analyze_material_aware_harmonic_response( &
        system, manager, frequencies, 293.15_dp, excitation, bindings)
    end if
  end subroutine exercise_invalid_case

  !> Fixed hub ve tek ring ataletli TVD fixture'ı kurar. Stored K'/K''
  !! nominal/frozen değerler, c frequency-independent viscous katsayıdır.
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

    call add_torsional_node(system, torsional_node_t(1, 0.5_dp, 0.0_dp, .false.))
    call add_torsional_node(system, torsional_node_t(2, inertia, 0.0_dp, .false.))
    call add_torsional_element(system, torsional_element_t( &
      element_id, 1, 2, stiffness, damping, loss_stiffness))
    call initialize_constraint_manager(manager)
    fixed_constraint_record = constraint_t( &
      1, 1, TORSIONAL_ROTATION, 0.0_dp, FIXED_CONSTRAINT)
    call add_constraint(manager, fixed_constraint_record, system)
  end subroutine build_fixed_system

  !> Fixed node 1, free J2/J3 ve iki seri elemanlı mixed-system fixture'ı.
  subroutine build_three_node_system(system, manager)
    type(torsional_system_t), intent(out) :: system
    type(constraint_manager_t), intent(out) :: manager

    type(constraint_t) :: fixed_constraint_record

    call add_torsional_node(system, torsional_node_t(1, 0.5_dp, 0.0_dp, .false.))
    call add_torsional_node(system, torsional_node_t(2, 0.2_dp, 0.0_dp, .false.))
    call add_torsional_node(system, torsional_node_t(3, 0.3_dp, 0.0_dp, .false.))
    call add_torsional_element(system, torsional_element_t( &
      10, 1, 2, 7000.0_dp, 1.0_dp, 700.0_dp))
    call add_torsional_element(system, torsional_element_t( &
      20, 2, 3, 300.0_dp, 0.5_dp, 30.0_dp))
    call initialize_constraint_manager(manager)
    fixed_constraint_record = constraint_t( &
      1, 1, TORSIONAL_ROTATION, 0.0_dp, FIXED_CONSTRAINT)
    call add_constraint(manager, fixed_constraint_record, system)
  end subroutine build_three_node_system

  !> İki serbest koordinatlı seri torsional sistemin complex 2x2 inverse
  !! analitik çözümünü verir. K'/K'' [N*m/rad], c [N*m*s/rad], J [kg*m^2],
  !! f [Hz], torque [N*m], output theta [rad] birimindedir.
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

  !> K=omega^2*J işleminin aynı IEEE arithmetic sırasıyla exact olacağı
  !! pozitif J [kg*m^2] değerini nearest-neighbor aramayla bulur. Bu yalnız
  !! singular status regresyon fixture'ı üretir; production fizik modeli değildir.
  pure function find_exact_resonant_inertia(stiffness, omega_squared) &
      result(inertia)
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
    error stop "Exact resonance inertia fixture oluşturulamadı."
  end function find_exact_resonant_inertia

  pure function make_provider( &
      dataset_id, storage_1, storage_2, loss_1, loss_2) result(provider)
    character(len=*), intent(in) :: dataset_id
    real(dp), intent(in) :: storage_1
    real(dp), intent(in) :: storage_2
    real(dp), intent(in) :: loss_1
    real(dp), intent(in) :: loss_2
    type(tabulated_dynamic_modulus_provider_t) :: provider

    type(material_frequency_point) :: points(2)

    points(1) = make_point(10.0_dp, storage_1, loss_1)
    points(2) = make_point(20.0_dp, storage_2, loss_2)
    provider = create_tabulated_dynamic_modulus_provider( &
      points, make_metadata(dataset_id), LINEAR_FREQUENCY)
  end function make_provider

  pure function make_point(frequency, storage, loss) result(point)
    real(dp), intent(in) :: frequency
    real(dp), intent(in) :: storage
    real(dp), intent(in) :: loss
    type(material_frequency_point) :: point

    point%frequency = frequency
    point%temperature = 293.15_dp
    point%storage_modulus = storage
    point%loss_modulus = loss
  end function make_point

  pure function make_metadata(dataset_id) result(metadata)
    character(len=*), intent(in) :: dataset_id
    type(dynamic_material_metadata_t) :: metadata

    metadata%dataset_identifier = dataset_id
    metadata%material_identifier = "EPDM"
    metadata%dataset_temperature_k = 293.15_dp
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

  subroutine assert_complex_close(actual, expected, tolerance, message)
    complex(dp), intent(in) :: actual
    complex(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (abs(actual-expected) > tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_complex_close

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (abs(actual-expected) > tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_close

end program test_material_aware_harmonic_analysis
