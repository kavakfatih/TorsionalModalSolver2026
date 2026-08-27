program test_harmonic_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
    ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_dof_types, only : TORSIONAL_ROTATION
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element, get_torsional_element
  use tms_constraint_types, only : constraint_t, FIXED_CONSTRAINT
  use tms_constraint_manager, only : constraint_manager_t, &
    active_dof_map_t, initialize_constraint_manager, &
    initialize_constraint_manager_from_system, add_constraint
  use tms_torsional_system, only : two_inertia_tvd_system_t, &
    build_generalized_two_inertia_system
  use tms_reduced_system, only : reduced_torsional_system_t, &
    build_reduced_torsional_system
  use tms_modal_analysis, only : analyze_reduced_torsional_system
  use tms_modal_result, only : modal_result_t, get_modal_frequencies_hz
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    get_stiffness_matrix_values
  use tms_loss_stiffness_matrix, only : loss_stiffness_matrix_t, &
    get_loss_stiffness_matrix_values
  use tms_damping_matrix, only : damping_matrix_t, &
    get_damping_matrix_values
  use tms_mass_matrix, only : mass_matrix_t, get_mass_matrix_values
  use tms_reduced_dynamic_system, only : &
    reduced_dynamic_torsional_system_t, &
    build_reduced_dynamic_torsional_system, &
    get_reduced_dynamic_stiffness, &
    get_reduced_dynamic_loss_stiffness, &
    get_reduced_dynamic_damping, get_reduced_dynamic_mass, &
    get_reduced_dynamic_active_dof_map
  use tms_dynamic_stiffness, only : dynamic_stiffness_matrix_t, &
    build_dynamic_stiffness, get_dynamic_stiffness_values
  use tms_harmonic_excitation, only : harmonic_excitation_t, &
    assemble_harmonic_load_vector
  use tms_harmonic_analysis, only : analyze_harmonic_response
  use tms_harmonic_response, only : harmonic_response_t, &
    get_harmonic_solution_statuses, get_harmonic_response_availability, &
    get_harmonic_reciprocal_condition_numbers, &
    get_harmonic_relative_residuals, &
    get_harmonic_forward_error_bounds, get_harmonic_backward_errors, &
    get_reduced_complex_response, get_physical_complex_response, &
    get_physical_response_magnitudes, get_physical_response_phases_rad, &
    get_physical_angular_velocities, &
    get_physical_angular_accelerations, get_harmonic_backend_identity, &
    calculate_element_relative_angle, calculate_element_dynamic_torque, &
    calculate_element_average_dissipated_power, &
    calculate_element_dissipated_energy_per_cycle, &
    get_element_relative_angle_response, &
    get_element_dynamic_torque_response, &
    get_element_dynamic_torque_magnitudes, &
    get_element_average_dissipated_power_response, &
    get_element_dissipated_energy_response
  use tms_frf, only : calculate_rotational_receptance, &
    calculate_rotational_mobility, calculate_rotational_accelerance
  use tms_complex_linear_solution, only : COMPLEX_SOLVE_SOLVED, &
    COMPLEX_SOLVE_SINGULAR
  implicit none

  real(dp), parameter :: tight_tolerance = 5.0e-12_dp
  real(dp), parameter :: engineering_tolerance = 1.0e-9_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_fixed_one_dof_viscous_response()
  call test_loss_and_combined_damping_response()
  call test_phase_convention()
  call test_fixed_hub_two_inertia_bridge_response()
  call test_free_free_balanced_response()
  call test_low_frequency_limit()
  call test_singular_sweep_status()
  call test_dynamic_stiffness_contract()
  call test_excitation_scatter_add()
  call test_passivity_and_element_results()
  call test_element_sweep_results()
  call test_rotational_frf_helpers()
  call test_modal_harmonic_cross_validation()

  print *, "V0.6 direct frequency-domain harmonic response doğrulandı."

contains

  !> Fixed 1-DOF viskoz sistemin complex response'unu analitik çözümle sınar.
  !! Fiziksel model: m*theta''+c*theta'+k*theta=T ve K''=0.
  !! Matematiksel referans: theta_hat=T_hat/(k-m*omega^2+i*omega*c).
  !! Kütle [kg*m^2], k [N*m/rad], c [N*m*s/rad], f [Hz], response [rad].
  subroutine test_fixed_one_dof_viscous_response()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: response
    type(torsional_element_t) :: element_before
    type(torsional_element_t) :: element_after
    complex(dp), allocatable :: physical(:, :)
    complex(dp), allocatable :: reduced(:, :)
    complex(dp), allocatable :: velocities(:, :)
    complex(dp), allocatable :: accelerations(:, :)
    real(dp), allocatable :: magnitudes(:, :)
    real(dp), allocatable :: phases(:, :)
    real(dp), allocatable :: residuals(:)
    real(dp), allocatable :: rcond(:)
    real(dp), allocatable :: ferr(:)
    real(dp), allocatable :: berr(:)
    complex(dp) :: expected
    complex(dp) :: torque_before
    real(dp) :: frequency_before
    real(dp) :: frequency_hz(1)
    real(dp) :: angular_frequency
    character(len=:), allocatable :: backend

    call build_fixed_one_dof_system( &
      system, manager, 1000.0_dp, 0.0_dp, 2.0_dp, 0.25_dp)
    excitation(1) = make_excitation(2, cmplx(3.0_dp, 2.0_dp, kind=dp))
    frequency_hz = [5.0_dp]
    element_before = get_torsional_element(system, 1)
    torque_before = excitation(1)%torque_amplitude_nm
    frequency_before = frequency_hz(1)

    response = analyze_harmonic_response( &
      system, manager, frequency_hz, excitation)
    angular_frequency = 2.0_dp*pi*frequency_hz(1)
    expected = torque_before / cmplx( &
      1000.0_dp-0.25_dp*angular_frequency**2, &
      angular_frequency*2.0_dp, kind=dp)

    physical = get_physical_complex_response(response)
    reduced = get_reduced_complex_response(response)
    call assert_complex_close( &
      physical(1, 1), cmplx(0.0_dp, 0.0_dp, kind=dp), &
      tight_tolerance, "Fixed DOF harmonic phasor sıfır değil.")
    call assert_complex_close(physical(2, 1), expected, tight_tolerance, &
      "1-DOF viskoz physical response analitik değerle uyuşmuyor.")
    call assert_complex_close(reduced(1, 1), expected, tight_tolerance, &
      "1-DOF viskoz reduced response analitik değerle uyuşmuyor.")

    magnitudes = get_physical_response_magnitudes(response)
    phases = get_physical_response_phases_rad(response)
    velocities = get_physical_angular_velocities(response)
    accelerations = get_physical_angular_accelerations(response)
    call assert_real_close(magnitudes(2, 1), abs(expected), tight_tolerance, &
      "Response magnitude abs(theta) ile uyuşmuyor.")
    call assert_real_close(phases(2, 1), &
      atan2(aimag(expected), real(expected, dp)), tight_tolerance, &
      "Response phase atan2 convention ile uyuşmuyor.")
    call assert_complex_close(velocities(2, 1), &
      cmplx(0.0_dp, angular_frequency, kind=dp)*expected, &
      tight_tolerance, "Angular velocity i*omega*theta değil.")
    call assert_complex_close(accelerations(2, 1), &
      -angular_frequency**2*expected, tight_tolerance, &
      "Angular acceleration -omega^2*theta değil.")

    residuals = get_harmonic_relative_residuals(response)
    rcond = get_harmonic_reciprocal_condition_numbers(response)
    ferr = get_harmonic_forward_error_bounds(response)
    berr = get_harmonic_backward_errors(response)
    call assert_finite_nonnegative_diagnostics(residuals, "Residual")
    call assert_finite_nonnegative_diagnostics(rcond, "RCOND")
    call assert_finite_nonnegative_diagnostics(ferr, "FERR")
    call assert_finite_nonnegative_diagnostics(berr, "BERR")
    if (residuals(1) > 1.0e-13_dp) then
      error stop "1-DOF complex relative residual beklenenden büyük."
    end if
    backend = get_harmonic_backend_identity(response)
    if (index(backend, "ZSYSVX") == 0) then
      error stop "Harmonic result ZSYSVX backend kimliğini korumadı."
    end if

    ! Solver authoritative system/load/frequency girdilerini değiştirmemelidir.
    element_after = get_torsional_element(system, 1)
    if (abs(element_after%stiffness_nm_per_rad - &
        element_before%stiffness_nm_per_rad) > 0.0_dp .or. &
        abs(element_after%loss_stiffness_nm_per_rad - &
        element_before%loss_stiffness_nm_per_rad) > 0.0_dp .or. &
        abs(element_after%damping_nms_per_rad - &
        element_before%damping_nms_per_rad) > 0.0_dp .or. &
        abs(excitation(1)%torque_amplitude_nm-torque_before) > 0.0_dp .or. &
        abs(frequency_hz(1)-frequency_before) > 0.0_dp) then
      error stop "Harmonic analysis authoritative girdileri değiştirdi."
    end if
  end subroutine test_fixed_one_dof_viscous_response

  !> Kayıp rijitliği ile viskoz sönümün ayrı ve birlikte doğru etki ettiğini
  !! analitik 1-DOF response üzerinden doğrular.
  !! Model: theta=T/[K'-omega^2*M+i(K''+omega*C)]. K'' [N*m/rad] ile
  !! C [N*m*s/rad] ayrı kanallar olarak girer; aralarında dönüşüm yapılmaz.
  subroutine test_loss_and_combined_damping_response()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: response
    complex(dp), allocatable :: physical(:, :)
    complex(dp) :: expected
    complex(dp) :: torque
    real(dp) :: angular_frequency
    real(dp) :: frequencies(1)

    frequencies = [7.0_dp]
    angular_frequency = 2.0_dp*pi*frequencies(1)
    torque = cmplx(1.2_dp, -0.4_dp, kind=dp)
    excitation(1) = make_excitation(2, torque)

    ! Yalnız K'' kanalı.
    call build_fixed_one_dof_system( &
      system, manager, 500.0_dp, 50.0_dp, 0.0_dp, 0.1_dp)
    response = analyze_harmonic_response( &
      system, manager, frequencies, excitation)
    physical = get_physical_complex_response(response)
    expected = torque/cmplx( &
      500.0_dp-0.1_dp*angular_frequency**2, 50.0_dp, kind=dp)
    call assert_complex_close(physical(2, 1), expected, tight_tolerance, &
      "K''-only harmonic response analitik değerle uyuşmuyor.")

    ! K'' ve C aynı anda, fakat ayrı boyutsal terimler olarak etkindir.
    call build_fixed_one_dof_system( &
      system, manager, 500.0_dp, 50.0_dp, 3.0_dp, 0.1_dp)
    response = analyze_harmonic_response( &
      system, manager, frequencies, excitation)
    physical = get_physical_complex_response(response)
    expected = torque/cmplx( &
      500.0_dp-0.1_dp*angular_frequency**2, &
      50.0_dp+angular_frequency*3.0_dp, kind=dp)
    call assert_complex_close(physical(2, 1), expected, tight_tolerance, &
      "Combined K''+C harmonic response analitik değerle uyuşmuyor.")
  end subroutine test_loss_and_combined_damping_response

  !> exp(+i*omega*t) faz işaretini passive 1-DOF sistem üzerinde kilitler.
  !! Pozitif real torque için response fazı resonance altında (0,-pi/2),
  !! resonance'da -pi/2 ve üstünde (-pi,-pi/2) aralığında olmalıdır.
  subroutine test_phase_convention()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: response
    real(dp), allocatable :: phases(:, :)
    real(dp) :: natural_frequency
    real(dp) :: frequencies(3)

    natural_frequency = sqrt(1000.0_dp/0.25_dp)/(2.0_dp*pi)
    frequencies = [ &
      0.5_dp*natural_frequency, natural_frequency, &
      2.0_dp*natural_frequency]
    call build_fixed_one_dof_system( &
      system, manager, 1000.0_dp, 0.0_dp, 1.0_dp, 0.25_dp)
    excitation(1) = make_excitation( &
      2, cmplx(1.0_dp, 0.0_dp, kind=dp))
    response = analyze_harmonic_response( &
      system, manager, frequencies, excitation)
    phases = get_physical_response_phases_rad(response)

    if (.not. (phases(2, 1) < 0.0_dp .and. &
        phases(2, 1) > -0.5_dp*pi)) then
      error stop "Resonance altı faz exp(+iwt) convention ile uyumsuz."
    end if
    call assert_real_close( &
      phases(2, 2), -0.5_dp*pi, 1.0e-10_dp, &
      "Resonance fazı -pi/2 değil.")
    if (.not. (phases(2, 3) < -0.5_dp*pi .and. &
        phases(2, 3) > -pi)) then
      error stop "Resonance üstü faz exp(+iwt) convention ile uyumsuz."
    end if
  end subroutine test_phase_convention

  !> Mevcut two-inertia TVD köprüsünün K'' değerini viskoz c'ye çevirmeden
  !! fixed-hub harmonic response zincirine taşıdığını doğrular.
  !! Matematiksel referans: theta_ring=T/[K'-omega^2*J_r+i*K''].
  !! J_r [kg*m^2], K'/K'' [N*m/rad], torque [N*m], response [rad].
  subroutine test_fixed_hub_two_inertia_bridge_response()
    type(two_inertia_tvd_system_t) :: source_system
    type(torsional_system_t) :: generalized_system
    type(constraint_manager_t) :: manager
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: response
    type(torsional_element_t) :: element
    complex(dp), allocatable :: physical(:, :)
    complex(dp) :: expected
    complex(dp) :: torque
    real(dp) :: angular_frequency
    real(dp) :: frequencies(1)

    source_system = two_inertia_tvd_system_t( &
      hub_polar_inertia_kg_m2=0.12_dp, &
      ring_polar_inertia_kg_m2=0.08_dp, &
      storage_stiffness_nm_per_rad=750.0_dp, &
      loss_stiffness_nm_per_rad=75.0_dp, loss_factor=0.1_dp, &
      material_reference_frequency_hz=20.0_dp, &
      material_temperature_k=293.15_dp)
    generalized_system = &
      build_generalized_two_inertia_system(source_system, .true.)
    call initialize_constraint_manager_from_system( &
      manager, generalized_system)
    element = get_torsional_element(generalized_system, 1)
    if (abs(element%loss_stiffness_nm_per_rad - &
        source_system%loss_stiffness_nm_per_rad) > 0.0_dp .or. &
        abs(element%damping_nms_per_rad) > 0.0_dp) then
      error stop "Two-inertia bridge K'' değerini ayrı kanalda korumadı."
    end if

    frequencies = [9.0_dp]
    torque = cmplx(2.0_dp, 0.5_dp, kind=dp)
    excitation(1) = make_excitation(2, torque)
    angular_frequency = 2.0_dp*pi*frequencies(1)
    expected = torque/cmplx( &
      source_system%storage_stiffness_nm_per_rad - &
      angular_frequency**2*source_system%ring_polar_inertia_kg_m2, &
      source_system%loss_stiffness_nm_per_rad, kind=dp)
    response = analyze_harmonic_response( &
      generalized_system, manager, frequencies, excitation)
    physical = get_physical_complex_response(response)
    call assert_complex_close(physical(1, 1), &
      cmplx(0.0_dp, 0.0_dp, kind=dp), 0.0_dp, &
      "Fixed hub physical harmonic response sıfır değil.")
    call assert_complex_close(physical(2, 1), expected, tight_tolerance, &
      "Fixed-hub two-inertia harmonic response analitik değerle uyuşmuyor.")
  end subroutine test_fixed_hub_two_inertia_bridge_response

  !> Free-free iki ataletli sistemin finite-frequency balanced-torque
  !! relative response'unu analitik reduced-inertia bağıntısıyla doğrular.
  !! Matematiksel model: Delta=T/[K'+i(K''+omega*c)-omega^2*J_eq],
  !! J_eq=J1*J2/(J1+J2). Rigid mode finite f çözümünü engellememelidir.
  subroutine test_free_free_balanced_response()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(harmonic_excitation_t) :: excitations(2)
    type(harmonic_response_t) :: response
    complex(dp), allocatable :: physical(:, :)
    complex(dp) :: expected_relative_angle
    complex(dp) :: torque
    real(dp) :: angular_frequency
    real(dp) :: equivalent_inertia
    real(dp) :: frequencies(1)

    call build_free_two_inertia_system( &
      system, manager, 0.1_dp, 0.2_dp, 1000.0_dp, 50.0_dp, 1.5_dp)
    torque = cmplx(2.0_dp, -0.5_dp, kind=dp)
    excitations(1) = make_excitation(1, torque)
    excitations(2) = make_excitation(2, -torque)
    frequencies = [12.0_dp]
    angular_frequency = 2.0_dp*pi*frequencies(1)
    equivalent_inertia = 0.1_dp*0.2_dp/(0.1_dp+0.2_dp)
    expected_relative_angle = torque/cmplx( &
      1000.0_dp-angular_frequency**2*equivalent_inertia, &
      50.0_dp+angular_frequency*1.5_dp, kind=dp)

    response = analyze_harmonic_response( &
      system, manager, frequencies, excitations)
    physical = get_physical_complex_response(response)
    call assert_complex_close( &
      physical(1, 1)-physical(2, 1), expected_relative_angle, &
      engineering_tolerance, &
      "Free-free balanced relative response analitik değerle uyuşmuyor.")
  end subroutine test_free_free_balanced_response

  !> Pozitif çok küçük frekansta constrained undamped 1-DOF response'un
  !! quasi-static compliance T/K limitine yaklaştığını doğrular. 0 Hz
  !! kullanılmaz; static solver davranışı harmonic API'ye eklenmez.
  subroutine test_low_frequency_limit()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: response
    complex(dp), allocatable :: physical(:, :)
    real(dp) :: frequencies(1)

    call build_fixed_one_dof_system( &
      system, manager, 2000.0_dp, 0.0_dp, 0.0_dp, 0.5_dp)
    excitation(1) = make_excitation( &
      2, cmplx(4.0_dp, 0.0_dp, kind=dp))
    frequencies = [1.0e-5_dp]
    response = analyze_harmonic_response( &
      system, manager, frequencies, excitation)
    physical = get_physical_complex_response(response)
    call assert_complex_close(physical(2, 1), &
      cmplx(4.0_dp/2000.0_dp, 0.0_dp, kind=dp), 1.0e-10_dp, &
      "Low-frequency response T/K limitine yaklaşmadı.")
  end subroutine test_low_frequency_limit

  !> Undamped exact resonance içeren sweep'in solved-singular-solved status
  !! dizisini döndürdüğünü ve singular sütunda uydurma X üretmediğini doğrular.
  !! K=(2*pi)^2, M=1 ve f=1 Hz seçimi Z=0 değerini deterministic kurar.
  subroutine test_singular_sweep_status()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: response
    complex(dp), allocatable :: physical(:, :)
    integer, allocatable :: statuses(:)
    logical, allocatable :: availability(:)
    real(dp), allocatable :: residuals(:)
    real(dp) :: frequencies(3)
    real(dp) :: exact_stiffness

    exact_stiffness = (2.0_dp*pi)*(2.0_dp*pi)
    call build_fixed_one_dof_system( &
      system, manager, exact_stiffness, 0.0_dp, 0.0_dp, 1.0_dp)
    excitation(1) = make_excitation( &
      2, cmplx(1.0_dp, 0.0_dp, kind=dp))
    frequencies = [0.5_dp, 1.0_dp, 1.5_dp]
    response = analyze_harmonic_response( &
      system, manager, frequencies, excitation)
    statuses = get_harmonic_solution_statuses(response)
    availability = get_harmonic_response_availability(response)
    physical = get_physical_complex_response(response)
    residuals = get_harmonic_relative_residuals(response)

    if (any(statuses /= [ &
        COMPLEX_SOLVE_SOLVED, COMPLEX_SOLVE_SINGULAR, &
        COMPLEX_SOLVE_SOLVED])) then
      error stop "Exact resonance sweep status dizisi yanlış."
    end if
    if (any(availability .neqv. [.true., .false., .true.])) then
      error stop "Singular sweep response availability dizisi yanlış."
    end if
    if (.not. ieee_is_nan(real(physical(1, 2), dp)) .or. &
        .not. ieee_is_nan(aimag(physical(1, 2))) .or. &
        .not. ieee_is_nan(real(physical(2, 2), dp)) .or. &
        .not. ieee_is_nan(aimag(physical(2, 2))) .or. &
        .not. ieee_is_nan(residuals(2))) then
      error stop "Singular frequency point unavailable sentinel taşımıyor."
    end if
    if (.not. all(ieee_is_finite(real(physical(:, 1), dp))) .or. &
        .not. all(ieee_is_finite(real(physical(:, 3), dp)))) then
      error stop "Singular nokta komşu solved response'ları bozdu."
    end if
  end subroutine test_singular_sweep_status

  !> K', K'', C ve M'den kurulan Z matrisinin exact real/imag katsayılarını,
  !! complex symmetry/non-Hermitian ayrımını ve input immutability'yi sınar.
  !! Z=K'-omega^2*M+i(K''+omega*C) [N*m/rad].
  subroutine test_dynamic_stiffness_contract()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_dynamic_torsional_system_t) :: reduced_system
    type(stiffness_matrix_t) :: stiffness
    type(loss_stiffness_matrix_t) :: loss_stiffness
    type(damping_matrix_t) :: damping
    type(mass_matrix_t) :: mass
    type(dynamic_stiffness_matrix_t) :: dynamic_stiffness
    complex(dp), allocatable :: z(:, :)
    real(dp), allocatable :: stiffness_before(:, :)
    real(dp), allocatable :: loss_before(:, :)
    real(dp), allocatable :: damping_before(:, :)
    real(dp), allocatable :: mass_before(:, :)
    real(dp) :: angular_frequency

    call build_free_two_inertia_system( &
      system, manager, 0.1_dp, 0.2_dp, 100.0_dp, 10.0_dp, 1.0_dp)
    reduced_system = &
      build_reduced_dynamic_torsional_system(system, manager)
    stiffness = get_reduced_dynamic_stiffness(reduced_system)
    loss_stiffness = &
      get_reduced_dynamic_loss_stiffness(reduced_system)
    damping = get_reduced_dynamic_damping(reduced_system)
    mass = get_reduced_dynamic_mass(reduced_system)
    stiffness_before = get_stiffness_matrix_values(stiffness)
    loss_before = get_loss_stiffness_matrix_values(loss_stiffness)
    damping_before = get_damping_matrix_values(damping)
    mass_before = get_mass_matrix_values(mass)

    dynamic_stiffness = build_dynamic_stiffness( &
      stiffness, loss_stiffness, damping, mass, 3.0_dp)
    z = get_dynamic_stiffness_values(dynamic_stiffness)
    angular_frequency = 2.0_dp*pi*3.0_dp
    call assert_real_matrix_close(real(z, dp), &
      stiffness_before-angular_frequency**2*mass_before, &
      tight_tolerance, "Dynamic stiffness real(Z) yanlış.")
    call assert_real_matrix_close(aimag(z), &
      loss_before+angular_frequency*damping_before, &
      tight_tolerance, "Dynamic stiffness imag(Z) yanlış.")
    call assert_complex_matrix_close(z, transpose(z), tight_tolerance, &
      "Dynamic stiffness complex symmetric değil.")
    if (maxval(abs(z-conjg(transpose(z)))) <= 1.0e-12_dp) then
      error stop "Loss içeren dynamic stiffness yanlışlıkla Hermitian."
    end if

    call assert_real_matrix_close( &
      get_stiffness_matrix_values(stiffness), stiffness_before, 0.0_dp, &
      "Dynamic stiffness builder K' girdisini değiştirdi.")
    call assert_real_matrix_close( &
      get_loss_stiffness_matrix_values(loss_stiffness), loss_before, &
      0.0_dp, "Dynamic stiffness builder K'' girdisini değiştirdi.")
    call assert_real_matrix_close( &
      get_damping_matrix_values(damping), damping_before, 0.0_dp, &
      "Dynamic stiffness builder C girdisini değiştirdi.")
    call assert_real_matrix_close( &
      get_mass_matrix_values(mass), mass_before, 0.0_dp, &
      "Dynamic stiffness builder M girdisini değiştirdi.")
  end subroutine test_dynamic_stiffness_contract

  !> Aynı active DOF'a gelen complex torque katkılarının scatter-add ile
  !! toplandığını doğrular. Girdi ve çıktı peak torque [N*m] birimindedir.
  subroutine test_excitation_scatter_add()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_dynamic_torsional_system_t) :: reduced_system
    type(active_dof_map_t) :: mapping
    type(harmonic_excitation_t) :: excitations(2)
    complex(dp), allocatable :: load_vector(:)

    call build_fixed_one_dof_system( &
      system, manager, 100.0_dp, 0.0_dp, 0.0_dp, 1.0_dp)
    reduced_system = &
      build_reduced_dynamic_torsional_system(system, manager)
    mapping = get_reduced_dynamic_active_dof_map(reduced_system)
    excitations(1) = make_excitation( &
      2, cmplx(1.0_dp, 2.0_dp, kind=dp))
    excitations(2) = make_excitation( &
      2, cmplx(3.0_dp, -1.0_dp, kind=dp))
    load_vector = assemble_harmonic_load_vector(mapping, excitations)
    call assert_complex_close(load_vector(1), &
      cmplx(4.0_dp, 1.0_dp, kind=dp), 0.0_dp, &
      "Harmonic excitation scatter-add sonucu yanlış.")
  end subroutine test_excitation_scatter_add

  !> Element orientation, dynamic torque ve passive enerji bağıntılarını
  !! K''-only, c-only, combined ve lossless durumlarda sınar.
  !! P=(omega/2)(K''+omega*c)|Delta|^2 [W] ve
  !! E=pi(K''+omega*c)|Delta|^2 [J/cycle].
  subroutine test_passivity_and_element_results()
    type(torsional_element_t) :: element
    complex(dp) :: relative_angle
    complex(dp) :: expected_torque
    complex(dp) :: dynamic_torque
    real(dp) :: angular_frequency
    real(dp) :: expected_power
    real(dp) :: expected_energy
    real(dp) :: frequency_hz
    real(dp) :: power
    real(dp) :: energy

    frequency_hz = 8.0_dp
    angular_frequency = 2.0_dp*pi*frequency_hz
    relative_angle = calculate_element_relative_angle( &
      cmplx(0.2_dp, 0.1_dp, kind=dp), &
      cmplx(-0.1_dp, 0.3_dp, kind=dp))
    call assert_complex_close(relative_angle, &
      cmplx(0.3_dp, -0.2_dp, kind=dp), tight_tolerance, &
      "Element relative angle node_i-node_j orientation ile uyuşmuyor.")

    element = torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=500.0_dp, &
      damping_nms_per_rad=2.0_dp, &
      loss_stiffness_nm_per_rad=25.0_dp)
    dynamic_torque = calculate_element_dynamic_torque( &
      element, frequency_hz, relative_angle)
    expected_torque = cmplx(500.0_dp, &
      25.0_dp+angular_frequency*2.0_dp, kind=dp)*relative_angle
    call assert_complex_close(dynamic_torque, expected_torque, &
      tight_tolerance, "Element dynamic torque bağıntısı yanlış.")

    expected_power = 0.5_dp*angular_frequency * &
      (25.0_dp+angular_frequency*2.0_dp)*abs(relative_angle)**2
    expected_energy = pi*(25.0_dp+angular_frequency*2.0_dp) * &
      abs(relative_angle)**2
    power = calculate_element_average_dissipated_power( &
      element, frequency_hz, relative_angle)
    energy = calculate_element_dissipated_energy_per_cycle( &
      element, frequency_hz, relative_angle)
    call assert_real_close(power, expected_power, tight_tolerance, &
      "Element average dissipated power yanlış.")
    call assert_real_close(energy, expected_energy, tight_tolerance, &
      "Element dissipated energy-per-cycle yanlış.")
    if (power < 0.0_dp .or. energy < 0.0_dp) then
      error stop "Passive element negatif enerji üretti."
    end if

    ! Lossless eleman tam olarak sıfır dissipated power/energy üretmelidir.
    element%loss_stiffness_nm_per_rad = 0.0_dp
    element%damping_nms_per_rad = 0.0_dp
    if (calculate_element_average_dissipated_power( &
        element, frequency_hz, relative_angle) > 0.0_dp .or. &
        calculate_element_dissipated_energy_per_cycle( &
        element, frequency_hz, relative_angle) > 0.0_dp) then
      error stop "Lossless element dissipated enerji üretti."
    end if

    ! K''-only ve c-only kanalları ayrı ayrı passive kalmalıdır.
    element%loss_stiffness_nm_per_rad = 25.0_dp
    if (calculate_element_average_dissipated_power( &
        element, frequency_hz, relative_angle) <= 0.0_dp) then
      error stop "K''-only passive power pozitif değil."
    end if
    element%loss_stiffness_nm_per_rad = 0.0_dp
    element%damping_nms_per_rad = 2.0_dp
    if (calculate_element_average_dissipated_power( &
        element, frequency_hz, relative_angle) <= 0.0_dp) then
      error stop "c-only passive power pozitif değil."
    end if
  end subroutine test_passivity_and_element_results

  !> TVD element derived-result sweep yordamlarını gerçek harmonic çözümün
  !! physical response'u üzerinden doğrular. Relative angle node_i-node_j,
  !! internal torque complex peak ve transmitted magnitude abs(T_hat_e)'dir;
  !! power [W] ile energy [J/cycle] aynı production scalar yordamlarıyla
  !! karşılaştırılır. Bu test solver denklemini yeniden kopyalamaz.
  subroutine test_element_sweep_results()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_dynamic_torsional_system_t) :: reduced_system
    type(active_dof_map_t) :: mapping
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: response
    type(torsional_element_t) :: element
    complex(dp), allocatable :: physical(:, :)
    complex(dp), allocatable :: relative_angles(:)
    complex(dp), allocatable :: dynamic_torques(:)
    real(dp), allocatable :: torque_magnitudes(:)
    real(dp), allocatable :: average_power(:)
    real(dp), allocatable :: dissipated_energy(:)
    complex(dp) :: expected_relative_angle
    complex(dp) :: expected_dynamic_torque
    real(dp) :: frequencies(2)
    integer :: frequency_index

    call build_fixed_one_dof_system( &
      system, manager, 600.0_dp, 30.0_dp, 1.5_dp, 0.2_dp)
    reduced_system = &
      build_reduced_dynamic_torsional_system(system, manager)
    mapping = get_reduced_dynamic_active_dof_map(reduced_system)
    element = get_torsional_element(system, 1)
    excitation(1) = make_excitation( &
      2, cmplx(1.5_dp, -0.25_dp, kind=dp))
    frequencies = [4.0_dp, 6.0_dp]
    response = analyze_harmonic_response( &
      system, manager, frequencies, excitation)

    physical = get_physical_complex_response(response)
    relative_angles = get_element_relative_angle_response( &
      response, mapping, element)
    dynamic_torques = get_element_dynamic_torque_response( &
      response, mapping, element)
    torque_magnitudes = get_element_dynamic_torque_magnitudes( &
      response, mapping, element)
    average_power = get_element_average_dissipated_power_response( &
      response, mapping, element)
    dissipated_energy = get_element_dissipated_energy_response( &
      response, mapping, element)

    do frequency_index = 1, size(frequencies)
      expected_relative_angle = calculate_element_relative_angle( &
        physical(1, frequency_index), physical(2, frequency_index))
      expected_dynamic_torque = calculate_element_dynamic_torque( &
        element, frequencies(frequency_index), expected_relative_angle)
      call assert_complex_close( &
        relative_angles(frequency_index), expected_relative_angle, &
        tight_tolerance, "Element relative-angle sweep sonucu yanlış.")
      call assert_complex_close( &
        dynamic_torques(frequency_index), expected_dynamic_torque, &
        tight_tolerance, "Element dynamic-torque sweep sonucu yanlış.")
      call assert_real_close( &
        torque_magnitudes(frequency_index), abs(expected_dynamic_torque), &
        tight_tolerance, "Transmitted torque magnitude sonucu yanlış.")
      call assert_real_close(average_power(frequency_index), &
        calculate_element_average_dissipated_power( &
          element, frequencies(frequency_index), expected_relative_angle), &
        tight_tolerance, "Element power sweep sonucu yanlış.")
      call assert_real_close(dissipated_energy(frequency_index), &
        calculate_element_dissipated_energy_per_cycle( &
          element, frequencies(frequency_index), expected_relative_angle), &
        tight_tolerance, "Element energy sweep sonucu yanlış.")
    end do
  end subroutine test_element_sweep_results

  !> Tanımlı tek torque input channel için receptance, mobility ve
  !! accelerance helper'larını exp(+iwt) türevleriyle doğrular.
  subroutine test_rotational_frf_helpers()
    complex(dp) :: angular_response
    complex(dp) :: input_torque
    complex(dp) :: receptance
    real(dp) :: angular_frequency
    real(dp) :: frequency_hz

    angular_response = cmplx(0.02_dp, -0.01_dp, kind=dp)
    input_torque = cmplx(2.0_dp, 1.0_dp, kind=dp)
    frequency_hz = 6.0_dp
    angular_frequency = 2.0_dp*pi*frequency_hz
    receptance = angular_response/input_torque
    call assert_complex_close( &
      calculate_rotational_receptance(angular_response, input_torque), &
      receptance, tight_tolerance, "Rotational receptance yanlış.")
    call assert_complex_close( &
      calculate_rotational_mobility( &
        angular_response, input_torque, frequency_hz), &
      cmplx(0.0_dp, angular_frequency, kind=dp)*receptance, &
      tight_tolerance, "Rotational mobility yanlış.")
    call assert_complex_close( &
      calculate_rotational_accelerance( &
        angular_response, input_torque, frequency_hz), &
      -angular_frequency**2*receptance, tight_tolerance, &
      "Rotational accelerance yanlış.")
  end subroutine test_rotational_frf_helpers

  !> Düşük loss ile harmonic response peak'inin V0.5 DSYGV doğal frekansına
  !! yaklaşmasını sınar. Coarse grid peak hatası, refined grid ile azalmalıdır;
  !! exact equality veya mode-superposition kullanılmaz.
  subroutine test_modal_harmonic_cross_validation()
    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(reduced_torsional_system_t) :: modal_system
    type(modal_result_t) :: modal_result
    type(harmonic_excitation_t) :: excitation(1)
    type(harmonic_response_t) :: coarse_response
    type(harmonic_response_t) :: refined_response
    real(dp), allocatable :: modal_frequencies(:)
    real(dp), allocatable :: coarse_magnitudes(:, :)
    real(dp), allocatable :: refined_magnitudes(:, :)
    real(dp) :: coarse_frequencies(4)
    real(dp) :: refined_frequencies(5)
    real(dp) :: coarse_peak_frequency
    real(dp) :: refined_peak_frequency
    real(dp) :: natural_frequency

    call build_fixed_one_dof_system( &
      system, manager, 400.0_dp, 0.1_dp, 0.0_dp, 1.0_dp)
    modal_system = build_reduced_torsional_system(system, manager)
    modal_result = analyze_reduced_torsional_system(modal_system)
    modal_frequencies = get_modal_frequencies_hz(modal_result)
    if (size(modal_frequencies) /= 1) then
      error stop "Fixed 1-DOF V0.5 modal oracle tek mod döndürmedi."
    end if
    natural_frequency = modal_frequencies(1)

    coarse_frequencies = natural_frequency * &
      [0.75_dp, 0.90_dp, 1.05_dp, 1.20_dp]
    refined_frequencies = natural_frequency * &
      [0.94_dp, 0.97_dp, 1.00_dp, 1.03_dp, 1.06_dp]
    excitation(1) = make_excitation( &
      2, cmplx(1.0_dp, 0.0_dp, kind=dp))
    coarse_response = analyze_harmonic_response( &
      system, manager, coarse_frequencies, excitation)
    refined_response = analyze_harmonic_response( &
      system, manager, refined_frequencies, excitation)
    coarse_magnitudes = get_physical_response_magnitudes(coarse_response)
    refined_magnitudes = get_physical_response_magnitudes(refined_response)
    coarse_peak_frequency = coarse_frequencies( &
      maxloc(coarse_magnitudes(2, :), dim=1))
    refined_peak_frequency = refined_frequencies( &
      maxloc(refined_magnitudes(2, :), dim=1))

    if (abs(refined_peak_frequency-natural_frequency) > &
        abs(coarse_peak_frequency-natural_frequency)) then
      error stop "Grid refinement harmonic peak'i modal frekansa yaklaştırmadı."
    end if
    if (abs(refined_peak_frequency-natural_frequency) / natural_frequency > &
        0.031_dp) then
      error stop "Harmonic peak V0.5 doğal frekans çevresinde değil."
    end if
  end subroutine test_modal_harmonic_cross_validation

  !> Selector tabanlı invalid-input testlerini production public API üzerinden
  !! çalıştırır. Geçersiz frekans [Hz], excitation [N*m] veya fully constrained
  !! sistem error stop üretmelidir; CTest WILL_FAIL bu reddi başarı sayar.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_system_t) :: system
    type(constraint_manager_t) :: manager
    type(constraint_t) :: second_constraint
    type(harmonic_excitation_t), allocatable :: excitations(:)
    type(harmonic_response_t) :: rejected_response
    complex(dp) :: rejected_frf
    real(dp), allocatable :: frequencies(:)
    real(dp) :: invalid_value

    call build_fixed_one_dof_system( &
      system, manager, 100.0_dp, 1.0_dp, 0.5_dp, 1.0_dp)
    allocate(excitations(1))
    excitations(1) = make_excitation( &
      2, cmplx(1.0_dp, 0.0_dp, kind=dp))
    frequencies = [1.0_dp]

    select case (case_name)
      case ("empty_frequency")
        deallocate(frequencies)
        allocate(frequencies(0))
      case ("zero_frequency")
        frequencies = [0.0_dp]
      case ("negative_frequency")
        frequencies = [-1.0_dp]
      case ("nan_frequency")
        invalid_value = ieee_value(0.0_dp, ieee_quiet_nan)
        frequencies = [invalid_value]
      case ("positive_infinity_frequency")
        invalid_value = ieee_value(0.0_dp, ieee_positive_inf)
        frequencies = [invalid_value]
      case ("negative_infinity_frequency")
        invalid_value = ieee_value(0.0_dp, ieee_negative_inf)
        frequencies = [invalid_value]
      case ("unordered_frequency")
        frequencies = [2.0_dp, 1.0_dp]
      case ("duplicate_frequency")
        frequencies = [1.0_dp, 1.0_dp]
      case ("empty_excitation")
        deallocate(excitations)
        allocate(excitations(0))
      case ("unknown_node")
        excitations(1)%node_id = 999
      case ("unsupported_dof")
        excitations(1)%dof_type = 999
      case ("constrained_target")
        excitations(1)%node_id = 1
      case ("nan_torque_real")
        invalid_value = ieee_value(0.0_dp, ieee_quiet_nan)
        excitations(1)%torque_amplitude_nm = &
          cmplx(invalid_value, 0.0_dp, kind=dp)
      case ("nan_torque_imaginary")
        invalid_value = ieee_value(0.0_dp, ieee_quiet_nan)
        excitations(1)%torque_amplitude_nm = &
          cmplx(0.0_dp, invalid_value, kind=dp)
      case ("infinite_torque_real")
        invalid_value = ieee_value(0.0_dp, ieee_positive_inf)
        excitations(1)%torque_amplitude_nm = &
          cmplx(invalid_value, 0.0_dp, kind=dp)
      case ("infinite_torque_imaginary")
        invalid_value = ieee_value(0.0_dp, ieee_negative_inf)
        excitations(1)%torque_amplitude_nm = &
          cmplx(0.0_dp, invalid_value, kind=dp)
      case ("overflow_torque_sum")
        deallocate(excitations)
        allocate(excitations(2))
        excitations(1) = make_excitation( &
          2, cmplx(0.75_dp*huge(1.0_dp), 0.0_dp, kind=dp))
        excitations(2) = excitations(1)
      case ("fully_constrained")
        second_constraint = constraint_t( &
          constraint_id=2, node_id=2, dof_type=TORSIONAL_ROTATION, &
          value=0.0_dp, constraint_type=FIXED_CONSTRAINT)
        call add_constraint(manager, second_constraint, system)
      case ("zero_frf_torque")
        rejected_frf = calculate_rotational_receptance( &
          cmplx(1.0_dp, 0.0_dp, kind=dp), &
          cmplx(0.0_dp, 0.0_dp, kind=dp))
        print *, "Geçersiz FRF beklenmedik biçimde kabul edildi:", rejected_frf
        return
      case default
        error stop "Bilinmeyen harmonic invalid test selector'ı."
    end select

    rejected_response = analyze_harmonic_response( &
      system, manager, frequencies, excitations)
    print *, "Geçersiz harmonic girdi kabul edildi:", case_name, &
      get_harmonic_backend_identity(rejected_response)
  end subroutine exercise_invalid_case

  !> Fixed 1-DOF test fixture'ını üretir. Node 1 explicit fixed, node 2 polar
  !! atalet M [kg*m^2] taşır; eleman K'/K'' [N*m/rad] ve c [N*m*s/rad]
  !! değerleriyle iki düğümü bağlar.
  subroutine build_fixed_one_dof_system( &
      system, manager, stiffness, loss_stiffness, damping, inertia)
    type(torsional_system_t), intent(out) :: system
    type(constraint_manager_t), intent(out) :: manager
    real(dp), intent(in) :: stiffness
    real(dp), intent(in) :: loss_stiffness
    real(dp), intent(in) :: damping
    real(dp), intent(in) :: inertia

    type(constraint_t) :: fixed_constraint_record

    call add_torsional_node(system, torsional_node_t( &
      id=1, polar_inertia_kg_m2=0.5_dp, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=2, polar_inertia_kg_m2=inertia, constrained=.false.))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=stiffness, &
      damping_nms_per_rad=damping, &
      loss_stiffness_nm_per_rad=loss_stiffness))
    call initialize_constraint_manager(manager)
    fixed_constraint_record = constraint_t( &
      constraint_id=1, node_id=1, dof_type=TORSIONAL_ROTATION, &
      value=0.0_dp, constraint_type=FIXED_CONSTRAINT)
    call add_constraint(manager, fixed_constraint_record, system)
  end subroutine build_fixed_one_dof_system

  !> Serbest-serbest iki ataletli harmonic test fixture'ını kurar. J1/J2
  !! [kg*m^2], K'/K'' [N*m/rad] ve c [N*m*s/rad] SI birimlerindedir.
  subroutine build_free_two_inertia_system( &
      system, manager, inertia_1, inertia_2, stiffness, loss_stiffness, damping)
    type(torsional_system_t), intent(out) :: system
    type(constraint_manager_t), intent(out) :: manager
    real(dp), intent(in) :: inertia_1
    real(dp), intent(in) :: inertia_2
    real(dp), intent(in) :: stiffness
    real(dp), intent(in) :: loss_stiffness
    real(dp), intent(in) :: damping

    call add_torsional_node(system, torsional_node_t( &
      id=1, polar_inertia_kg_m2=inertia_1, constrained=.false.))
    call add_torsional_node(system, torsional_node_t( &
      id=2, polar_inertia_kg_m2=inertia_2, constrained=.false.))
    call add_torsional_element(system, torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=stiffness, &
      damping_nms_per_rad=damping, &
      loss_stiffness_nm_per_rad=loss_stiffness))
    call initialize_constraint_manager(manager)
  end subroutine build_free_two_inertia_system

  !> Complex peak nodal torque fixture'ını [N*m] üretir.
  pure function make_excitation(node_id, torque) result(excitation)
    integer, intent(in) :: node_id
    complex(dp), intent(in) :: torque
    type(harmonic_excitation_t) :: excitation

    excitation = harmonic_excitation_t( &
      node_id=node_id, dof_type=TORSIONAL_ROTATION, &
      torque_amplitude_nm=torque)
  end function make_excitation

  !> Complex scalar değerleri scale-aware mutlak/bağıl toleransla kıyaslar.
  subroutine assert_complex_close(actual, expected, tolerance, message)
    complex(dp), intent(in) :: actual
    complex(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(real(actual, dp)) .or. &
        .not. ieee_is_finite(aimag(actual)) .or. &
        abs(actual-expected) > tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_complex_close

  !> Gerçek scalar değerleri scale-aware mutlak/bağıl toleransla kıyaslar.
  subroutine assert_real_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        abs(actual-expected) > tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_real_close

  !> Gerçek matrisleri aynı boyutta scale-aware toleransla kıyaslar.
  subroutine assert_real_matrix_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:, :)
    real(dp), intent(in) :: expected(:, :)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (any(shape(actual) /= shape(expected))) error stop message
    if (.not. all(ieee_is_finite(actual)) .or. &
        maxval(abs(actual-expected)) > &
        tolerance*max(1.0_dp, maxval(abs(expected)))) then
      error stop message
    end if
  end subroutine assert_real_matrix_close

  !> Complex matrisleri aynı boyutta scale-aware toleransla kıyaslar.
  subroutine assert_complex_matrix_close(actual, expected, tolerance, message)
    complex(dp), intent(in) :: actual(:, :)
    complex(dp), intent(in) :: expected(:, :)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (any(shape(actual) /= shape(expected))) error stop message
    if (.not. all(ieee_is_finite(real(actual, dp))) .or. &
        .not. all(ieee_is_finite(aimag(actual))) .or. &
        maxval(abs(actual-expected)) > &
        tolerance*max(1.0_dp, maxval(abs(expected)))) then
      error stop message
    end if
  end subroutine assert_complex_matrix_close

  !> Sayısal tanı dizisinin sonlu ve negatif olmayan olduğunu doğrular.
  subroutine assert_finite_nonnegative_diagnostics(values, label)
    real(dp), intent(in) :: values(:)
    character(len=*), intent(in) :: label

    if (.not. all(ieee_is_finite(values)) .or. any(values < 0.0_dp)) then
      error stop label//" değerleri sonlu ve negatif olmayan olmalıdır."
    end if
  end subroutine assert_finite_nonnegative_diagnostics

end program test_harmonic_analysis
