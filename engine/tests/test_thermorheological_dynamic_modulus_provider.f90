program test_thermorheological_dynamic_modulus_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan, ieee_positive_inf
  use tms_kinds, only : dp
  use tms_material_frequency, only : material_frequency_point
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    DYNAMIC_DEFORMATION_MODE_SHEAR
  use tms_dynamic_modulus_provider, only : dynamic_modulus_evaluation_t, &
    LINEAR_LOG_FREQUENCY, evaluate_dynamic_shear_modulus, &
    get_dynamic_modulus_provider_metadata
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t, &
    create_tabulated_dynamic_modulus_provider
  use tms_tabulated_temperature_shift, only : &
    tabulated_log10_shift_provider_t, &
    create_tabulated_temperature_shift_provider
  use tms_thermorheological_dynamic_modulus_provider, only : &
    thermorheological_dynamic_modulus_provider_t, &
    create_thermorheological_dynamic_modulus_provider
  use tms_temperature_shift_types, only : TABULATED_LOG10_SHIFT
  implicit none

  real(dp), parameter :: tolerance = 3.0e-12_dp
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_reference_temperature_identity()
  call test_physical_reduced_coordinates_and_interpolation()
  call test_passivity_and_independent_storage()

  print *, "V0.8 thermorheological dynamic-modulus provider doğrulandı."

contains

  !> T=T_ref için a_T=1 ve f_r=f olduğundan thermorheological provider'ın
  !! V0.7 master provider ile aynı G'/G'' [Pa] değerini döndürdüğünü doğrular.
  !! Returned modulus physical f [Hz], operating T [K] semantiğini korur.
  subroutine test_reference_temperature_identity()
    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    type(tabulated_log10_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: master_evaluation
    type(dynamic_modulus_evaluation_t) :: evaluation
    type(dynamic_material_metadata_t) :: metadata

    master_curve = make_master_curve("IDENTITY", reference_temperature_k)
    shift_provider = make_shift_provider()
    provider = create_thermorheological_dynamic_modulus_provider( &
      master_curve, shift_provider)
    master_evaluation = evaluate_dynamic_shear_modulus( &
      master_curve, 10.0_dp, reference_temperature_k)
    evaluation = evaluate_dynamic_shear_modulus( &
      provider, 10.0_dp, reference_temperature_k)

    call assert_close(evaluation%modulus%storage_modulus, &
      master_evaluation%modulus%storage_modulus, 0.0_dp, &
      "T_ref identity G' regression sonucu değişti.")
    call assert_close(evaluation%modulus%loss_modulus, &
      master_evaluation%modulus%loss_modulus, 0.0_dp, &
      "T_ref identity G'' regression sonucu değişti.")
    call assert_close(evaluation%physical_frequency_hz, &
      10.0_dp, 0.0_dp, "T_ref physical frequency hatalı.")
    call assert_close(evaluation%lookup_frequency_hz, &
      10.0_dp, 0.0_dp, "T_ref lookup frequency physical f'ye eşit değil.")
    call assert_close(evaluation%modulus%frequency, &
      10.0_dp, 0.0_dp, "Returned modulus physical frequency taşımıyor.")
    call assert_close(evaluation%modulus%temperature, &
      reference_temperature_k, 0.0_dp, &
      "Returned modulus operating sıcaklığı taşımıyor.")
    call assert_close(evaluation%a_t, 1.0_dp, 0.0_dp, &
      "T_ref thermorheological a_T değeri bir değil.")
    if (.not. evaluation%temperature_shift_applied) then
      error stop "Thermorheological evaluation shift context'i taşımıyor."
    end if

    metadata = get_dynamic_modulus_provider_metadata(provider)
    call assert_close(metadata%dataset_temperature_k, &
      reference_temperature_k, 0.0_dp, &
      "Thermorheological metadata operating T ile yeniden tanımlandı.")
  end subroutine test_reference_temperature_identity

  !> Physical f=100 Hz ve tabulated shift'te a_T=0.1 için reduced lookup
  !! f_r=10 Hz olduğunu doğrular. İkinci sorgu sıcaklık ve log-frequency
  !! midpoint'lerinde alpha=0.5 üretir; bracket physical 100 Hz'e göre değil,
  !! master curve reduced-frequency eksenine göre tutulur.
  subroutine test_physical_reduced_coordinates_and_interpolation()
    type(thermorheological_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: exact_evaluation
    type(dynamic_modulus_evaluation_t) :: interpolated_evaluation

    provider = make_thermorheological_provider("TRACE")
    exact_evaluation = evaluate_dynamic_shear_modulus( &
      provider, 100.0_dp, 313.15_dp)

    call assert_close(exact_evaluation%physical_frequency_hz, &
      100.0_dp, tolerance, "Physical frequency trace hatalı.")
    call assert_close(exact_evaluation%lookup_frequency_hz, &
      10.0_dp, tolerance, "Reduced lookup frequency f_r=a_T*f değil.")
    call assert_close(exact_evaluation%modulus%frequency, &
      100.0_dp, tolerance, "Returned modulus reduced frequency taşıyor.")
    call assert_close(exact_evaluation%modulus%temperature, &
      313.15_dp, tolerance, "Returned modulus reference T taşıyor.")
    call assert_close(exact_evaluation%log10_a_t, &
      -1.0_dp, tolerance, "Shift trace log10(a_T) hatalı.")
    call assert_close(exact_evaluation%a_t, &
      0.1_dp, tolerance, "Shift trace a_T hatalı.")
    if (exact_evaluation%shift_model_kind /= TABULATED_LOG10_SHIFT .or. &
        .not. exact_evaluation%shift_exact_temperature_point .or. &
        .not. exact_evaluation%exact_table_point) then
      error stop "Exact temperature/master point trace semantiği hatalı."
    end if
    call assert_close(exact_evaluation%lower_frequency_hz, &
      10.0_dp, tolerance, "Master exact lower bracket physical f kullanıyor.")
    call assert_close(exact_evaluation%upper_frequency_hz, &
      10.0_dp, tolerance, "Master exact upper bracket physical f kullanıyor.")

    interpolated_evaluation = evaluate_dynamic_shear_modulus( &
      provider, 100.0_dp, 303.15_dp)
    call assert_close(interpolated_evaluation%log10_a_t, &
      -0.5_dp, tolerance, "Temperature shift midpoint s değeri hatalı.")
    call assert_close(interpolated_evaluation%lookup_frequency_hz, &
      sqrt(1000.0_dp), tolerance, "Midpoint reduced frequency hatalı.")
    call assert_close(interpolated_evaluation%interpolation_alpha, &
      0.5_dp, tolerance, "Master log-frequency alpha hatalı.")
    call assert_close( &
      interpolated_evaluation%temperature_interpolation_alpha, &
      0.5_dp, tolerance, "Temperature interpolation alpha hatalı.")
    call assert_close(interpolated_evaluation%lower_frequency_hz, &
      10.0_dp, tolerance, "Reduced lower master bracket hatalı.")
    call assert_close(interpolated_evaluation%upper_frequency_hz, &
      100.0_dp, tolerance, "Reduced upper master bracket hatalı.")
    call assert_close(interpolated_evaluation%modulus%storage_modulus, &
      3.0e6_dp, tolerance, "Shifted master-curve G' interpolation'ı hatalı.")
    call assert_close(interpolated_evaluation%modulus%loss_modulus, &
      0.3e6_dp, tolerance, "Shifted master-curve G'' interpolation'ı hatalı.")
  end subroutine test_physical_reduced_coordinates_and_interpolation

  !> Horizontal shifting'in yalnız validated master curve'den değer seçtiğini,
  !! G'>0 ve G''>=0 passivity koşullarını koruduğunu sınar. Constructor'ın
  !! master/shift inputlarını independent copy olarak sakladığını da doğrular.
  subroutine test_passivity_and_independent_storage()
    type(material_frequency_point) :: points(3)
    type(dynamic_material_metadata_t) :: metadata
    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    type(tabulated_log10_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: evaluation
    real(dp) :: temperatures(3)
    real(dp) :: shifts(3)

    points = make_master_points(reference_temperature_k)
    metadata = make_metadata("PRIVATE-COPY", reference_temperature_k)
    master_curve = create_tabulated_dynamic_modulus_provider( &
      points, metadata, LINEAR_LOG_FREQUENCY)
    temperatures = [273.15_dp, 293.15_dp, 313.15_dp]
    shifts = [1.0_dp, 0.0_dp, -1.0_dp]
    shift_provider = create_tabulated_temperature_shift_provider( &
      temperatures, shifts, reference_temperature_k)
    provider = create_thermorheological_dynamic_modulus_provider( &
      master_curve, shift_provider)

    points(2)%storage_modulus = 99.0e6_dp
    metadata%dataset_identifier = "caller-mutated"
    shifts(3) = -9.0_dp
    evaluation = evaluate_dynamic_shear_modulus( &
      provider, 100.0_dp, 313.15_dp)

    call assert_close(evaluation%modulus%storage_modulus, &
      2.0e6_dp, tolerance, "Thermorheological private master kopyası değişti.")
    if (evaluation%modulus%storage_modulus <= 0.0_dp .or. &
        evaluation%modulus%loss_modulus < 0.0_dp) then
      error stop "Horizontal shift master-curve passivity'sini bozdu."
    end if
    metadata = get_dynamic_modulus_provider_metadata(provider)
    if (trim(metadata%dataset_identifier) /= "PRIVATE-COPY") then
      error stop "Thermorheological private metadata kopyası değişti."
    end if
  end subroutine test_passivity_and_independent_storage

  !> Reference-temperature uyumu, dual domain, finite physical input ve
  !! log-space reduced-frequency prevalidation error-stop yollarını tetikler.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    type(tabulated_log10_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: evaluation
    real(dp) :: temperatures(3)
    real(dp) :: shifts(3)
    real(dp) :: nan_value
    real(dp) :: infinity

    master_curve = make_master_curve("INVALID-THERMO", &
      reference_temperature_k)
    shift_provider = make_shift_provider()
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    select case (case_name)
    case ("reference_temperature_mismatch")
      temperatures = [278.15_dp, 298.15_dp, 318.15_dp]
      shifts = [1.0_dp, 0.0_dp, -1.0_dp]
      shift_provider = create_tabulated_temperature_shift_provider( &
        temperatures, shifts, 298.15_dp)
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
    case ("temperature_below_domain")
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, 10.0_dp, 273.0_dp)
    case ("temperature_above_domain")
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, 10.0_dp, 313.2_dp)
    case ("reduced_frequency_below_domain")
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, 1.0_dp, 313.15_dp)
    case ("reduced_frequency_above_domain")
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, 100.0_dp, 273.15_dp)
    case ("zero_physical_frequency")
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, 0.0_dp, reference_temperature_k)
    case ("nan_physical_frequency")
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, nan_value, reference_temperature_k)
    case ("infinite_physical_frequency")
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, infinity, reference_temperature_k)
    case ("logspace_overflow_outside_domain")
      temperatures = [273.15_dp, 293.15_dp, 313.15_dp]
      shifts = [-300.0_dp, 0.0_dp, 300.0_dp]
      shift_provider = create_tabulated_temperature_shift_provider( &
        temperatures, shifts, 293.15_dp)
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, 1.0e100_dp, 313.15_dp)
    case ("logspace_underflow_outside_domain")
      ! Reduced frequency doğrudan 10**s*f çarpımıyla kurulmadan önce
      ! log10(f_r)=log10(f)+s domain dışında güvenle reddedilmelidir.
      temperatures = [273.15_dp, 293.15_dp, 313.15_dp]
      shifts = [300.0_dp, 0.0_dp, -300.0_dp]
      shift_provider = create_tabulated_temperature_shift_provider( &
        temperatures, shifts, 293.15_dp)
      provider = create_thermorheological_dynamic_modulus_provider( &
        master_curve, shift_provider)
      evaluation = evaluate_dynamic_shear_modulus( &
        provider, 1.0e-100_dp, 313.15_dp)
    case default
      error stop "Bilinmeyen thermorheological provider validation selector."
    end select
  end subroutine exercise_invalid_case

  function make_thermorheological_provider(dataset_id) result(provider)
    character(len=*), intent(in) :: dataset_id
    type(thermorheological_dynamic_modulus_provider_t) :: provider

    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    type(tabulated_log10_shift_provider_t) :: shift_provider

    master_curve = make_master_curve(dataset_id, reference_temperature_k)
    shift_provider = make_shift_provider()
    provider = create_thermorheological_dynamic_modulus_provider( &
      master_curve, shift_provider)
  end function make_thermorheological_provider

  pure function make_master_curve(dataset_id, temperature_k) result(provider)
    character(len=*), intent(in) :: dataset_id
    real(dp), intent(in) :: temperature_k
    type(tabulated_dynamic_modulus_provider_t) :: provider

    provider = create_tabulated_dynamic_modulus_provider( &
      make_master_points(temperature_k), &
      make_metadata(dataset_id, temperature_k), LINEAR_LOG_FREQUENCY)
  end function make_master_curve

  pure function make_shift_provider() result(provider)
    type(tabulated_log10_shift_provider_t) :: provider
    real(dp) :: temperature_k(3)
    real(dp) :: log10_a_t(3)

    temperature_k = [273.15_dp, 293.15_dp, 313.15_dp]
    log10_a_t = [1.0_dp, 0.0_dp, -1.0_dp]
    provider = create_tabulated_temperature_shift_provider( &
      temperature_k, log10_a_t, reference_temperature_k)
  end function make_shift_provider

  pure function make_master_points(temperature_k) result(points)
    real(dp), intent(in) :: temperature_k
    type(material_frequency_point) :: points(3)

    points(1) = make_point(1.0_dp, 1.0e6_dp, 0.1e6_dp, temperature_k)
    points(2) = make_point(10.0_dp, 2.0e6_dp, 0.2e6_dp, temperature_k)
    points(3) = make_point(100.0_dp, 4.0e6_dp, 0.4e6_dp, temperature_k)
  end function make_master_points

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

end program test_thermorheological_dynamic_modulus_provider
