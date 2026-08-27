program test_tabulated_dynamic_modulus_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
    ieee_positive_inf
  use tms_kinds, only : dp
  use tms_material_frequency, only : material_frequency_point
  use tms_dynamic_modulus, only : calculate_loss_factor
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    DYNAMIC_DEFORMATION_MODE_SHEAR
  use tms_dynamic_modulus_provider, only : dynamic_modulus_evaluation_t, &
    LINEAR_FREQUENCY, LINEAR_LOG_FREQUENCY, &
    evaluate_dynamic_shear_modulus, get_dynamic_modulus_provider_metadata
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t, &
    create_tabulated_dynamic_modulus_provider, &
    get_tabulated_frequency_points, get_tabulated_interpolation_policy
  implicit none

  real(dp), parameter :: tolerance = 2.0e-13_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_exact_points_and_machine_match()
  call test_linear_frequency_interpolation()
  call test_linear_log_frequency_interpolation()
  call test_independent_storage_and_metadata()

  print *, "V0.7 tabulated dynamic shear-modulus provider doğrulandı."

contains

  !> Stored frequency noktasında her iki policy'nin primary G'/G'' değerini
  !! aynen döndürdüğünü ve bir ULP farkın yalnız representation tolerance
  !! kapsamında exact kabul edildiğini doğrular. G'/G'' [Pa], f [Hz], T [K].
  subroutine test_exact_points_and_machine_match()
    type(material_frequency_point) :: points(3)
    type(dynamic_material_metadata_t) :: metadata
    type(tabulated_dynamic_modulus_provider_t) :: linear_provider
    type(tabulated_dynamic_modulus_provider_t) :: log_provider
    type(dynamic_modulus_evaluation_t) :: evaluation

    points = make_three_points()
    metadata = make_metadata()
    linear_provider = create_tabulated_dynamic_modulus_provider( &
      points, metadata)
    log_provider = create_tabulated_dynamic_modulus_provider( &
      points, metadata, LINEAR_LOG_FREQUENCY)

    if (get_tabulated_interpolation_policy(linear_provider) /= &
        LINEAR_FREQUENCY) then
      error stop "Default interpolation policy LINEAR_FREQUENCY değil."
    end if
    evaluation = evaluate_dynamic_shear_modulus( &
      linear_provider, 20.0_dp, 293.15_dp)
    call assert_exact_evaluation(evaluation, 2.0e6_dp, 2.0e5_dp)
    evaluation = evaluate_dynamic_shear_modulus( &
      log_provider, 20.0_dp, 293.15_dp)
    call assert_exact_evaluation(evaluation, 2.0e6_dp, 2.0e5_dp)

    evaluation = evaluate_dynamic_shear_modulus( &
      linear_provider, nearest(20.0_dp, 1.0_dp), &
      nearest(293.15_dp, 1.0_dp))
    call assert_exact_evaluation(evaluation, 2.0e6_dp, 2.0e5_dp)

    ! Machine tolerance'dan belirgin büyük, fakat fiziksel olarak çok küçük
    ! frekans farkı exact-point semantiğini kullanmamalı; interpolation'dır.
    evaluation = evaluate_dynamic_shear_modulus( &
      linear_provider, 20.0_dp*(1.0_dp+1.0e-12_dp), 293.15_dp)
    if (evaluation%exact_table_point) then
      error stop "Physical-level frequency farkı exact match kabul edildi."
    end if
  end subroutine test_exact_points_and_machine_match

  !> LINEAR_FREQUENCY scaling invariant'ını üretim provider sonucu üzerinden
  !! sınar. 10-20 Hz arasında 15 Hz için alpha=0.5 ve G bileşenleri endpoint
  !! ortalamasıdır; tan(delta) ayrıca G''/G' olarak türetilir.
  subroutine test_linear_frequency_interpolation()
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: evaluation

    provider = create_tabulated_dynamic_modulus_provider( &
      make_three_points(), make_metadata(), LINEAR_FREQUENCY)
    evaluation = evaluate_dynamic_shear_modulus( &
      provider, 15.0_dp, 293.15_dp)

    if (evaluation%exact_table_point) then
      error stop "Ara frekans exact table point olarak işaretlendi."
    end if
    call assert_close(evaluation%interpolation_alpha, 0.5_dp, tolerance, &
      "Linear-frequency alpha hatalı.")
    call assert_close(evaluation%modulus%storage_modulus, &
      1.5e6_dp, tolerance, "Linear-frequency G' hatalı.")
    call assert_close(evaluation%modulus%loss_modulus, &
      1.5e5_dp, tolerance, "Linear-frequency G'' hatalı.")
    call assert_close(calculate_loss_factor(evaluation%modulus), &
      0.1_dp, tolerance, "Derived tan(delta) hatalı.")
  end subroutine test_linear_frequency_interpolation

  !> LINEAR_LOG_FREQUENCY modelinde yalnız frequency axis'in log10 alındığını
  !! doğrular. 10 ve 100 Hz arasındaki sqrt(1000) Hz log-midpoint'te alpha=0.5;
  !! modüllerin logaritması alınmadan G'/G'' lineer ortalanır.
  subroutine test_linear_log_frequency_interpolation()
    type(material_frequency_point) :: points(2)
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: evaluation

    points(1) = make_point(10.0_dp, 1.0e6_dp, 0.1e6_dp)
    points(2) = make_point(100.0_dp, 5.0e6_dp, 0.5e6_dp)
    provider = create_tabulated_dynamic_modulus_provider( &
      points, make_metadata(), LINEAR_LOG_FREQUENCY)
    evaluation = evaluate_dynamic_shear_modulus( &
      provider, sqrt(1000.0_dp), 293.15_dp)

    call assert_close(evaluation%interpolation_alpha, 0.5_dp, tolerance, &
      "Log-frequency alpha hatalı.")
    call assert_close(evaluation%modulus%storage_modulus, &
      3.0e6_dp, tolerance, "Log-frequency G' hatalı.")
    call assert_close(evaluation%modulus%loss_modulus, &
      0.3e6_dp, tolerance, "Log-frequency G'' hatalı.")
  end subroutine test_linear_log_frequency_interpolation

  !> Provider'ın authoritative input points/metadata için independent copy
  !! tuttuğunu doğrular. Caller ve getter kopyalarının mutation'ı provider
  !! sorgusundaki G'/G'' [Pa] veya dataset kimliğini değiştirmemelidir.
  subroutine test_independent_storage_and_metadata()
    type(material_frequency_point) :: points(3)
    type(material_frequency_point), allocatable :: copied_points(:)
    type(dynamic_material_metadata_t) :: metadata
    type(dynamic_material_metadata_t) :: copied_metadata
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: evaluation

    points = make_three_points()
    metadata = make_metadata()
    provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    points(2)%storage_modulus = 99.0e6_dp
    metadata%dataset_identifier = "caller-mutated"
    copied_points = get_tabulated_frequency_points(provider)
    copied_points(2)%storage_modulus = 88.0e6_dp
    copied_metadata = get_dynamic_modulus_provider_metadata(provider)
    copied_metadata%dataset_identifier = "getter-mutated"

    evaluation = evaluate_dynamic_shear_modulus( &
      provider, 20.0_dp, 293.15_dp)
    call assert_close(evaluation%modulus%storage_modulus, &
      2.0e6_dp, 0.0_dp, "Provider private point kopyası değişti.")
    copied_metadata = get_dynamic_modulus_provider_metadata(provider)
    if (trim(copied_metadata%dataset_identifier) /= "EPDM-ISO-293K") then
      error stop "Provider private metadata kopyası değişti."
    end if
  end subroutine test_independent_storage_and_metadata

  !> Error-stop regresyon selector'ları provider data-quality, passivity,
  !! isotherm ve no-extrapolation sözleşmelerini ayrı süreçlerde tetikler.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(material_frequency_point) :: points(3)
    type(material_frequency_point) :: one_point(1)
    type(dynamic_material_metadata_t) :: metadata
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_modulus_evaluation_t) :: evaluation
    real(dp) :: nan_value
    real(dp) :: infinity

    points = make_three_points()
    one_point(1) = points(1)
    metadata = make_metadata()
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    select case (case_name)
    case ("less_than_two_points")
      provider = create_tabulated_dynamic_modulus_provider( &
        one_point, metadata)
    case ("zero_frequency")
      points(1)%frequency = 0.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("negative_frequency")
      points(1)%frequency = -1.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("nan_frequency")
      points(1)%frequency = nan_value
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("infinite_frequency")
      points(3)%frequency = infinity
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("duplicate_frequency")
      points(2)%frequency = points(1)%frequency
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("numerical_duplicate_frequency")
      points(2)%frequency = nearest(points(1)%frequency, 1.0_dp)
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("unordered_frequency")
      points(2)%frequency = 5.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("zero_storage_modulus")
      points(2)%storage_modulus = 0.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("negative_storage_modulus")
      points(2)%storage_modulus = -1.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("nan_storage_modulus")
      points(2)%storage_modulus = nan_value
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("infinite_storage_modulus")
      points(2)%storage_modulus = infinity
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("negative_loss_modulus")
      points(2)%loss_modulus = -1.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("nan_loss_modulus")
      points(2)%loss_modulus = nan_value
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("infinite_loss_modulus")
      points(2)%loss_modulus = infinity
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("invalid_dataset_temperature")
      metadata%dataset_temperature_k = 0.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("nan_dataset_temperature")
      metadata%dataset_temperature_k = nan_value
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("infinite_dataset_temperature")
      metadata%dataset_temperature_k = infinity
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("mixed_point_temperature")
      points(2)%temperature = 313.15_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("invalid_policy")
      provider = create_tabulated_dynamic_modulus_provider( &
        points, metadata, 999)
    case ("below_linear_domain")
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
      evaluation = provider%evaluate(5.0_dp, 293.15_dp)
    case ("above_linear_domain")
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
      evaluation = provider%evaluate(101.0_dp, 293.15_dp)
    case ("below_log_domain")
      provider = create_tabulated_dynamic_modulus_provider( &
        points, metadata, LINEAR_LOG_FREQUENCY)
      evaluation = provider%evaluate(5.0_dp, 293.15_dp)
    case ("above_log_domain")
      provider = create_tabulated_dynamic_modulus_provider( &
        points, metadata, LINEAR_LOG_FREQUENCY)
      evaluation = provider%evaluate(101.0_dp, 293.15_dp)
    case ("temperature_mismatch")
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
      evaluation = provider%evaluate(20.0_dp, 303.15_dp)
    case ("invalid_strain_amplitude")
      metadata%dynamic_shear_strain_amplitude = 0.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("invalid_prestrain")
      metadata%static_shear_prestrain = -0.01_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("magic_unavailable_strain")
      metadata%has_dynamic_shear_strain_amplitude = .false.
      metadata%dynamic_shear_strain_amplitude = -999.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("magic_unavailable_prestrain")
      metadata%has_static_shear_prestrain = .false.
      metadata%static_shear_prestrain = -999.0_dp
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case ("invalid_mode")
      metadata%deformation_mode = 999
      provider = create_tabulated_dynamic_modulus_provider(points, metadata)
    case default
      error stop "Bilinmeyen tabulated provider validation selector."
    end select
  end subroutine exercise_invalid_case

  pure function make_three_points() result(points)
    type(material_frequency_point) :: points(3)

    points(1) = make_point(10.0_dp, 1.0e6_dp, 0.1e6_dp)
    points(2) = make_point(20.0_dp, 2.0e6_dp, 0.2e6_dp)
    points(3) = make_point(100.0_dp, 4.0e6_dp, 0.4e6_dp)
  end function make_three_points

  pure function make_point(frequency_hz, storage_pa, loss_pa) result(point)
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: storage_pa
    real(dp), intent(in) :: loss_pa
    type(material_frequency_point) :: point

    point%frequency = frequency_hz
    point%temperature = 293.15_dp
    point%storage_modulus = storage_pa
    point%loss_modulus = loss_pa
  end function make_point

  pure function make_metadata() result(metadata)
    type(dynamic_material_metadata_t) :: metadata

    metadata%dataset_identifier = "EPDM-ISO-293K"
    metadata%material_identifier = "EPDM-DEMO"
    metadata%has_specimen_identifier = .true.
    metadata%specimen_identifier = "S-01"
    metadata%dataset_temperature_k = 293.15_dp
    metadata%has_dynamic_shear_strain_amplitude = .true.
    metadata%dynamic_shear_strain_amplitude = 0.005_dp
    metadata%has_static_shear_prestrain = .true.
    metadata%static_shear_prestrain = 0.02_dp
    metadata%deformation_mode = DYNAMIC_DEFORMATION_MODE_SHEAR
    metadata%has_conditioning_state = .true.
    metadata%conditioning_state = "10 precycle"
    metadata%has_test_method_source = .true.
    metadata%test_method_source = "DMA torsional shear"
    metadata%has_standard_reference = .true.
    metadata%standard_reference = "ASTM D5992 / ISO 4664-1 trace"
  end function make_metadata

  subroutine assert_exact_evaluation(evaluation, expected_storage, expected_loss)
    type(dynamic_modulus_evaluation_t), intent(in) :: evaluation
    real(dp), intent(in) :: expected_storage
    real(dp), intent(in) :: expected_loss

    if (.not. evaluation%exact_table_point) then
      error stop "Stored frequency numerical exact point olarak bulunamadı."
    end if
    call assert_close(evaluation%modulus%storage_modulus, &
      expected_storage, 0.0_dp, "Exact G' stored değeri değişti.")
    call assert_close(evaluation%modulus%loss_modulus, &
      expected_loss, 0.0_dp, "Exact G'' stored değeri değişti.")
  end subroutine assert_exact_evaluation

  subroutine assert_close(actual, expected, relative_tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: relative_tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > &
        relative_tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_close

end program test_tabulated_dynamic_modulus_provider
