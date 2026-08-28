program test_arrhenius_temperature_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan, ieee_positive_inf
  use tms_kinds, only : dp
  use tms_temperature_shift_types, only : temperature_shift_evaluation_t, &
    ARRHENIUS_TEMPERATURE_SHIFT
  use tms_temperature_shift_provider, only : evaluate_temperature_shift
  use tms_arrhenius_temperature_shift, only : &
    arrhenius_temperature_shift_provider_t, &
    create_arrhenius_temperature_shift_provider
  implicit none

  real(dp), parameter :: tolerance = 3.0e-13_dp
  real(dp), parameter :: activation_energy = 50000.0_dp
  real(dp), parameter :: gas_constant_reference = 8.31446261815324_dp
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: ln_10 = log(10.0_dp)
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_reference_identity()
  call test_analytical_values_and_direction()

  print *, "V0.8 Arrhenius temperature-shift provider doğrulandı."

contains

  !> T=T_ref için Arrhenius reciprocal-temperature farkının sıfır olduğunu;
  !! authoritative s=0 ve derived a_T=1 identity koşullarını doğrular.
  subroutine test_reference_identity()
    type(arrhenius_temperature_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: evaluation

    provider = make_provider()
    evaluation = evaluate_temperature_shift( &
      provider, reference_temperature_k)

    if (evaluation%shift_model_kind /= ARRHENIUS_TEMPERATURE_SHIFT) then
      error stop "Arrhenius shift model kind trace'i hatalı."
    end if
    if (evaluation%has_temperature_bracket) then
      error stop "Analitik Arrhenius modeli için yapay table bracket üretildi."
    end if
    call assert_close(evaluation%log10_a_t, 0.0_dp, 0.0_dp, &
      "Arrhenius reference sıcaklığında log10(a_T) sıfır değil.")
    call assert_close(evaluation%a_t, 1.0_dp, 0.0_dp, &
      "Arrhenius reference sıcaklığında a_T bir değil.")
  end subroutine test_reference_identity

  !> ln(a_T)=Ea/R*(1/T-1/T_ref) eşitliğini bağımsız sabit ve arithmetic ile
  !! doğrular. Pozitif Ea için T>T_ref yönünde a_T<1; T<T_ref yönünde a_T>1
  !! olması canonical f_r=a_T*f davranışını kilitler.
  subroutine test_analytical_values_and_direction()
    type(arrhenius_temperature_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: hot
    type(temperature_shift_evaluation_t) :: cold
    real(dp) :: expected_hot_ln
    real(dp) :: expected_cold_ln
    real(dp), parameter :: physical_frequency_hz = 100.0_dp

    provider = make_provider()
    hot = evaluate_temperature_shift(provider, 313.15_dp)
    cold = evaluate_temperature_shift(provider, 273.15_dp)
    expected_hot_ln = activation_energy/gas_constant_reference * &
      (1.0_dp/313.15_dp-1.0_dp/reference_temperature_k)
    expected_cold_ln = activation_energy/gas_constant_reference * &
      (1.0_dp/273.15_dp-1.0_dp/reference_temperature_k)

    call assert_close(hot%log10_a_t, expected_hot_ln/ln_10, tolerance, &
      "Arrhenius sıcak analytical log10(a_T) sonucu hatalı.")
    call assert_close(cold%log10_a_t, expected_cold_ln/ln_10, tolerance, &
      "Arrhenius soğuk analytical log10(a_T) sonucu hatalı.")
    call assert_close(hot%a_t, exp(expected_hot_ln), tolerance, &
      "Arrhenius sıcak derived a_T sonucu hatalı.")
    call assert_close(cold%a_t, exp(expected_cold_ln), tolerance, &
      "Arrhenius soğuk derived a_T sonucu hatalı.")
    if (hot%a_t >= 1.0_dp .or. &
        hot%a_t*physical_frequency_hz >= physical_frequency_hz) then
      error stop "Arrhenius sıcaklık artışı reduced frequency'yi azaltmadı."
    end if
    if (cold%a_t <= 1.0_dp .or. &
        cold%a_t*physical_frequency_hz <= physical_frequency_hz) then
      error stop "Arrhenius sıcaklık düşüşü reduced frequency'yi artırmadı."
    end if
  end subroutine test_analytical_values_and_direction

  !> Ea>0 ve sonlu pozitif Kelvin domain/reference/query sözleşmesinin
  !! constructor ile runtime error-stop yollarını ayrı süreçlerde tetikler.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(arrhenius_temperature_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: evaluation
    real(dp) :: nan_value
    real(dp) :: infinity

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    select case (case_name)
    case ("zero_activation_energy")
      provider = create_arrhenius_temperature_shift_provider( &
        0.0_dp, 293.15_dp, 253.15_dp, 353.15_dp)
    case ("negative_activation_energy")
      provider = create_arrhenius_temperature_shift_provider( &
        -1.0_dp, 293.15_dp, 253.15_dp, 353.15_dp)
    case ("nan_activation_energy")
      provider = create_arrhenius_temperature_shift_provider( &
        nan_value, 293.15_dp, 253.15_dp, 353.15_dp)
    case ("infinite_activation_energy")
      provider = create_arrhenius_temperature_shift_provider( &
        infinity, 293.15_dp, 253.15_dp, 353.15_dp)
    case ("invalid_reference_temperature")
      provider = create_arrhenius_temperature_shift_provider( &
        activation_energy, 0.0_dp, 253.15_dp, 353.15_dp)
    case ("invalid_minimum_temperature")
      provider = create_arrhenius_temperature_shift_provider( &
        activation_energy, 293.15_dp, 0.0_dp, 353.15_dp)
    case ("reversed_domain")
      provider = create_arrhenius_temperature_shift_provider( &
        activation_energy, 293.15_dp, 353.15_dp, 253.15_dp)
    case ("reference_outside_domain")
      provider = create_arrhenius_temperature_shift_provider( &
        activation_energy, 293.15_dp, 303.15_dp, 353.15_dp)
    case ("below_domain")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 253.0_dp)
    case ("above_domain")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 353.2_dp)
    case ("zero_query_temperature")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 0.0_dp)
    case ("nan_query_temperature")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, nan_value)
    case default
      error stop "Bilinmeyen Arrhenius validation selector."
    end select
  end subroutine exercise_invalid_case

  pure function make_provider() result(provider)
    type(arrhenius_temperature_shift_provider_t) :: provider

    provider = create_arrhenius_temperature_shift_provider( &
      activation_energy, reference_temperature_k, 253.15_dp, 353.15_dp)
  end function make_provider

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

end program test_arrhenius_temperature_shift
