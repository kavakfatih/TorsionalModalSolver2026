program test_wlf_temperature_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan, ieee_positive_inf
  use tms_kinds, only : dp
  use tms_temperature_shift_types, only : temperature_shift_evaluation_t, &
    WLF_TEMPERATURE_SHIFT
  use tms_temperature_shift_provider, only : evaluate_temperature_shift
  use tms_wlf_temperature_shift, only : wlf_temperature_shift_provider_t, &
    create_wlf_temperature_shift_provider
  implicit none

  real(dp), parameter :: tolerance = 2.0e-13_dp
  real(dp), parameter :: c1 = 8.86_dp
  real(dp), parameter :: c2_k = 101.6_dp
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_reference_identity()
  call test_temperature_direction_and_analytical_values()

  print *, "V0.8 WLF temperature-shift provider doğrulandı."

contains

  !> T=T_ref için WLF payının sıfır olduğunu; s=0 ve a_T=1 identity
  !! koşullarını doğrular. Canonical TMS26 convention f_r=a_T*f'tir.
  subroutine test_reference_identity()
    type(wlf_temperature_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: evaluation

    provider = make_provider()
    evaluation = evaluate_temperature_shift( &
      provider, reference_temperature_k)

    if (evaluation%shift_model_kind /= WLF_TEMPERATURE_SHIFT) then
      error stop "WLF shift model kind trace'i hatalı."
    end if
    if (evaluation%has_temperature_bracket) then
      error stop "Analitik WLF modeli için yapay table bracket üretildi."
    end if
    call assert_close(evaluation%log10_a_t, 0.0_dp, 0.0_dp, &
      "WLF reference sıcaklığında log10(a_T) sıfır değil.")
    call assert_close(evaluation%a_t, 1.0_dp, 0.0_dp, &
      "WLF reference sıcaklığında a_T bir değil.")
  end subroutine test_reference_identity

  !> WLF modelini bağımsız s=-C1*dT/(C2+dT) hesabıyla doğrular.
  !! T>T_ref için a_T<1 ve f_r<f; T<T_ref için a_T>1 ve f_r>f
  !! olması, TMS26 canonical thermally activated direction testidir.
  subroutine test_temperature_direction_and_analytical_values()
    type(wlf_temperature_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: hot
    type(temperature_shift_evaluation_t) :: cold
    real(dp) :: expected_hot_s
    real(dp) :: expected_cold_s
    real(dp), parameter :: physical_frequency_hz = 100.0_dp

    provider = make_provider()
    hot = evaluate_temperature_shift(provider, 303.15_dp)
    cold = evaluate_temperature_shift(provider, 283.15_dp)
    expected_hot_s = -c1*10.0_dp/(c2_k+10.0_dp)
    expected_cold_s = -c1*(-10.0_dp)/(c2_k-10.0_dp)

    call assert_close(hot%log10_a_t, expected_hot_s, tolerance, &
      "WLF sıcak analytical log10(a_T) sonucu hatalı.")
    call assert_close(cold%log10_a_t, expected_cold_s, tolerance, &
      "WLF soğuk analytical log10(a_T) sonucu hatalı.")
    call assert_close(hot%a_t, 10.0_dp**expected_hot_s, tolerance, &
      "WLF sıcak derived a_T sonucu hatalı.")
    call assert_close(cold%a_t, 10.0_dp**expected_cold_s, tolerance, &
      "WLF soğuk derived a_T sonucu hatalı.")
    if (hot%a_t >= 1.0_dp .or. &
        hot%a_t*physical_frequency_hz >= physical_frequency_hz) then
      error stop "WLF sıcaklık artışı reduced frequency'yi azaltmadı."
    end if
    if (cold%a_t <= 1.0_dp .or. &
        cold%a_t*physical_frequency_hz <= physical_frequency_hz) then
      error stop "WLF sıcaklık düşüşü reduced frequency'yi artırmadı."
    end if
  end subroutine test_temperature_direction_and_analytical_values

  !> Pozitif C1/C2, sonlu Kelvin domain'i, reference containment, domain
  !! sorgusu ve WLF denominator D(T)>0 koşullarının hata yollarını tetikler.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(wlf_temperature_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: evaluation
    real(dp) :: nan_value
    real(dp) :: infinity

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    select case (case_name)
    case ("zero_c1")
      provider = create_wlf_temperature_shift_provider( &
        0.0_dp, c2_k, 293.15_dp, 273.15_dp, 333.15_dp)
    case ("negative_c1")
      provider = create_wlf_temperature_shift_provider( &
        -1.0_dp, c2_k, 293.15_dp, 273.15_dp, 333.15_dp)
    case ("nan_c1")
      provider = create_wlf_temperature_shift_provider( &
        nan_value, c2_k, 293.15_dp, 273.15_dp, 333.15_dp)
    case ("zero_c2")
      provider = create_wlf_temperature_shift_provider( &
        c1, 0.0_dp, 293.15_dp, 273.15_dp, 333.15_dp)
    case ("negative_c2")
      provider = create_wlf_temperature_shift_provider( &
        c1, -1.0_dp, 293.15_dp, 273.15_dp, 333.15_dp)
    case ("infinite_c2")
      provider = create_wlf_temperature_shift_provider( &
        c1, infinity, 293.15_dp, 273.15_dp, 333.15_dp)
    case ("invalid_reference_temperature")
      provider = create_wlf_temperature_shift_provider( &
        c1, c2_k, 0.0_dp, 273.15_dp, 333.15_dp)
    case ("invalid_minimum_temperature")
      provider = create_wlf_temperature_shift_provider( &
        c1, c2_k, 293.15_dp, 0.0_dp, 333.15_dp)
    case ("reversed_domain")
      provider = create_wlf_temperature_shift_provider( &
        c1, c2_k, 293.15_dp, 333.15_dp, 273.15_dp)
    case ("reference_outside_domain")
      provider = create_wlf_temperature_shift_provider( &
        c1, c2_k, 293.15_dp, 303.15_dp, 333.15_dp)
    case ("pole_in_domain")
      provider = create_wlf_temperature_shift_provider( &
        c1, 20.0_dp, 293.15_dp, 273.15_dp, 313.15_dp)
    case ("below_domain")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 273.0_dp)
    case ("above_domain")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 333.2_dp)
    case ("invalid_query_temperature")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, nan_value)
    case default
      error stop "Bilinmeyen WLF validation selector."
    end select
  end subroutine exercise_invalid_case

  pure function make_provider() result(provider)
    type(wlf_temperature_shift_provider_t) :: provider

    provider = create_wlf_temperature_shift_provider( &
      c1, c2_k, reference_temperature_k, 273.15_dp, 333.15_dp)
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

end program test_wlf_temperature_shift
