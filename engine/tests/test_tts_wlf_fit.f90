program test_tts_wlf_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tms_kinds, only : dp
  use tms_tts_shift_law_types, only : tts_pair_shift_observation_t, &
    tts_wlf_fit_result_t, SHIFT_LAW_FIT_SUCCESS, &
    SHIFT_LAW_FIT_INSUFFICIENT_DATA, SHIFT_LAW_FIT_INVALID_INPUT, &
    WLF_FIT_NO_INTERIOR_BRACKET, WLF_FIT_POORLY_IDENTIFIED
  use tms_tts_shift_law_fit, only : fit_tts_wlf_shift_law
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures_k(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: known_c1 = 8.86_dp
  real(dp), parameter :: known_c2_k = 101.6_dp
  type(tts_pair_shift_observation_t), allocatable :: observations(:)
  type(tts_pair_shift_observation_t), allocatable :: linear_observations(:)
  type(tts_wlf_fit_result_t) :: fit
  type(tts_wlf_fit_result_t) :: reparameterized_fit
  type(tts_wlf_fit_result_t) :: weak_fit
  real(dp) :: transformed_c1
  real(dp) :: transformed_c2_k

  observations = make_wlf_observations(temperatures_k, &
    reference_temperature_k, known_c1, known_c2_k)
  fit = fit_tts_wlf_shift_law(observations, reference_temperature_k)
  call assert_true(fit%status == SHIFT_LAW_FIT_SUCCESS .and. &
    fit%fit_available .and. fit%parameter_identifiable, &
    "Exact WLF profiled fit identifiable sonuç üretmedi.")
  call assert_close(fit%c1, known_c1, 2.0e-7_dp, &
    "WLF C1 recovery hatalı.")
  call assert_close(fit%c2_k, known_c2_k, 2.0e-7_dp, &
    "WLF C2 recovery hatalı.")
  call assert_close(fit%p_c1_over_c2_per_k, known_c1/known_c2_k, &
    2.0e-7_dp, "WLF p=C1/C2 diagnostic hatalı.")
  call assert_close(fit%q_inverse_c2_per_k, 1.0_dp/known_c2_k, &
    2.0e-7_dp, "WLF q=1/C2 diagnostic hatalı.")
  call assert_true(fit%pair_rmse < 3.0e-8_dp .and. &
    fit%pair_max_abs_residual < 3.0e-8_dp, &
    "Exact WLF pair residual yeterince küçük değil.")

  ! Aynı physical WLF law yeni measured reference'a taşındığında
  ! C2'=C2+d, C1'=C1*C2/(C2+d) ve C1'*C2'=C1*C2 olmalıdır.
  transformed_c2_k = known_c2_k+20.0_dp
  transformed_c1 = known_c1*known_c2_k/transformed_c2_k
  reparameterized_fit = fit_tts_wlf_shift_law(observations, 313.15_dp)
  call assert_true(reparameterized_fit%status == SHIFT_LAW_FIT_SUCCESS .and. &
    reparameterized_fit%parameter_identifiable, &
    "Yeni measured reference ile WLF fit identifiable değil.")
  call assert_close(reparameterized_fit%c2_k, transformed_c2_k, &
    2.0e-7_dp, "WLF reference dönüşümünde C2' hatalı.")
  call assert_close(reparameterized_fit%c1, transformed_c1, &
    2.0e-7_dp, "WLF reference dönüşümünde C1' hatalı.")
  call assert_close(reparameterized_fit%c1*reparameterized_fit%c2_k, &
    known_c1*known_c2_k, 4.0e-7_dp, &
    "WLF C1*C2 reference invariant'ı bozuldu.")
  call assert_true(maxval(abs(reparameterized_fit%predicted_pair_shifts - &
    fit%predicted_pair_shifts)) < 1.0e-8_dp, &
    "WLF relative physical predictions reference ile değişti.")

  ! Exact linear-in-temperature shift, WLF'nin C2->infinity limitidir.
  ! Residual küçük olsa bile ayrı finite C1/C2 tanımlanamaz; yalnız p belirlenir.
  linear_observations = make_linear_observations(temperatures_k, 0.05_dp)
  weak_fit = fit_tts_wlf_shift_law(linear_observations, &
    reference_temperature_k)
  call assert_true(weak_fit%status == WLF_FIT_POORLY_IDENTIFIED .and. &
    weak_fit%fit_available .and. .not. weak_fit%parameter_identifiable, &
    "Large-C2 WLF degeneracy explicit tanınmadı.")
  call assert_close(weak_fit%p_c1_over_c2_per_k, 0.05_dp, 2.0e-12_dp, &
    "Poorly identified WLF linear-limit slope p hatalı.")
  call assert_true(weak_fit%c1 == 0.0_dp .and. weak_fit%c2_k == 0.0_dp .and. &
    weak_fit%q_inverse_c2_per_k == 0.0_dp, &
    "Poorly identified WLF absurd finite C1/C2 olarak sunuldu.")

  linear_observations(1)%delta_s_j_minus_i = -100.0_dp
  linear_observations(2)%delta_s_j_minus_i = -1.0_dp
  linear_observations(3)%delta_s_j_minus_i = -0.01_dp
  linear_observations(4)%delta_s_j_minus_i = -0.0001_dp
  weak_fit = fit_tts_wlf_shift_law(linear_observations, &
    reference_temperature_k)
  call assert_true(weak_fit%status == WLF_FIT_NO_INTERIOR_BRACKET .and. &
    .not. weak_fit%fit_available, &
    "Pole sınırına kaçan WLF objective interior bracket başarı sayıldı.")

  call test_invalid_inputs()

  print *, "V0.8.2 profiled WLF fit ve identifiability doğrulandı."

contains

  pure function make_wlf_observations( &
      temperatures, reference_temperature, c1, c2_k) result(pair_observations)
    real(dp), intent(in) :: temperatures(:)
    real(dp), intent(in) :: reference_temperature
    real(dp), intent(in) :: c1
    real(dp), intent(in) :: c2_k
    type(tts_pair_shift_observation_t), allocatable :: pair_observations(:)
    real(dp) :: shift_i
    real(dp) :: shift_j
    integer :: i

    allocate(pair_observations(size(temperatures)-1))
    do i = 1, size(pair_observations)
      shift_i = -c1*(temperatures(i)-reference_temperature) / &
        (c2_k+temperatures(i)-reference_temperature)
      shift_j = -c1*(temperatures(i+1)-reference_temperature) / &
        (c2_k+temperatures(i+1)-reference_temperature)
      pair_observations(i)%isotherm_i_index = i
      pair_observations(i)%isotherm_j_index = i+1
      pair_observations(i)%temperature_i_k = temperatures(i)
      pair_observations(i)%temperature_j_k = temperatures(i+1)
      pair_observations(i)%delta_s_j_minus_i = shift_j-shift_i
    end do
  end function make_wlf_observations

  pure function make_linear_observations(temperatures, p) &
      result(pair_observations)
    real(dp), intent(in) :: temperatures(:)
    real(dp), intent(in) :: p
    type(tts_pair_shift_observation_t), allocatable :: pair_observations(:)
    integer :: i

    allocate(pair_observations(size(temperatures)-1))
    do i = 1, size(pair_observations)
      pair_observations(i)%temperature_i_k = temperatures(i)
      pair_observations(i)%temperature_j_k = temperatures(i+1)
      pair_observations(i)%delta_s_j_minus_i = &
        -p*(temperatures(i+1)-temperatures(i))
    end do
  end function make_linear_observations

  subroutine test_invalid_inputs()
    type(tts_pair_shift_observation_t), allocatable :: invalid(:)
    type(tts_wlf_fit_result_t) :: invalid_fit
    real(dp) :: nan_value

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    allocate(invalid(1))
    invalid(1)%temperature_i_k = 273.15_dp
    invalid(1)%temperature_j_k = 293.15_dp
    invalid(1)%delta_s_j_minus_i = -1.0_dp
    invalid_fit = fit_tts_wlf_shift_law(invalid, &
      reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INSUFFICIENT_DATA, &
      "Tek pair WLF fit için yetersiz sayılmadı.")

    deallocate(invalid)
    allocate(invalid(2))
    invalid(1)%temperature_i_k = 273.15_dp
    invalid(1)%temperature_j_k = 273.15_dp
    invalid(1)%delta_s_j_minus_i = 0.0_dp
    invalid(2)%temperature_i_k = 293.15_dp
    invalid(2)%temperature_j_k = 313.15_dp
    invalid(2)%delta_s_j_minus_i = -1.0_dp
    invalid_fit = fit_tts_wlf_shift_law(invalid, reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Duplicate WLF pair temperatures reddedilmedi.")

    invalid(1)%temperature_j_k = 293.15_dp
    invalid(1)%temperature_i_k = nan_value
    invalid_fit = fit_tts_wlf_shift_law(invalid, reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Nonfinite WLF temperature reddedilmedi.")

    invalid(1)%temperature_i_k = 273.15_dp
    invalid(1)%delta_s_j_minus_i = nan_value
    invalid_fit = fit_tts_wlf_shift_law(invalid, reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Nonfinite WLF pair shift reddedilmedi.")
  end subroutine test_invalid_inputs

end program test_tts_wlf_fit
