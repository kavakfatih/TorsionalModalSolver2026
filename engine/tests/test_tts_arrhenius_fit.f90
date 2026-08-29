program test_tts_arrhenius_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tms_kinds, only : dp
  use tms_constants, only : universal_gas_constant_j_per_mol_k
  use tms_tts_shift_law_types, only : tts_pair_shift_observation_t, &
    tts_arrhenius_fit_result_t, SHIFT_LAW_FIT_SUCCESS, &
    SHIFT_LAW_FIT_INSUFFICIENT_DATA, SHIFT_LAW_FIT_INVALID_INPUT, &
    ARRHENIUS_FIT_INVALID_SLOPE
  use tms_tts_shift_law_fit, only : fit_tts_arrhenius_shift_law
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures_k(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: apparent_activation_energy_j_per_mol = 48000.0_dp
  type(tts_pair_shift_observation_t), allocatable :: observations(:)
  type(tts_arrhenius_fit_result_t) :: fit
  type(tts_arrhenius_fit_result_t) :: reparameterized_fit
  type(tts_arrhenius_fit_result_t) :: single_pair_fit
  real(dp) :: beta_k
  real(dp) :: old_shift
  real(dp) :: new_shift
  real(dp) :: shifted_old_reference

  beta_k = apparent_activation_energy_j_per_mol / &
    (universal_gas_constant_j_per_mol_k*log(10.0_dp))
  observations = make_arrhenius_observations(temperatures_k, beta_k)
  fit = fit_tts_arrhenius_shift_law(observations, reference_temperature_k)

  call assert_true(fit%status == SHIFT_LAW_FIT_SUCCESS .and. &
    fit%fit_available, "Exact Arrhenius adjacent-pair fit hazır değil.")
  call assert_close(fit%beta_k, beta_k, 2.0e-13_dp, &
    "Arrhenius analytical beta recovery hatalı.")
  call assert_close(fit%apparent_activation_energy_j_per_mol, &
    apparent_activation_energy_j_per_mol, 2.0e-13_dp, &
    "Arrhenius apparent activation energy recovery hatalı.")
  call assert_true(fit%residual_validation_available .and. &
    fit%validation_degree_of_freedom == 3, &
    "Arrhenius residual degree-of-freedom ayrımı hatalı.")
  call assert_true(fit%pair_rmse < 1.0e-13_dp .and. &
    fit%pair_max_abs_residual < 1.0e-13_dp, &
    "Exact Arrhenius pair residual sıfıra yakın değil.")

  ! Ea_app relative pair denkleminden geldiği için measured reference değişse
  ! de aynı kalmalıdır. Absolute shift yalnız yeni reference'a göre yeniden
  ! sıfırlanır: s_new(T)=s_old(T)-s_old(T_ref,new).
  reparameterized_fit = fit_tts_arrhenius_shift_law( &
    observations, 313.15_dp)
  call assert_true(reparameterized_fit%status == SHIFT_LAW_FIT_SUCCESS, &
    "İkinci measured reference ile Arrhenius fit başarısız.")
  call assert_close(reparameterized_fit%apparent_activation_energy_j_per_mol, &
    fit%apparent_activation_energy_j_per_mol, 2.0e-13_dp, &
    "Ea_app reference-temperature invariant değil.")
  old_shift = beta_k*(1.0_dp/253.15_dp-1.0_dp/reference_temperature_k)
  shifted_old_reference = beta_k*(1.0_dp/313.15_dp - &
    1.0_dp/reference_temperature_k)
  new_shift = beta_k*(1.0_dp/253.15_dp-1.0_dp/313.15_dp)
  call assert_close(new_shift, old_shift-shifted_old_reference, &
    2.0e-13_dp, "Arrhenius reference shift reconstruction hatalı.")

  ! Tek nonzero pair slope'u matematiksel olarak belirler; fakat residual
  ! degree-of-freedom ve predictive validation sağladığı iddia edilmez.
  single_pair_fit = fit_tts_arrhenius_shift_law( &
    observations(1:1), temperatures_k(1))
  call assert_true(single_pair_fit%fit_available .and. &
    .not. single_pair_fit%residual_validation_available .and. &
    .not. single_pair_fit%predictive_validation_available, &
    "Tek-pair Arrhenius fit/validation semantics ayrılmadı.")

  call test_invalid_inputs(beta_k)

  print *, "V0.8.2 analytical adjacent-pair Arrhenius fit doğrulandı."

contains

  pure function make_arrhenius_observations(temperatures, beta) &
      result(pair_observations)
    real(dp), intent(in) :: temperatures(:)
    real(dp), intent(in) :: beta
    type(tts_pair_shift_observation_t), allocatable :: pair_observations(:)
    integer :: i

    allocate(pair_observations(size(temperatures)-1))
    do i = 1, size(pair_observations)
      pair_observations(i)%isotherm_i_index = i
      pair_observations(i)%isotherm_j_index = i+1
      pair_observations(i)%temperature_i_k = temperatures(i)
      pair_observations(i)%temperature_j_k = temperatures(i+1)
      pair_observations(i)%delta_s_j_minus_i = beta * &
        (1.0_dp/temperatures(i+1)-1.0_dp/temperatures(i))
    end do
  end function make_arrhenius_observations

  subroutine test_invalid_inputs(beta)
    real(dp), intent(in) :: beta
    type(tts_pair_shift_observation_t), allocatable :: invalid(:)
    type(tts_arrhenius_fit_result_t) :: invalid_fit
    real(dp) :: nan_value

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    allocate(invalid(0))
    invalid_fit = fit_tts_arrhenius_shift_law(invalid, &
      reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INSUFFICIENT_DATA, &
      "Sıfır pair observation Arrhenius fit tarafından reddedilmedi.")

    deallocate(invalid)
    allocate(invalid(1))
    invalid(1)%temperature_i_k = 293.15_dp
    invalid(1)%temperature_j_k = 293.15_dp
    invalid(1)%delta_s_j_minus_i = 0.0_dp
    invalid_fit = fit_tts_arrhenius_shift_law(invalid, &
      reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Duplicate pair temperatures reddedilmedi.")

    invalid(1)%temperature_j_k = 313.15_dp
    invalid(1)%temperature_i_k = nan_value
    invalid_fit = fit_tts_arrhenius_shift_law(invalid, &
      reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Nonfinite Arrhenius temperature reddedilmedi.")

    invalid(1)%temperature_i_k = 293.15_dp
    invalid(1)%delta_s_j_minus_i = nan_value
    invalid_fit = fit_tts_arrhenius_shift_law(invalid, &
      reference_temperature_k)
    call assert_true(invalid_fit%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Nonfinite Arrhenius pair shift reddedilmedi.")

    invalid(1)%delta_s_j_minus_i = -beta * &
      (1.0_dp/313.15_dp-1.0_dp/293.15_dp)
    invalid_fit = fit_tts_arrhenius_shift_law(invalid, &
      reference_temperature_k)
    call assert_true(invalid_fit%status == ARRHENIUS_FIT_INVALID_SLOPE .and. &
      .not. invalid_fit%fit_available, &
      "Negatif/nonphysical Ea_app slope başarı sayıldı.")
  end subroutine test_invalid_inputs

end program test_tts_arrhenius_fit
