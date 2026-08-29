program test_tts_shift_law_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tms_kinds, only : dp
  use tms_constants, only : universal_gas_constant_j_per_mol_k
  use tms_tts_types, only : tts_identification_result_t
  use tms_tts_shift_law_types, only : &
    tts_shift_law_identification_result_t, SHIFT_LAW_FIT_SUCCESS, &
    SHIFT_LAW_FIT_INVALID_INPUT, SHIFT_LAW_FIT_INSUFFICIENT_DATA
  use tms_tts_shift_law_validation, only : fit_tts_shift_laws
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures_k(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: apparent_activation_energy_j_per_mol = 12000.0_dp
  real(dp), parameter :: wlf_c1 = 2.0_dp
  real(dp), parameter :: wlf_c2_k = 120.0_dp
  type(tts_identification_result_t) :: arrhenius_identification
  type(tts_identification_result_t) :: wlf_identification
  type(tts_shift_law_identification_result_t) :: laws
  real(dp), allocatable :: empirical_snapshot(:)
  real(dp), allocatable :: shifts(:)
  real(dp) :: beta_k
  integer :: i

  beta_k = apparent_activation_energy_j_per_mol / &
    (universal_gas_constant_j_per_mol_k*log(10.0_dp))
  allocate(shifts(size(temperatures_k)))
  do i = 1, size(temperatures_k)
    shifts(i) = beta_k*(1.0_dp/temperatures_k(i) - &
      1.0_dp/reference_temperature_k)
  end do
  arrhenius_identification = make_stub_identification(temperatures_k, shifts)

  ! Cumulative empirical shift tablosuna kasıtlı farklı değerler yazılır.
  ! Fit doğru kalıyorsa primary observations'ın yalnız adjacent pair delta_s
  ! olduğu ve empirical authoritative tablonun mutate edilmediği kanıtlanır.
  do i = 1, size(arrhenius_identification%empirical_shifts)
    arrhenius_identification%empirical_shifts(i)%log10_a_t = &
      100.0_dp+real(i, dp)
  end do
  empirical_snapshot = &
    arrhenius_identification%empirical_shifts%log10_a_t
  laws = fit_tts_shift_laws(arrhenius_identification)
  call assert_true(laws%status == SHIFT_LAW_FIT_SUCCESS .and. &
    laws%pair_observations_available .and. laws%pair_observation_count == 4, &
    "V0.8.1 adjacent pair observations extraction başarısız.")
  call assert_close(laws%arrhenius%apparent_activation_energy_j_per_mol, &
    apparent_activation_energy_j_per_mol, 2.0e-13_dp, &
    "Shift-law layer cumulative absolute shift'i fit observation yaptı.")
  call assert_true(all(arrhenius_identification%empirical_shifts%log10_a_t == &
    empirical_snapshot), "Parametric fit empirical shift table'ı mutate etti.")
  call assert_true(laws%arrhenius%predictive_validation_available .and. &
    allocated(laws%arrhenius%loto_diagnostics) .and. &
    size(laws%arrhenius%loto_diagnostics) == 5 .and. &
    all(laws%arrhenius%loto_diagnostics%available), &
    "Beş-temperature Arrhenius LOTO diagnostics tamamlanmadı.")
  call assert_true(laws%arrhenius%loto_rmse < 1.0e-12_dp, &
    "Exact Arrhenius LOTO prediction residual yüksek.")

  do i = 1, size(temperatures_k)
    shifts(i) = -wlf_c1*(temperatures_k(i)-reference_temperature_k) / &
      (wlf_c2_k+temperatures_k(i)-reference_temperature_k)
  end do
  wlf_identification = make_stub_identification(temperatures_k, shifts)
  laws = fit_tts_shift_laws(wlf_identification)
  call assert_true(laws%wlf%fit_available .and. &
    laws%wlf%parameter_identifiable, &
    "Exact WLF identification result profiled fit'e taşınmadı.")
  call assert_true(laws%wlf%predictive_validation_available .and. &
    allocated(laws%wlf%loto_diagnostics) .and. &
    size(laws%wlf%loto_diagnostics) == 5 .and. &
    all(laws%wlf%loto_diagnostics%available), &
    "Beş-temperature WLF LOTO diagnostics tamamlanmadı.")
  call assert_true(laws%wlf%loto_rmse < 5.0e-7_dp, &
    "Exact WLF LOTO prediction residual yüksek.")

  call test_invalid_identification(arrhenius_identification)

  print *, "V0.8.2 pair extraction, immutability ve LOTO doğrulandı."

contains

  function make_stub_identification(temperatures, absolute_shifts) &
      result(identification)
    real(dp), intent(in) :: temperatures(:)
    real(dp), intent(in) :: absolute_shifts(:)
    type(tts_identification_result_t) :: identification
    integer :: pair_index

    if (size(temperatures) /= size(absolute_shifts)) then
      error stop "Stub identification temperature/shift boyutları eşit olmalı."
    end if
    identification%reference_isotherm_index = 3
    identification%reference_temperature_k = temperatures(3)
    allocate(identification%source_family%isotherms(size(temperatures)))
    allocate(identification%empirical_shifts(size(temperatures)))
    do pair_index = 1, size(temperatures)
      identification%source_family%isotherms(pair_index)%temperature_k = &
        temperatures(pair_index)
      identification%empirical_shifts(pair_index)%temperature_k = &
        temperatures(pair_index)
      identification%empirical_shifts(pair_index)%log10_a_t = &
        absolute_shifts(pair_index)
    end do
    allocate(identification%pair_shift_results(size(temperatures)-1))
    do pair_index = 1, size(identification%pair_shift_results)
      identification%pair_shift_results(pair_index) &
        %reference_isotherm_index = pair_index
      identification%pair_shift_results(pair_index)%moving_isotherm_index = &
        pair_index+1
      identification%pair_shift_results(pair_index)%shift_available = .true.
      identification%pair_shift_results(pair_index)%delta_s = &
        absolute_shifts(pair_index+1)-absolute_shifts(pair_index)
    end do
  end function make_stub_identification

  subroutine test_invalid_identification(valid_identification)
    type(tts_identification_result_t), intent(in) :: valid_identification
    type(tts_identification_result_t) :: invalid
    type(tts_shift_law_identification_result_t) :: invalid_laws

    invalid%reference_isotherm_index = 1
    invalid%reference_temperature_k = reference_temperature_k
    invalid_laws = fit_tts_shift_laws(invalid)
    call assert_true(invalid_laws%status == &
      SHIFT_LAW_FIT_INSUFFICIENT_DATA, &
      "Pair results içermeyen identification reddedilmedi.")

    invalid = valid_identification
    invalid%pair_shift_results(1)%shift_available = .false.
    invalid_laws = fit_tts_shift_laws(invalid)
    call assert_true(invalid_laws%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Unavailable pair shift parametric fit'e girdi.")

    invalid = valid_identification
    invalid%pair_shift_results(1)%delta_s = &
      ieee_value(0.0_dp, ieee_quiet_nan)
    invalid_laws = fit_tts_shift_laws(invalid)
    call assert_true(invalid_laws%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Nonfinite pair shift parametric fit tarafından reddedilmedi.")

    invalid = valid_identification
    invalid%reference_temperature_k = 300.0_dp
    invalid_laws = fit_tts_shift_laws(invalid)
    call assert_true(invalid_laws%status == SHIFT_LAW_FIT_INVALID_INPUT, &
      "Reference temperature provenance mismatch reddedilmedi.")
  end subroutine test_invalid_identification

end program test_tts_shift_law_validation
