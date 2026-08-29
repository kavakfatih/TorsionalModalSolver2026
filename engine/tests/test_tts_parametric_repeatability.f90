program test_tts_parametric_repeatability
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  use tms_constants, only : universal_gas_constant_j_per_mol_k
  use tms_tts_repeatability_types, only : tts_repeatability_campaign_t, &
    tts_repeatability_study_configuration_t, &
    tts_repeatability_study_result_t, INDEPENDENT_SPECIMEN_CAMPAIGN, &
    TTS_REPEATABILITY_SUCCESS
  use tms_tts_repeatability_analysis, only : analyze_tts_repeatability
  use tms_tts_repeatability_test_support, only : &
    make_repeatability_campaign, make_arrhenius_shifts, make_wlf_shifts, &
    make_linear_temperature_shifts, replace_campaign_shift_evidence
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures_k(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: arrhenius_energies(3) = &
    [45000.0_dp, 48000.0_dp, 51000.0_dp]
  real(dp), parameter :: wlf_c1_values(3) = [1.9_dp, 2.0_dp, 2.1_dp]
  real(dp), parameter :: wlf_c2_k = 120.0_dp
  real(dp), parameter :: steep_absolute_shifts(5) = &
    [0.0_dp, -100.0_dp, -101.0_dp, -101.01_dp, -101.0101_dp]
  type(tts_repeatability_campaign_t) :: arrhenius_campaigns(3)
  type(tts_repeatability_campaign_t) :: wlf_campaigns(5)
  type(tts_repeatability_study_configuration_t) :: configuration
  type(tts_repeatability_study_result_t) :: study
  real(dp) :: expected_beta_mean
  integer :: i
  character(len=32) :: identifier

  do i = 1, size(arrhenius_campaigns)
    write(identifier, '("ARR-",I0)') i
    arrhenius_campaigns(i) = make_repeatability_campaign(temperatures_k, &
      make_arrhenius_shifts(temperatures_k, reference_temperature_k, &
        arrhenius_energies(i)), reference_temperature_k, trim(identifier), &
      INDEPENDENT_SPECIMEN_CAMPAIGN)
  end do
  configuration%study_identifier = "ARRHENIUS-COHORT"
  configuration%bootstrap%draw_count = 200
  configuration%bootstrap%seed = 811_int64
  study = analyze_tts_repeatability(arrhenius_campaigns, &
    reference_temperature_k, configuration)
  call assert_true(study%status == TTS_REPEATABILITY_SUCCESS .and. &
    study%arrhenius%fit_available_count == 3 .and. &
    study%arrhenius%fit_unavailable_count == 0, &
    "Valid Arrhenius campaign cohort fits korunmadı.")
  call assert_close(study%arrhenius &
    %apparent_activation_energy_statistics%mean, 48000.0_dp, 3.0e-7_dp, &
    "Ea_app cohort mean hatalı.")
  call assert_close(study%arrhenius &
    %apparent_activation_energy_statistics%sample_standard_deviation, &
    3000.0_dp, 3.0e-7_dp, "Ea_app cohort sample SD hatalı.")
  call assert_close(study%arrhenius &
    %apparent_activation_energy_statistics%standard_error, &
    3000.0_dp/sqrt(3.0_dp), 3.0e-7_dp, "Ea_app cohort SE hatalı.")
  call assert_close(study%arrhenius &
    %apparent_activation_energy_statistics%median, 48000.0_dp, 3.0e-7_dp, &
    "Ea_app cohort median hatalı.")
  call assert_close(study%arrhenius &
    %apparent_activation_energy_statistics%median_absolute_deviation, &
    3000.0_dp, 3.0e-7_dp, "Ea_app cohort MAD hatalı.")
  expected_beta_mean = 48000.0_dp / &
    (universal_gas_constant_j_per_mol_k*log(10.0_dp))
  call assert_close(study%arrhenius%beta_k_statistics%mean, &
    expected_beta_mean, 3.0e-7_dp, "Arrhenius beta cohort mean hatalı.")

  do i = 1, 3
    write(identifier, '("WLF-",I0)') i
    wlf_campaigns(i) = make_repeatability_campaign(temperatures_k, &
      make_wlf_shifts(temperatures_k, reference_temperature_k, &
        wlf_c1_values(i), wlf_c2_k), reference_temperature_k, &
      trim(identifier), INDEPENDENT_SPECIMEN_CAMPAIGN)
  end do
  wlf_campaigns(4) = make_repeatability_campaign(temperatures_k, &
    make_linear_temperature_shifts(temperatures_k, reference_temperature_k, &
      0.01_dp), reference_temperature_k, "WLF-POOR", &
    INDEPENDENT_SPECIMEN_CAMPAIGN)
  call replace_campaign_shift_evidence(wlf_campaigns(4), temperatures_k, &
    make_linear_temperature_shifts(temperatures_k, &
      reference_temperature_k, 0.01_dp))
  wlf_campaigns(5) = make_repeatability_campaign(temperatures_k, &
    make_wlf_shifts(temperatures_k, reference_temperature_k, 2.0_dp, &
      wlf_c2_k), reference_temperature_k, "WLF-INVALID", &
    INDEPENDENT_SPECIMEN_CAMPAIGN)
  call replace_campaign_shift_evidence(wlf_campaigns(5), temperatures_k, &
    steep_absolute_shifts)

  configuration%study_identifier = "WLF-PARTIAL-COHORT"
  configuration%bootstrap%draw_count = 400
  configuration%bootstrap%seed = 991_int64
  study = analyze_tts_repeatability(wlf_campaigns, &
    reference_temperature_k, configuration)
  call assert_true(study%status == TTS_REPEATABILITY_SUCCESS .and. &
    study%wlf%parameter_identifiable_count == 3 .and. &
    study%wlf%poorly_identified_count == 1 .and. &
    study%wlf%invalid_fit_count == 1, &
    "WLF identifiable/poor/invalid cohort counts hatalı.")
  call assert_true(study%wlf%c1_statistics%sample_count == 3 .and. &
    study%wlf%c2_k_statistics%sample_count == 3 .and. &
    study%wlf%p_c1_over_c2_per_k_statistics%sample_count == 3 .and. &
    study%wlf%q_inverse_c2_per_k_statistics%sample_count == 3, &
    "Poor/invalid WLF values parameter statistics'e placeholder olarak girdi.")
  call assert_close(study%wlf%c1_statistics%mean, 2.0_dp, 5.0e-6_dp, &
    "Identifiable WLF C1 cohort mean hatalı.")
  call assert_close(study%wlf%c2_k_statistics%mean, wlf_c2_k, 5.0e-6_dp, &
    "Identifiable WLF C2 cohort mean hatalı.")
  call assert_close(study%wlf%p_c1_over_c2_per_k_statistics%mean, &
    2.0_dp/wlf_c2_k, 5.0e-6_dp, "Identifiable WLF p mean hatalı.")
  call assert_close(study%wlf%q_inverse_c2_per_k_statistics%mean, &
    1.0_dp/wlf_c2_k, 5.0e-6_dp, "Identifiable WLF q mean hatalı.")
  call assert_true(study%wlf%c1_mean_bootstrap_interval &
    %requested_bootstrap_draw_count == 400 .and. &
    study%wlf%c1_mean_bootstrap_interval%valid_bootstrap_draw_count > 0 .and. &
    study%wlf%c1_mean_bootstrap_interval &
      %unavailable_bootstrap_draw_count > 0 .and. &
    study%wlf%c1_mean_bootstrap_interval%valid_bootstrap_draw_count + &
      study%wlf%c1_mean_bootstrap_interval &
        %unavailable_bootstrap_draw_count == 400, &
    "Partial WLF bootstrap valid/unavailable draw accounting hatalı.")
  call assert_true(study%wlf%c1_mean_bootstrap_interval &
      %valid_bootstrap_draw_count == &
        study%wlf%c2_k_mean_bootstrap_interval%valid_bootstrap_draw_count .and. &
    study%wlf%c1_mean_bootstrap_interval%valid_bootstrap_draw_count == &
      study%wlf%p_mean_bootstrap_interval%valid_bootstrap_draw_count .and. &
    study%wlf%c1_mean_bootstrap_interval%valid_bootstrap_draw_count == &
      study%wlf%q_mean_bootstrap_interval%valid_bootstrap_draw_count, &
    "WLF quantities aynı campaign draw/usability planını paylaşmıyor.")

  print *, "V0.8.3 Arrhenius/WLF cohort variability doğrulandı."
end program test_tts_parametric_repeatability
