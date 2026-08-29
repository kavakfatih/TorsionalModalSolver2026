program test_tts_repeatability_bootstrap
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  use tms_tts_repeatability_types, only : tts_repeatability_campaign_t, &
    tts_repeatability_study_configuration_t, &
    tts_repeatability_study_result_t, INDEPENDENT_SPECIMEN_CAMPAIGN, &
    SAME_SPECIMEN_RERUN, TTS_REPEATABILITY_SUCCESS
  use tms_tts_repeatability_analysis, only : analyze_tts_repeatability
  use tms_tts_repeatability_test_support, only : &
    make_repeatability_campaign, make_linear_temperature_shifts
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures_k(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: slopes_per_k(5) = &
    [0.008_dp, 0.009_dp, 0.010_dp, 0.011_dp, 0.012_dp]
  type(tts_repeatability_campaign_t) :: campaigns(5)
  type(tts_repeatability_study_configuration_t) :: configuration
  type(tts_repeatability_study_result_t) :: first_study
  type(tts_repeatability_study_result_t) :: repeated_study
  type(tts_repeatability_study_result_t) :: independent_only_descriptive
  type(tts_repeatability_study_result_t) :: one_independent_study
  integer :: i
  character(len=32) :: identifier

  do i = 1, size(campaigns)
    write(identifier, '("MIXED-",I0)') i
    campaigns(i) = make_repeatability_campaign(temperatures_k, &
      make_linear_temperature_shifts(temperatures_k, &
        reference_temperature_k, slopes_per_k(i)), &
      reference_temperature_k, trim(identifier), &
      INDEPENDENT_SPECIMEN_CAMPAIGN)
  end do
  campaigns(4:5)%replicate_basis = SAME_SPECIMEN_RERUN

  configuration%study_identifier = "MIXED-INDEPENDENCE-STUDY"
  configuration%bootstrap%draw_count = 128
  configuration%bootstrap%confidence_level = 0.90_dp
  configuration%bootstrap%seed = 4242_int64
  first_study = analyze_tts_repeatability(campaigns, &
    reference_temperature_k, configuration)
  repeated_study = analyze_tts_repeatability(campaigns, &
    reference_temperature_k, configuration)

  call assert_true(first_study%status == TTS_REPEATABILITY_SUCCESS .and. &
    first_study%total_campaign_count == 5 .and. &
    first_study%descriptive_campaign_count == 5 .and. &
    first_study%independent_campaign_count == 3 .and. &
    first_study%same_specimen_rerun_count == 2, &
    "Mixed independent/rerun campaign counts hatalı.")
  call assert_true(first_study%independent_cluster_bootstrap_available .and. &
    size(first_study%bootstrap_population_campaign_identifiers) == 3 .and. &
    all(first_study%bootstrap_population_campaign_identifiers == &
      [character(len=256) :: "MIXED-1", "MIXED-2", "MIXED-3"]), &
    "Same-specimen rerun independent bootstrap population'a girdi.")
  call assert_true(first_study%pair_results(1)%delta_s_statistics &
    %sample_count == 5, &
    "Default policy same-specimen rerun descriptive evidence'ı gizledi.")
  call assert_true(first_study%pair_results(1)%mean_bootstrap_interval &
    %requested_bootstrap_draw_count == 128 .and. &
    first_study%pair_results(1)%mean_bootstrap_interval%seed == 4242_int64, &
    "Bootstrap configuration/provenance result'a taşınmadı.")
  call assert_close(first_study%pair_results(1)%mean_bootstrap_interval%lower, &
    repeated_study%pair_results(1)%mean_bootstrap_interval%lower, 0.0_dp, &
    "Aynı seed repeatability CI lower değerini değiştirdi.")
  call assert_close(first_study%pair_results(1)%mean_bootstrap_interval%upper, &
    repeated_study%pair_results(1)%mean_bootstrap_interval%upper, 0.0_dp, &
    "Aynı seed repeatability CI upper değerini değiştirdi.")
  call assert_true(first_study%pair_results(1)%mean_bootstrap_interval &
      %valid_bootstrap_draw_count == repeated_study%pair_results(1) &
        %mean_bootstrap_interval%valid_bootstrap_draw_count .and. &
    first_study%pair_results(1)%mean_bootstrap_interval &
      %unavailable_bootstrap_draw_count == repeated_study%pair_results(1) &
        %mean_bootstrap_interval%unavailable_bootstrap_draw_count, &
    "Deterministic bootstrap draw accounting değişti.")

  configuration%include_same_specimen_reruns_in_descriptive = .false.
  independent_only_descriptive = analyze_tts_repeatability(campaigns, &
    reference_temperature_k, configuration)
  call assert_true(independent_only_descriptive%status == &
      TTS_REPEATABILITY_SUCCESS .and. &
    independent_only_descriptive%descriptive_campaign_count == 3 .and. &
    independent_only_descriptive%pair_results(1)%delta_s_statistics &
      %sample_count == 3, &
    "Study policy rerun descriptive inclusion'ını kontrol etmiyor.")

  campaigns(2:5)%replicate_basis = SAME_SPECIMEN_RERUN
  one_independent_study = analyze_tts_repeatability(campaigns, &
    reference_temperature_k, configuration)
  call assert_true(one_independent_study%status == &
      TTS_REPEATABILITY_SUCCESS .and. &
    one_independent_study%independent_campaign_count == 1 .and. &
    .not. one_independent_study%independent_cluster_bootstrap_available, &
    "Tek independent campaign cluster uncertainty üretti.")

  print *, "V0.8.3 independence policy ve deterministic CI doğrulandı."
end program test_tts_repeatability_bootstrap
