program test_tts_repeatability
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  use tms_tts_repeatability_types, only : tts_repeatability_campaign_t, &
    tts_repeatability_study_configuration_t, &
    tts_repeatability_study_result_t, INDEPENDENT_SPECIMEN_CAMPAIGN, &
    SAME_SPECIMEN_RERUN, TTS_REPEATABILITY_SUCCESS, &
    TTS_REPEATABILITY_INVALID_INPUT, &
    TTS_REPEATABILITY_INCOMPATIBLE_STATE, &
    TTS_REPEATABILITY_TEMPERATURE_SET_MISMATCH, &
    TTS_REPEATABILITY_REFERENCE_NOT_FOUND, &
    TTS_REPEATABILITY_NO_INDEPENDENT_CAMPAIGNS, &
    TTS_REPEATABILITY_NONFINITE_DATA
  use tms_tts_repeatability_analysis, only : analyze_tts_repeatability
  use tms_tts_repeatability_test_support, only : &
    make_repeatability_campaign, make_linear_temperature_shifts, &
    reverse_pair_order_and_orientation
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: temperatures_a(5) = &
    [293.15_dp, 273.15_dp, 313.15_dp, 333.15_dp, 253.15_dp]
  real(dp), parameter :: temperatures_b(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: temperatures_c(5) = &
    [333.15_dp, 313.15_dp, 293.15_dp, 273.15_dp, 253.15_dp]
  type(tts_repeatability_campaign_t) :: campaigns(3)
  type(tts_repeatability_campaign_t) :: altered(3)
  type(tts_repeatability_campaign_t) :: reference_campaigns(3)
  type(tts_repeatability_study_configuration_t) :: configuration
  type(tts_repeatability_study_result_t) :: study
  real(dp), allocatable :: empirical_snapshot(:)
  real(dp), allocatable :: pair_snapshot(:)

  campaigns(1) = make_repeatability_campaign(temperatures_a, &
    make_linear_temperature_shifts(temperatures_a, reference_temperature_k, &
      0.009_dp), reference_temperature_k, "CAMPAIGN-A", &
    INDEPENDENT_SPECIMEN_CAMPAIGN)
  campaigns(2) = make_repeatability_campaign(temperatures_b, &
    make_linear_temperature_shifts(temperatures_b, reference_temperature_k, &
      0.010_dp), reference_temperature_k, "CAMPAIGN-B", &
    INDEPENDENT_SPECIMEN_CAMPAIGN)
  campaigns(3) = make_repeatability_campaign(temperatures_c, &
    make_linear_temperature_shifts(temperatures_c, reference_temperature_k, &
      0.011_dp), reference_temperature_k, "CAMPAIGN-C", &
    INDEPENDENT_SPECIMEN_CAMPAIGN)
  call reverse_pair_order_and_orientation(campaigns(1))
  call reverse_pair_order_and_orientation(campaigns(3))
  empirical_snapshot = campaigns(1)%identification%empirical_shifts%log10_a_t
  pair_snapshot = campaigns(1)%identification%pair_shift_results%delta_s

  configuration%study_identifier = "ORDER-AND-SCATTER-STUDY"
  configuration%bootstrap%draw_count = 200
  configuration%bootstrap%seed = 1234_int64
  study = analyze_tts_repeatability(campaigns, reference_temperature_k, &
    configuration)
  call assert_true(study%status == TTS_REPEATABILITY_SUCCESS .and. &
    study%descriptive_statistics_available .and. &
    study%independent_cluster_bootstrap_available, &
    "Üç independent complete campaign repeatability study başarısız.")
  call assert_true(study%total_campaign_count == 3 .and. &
    study%descriptive_campaign_count == 3 .and. &
    study%independent_campaign_count == 3 .and. &
    study%same_specimen_rerun_count == 0, &
    "Campaign independence counts hatalı.")
  call assert_true(study%intralaboratory_context_explicit, &
    "Explicit synthetic intralaboratory context tanınmadı.")
  call assert_true(size(study%pair_results) == 4 .and. &
    size(study%absolute_shift_results) == 5, &
    "Canonical temperature/pair result boyutları hatalı.")
  call assert_close(study%pair_results(1)%lower_temperature_k, 253.15_dp, &
    1.0e-14_dp, "Canonical pair lower temperature hatalı.")
  call assert_close(study%pair_results(1)%upper_temperature_k, 273.15_dp, &
    1.0e-14_dp, "Canonical pair upper temperature hatalı.")
  call assert_close(study%pair_results(1)%delta_s_statistics%mean, &
    -0.2_dp, 2.0e-8_dp, "Canonical delta_s cohort mean hatalı.")
  call assert_close( &
    study%pair_results(1)%delta_s_statistics%sample_standard_deviation, &
    0.02_dp, 2.0e-8_dp, "Canonical delta_s sample SD hatalı.")
  call assert_close(study%pair_results(1)%delta_s_statistics%standard_error, &
    0.02_dp/sqrt(3.0_dp), 2.0e-8_dp, "Canonical delta_s SE hatalı.")
  call assert_close( &
    study%pair_results(1)%delta_s_statistics%median_absolute_deviation, &
    0.02_dp, 2.0e-8_dp, "Canonical delta_s MAD hatalı.")
  call assert_close( &
    study%absolute_shift_results(1)%log10_a_t_statistics%mean, &
    0.4_dp, 2.0e-8_dp, "Common-reference absolute shift mean hatalı.")
  call assert_close(study%absolute_shift_results(1) &
    %log10_a_t_statistics%sample_standard_deviation, 0.04_dp, &
    2.0e-8_dp, "Common-reference absolute shift SD hatalı.")

  call assert_true(study%absolute_shift_results(3)%is_reference_anchor .and. &
    .not. study%absolute_shift_results(3)%uncertainty_informative, &
    "Common reference structural anchor semantics işaretlenmedi.")
  call assert_close(study%absolute_shift_results(3) &
    %log10_a_t_statistics%mean, 0.0_dp, 0.0_dp, &
    "Structural reference anchor mean sıfır değil.")
  call assert_close(study%absolute_shift_results(3) &
    %log10_a_t_statistics%sample_standard_deviation, 0.0_dp, 0.0_dp, &
    "Structural reference anchor descriptive spread sıfır değil.")
  call assert_true(study%absolute_shift_results(3) &
    %mean_bootstrap_interval%available, &
    "Reference anchor bootstrap numerical interval unavailable.")
  call assert_close(study%absolute_shift_results(3) &
    %mean_bootstrap_interval%lower, 0.0_dp, 0.0_dp, &
    "Reference anchor bootstrap lower sıfır değil.")
  call assert_close(study%absolute_shift_results(3) &
    %mean_bootstrap_interval%upper, 0.0_dp, 0.0_dp, &
    "Reference anchor bootstrap upper sıfır değil.")

  call assert_true(maxval(abs(campaigns(1)%identification%empirical_shifts &
    %log10_a_t - empirical_snapshot)) <= 0.0_dp, &
    "Repeatability analysis authoritative empirical shifts'i mutate etti.")
  call assert_true(maxval(abs(campaigns(1)%identification%pair_shift_results &
    %delta_s - pair_snapshot)) <= 0.0_dp, &
    "Repeatability analysis authoritative pair shifts'i mutate etti.")
  call assert_true(allocated(study%campaign_provenance(1) &
    %specimen_identifiers) .and. &
    trim(study%campaign_provenance(1)%specimen_identifiers(1)) /= &
      trim(study%campaign_provenance(2)%specimen_identifiers(1)), &
    "Specimen provenance korunmadı veya farklı specimen false gate oldu.")

  call test_reference_reparameterization(reference_campaigns, configuration)
  call test_validation_failures(campaigns, altered, configuration)
  call test_single_and_same_specimen_semantics(campaigns, configuration)

  print *, "V0.8.3 canonical repeatability ve anchor semantics doğrulandı."

contains

  subroutine test_reference_reparameterization( &
      local_campaigns, local_configuration)
    type(tts_repeatability_campaign_t), intent(out) :: local_campaigns(3)
    type(tts_repeatability_study_configuration_t), intent(in) :: &
      local_configuration
    type(tts_repeatability_study_result_t) :: local_study
    integer :: i

    local_campaigns(1) = make_repeatability_campaign(temperatures_a, &
      make_linear_temperature_shifts(temperatures_a, &
        reference_temperature_k, 0.01_dp), 293.15_dp, "REF-293", &
      INDEPENDENT_SPECIMEN_CAMPAIGN)
    local_campaigns(2) = make_repeatability_campaign(temperatures_b, &
      make_linear_temperature_shifts(temperatures_b, &
        reference_temperature_k, 0.01_dp), 313.15_dp, "REF-313", &
      INDEPENDENT_SPECIMEN_CAMPAIGN)
    local_campaigns(3) = make_repeatability_campaign(temperatures_c, &
      make_linear_temperature_shifts(temperatures_c, &
        reference_temperature_k, 0.01_dp), 273.15_dp, "REF-273", &
      INDEPENDENT_SPECIMEN_CAMPAIGN)
    call reverse_pair_order_and_orientation(local_campaigns(3))
    local_study = analyze_tts_repeatability(local_campaigns, &
      reference_temperature_k, local_configuration)
    call assert_true(local_study%status == TTS_REPEATABILITY_SUCCESS, &
      "Farklı original reference campaign normalizasyonu başarısız.")
    do i = 1, size(local_study%absolute_shift_results)
      call assert_true(local_study%absolute_shift_results(i) &
        %log10_a_t_statistics%spread_statistics_available, &
        "Reference-normalized cohort spread availability eksik.")
      call assert_true(local_study%absolute_shift_results(i) &
        %log10_a_t_statistics%sample_standard_deviation < 2.0e-8_dp, &
        "Aynı physical shifts farklı original reference ile eşleşmedi.")
    end do
  end subroutine test_reference_reparameterization

  subroutine test_validation_failures( &
      valid_campaigns, local_altered, local_configuration)
    type(tts_repeatability_campaign_t), intent(in) :: valid_campaigns(3)
    type(tts_repeatability_campaign_t), intent(out) :: local_altered(3)
    type(tts_repeatability_study_configuration_t), intent(in) :: &
      local_configuration
    type(tts_repeatability_study_result_t) :: invalid_study
    real(dp), parameter :: mismatch_temperatures(5) = &
      [253.15_dp, 273.15_dp, 293.15_dp, 318.15_dp, 333.15_dp]

    local_altered = valid_campaigns
    local_altered(2) = make_repeatability_campaign(mismatch_temperatures, &
      make_linear_temperature_shifts(mismatch_temperatures, &
        reference_temperature_k, 0.01_dp), reference_temperature_k, &
      "GRID-MISMATCH", INDEPENDENT_SPECIMEN_CAMPAIGN)
    invalid_study = analyze_tts_repeatability(local_altered, &
      reference_temperature_k, local_configuration)
    call assert_true(invalid_study%status == &
      TTS_REPEATABILITY_TEMPERATURE_SET_MISMATCH, &
      "Mismatched measured temperature grid reddedilmedi.")

    local_altered = valid_campaigns
    local_altered(2)%identification%source_family%common_state &
      %material_identifier = "OTHER-MATERIAL"
    call assert_incompatible(local_altered, local_configuration, &
      "Material identifier mismatch incompatible-state olmadı.")
    local_altered = valid_campaigns
    local_altered(2)%identification%source_family%common_state &
      %batch_state_identifier = "OTHER-BATCH"
    call assert_incompatible(local_altered, local_configuration, &
      "Batch-state mismatch incompatible-state olmadı.")
    local_altered = valid_campaigns
    local_altered(2)%identification%source_family%common_state &
      %dynamic_strain_amplitude_ratio = 0.01_dp
    call assert_incompatible(local_altered, local_configuration, &
      "Strain-amplitude mismatch incompatible-state olmadı.")
    local_altered = valid_campaigns
    local_altered(2)%identification%source_family%common_state &
      %static_prestrain_ratio = 0.03_dp
    call assert_incompatible(local_altered, local_configuration, &
      "Static-prestrain mismatch incompatible-state olmadı.")
    local_altered = valid_campaigns
    local_altered(2)%identification%source_family%common_state &
      %deformation_mode = 99
    call assert_incompatible(local_altered, local_configuration, &
      "Deformation-mode mismatch incompatible-state olmadı.")

    invalid_study = analyze_tts_repeatability(valid_campaigns, 300.0_dp, &
      local_configuration)
    call assert_true(invalid_study%status == &
      TTS_REPEATABILITY_REFERENCE_NOT_FOUND, &
      "Measured olmayan common reference reddedilmedi.")

    local_altered = valid_campaigns
    local_altered(2)%campaign_identifier = &
      local_altered(1)%campaign_identifier
    invalid_study = analyze_tts_repeatability(local_altered, &
      reference_temperature_k, local_configuration)
    call assert_true(invalid_study%status == TTS_REPEATABILITY_INVALID_INPUT, &
      "Duplicate campaign identifier reddedilmedi.")

    local_altered = valid_campaigns
    local_altered(1)%identification%pair_shift_results(1)%delta_s = &
      ieee_value(0.0_dp, ieee_quiet_nan)
    invalid_study = analyze_tts_repeatability(local_altered, &
      reference_temperature_k, local_configuration)
    call assert_true(invalid_study%status == TTS_REPEATABILITY_NONFINITE_DATA, &
      "Nonfinite authoritative pair shift reddedilmedi.")
  end subroutine test_validation_failures

  subroutine assert_incompatible( &
      local_campaigns, local_configuration, message)
    type(tts_repeatability_campaign_t), intent(in) :: local_campaigns(3)
    type(tts_repeatability_study_configuration_t), intent(in) :: &
      local_configuration
    character(len=*), intent(in) :: message
    type(tts_repeatability_study_result_t) :: local_study

    local_study = analyze_tts_repeatability(local_campaigns, &
      reference_temperature_k, local_configuration)
    call assert_true(local_study%status == &
      TTS_REPEATABILITY_INCOMPATIBLE_STATE, message)
  end subroutine assert_incompatible

  subroutine test_single_and_same_specimen_semantics( &
      valid_campaigns, local_configuration)
    type(tts_repeatability_campaign_t), intent(in) :: valid_campaigns(3)
    type(tts_repeatability_study_configuration_t), intent(in) :: &
      local_configuration
    type(tts_repeatability_campaign_t) :: local_campaigns(3)
    type(tts_repeatability_study_result_t) :: local_study

    local_study = analyze_tts_repeatability(valid_campaigns(1:1), &
      reference_temperature_k, local_configuration)
    call assert_true(local_study%status == TTS_REPEATABILITY_SUCCESS .and. &
      local_study%descriptive_statistics_available .and. &
      .not. local_study%pair_results(1)%delta_s_statistics &
        %spread_statistics_available .and. &
      .not. local_study%independent_cluster_bootstrap_available, &
      "Tek campaign fake repeatability uncertainty üretti.")

    local_campaigns = valid_campaigns
    local_campaigns%replicate_basis = SAME_SPECIMEN_RERUN
    local_study = analyze_tts_repeatability(local_campaigns, &
      reference_temperature_k, local_configuration)
    call assert_true(local_study%status == TTS_REPEATABILITY_SUCCESS .and. &
      local_study%descriptive_statistics_available .and. &
      local_study%same_specimen_rerun_count == 3 .and. &
      local_study%independent_campaign_count == 0 .and. &
      .not. local_study%independent_cluster_bootstrap_available .and. &
      local_study%bootstrap_status == &
        TTS_REPEATABILITY_NO_INDEPENDENT_CAMPAIGNS, &
      "Same-specimen-only study independent bootstrap iddia etti.")
  end subroutine test_single_and_same_specimen_semantics

end program test_tts_repeatability
