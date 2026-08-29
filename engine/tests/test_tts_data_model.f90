program test_tts_data_model
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_validation_result_t, &
    TTS_IDENTIFICATION_SUCCESS, TTS_IDENTIFICATION_INVALID_INPUT, &
    TTS_IDENTIFICATION_INCONSISTENT_STATE, &
    TTS_IDENTIFICATION_INSUFFICIENT_ISOTHERMS, MEASUREMENT_VALID, &
    BELOW_RELIABLE_FLOOR, MEASUREMENT_UNAVAILABLE, MEASUREMENT_REJECTED, &
    validate_tts_material_family, is_measurement_quality_known
  use tms_tts_test_support, only : make_exact_trs_family, assert_true
  implicit none

  type(tts_material_family_t) :: family
  type(tts_validation_result_t) :: validation

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms(2)%specimen_identifier = "DIFFERENT-SPECIMEN"
  validation = validate_tts_material_family(family)
  call assert_true(validation%status == TTS_IDENTIFICATION_SUCCESS .and. &
    validation%valid, &
    "Aynı common state içindeki farklı specimen kimliği reddedildi.")
  call assert_true(is_measurement_quality_known(MEASUREMENT_VALID) .and. &
    is_measurement_quality_known(BELOW_RELIABLE_FLOOR) .and. &
    is_measurement_quality_known(MEASUREMENT_UNAVAILABLE) .and. &
    is_measurement_quality_known(MEASUREMENT_REJECTED) .and. &
    .not. is_measurement_quality_known(99), &
    "Explicit measurement quality enum sözleşmesi hatalı.")

  family%common_state%batch_state_identifier = ""
  validation = validate_tts_material_family(family)
  call assert_true(validation%status == TTS_IDENTIFICATION_INVALID_INPUT, &
    "Eksik authoritative batch/material-state reddedilmedi.")

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%common_state%dynamic_strain_amplitude_ratio = -0.01_dp
  validation = validate_tts_material_family(family)
  call assert_true(validation%status == TTS_IDENTIFICATION_INCONSISTENT_STATE, &
    "Geçersiz common strain state reddedilmedi.")

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms(1)%temperature_k = 0.0_dp
  validation = validate_tts_material_family(family)
  call assert_true(validation%status == TTS_IDENTIFICATION_INVALID_INPUT, &
    "Pozitif olmayan Kelvin sıcaklığı reddedilmedi.")

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms(1)%points(1)%storage_modulus_pa = 0.0_dp
  family%isotherms(1)%points(1)%storage_quality = MEASUREMENT_VALID
  validation = validate_tts_material_family(family)
  call assert_true(.not. validation%valid, &
    "VALID G'<=0 experimental point reddedilmedi.")

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms(1)%points(1)%loss_modulus_pa = -1.0_dp
  family%isotherms(1)%points(1)%loss_quality = MEASUREMENT_VALID
  validation = validate_tts_material_family(family)
  call assert_true(.not. validation%valid, &
    "VALID G''<0 experimental point reddedilmedi.")

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms(1)%points(1)%loss_quality = 99
  validation = validate_tts_material_family(family)
  call assert_true(.not. validation%valid, &
    "Magic numeric measurement quality reddedilmedi.")

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms = family%isotherms(1:1)
  validation = validate_tts_material_family(family)
  call assert_true(validation%status == &
    TTS_IDENTIFICATION_INSUFFICIENT_ISOTHERMS, &
    "Tek isotherm insufficient-isotherms status vermedi.")

  print *, "V0.8.1 experimental data/quality model doğrulandı."
end program test_tts_data_model
