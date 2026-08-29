program test_tts_runtime_roundtrip
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_identification_result_t, TTS_IDENTIFICATION_SUCCESS
  use tms_tts_identification, only : identify_tts_master_curve
  use tms_tts_runtime_export, only : tts_runtime_export_t, &
    create_tts_runtime_export
  use tms_dynamic_modulus_provider, only : dynamic_modulus_evaluation_t, &
    evaluate_dynamic_shear_modulus
  use tms_tts_test_support, only : make_exact_trs_family, &
    truth_storage_modulus, truth_loss_modulus, assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures(3) = &
    [273.15_dp, 293.15_dp, 313.15_dp]
  real(dp), parameter :: shifts(3) = [1.0_dp, 0.0_dp, -1.0_dp]
  real(dp), parameter :: physical_frequencies(3) = &
    [0.1_dp, 1.0_dp, 10.0_dp]
  type(tts_material_family_t) :: family
  type(tts_identification_result_t) :: identification
  type(tts_runtime_export_t) :: runtime_export
  type(dynamic_modulus_evaluation_t) :: evaluation
  real(dp) :: expected_reduced_frequency
  integer :: i

  family = make_exact_trs_family(temperatures, shifts)
  identification = identify_tts_master_curve(family, "ISO-2")
  call assert_true(identification%status == TTS_IDENTIFICATION_SUCCESS .and. &
    identification%runtime_export_ready, &
    "Runtime round-trip identification hazır değil.")
  runtime_export = create_tts_runtime_export(identification)
  call assert_true(runtime_export%available, &
    "V0.8.0 provider nesneleri V0.8.1 output'tan oluşturulamadı.")

  ! Üç farklı physical (f,T) sorgusu aynı reduced master point f_r=1 Hz'e
  ! düşer. Expected G'/G'' değerleri provider'dan değil independent synthetic
  ! truth fonksiyonundan alınır; canonical f_r=a_T*f convention'ı korunur.
  do i = 1, 3
    evaluation = evaluate_dynamic_shear_modulus( &
      runtime_export%thermorheological_provider, physical_frequencies(i), &
      temperatures(i))
    expected_reduced_frequency = &
      physical_frequencies(i)*10.0_dp**shifts(i)
    call assert_close(evaluation%lookup_frequency_hz, &
      expected_reduced_frequency, 2.0e-6_dp, &
      "Round-trip f_r=a_T*f lookup frequency hatalı.")
    call assert_close(evaluation%physical_frequency_hz, &
      physical_frequencies(i), 0.0_dp, &
      "Round-trip physical frequency trace değişti.")
    call assert_close(evaluation%modulus%temperature, temperatures(i), &
      0.0_dp, "Round-trip operating temperature trace değişti.")
    call assert_close(evaluation%modulus%storage_modulus, &
      truth_storage_modulus(expected_reduced_frequency), 2.0e-6_dp, &
      "Round-trip G'(f,T) independent truth ile uyuşmuyor.")
    call assert_close(evaluation%modulus%loss_modulus, &
      truth_loss_modulus(expected_reduced_frequency), 2.0e-6_dp, &
      "Round-trip G''(f,T) independent truth ile uyuşmuyor.")
    call assert_close(evaluation%log10_a_t, shifts(i), 2.0e-6_dp, &
      "Empirical tabulated log10(a_T) round-trip hatalı.")
  end do

  print *, "V0.8.1 identification -> V0.8.0 provider round-trip doğrulandı."
end program test_tts_runtime_roundtrip
