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
  real(dp) :: expected_loss_modulus
  real(dp) :: expected_reduced_frequency
  real(dp) :: expected_storage_modulus
  real(dp) :: intermediate_reduced_frequency
  real(dp) :: off_reference_physical_frequency
  integer :: lower_index
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

  ! Exact master point dışında iki adjacent runtime noktası arasındaki
  ! geometric midpoint seçilir. LINEAR_LOG_FREQUENCY alpha=0.5 üretirken G'
  ! ve G'' değerlerinin kendileri linear interpolate edilir; log(G) alınmaz.
  lower_index = size(identification%runtime_master_table)/2
  intermediate_reduced_frequency = 10.0_dp**(0.5_dp*( &
    log10(identification%runtime_master_table(lower_index) &
      %reduced_frequency_hz) + &
    log10(identification%runtime_master_table(lower_index + 1) &
      %reduced_frequency_hz)))
  expected_storage_modulus = 0.5_dp*( &
    identification%runtime_master_table(lower_index)%storage_modulus_pa + &
    identification%runtime_master_table(lower_index + 1)%storage_modulus_pa)
  expected_loss_modulus = 0.5_dp*( &
    identification%runtime_master_table(lower_index)%loss_modulus_pa + &
    identification%runtime_master_table(lower_index + 1)%loss_modulus_pa)

  evaluation = evaluate_dynamic_shear_modulus( &
    runtime_export%thermorheological_provider, &
    intermediate_reduced_frequency, temperatures(2))
  call assert_true(.not. evaluation%exact_table_point, &
    "Intermediate reference query exact master point'e düştü.")
  call assert_close(evaluation%interpolation_alpha, 0.5_dp, 2.0e-12_dp, &
    "LINEAR_LOG_FREQUENCY midpoint alpha değeri hatalı.")
  call assert_close(evaluation%lookup_frequency_hz, &
    intermediate_reduced_frequency, 2.0e-12_dp, &
    "Intermediate reference lookup frequency trace hatalı.")
  call assert_close(evaluation%modulus%storage_modulus, &
    expected_storage_modulus, 2.0e-12_dp, &
    "Intermediate reference G' linear interpolation hatalı.")
  call assert_close(evaluation%modulus%loss_modulus, &
    expected_loss_modulus, 2.0e-12_dp, &
    "Intermediate reference G'' linear interpolation hatalı.")

  ! Aynı reduced midpoint measured off-reference sıcaklıkta physical
  ! f=f_r/a_T ile sorgulanır. Expected modulus yine provider'dan bağımsız
  ! arithmetic interpolation denkleminden gelir ve bütün trace doğrulanır.
  off_reference_physical_frequency = &
    intermediate_reduced_frequency/(10.0_dp**shifts(1))
  evaluation = evaluate_dynamic_shear_modulus( &
    runtime_export%thermorheological_provider, &
    off_reference_physical_frequency, temperatures(1))
  call assert_true(.not. evaluation%exact_table_point, &
    "Off-reference intermediate query exact master point'e düştü.")
  call assert_close(evaluation%physical_frequency_hz, &
    off_reference_physical_frequency, 0.0_dp, &
    "Off-reference physical frequency trace değişti.")
  call assert_close(evaluation%lookup_frequency_hz, &
    intermediate_reduced_frequency, 2.0e-12_dp, &
    "Off-reference reduced lookup frequency hatalı.")
  call assert_close(evaluation%modulus%temperature, temperatures(1), &
    0.0_dp, "Off-reference operating temperature trace hatalı.")
  call assert_close(evaluation%interpolation_alpha, 0.5_dp, 2.0e-12_dp, &
    "Off-reference interpolation alpha değeri hatalı.")
  call assert_close(evaluation%modulus%storage_modulus, &
    expected_storage_modulus, 2.0e-12_dp, &
    "Off-reference intermediate G' interpolation hatalı.")
  call assert_close(evaluation%modulus%loss_modulus, &
    expected_loss_modulus, 2.0e-12_dp, &
    "Off-reference intermediate G'' interpolation hatalı.")
  call assert_close(evaluation%log10_a_t, shifts(1), 2.0e-12_dp, &
    "Off-reference empirical log10(a_T) trace hatalı.")

  print *, "V0.8.1 identification -> V0.8.0 provider round-trip doğrulandı."
end program test_tts_runtime_roundtrip
