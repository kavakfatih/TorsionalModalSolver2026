program test_tts_identification
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_identification_result_t, tts_pair_shift_result_t, &
    TTS_IDENTIFICATION_SUCCESS, TTS_IDENTIFICATION_INVALID_INPUT, &
    TTS_IDENTIFICATION_CHAIN_BROKEN, &
    TTS_IDENTIFICATION_RUNTIME_DOMAIN_GAP, BELOW_RELIABLE_FLOOR, &
    MEASUREMENT_UNAVAILABLE
  use tms_tts_identification, only : identify_tts_master_curve
  use tms_tts_runtime_export, only : tts_runtime_export_t, &
    create_tts_runtime_export
  use tms_tts_pair_shift, only : identify_tts_pair_shift
  use tms_tts_test_support, only : make_exact_trs_family, &
    replace_loss_truth_shift, find_empirical_shift, assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: exact_shifts(5) = &
    [2.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, -2.0_dp]
  real(dp), parameter :: nontrs_loss_shifts(5) = &
    [2.6_dp, 1.25_dp, 0.0_dp, -1.25_dp, -2.6_dp]
  type(tts_material_family_t) :: family
  type(tts_identification_result_t) :: exact_result
  type(tts_identification_result_t) :: nontrs_result
  type(tts_identification_result_t) :: result
  type(tts_pair_shift_result_t) :: strong_pair
  type(tts_pair_shift_result_t) :: weak_pair
  type(tts_runtime_export_t) :: unsafe_export
  real(dp) :: exact_objective_sum
  real(dp) :: nontrs_objective_sum
  real(dp) :: stored_original_value
  integer :: i

  family = make_exact_trs_family(temperatures, exact_shifts)
  stored_original_value = family%isotherms(1)%points(1)%storage_modulus_pa
  exact_result = identify_tts_master_curve(family, "ISO-3")
  call assert_true(exact_result%status == TTS_IDENTIFICATION_SUCCESS .and. &
    exact_result%shift_chain_available .and. &
    exact_result%master_cloud_available .and. &
    exact_result%runtime_export_ready, &
    "Exact synthetic TRS top-level identification başarısız.")
  do i = 1, 5
    call assert_close(find_empirical_shift(exact_result%empirical_shifts, i), &
      exact_shifts(i), 1.5e-6_dp, &
      "Exact synthetic absolute empirical shift hatalı.")
  end do
  call assert_true(exact_result%diagnostics%full_complex_pair_support .and. &
    size(exact_result%master_cloud) == 35 .and. &
    size(exact_result%runtime_master_table) >= 7 .and. &
    size(exact_result%diagnostics%vgp_points) == 35 .and. &
    size(exact_result%diagnostics%cole_cole_points) == 35, &
    "Exact TRS result model/evidence alanları eksik.")
  exact_objective_sum = sum(exact_result%pair_shift_results%objective_minimum)

  ! Input family daha sonra değiştirilse bile result authoritative deep copy ve
  ! master cloud değerleri değişmemelidir.
  family%isotherms(1)%points(1)%storage_modulus_pa = 99.0e9_dp
  call assert_close(exact_result%source_family%isotherms(1)%points(1) &
    %storage_modulus_pa, stored_original_value, 0.0_dp, &
    "Identification result input mutation'dan etkilendi.")
  call assert_close(exact_result%master_cloud(1)%storage_modulus_pa, &
    stored_original_value, 0.0_dp, &
    "Master cloud original measurement overwrite edildi.")

  ! Deterministic non-TRS data G' ve G'' için farklı horizontal shifts taşır.
  ! Algorithm bunu zorla kusursuz collapse yapmamalı; residual ve channel-shift
  ! discrepancy exact TRS'den belirgin yüksek kalmalıdır.
  family = make_exact_trs_family(temperatures, exact_shifts)
  do i = 1, 5
    call replace_loss_truth_shift(family%isotherms(i), nontrs_loss_shifts(i))
  end do
  nontrs_result = identify_tts_master_curve(family, "ISO-3")
  call assert_true(nontrs_result%status == TTS_IDENTIFICATION_SUCCESS, &
    "Deterministic non-TRS dataset numerical result üretemedi.")
  nontrs_objective_sum = &
    sum(nontrs_result%pair_shift_results%objective_minimum)
  call assert_true(nontrs_objective_sum > exact_objective_sum + 1.0e-6_dp, &
    "Non-TRS joint residual exact TRS'den ayrışmadı.")
  call assert_true(maxval(nontrs_result%pair_shift_results &
    %storage_loss_shift_discrepancy) > 0.1_dp, &
    "Non-TRS storage-loss shift discrepancy görünür değil.")

  ! Nearly flat frequency response numerical minimum üretebilir fakat curvature
  ! güçlü tanımlanmış response'tan çok daha düşük olmalıdır.
  family = make_exact_trs_family([293.15_dp, 313.15_dp], &
    [0.0_dp, 0.75_dp])
  strong_pair = identify_tts_pair_shift(family%isotherms(1), &
    family%isotherms(2))
  family = make_exact_trs_family([293.15_dp, 313.15_dp], &
    [0.0_dp, 0.75_dp], 1.0e-4_dp, 1.0e-4_dp)
  weak_pair = identify_tts_pair_shift(family%isotherms(1), &
    family%isotherms(2))
  call assert_true(weak_pair%shift_available .and. &
    weak_pair%objective_curvature < &
      1.0e-4_dp*strong_pair%objective_curvature, &
    "Weak-identifiability plateau curvature diagnostic ile ayrışmadı.")

  ! Invalid experimental ordering status döndürmeli; process error-stop ile
  ! sonlanmamalı ve silent sort uygulanmamalıdır.
  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms(1)%points(2)%frequency_hz = &
    family%isotherms(1)%points(1)%frequency_hz
  result = identify_tts_master_curve(family, "ISO-1")
  call assert_true(result%status == TTS_IDENTIFICATION_INVALID_INPUT .and. &
    .not. result%runtime_export_ready, &
    "Duplicate/unordered frequency clean invalid-input status vermedi.")

  ! Bir zorunlu adjacent link support taşımıyorsa complete empirical shift ve
  ! runtime export unavailable kalmalıdır.
  family = make_exact_trs_family(temperatures, exact_shifts)
  do i = 1, size(family%isotherms(5)%points)
    family%isotherms(5)%points(i)%storage_quality = MEASUREMENT_UNAVAILABLE
    family%isotherms(5)%points(i)%loss_quality = MEASUREMENT_UNAVAILABLE
  end do
  result = identify_tts_master_curve(family, "ISO-3")
  call assert_true(result%status == TTS_IDENTIFICATION_CHAIN_BROKEN .and. &
    .not. result%shift_chain_available .and. &
    .not. result%runtime_export_ready, &
    "Broken chain partial runtime export üretti.")

  ! Pair alignment mümkün olsa da bütün isotherm'lerde aynı internal loss
  ! quality hole kalırsa master cloud korunmalı, continuous runtime export ise
  ! explicit domain-gap status'u ile reddedilmelidir.
  family = make_exact_trs_family( &
    [293.15_dp, 313.15_dp], [0.0_dp, 0.0_dp])
  family%isotherms(1)%points(4)%loss_quality = BELOW_RELIABLE_FLOOR
  family%isotherms(2)%points(4)%loss_quality = BELOW_RELIABLE_FLOOR
  result = identify_tts_master_curve(family, "ISO-1")
  call assert_true( &
    result%status == TTS_IDENTIFICATION_RUNTIME_DOMAIN_GAP .and. &
    result%master_cloud_available .and. .not. result%runtime_export_ready, &
    "Internal quality hole top-level runtime readiness'i kapatmadı.")
  unsafe_export = create_tts_runtime_export(result)
  call assert_true(.not. unsafe_export%available, &
    "Unsafe identification sonucundan runtime provider oluşturuldu.")

  print *, "V0.8.1 exact/non-TRS/weak identification gates doğrulandı."
end program test_tts_identification
