program test_tts_pair_shift
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_log_segment_t, &
    tts_pair_shift_result_t, tts_pair_objective_evaluation_t, &
    TTS_CHANNEL_LOSS, TTS_CHANNEL_JOINT, SHIFT_FROM_JOINT, &
    SHIFT_FROM_STORAGE_ONLY, PAIR_SHIFT_SUCCESS, &
    PAIR_SHIFT_STORAGE_ONLY, PAIR_SHIFT_INSUFFICIENT_SUPPORT, &
    PAIR_SHIFT_NO_INTERIOR_MINIMUM, BELOW_RELIABLE_FLOOR, &
    MEASUREMENT_VALID, MEASUREMENT_UNAVAILABLE, &
    is_loss_log_usable, is_runtime_export_usable
  use tms_tts_pair_shift, only : identify_tts_pair_shift, &
    evaluate_tts_pair_objective, build_tts_valid_log_segments, &
    get_tts_pair_feasible_shift_domain, integrate_linear_squared_residual
  use tms_tts_test_support, only : make_exact_trs_family, &
    populate_isotherm_from_log_grid, assert_true, assert_close
  implicit none

  type(tts_material_family_t) :: family
  type(tts_pair_shift_result_t) :: pair
  type(tts_pair_objective_evaluation_t) :: objective
  type(tts_log_segment_t), allocatable :: segments(:)
  real(dp) :: lower_shift
  real(dp) :: upper_shift
  logical :: domain_valid
  real(dp), parameter :: temperatures(2) = [293.15_dp, 313.15_dp]
  real(dp), parameter :: truth_shifts(2) = [0.0_dp, 0.75_dp]
  real(dp), parameter :: irregular_grid(8) = &
    [-2.0_dp, -1.7_dp, -1.1_dp, -0.2_dp, 0.35_dp, 0.9_dp, 1.55_dp, 2.0_dp]
  integer :: i

  family = make_exact_trs_family(temperatures, truth_shifts)
  pair = identify_tts_pair_shift(family%isotherms(1), family%isotherms(2))
  call assert_true(pair%status == PAIR_SHIFT_SUCCESS .and. &
    pair%production_channel == SHIFT_FROM_JOINT, &
    "Exact TRS pair joint production shift üretmedi.")
  call assert_close(pair%delta_s, 0.75_dp, 8.0e-7_dp, &
    "Exact adjacent relative shift geri kazanılamadı.")
  call assert_close(pair%delta_s_storage, 0.75_dp, 8.0e-7_dp, &
    "Storage diagnostic shift hatalı.")
  call assert_close(pair%delta_s_loss, 0.75_dp, 8.0e-7_dp, &
    "Loss diagnostic shift hatalı.")
  call assert_true(pair%objective_minimum < 1.0e-12_dp .and. &
    pair%overlap_width_decades > 0.0_dp .and. &
    pair%overlap_fraction > 0.0_dp .and. &
    pair%objective_curvature > 0.0_dp, &
    "Pair residual/overlap/curvature tanıları eksik.")

  call get_tts_pair_feasible_shift_domain(family%isotherms(1), &
    family%isotherms(2), lower_shift, upper_shift, domain_valid)
  call assert_true(domain_valid, "Measured feasible shift domain kurulamadı.")
  call assert_close(lower_shift, -4.0_dp, 1.0e-14_dp, &
    "Feasible lower shift measured domain'den türetilmedi.")
  call assert_close(upper_shift, 4.0_dp, 1.0e-14_dp, &
    "Feasible upper shift measured domain'den türetilmedi.")
  call assert_close(integrate_linear_squared_residual(2.0_dp, 1.0_dp, &
    3.0_dp), 26.0_dp/3.0_dp, 1.0e-14_dp, &
    "Exact linear residual integral formülü hatalı.")
  objective = evaluate_tts_pair_objective(family%isotherms(1), &
    family%isotherms(2), 0.75_dp, TTS_CHANNEL_JOINT)
  call assert_true(objective%valid .and. objective%objective < 1.0e-25_dp, &
    "Known shift exact joint objective sıfıra yaklaşmadı.")

  ! Irregular frequency spacing ve farklı sample density, grid-sampling'e
  ! bağımlı olmadan aynı piecewise-linear truth shift'ini vermelidir.
  call populate_isotherm_from_log_grid(family%isotherms(2), irregular_grid, &
    0.75_dp, 0.25_dp, 0.15_dp)
  pair = identify_tts_pair_shift(family%isotherms(1), family%isotherms(2))
  call assert_true(pair%status == PAIR_SHIFT_SUCCESS, &
    "Irregular/sample-density pair shift başarısız.")
  call assert_close(pair%delta_s, 0.75_dp, 8.0e-7_dp, &
    "Irregular/sample-density shift değişti.")

  ! Loss quality gap iki ayrı contiguous segment üretmeli; floor noktasının
  ! üzerinden log-G'' interpolation yapılmamalıdır.
  family = make_exact_trs_family(temperatures, truth_shifts)
  family%isotherms(1)%points(3)%loss_quality = BELOW_RELIABLE_FLOOR
  segments = build_tts_valid_log_segments(family%isotherms(1), TTS_CHANNEL_LOSS)
  call assert_true(size(segments) == 2, &
    "Quality gap contiguous loss segmentlerini ayırmadı.")
  objective = evaluate_tts_pair_objective(family%isotherms(1), &
    family%isotherms(2), 0.75_dp, TTS_CHANNEL_LOSS)
  call assert_true(objective%valid .and. &
    objective%interpolation_interval_count >= 2, &
    "Quality-gap objective yalnız gerçek segment overlap'larını kullanmadı.")

  ! VALID G''=0 passivity/runtime açısından geçerli, log-loss objective için
  ! kullanılamaz; epsilon substitution yapılmaz.
  family%isotherms(1)%points(2)%loss_modulus_pa = 0.0_dp
  family%isotherms(1)%points(2)%loss_quality = MEASUREMENT_VALID
  call assert_true(.not. is_loss_log_usable(family%isotherms(1)%points(2)), &
    "Zero loss log-objective'e yanlışlıkla alındı.")
  call assert_true(is_runtime_export_usable(family%isotherms(1)%points(2)), &
    "Zero valid loss runtime export'tan yanlışlıkla çıkarıldı.")

  ! G'' pair support yoksa production shift G' üzerinden açık storage-only
  ! status ile sağlanır; full complex TRS desteği iddia edilmez.
  family = make_exact_trs_family(temperatures, truth_shifts)
  do i = 1, size(family%isotherms(2)%points)
    family%isotherms(2)%points(i)%loss_quality = MEASUREMENT_UNAVAILABLE
  end do
  pair = identify_tts_pair_shift(family%isotherms(1), family%isotherms(2))
  call assert_true(pair%status == PAIR_SHIFT_STORAGE_ONLY .and. &
    pair%production_channel == SHIFT_FROM_STORAGE_ONLY .and. &
    pair%shift_available .and. .not. pair%joint_shift_available, &
    "Storage-only fallback açık diagnostic üretmedi.")
  call assert_close(pair%delta_s, 0.75_dp, 8.0e-7_dp, &
    "Storage-only relative shift hatalı.")

  ! Yalnız G'' support diagnostic loss shift üretebilir fakat production için
  ! gerekli G' bulunmadığından pair başarı sayılmamalıdır.
  family = make_exact_trs_family(temperatures, truth_shifts)
  do i = 1, size(family%isotherms(2)%points)
    family%isotherms(2)%points(i)%storage_quality = MEASUREMENT_UNAVAILABLE
  end do
  pair = identify_tts_pair_shift(family%isotherms(1), family%isotherms(2))
  call assert_true(.not. pair%shift_available .and. pair%loss_shift_available &
    .and. pair%status == PAIR_SHIFT_INSUFFICIENT_SUPPORT, &
    "G''-only diagnostic production shift gibi kabul edildi.")

  ! Hiçbir channel interval support taşımıyorsa clean status dönmelidir.
  family = make_exact_trs_family(temperatures, truth_shifts)
  do i = 1, size(family%isotherms(2)%points)
    family%isotherms(2)%points(i)%storage_quality = MEASUREMENT_UNAVAILABLE
    family%isotherms(2)%points(i)%loss_quality = MEASUREMENT_UNAVAILABLE
  end do
  pair = identify_tts_pair_shift(family%isotherms(1), family%isotherms(2))
  call assert_true(pair%status == PAIR_SHIFT_INSUFFICIENT_SUPPORT .and. &
    .not. pair%shift_available, "No-support pair clean failure vermedi.")

  ! Tam plateau objective'in scan boyunca boundary/equal minimum üretmesi
  ! Brent'i çağırmamalı ve no-interior-minimum status'u vermelidir.
  family = make_exact_trs_family(temperatures, truth_shifts, 0.0_dp, 0.0_dp)
  pair = identify_tts_pair_shift(family%isotherms(1), family%isotherms(2))
  call assert_true(pair%status == PAIR_SHIFT_NO_INTERIOR_MINIMUM .and. &
    .not. pair%shift_available, &
    "Plateau/boundary minimum Brent başarısı gibi raporlandı.")

  print *, "V0.8.1 exact pair shift ve quality semantics doğrulandı."
end program test_tts_pair_shift
