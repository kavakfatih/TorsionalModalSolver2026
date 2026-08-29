program test_tts_uncertainty_propagation
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
    ieee_positive_inf
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, TTS_CHANNEL_STORAGE, &
    TTS_CHANNEL_LOSS
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, &
    tts_dynamic_modulus_uncertainty_isotherm_t, &
    tts_dynamic_modulus_uncertainty_point_t, &
    tts_log_uncertainty_result_t, tts_uncertainty_validation_result_t, &
    tts_uncertainty_log_segment_t, TTS_UNCERTAINTY_SUCCESS, &
    TTS_UNCERTAINTY_INVALID_INPUT, TTS_UNCERTAINTY_NONFINITE_DATA
  use tms_tts_uncertainty_validation, only : &
    propagate_log10_standard_uncertainty, &
    validate_tts_uncertainty_family, build_tts_uncertainty_log_segments
  use tms_tts_test_support, only : make_exact_trs_family, &
    populate_isotherm_from_log_grid, assert_true, assert_close
  use tms_tts_uncertainty_test_support, only : &
    make_relative_uncertainty_family
  implicit none

  type(tts_material_family_t) :: family
  type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty
  type(tts_dynamic_modulus_uncertainty_family_t) :: shuffled
  type(tts_dynamic_modulus_uncertainty_isotherm_t) :: temporary_isotherm
  type(tts_dynamic_modulus_uncertainty_point_t) :: temporary_point
  type(tts_log_uncertainty_result_t) :: propagated
  type(tts_uncertainty_validation_result_t) :: validation
  type(tts_uncertainty_log_segment_t), allocatable :: segments(:)
  real(dp), parameter :: temperatures(2) = [293.15_dp, 313.15_dp]
  real(dp), parameter :: shifts(2) = [0.0_dp, 0.4_dp]
  real(dp), parameter :: frequency_grid(5) = &
    [1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 20.0_dp]
  real(dp) :: log_grid(5)
  real(dp) :: expected
  integer :: i
  integer :: j
  logical :: zero_loss_index_found

  ! y=log10(G) için production propagation, bağımsız sensitivity coefficient
  ! u_y=u_G/(G ln(10)) sonucu ve Pa ölçek invariance'ıyla doğrulanır.
  propagated = propagate_log10_standard_uncertainty(1.0e6_dp, 2.0e4_dp)
  expected = 2.0e4_dp/(1.0e6_dp*log(10.0_dp))
  call assert_true(propagated%valid .and. &
    propagated%status == TTS_UNCERTAINTY_SUCCESS, &
    "Storage log-uncertainty propagation kullanılabilir sonuç vermedi.")
  call assert_close(propagated%standard_uncertainty, expected, 2.0e-15_dp, &
    "Storage u_G -> u_log10 propagation hatalı.")
  call assert_close(propagated%variance, expected**2, 2.0e-15_dp, &
    "Log-modulus variance u_log10^2 değil.")
  propagated = propagate_log10_standard_uncertainty(2.5e5_dp, 1.25e4_dp)
  expected = 1.25e4_dp/(2.5e5_dp*log(10.0_dp))
  call assert_close(propagated%standard_uncertainty, expected, 2.0e-15_dp, &
    "Loss u_G -> u_log10 propagation hatalı.")
  propagated = propagate_log10_standard_uncertainty(2.5e8_dp, 1.25e7_dp)
  call assert_close(propagated%standard_uncertainty, expected, 2.0e-15_dp, &
    "G ve u_G ortak Pa ölçeğinde çarpılınca u_log10 değişti.")

  family = make_exact_trs_family(temperatures, shifts)
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  shuffled = uncertainty
  temporary_isotherm = shuffled%isotherms(1)
  shuffled%isotherms(1) = shuffled%isotherms(2)
  shuffled%isotherms(2) = temporary_isotherm
  do i = 1, size(shuffled%isotherms)
    do j = 1, size(shuffled%isotherms(i)%points)/2
      temporary_point = shuffled%isotherms(i)%points(j)
      shuffled%isotherms(i)%points(j) = shuffled%isotherms(i)%points( &
        size(shuffled%isotherms(i)%points) - j + 1)
      shuffled%isotherms(i)%points( &
        size(shuffled%isotherms(i)%points) - j + 1) = temporary_point
    end do
  end do
  validation = validate_tts_uncertainty_family(family, shuffled)
  call assert_true(validation%valid .and. &
    validation%matched_point_count == 14, &
    "Uncertainty array sırası yerine unique physical (T,f) eşlemesi yapılmadı.")

  ! 1,2,5,10,20 Hz fixture'ında 5 Hz uncertainty eksikliği weighted support'u
  ! [1,2] ve [10,20] olarak ikiye bölmelidir; 2->10 gap bridge edilmez.
  do i = 1, size(log_grid)
    log_grid(i) = log10(frequency_grid(i))
  end do
  do i = 1, size(family%isotherms)
    call populate_isotherm_from_log_grid( &
      family%isotherms(i), log_grid, shifts(i), 0.25_dp, 0.15_dp)
  end do
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  uncertainty%isotherms(1)%points(3)%storage_uncertainty_available = .false.
  segments = build_tts_uncertainty_log_segments( &
    family%isotherms(1), uncertainty, TTS_CHANNEL_STORAGE)
  call assert_true(size(segments) == 2, &
    "Missing uncertainty contiguous support gap üretmedi.")
  call assert_close(segments(1)%x(1), log10(1.0_dp), 1.0e-15_dp, &
    "İlk uncertainty segment lower endpoint'i hatalı.")
  call assert_close(segments(1)%x(2), log10(2.0_dp), 1.0e-15_dp, &
    "İlk uncertainty segment upper endpoint'i hatalı.")
  call assert_close(segments(2)%x(1), log10(10.0_dp), 1.0e-15_dp, &
    "İkinci uncertainty segment lower endpoint'i hatalı.")
  call assert_close(segments(2)%x(2), log10(20.0_dp), 1.0e-15_dp, &
    "İkinci uncertainty segment upper endpoint'i hatalı.")

  ! Available=true ile zero/negative/NaN/Inf standard uncertainty geçersizdir;
  ! epsilon weight üretilmez ve expected data failure status ile döner.
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  uncertainty%isotherms(1)%points(1)%storage_standard_uncertainty_pa = 0.0_dp
  validation = validate_tts_uncertainty_family(family, uncertainty)
  call assert_true(.not. validation%valid .and. &
    validation%status == TTS_UNCERTAINTY_INVALID_INPUT, &
    "Zero standard uncertainty reddedilmedi.")
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  uncertainty%isotherms(1)%points(1)%loss_standard_uncertainty_pa = -1.0_dp
  validation = validate_tts_uncertainty_family(family, uncertainty)
  call assert_true(.not. validation%valid .and. &
    validation%status == TTS_UNCERTAINTY_INVALID_INPUT, &
    "Negative standard uncertainty reddedilmedi.")
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  uncertainty%isotherms(1)%points(1)%storage_standard_uncertainty_pa = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  validation = validate_tts_uncertainty_family(family, uncertainty)
  call assert_true(.not. validation%valid .and. &
    validation%status == TTS_UNCERTAINTY_NONFINITE_DATA, &
    "NaN standard uncertainty reddedilmedi.")
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  uncertainty%isotherms(1)%points(1)%loss_standard_uncertainty_pa = &
    ieee_value(0.0_dp, ieee_positive_inf)
  validation = validate_tts_uncertainty_family(family, uncertainty)
  call assert_true(.not. validation%valid .and. &
    validation%status == TTS_UNCERTAINTY_NONFINITE_DATA, &
    "Inf standard uncertainty reddedilmedi.")

  ! VALID G''=0, pozitif loss uncertainty olsa bile log-loss support'a
  ! giremez. Passive runtime semantiği epsilon ile weighted log'a çevrilmez.
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  family%isotherms(1)%points(3)%loss_modulus_pa = 0.0_dp
  segments = build_tts_uncertainty_log_segments( &
    family%isotherms(1), uncertainty, TTS_CHANNEL_LOSS)
  zero_loss_index_found = .false.
  do i = 1, size(segments)
    zero_loss_index_found = zero_loss_index_found .or. &
      any(segments(i)%source_point_indices == 3)
  end do
  call assert_true(size(segments) == 2 .and. .not. zero_loss_index_found, &
    "G''=0 weighted log-loss support'tan çıkarılmadı.")

  print *, "V0.8.4 uncertainty propagation, key matching ve gap doğrulandı."
end program test_tts_uncertainty_propagation
