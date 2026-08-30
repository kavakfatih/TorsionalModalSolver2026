program test_tts_covariance_propagation
  use tms_kinds, only : dp
  use tms_tts_covariance_types, only : tts_covariance_matrix_2x2_t, &
    tts_polar_covariance_propagation_result_t, &
    tts_log_covariance_propagation_result_t, TTS_COVARIANCE_SUCCESS
  use tms_tts_covariance_propagation, only : &
    propagate_magnitude_phase_covariance, &
    propagate_log10_modulus_covariance
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  type(tts_covariance_matrix_2x2_t) :: covariance
  type(tts_covariance_matrix_2x2_t) :: covariance_mpa
  type(tts_polar_covariance_propagation_result_t) :: polar
  type(tts_log_covariance_propagation_result_t) :: logarithmic
  type(tts_log_covariance_propagation_result_t) :: logarithmic_mpa
  real(dp), parameter :: magnitude_pa = 2.0e6_dp
  real(dp), parameter :: phase_rad = 0.4_dp
  real(dp), parameter :: variance_magnitude = (2.0e4_dp)**2
  real(dp), parameter :: variance_phase = (8.0e-3_dp)**2
  real(dp) :: cosine
  real(dp) :: expected_covariance
  real(dp) :: expected_loss_variance
  real(dp) :: expected_storage_variance
  real(dp) :: linear_correlation
  real(dp) :: sine

  ! Polar helper Jacobian referansı testte bağımsız açık denklemlerle kurulur.
  ! Cov(M,delta)=0 olsa da ortak M/delta girdileri G' ile G'' arasında sıfır
  ! olmayan covariance üretebilir; bağımsız polar girdiler bağımsız G kanalları
  ! anlamına gelmez.
  polar = propagate_magnitude_phase_covariance(magnitude_pa, phase_rad, &
    variance_magnitude, variance_phase, 0.0_dp)
  cosine = cos(phase_rad)
  sine = sin(phase_rad)
  expected_storage_variance = cosine**2*variance_magnitude + &
    magnitude_pa**2*sine**2*variance_phase
  expected_loss_variance = sine**2*variance_magnitude + &
    magnitude_pa**2*cosine**2*variance_phase
  expected_covariance = sine*cosine*variance_magnitude - &
    magnitude_pa**2*sine*cosine*variance_phase
  call assert_true(polar%status == TTS_COVARIANCE_SUCCESS .and. &
    polar%valid, "Magnitude/phase covariance propagation başarısız.")
  call assert_close(polar%storage_modulus_pa, magnitude_pa*cosine, &
    2.0e-15_dp, "Polar G' dönüşümü hatalı.")
  call assert_close(polar%loss_modulus_pa, magnitude_pa*sine, &
    2.0e-15_dp, "Polar G'' dönüşümü hatalı.")
  call assert_close(polar%covariance%storage_variance_pa2, &
    expected_storage_variance, 3.0e-15_dp, "Polar Var(G') hatalı.")
  call assert_close(polar%covariance%loss_variance_pa2, &
    expected_loss_variance, 3.0e-15_dp, "Polar Var(G'') hatalı.")
  call assert_close(polar%covariance%storage_loss_covariance_pa2, &
    expected_covariance, 3.0e-15_dp, "Polar Cov(G',G'') hatalı.")
  call assert_true(abs(expected_covariance) > 0.0_dp .and. &
    polar%first_order_covariance_propagation, &
    "Independent M/delta nonzero G covariance veya provenance üretmedi.")

  ! Physical Pa^2 covariance, D=diag(1/(G ln10)) ile dimensionless log-space
  ! covariance'ya dönüşür. Positive diagonal scaling correlation'ı değiştirmez.
  covariance%storage_variance_pa2 = (2.0e4_dp)**2
  covariance%loss_variance_pa2 = (1.0e4_dp)**2
  covariance%storage_loss_covariance_pa2 = &
    0.35_dp*2.0e4_dp*1.0e4_dp
  logarithmic = propagate_log10_modulus_covariance( &
    2.0e6_dp, 5.0e5_dp, covariance)
  linear_correlation = covariance%storage_loss_covariance_pa2 / &
    sqrt(covariance%storage_variance_pa2*covariance%loss_variance_pa2)
  call assert_true(logarithmic%valid .and. &
    logarithmic%first_order_covariance_propagation, &
    "Log covariance first-order propagation sonucu geçersiz.")
  call assert_close(logarithmic%storage_variance, &
    covariance%storage_variance_pa2/(2.0e6_dp*log(10.0_dp))**2, &
    2.0e-15_dp, "Var(log10 G') propagation hatalı.")
  call assert_close(logarithmic%loss_variance, &
    covariance%loss_variance_pa2/(5.0e5_dp*log(10.0_dp))**2, &
    2.0e-15_dp, "Var(log10 G'') propagation hatalı.")
  call assert_close(logarithmic%storage_loss_covariance, &
    covariance%storage_loss_covariance_pa2 / &
      (2.0e6_dp*5.0e5_dp*log(10.0_dp)**2), &
    2.0e-15_dp, "Cov(log10 G',log10 G'') propagation hatalı.")
  call assert_close(logarithmic%correlation, linear_correlation, &
    3.0e-15_dp, "Linear/log correlation invariance bozuldu.")

  ! Pa -> MPa ile G ve uncertainty/covariance birlikte ölçeklendiğinde
  ! dimensionless log covariance aynıdır; helper unit-system bağımsızdır.
  covariance_mpa%storage_variance_pa2 = &
    covariance%storage_variance_pa2*1.0e-12_dp
  covariance_mpa%loss_variance_pa2 = &
    covariance%loss_variance_pa2*1.0e-12_dp
  covariance_mpa%storage_loss_covariance_pa2 = &
    covariance%storage_loss_covariance_pa2*1.0e-12_dp
  logarithmic_mpa = propagate_log10_modulus_covariance( &
    2.0_dp, 0.5_dp, covariance_mpa)
  call assert_true(logarithmic_mpa%valid, &
    "MPa-consistent covariance propagation geçersiz.")
  call assert_close(logarithmic_mpa%storage_variance, &
    logarithmic%storage_variance, 2.0e-15_dp, &
    "Pa/MPa Var(log G') unit invariance bozuldu.")
  call assert_close(logarithmic_mpa%loss_variance, &
    logarithmic%loss_variance, 2.0e-15_dp, &
    "Pa/MPa Var(log G'') unit invariance bozuldu.")
  call assert_close(logarithmic_mpa%storage_loss_covariance, &
    logarithmic%storage_loss_covariance, 2.0e-15_dp, &
    "Pa/MPa log covariance unit invariance bozuldu.")

  logarithmic = propagate_log10_modulus_covariance( &
    2.0e6_dp, 0.0_dp, covariance)
  call assert_true(.not. logarithmic%valid, &
    "G''=0 bivariate log objective için reddedilmedi.")
  logarithmic = propagate_log10_modulus_covariance( &
    -2.0e6_dp, 5.0e5_dp, covariance)
  call assert_true(.not. logarithmic%valid, &
    "Negatif G' bivariate log propagation'da reddedilmedi.")

  print *, "V0.8.5 covariance propagation ve invariants doğrulandı."
end program test_tts_covariance_propagation
