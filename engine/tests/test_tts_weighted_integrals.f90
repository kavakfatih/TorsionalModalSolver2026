program test_tts_weighted_integrals
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
    ieee_positive_inf
  use tms_kinds, only : dp
  use tms_weighted_piecewise_integrals, only : &
    weighted_interval_integral_t, integrate_weighted_linear_interval, &
    WEIGHTED_INTEGRAL_SUCCESS, WEIGHTED_INTEGRAL_INVALID_INPUT, &
    WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  type(weighted_interval_integral_t) :: result
  real(dp) :: a
  real(dp) :: alpha
  real(dp) :: b
  real(dp) :: beta
  real(dp) :: c
  real(dp) :: d
  real(dp) :: expected
  real(dp) :: h
  real(dp) :: r0
  real(dp) :: r1
  real(dp) :: variance_perturbation
  real(dp) :: v0
  real(dp) :: v1

  ! Nontrivial b/=0 ve d/=0 vakasının truth değeri test içinde production
  ! helper'dan bağımsız u=c+dt substitution antiderivative'ından kurulur.
  h = 1.2_dp
  a = 0.4_dp
  b = 0.7_dp
  c = 0.5_dp
  d = 0.3_dp
  r0 = a
  r1 = a + b*h
  v0 = c
  v1 = c + d*h
  alpha = b/d
  beta = a - b*c/d
  expected = (0.5_dp*alpha*alpha*(v1*v1 - v0*v0) + &
    2.0_dp*alpha*beta*(v1 - v0) + &
    beta*beta*log(v1/v0))/d
  result = integrate_weighted_linear_interval(h, r0, r1, v0, v1)
  call assert_true(result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_SUCCESS, &
    "Nonconstant weighted analytical integral üretilemedi.")
  call assert_close(result%integral, expected, 2.0e-13_dp, &
    "Analytical r^2/v interval integral hatalı.")

  ! d=0 limitinde variance sabittir ve exact linear-residual L2 integralinin
  ! variance'a bölünmüş hali geri kazanılmalıdır.
  h = 2.0_dp
  r0 = -0.3_dp
  r1 = 0.8_dp
  v0 = 0.04_dp
  v1 = v0
  expected = h*(r0*r0 + r0*r1 + r1*r1)/(3.0_dp*v0)
  result = integrate_weighted_linear_interval(h, r0, r1, v0, v1)
  call assert_close(result%integral, expected, 2.0e-14_dp, &
    "Constant-variance weighted integral limit'i hatalı.")

  ! Machine-small variance eğimi, arbitrary engineering threshold olmadan
  ! constant-d limitine sürekli yaklaşmalıdır.
  variance_perturbation = 0.5_dp*sqrt(sqrt(epsilon(1.0_dp)))
  v1 = v0*(1.0_dp + variance_perturbation)
  result = integrate_weighted_linear_interval(h, r0, r1, v0, v1)
  call assert_true(result%valid, &
    "Nearly-constant variance dalı numerical failure üretti.")
  call assert_close(result%integral, expected, 2.0e-4_dp, &
    "d->0 weighted integral sürekliliği bozuldu.")

  ! Tiny fakat pozitif variance ve buna uyumlu residual ölçeği overflow
  ! üretmeden sonlu kalmalıdır.
  result = integrate_weighted_linear_interval( &
    1.0_dp, 1.0e-160_dp, 2.0e-160_dp, 1.0e-300_dp, 2.0e-300_dp)
  call assert_true(result%valid .and. result%integral > 0.0_dp, &
    "Tiny positive variance ölçeği sonlu değerlendirilemedi.")

  ! Sıfır/negatif/nonfinite variance clean invalid-input döndürür; denominator
  ! epsilon ile düzeltilmez. Gerçek overflow ayrı numerical-failure status'udur.
  result = integrate_weighted_linear_interval(1.0_dp, 1.0_dp, 2.0_dp, &
    0.0_dp, 1.0_dp)
  call assert_true(.not. result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_INVALID_INPUT, &
    "Zero variance reddedilmedi.")
  result = integrate_weighted_linear_interval(1.0_dp, 1.0_dp, 2.0_dp, &
    -1.0_dp, 1.0_dp)
  call assert_true(.not. result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_INVALID_INPUT, &
    "Negative variance reddedilmedi.")
  result = integrate_weighted_linear_interval(1.0_dp, 1.0_dp, 2.0_dp, &
    ieee_value(0.0_dp, ieee_quiet_nan), 1.0_dp)
  call assert_true(.not. result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_INVALID_INPUT, &
    "NaN variance reddedilmedi.")
  result = integrate_weighted_linear_interval(1.0_dp, 1.0_dp, 2.0_dp, &
    ieee_value(0.0_dp, ieee_positive_inf), 1.0_dp)
  call assert_true(.not. result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_INVALID_INPUT, &
    "Infinite variance reddedilmedi.")
  result = integrate_weighted_linear_interval( &
    huge(1.0_dp), huge(1.0_dp), huge(1.0_dp), tiny(1.0_dp), tiny(1.0_dp))
  call assert_true(.not. result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_NUMERICAL_FAILURE, &
    "Overflow analytical integral numerical failure vermedi.")

  print *, "V0.8.4 analytical weighted interval integralleri doğrulandı."
end program test_tts_weighted_integrals
