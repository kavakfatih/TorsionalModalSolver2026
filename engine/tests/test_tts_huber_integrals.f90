program test_tts_huber_integrals
  use tms_kinds, only : dp
  use tms_weighted_piecewise_integrals, only : &
    huber_interval_integral_t, integrate_huber_linear_interval, &
    WEIGHTED_INTEGRAL_SUCCESS, WEIGHTED_INTEGRAL_INVALID_INPUT
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  type(huber_interval_integral_t) :: result
  real(dp) :: expected
  real(dp) :: variance_perturbation

  ! A) Bütün interval |z|<=c: Huber tam olarak weighted quadratic loss'un
  ! yarısıdır. Truth, sabit variance linear-residual polinom integralidir.
  expected = 0.5_dp*(0.2_dp**2 + 0.2_dp*0.8_dp + 0.8_dp**2)/3.0_dp
  result = integrate_huber_linear_interval( &
    1.0_dp, 0.2_dp, 0.8_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call assert_true(result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_SUCCESS, &
    "Huber quadratic-region integral üretilemedi.")
  call assert_close(result%integral, expected, 2.0e-14_dp, &
    "Huber quadratic-region analytical integral hatalı.")
  call assert_close(result%quadratic_width_decades, 1.0_dp, 1.0e-14_dp, &
    "Quadratic-region width diagnostic hatalı.")

  ! B/C) Sabit positive ve negative tail işaretleri aynı |z| loss'unu verir.
  expected = 1.5_dp
  result = integrate_huber_linear_interval( &
    1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call assert_close(result%integral, expected, 2.0e-14_dp, &
    "Positive Huber tail integral hatalı.")
  call assert_close(result%tail_width_decades, 1.0_dp, 1.0e-14_dp, &
    "Positive tail width diagnostic hatalı.")
  result = integrate_huber_linear_interval( &
    1.0_dp, -2.0_dp, -2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call assert_close(result%integral, expected, 2.0e-14_dp, &
    "Negative Huber tail integral hatalı.")

  ! D) r=t, t=[0,2] bir threshold crossing taşır. Bağımsız parçalar:
  ! 0..1: 0.5*t^2; 1..2: t-0.5. Toplam 1/6+1=7/6.
  result = integrate_huber_linear_interval( &
    2.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call assert_close(result%integral, 7.0_dp/6.0_dp, 3.0e-14_dp, &
    "Tek Huber threshold crossing integral'i hatalı.")
  call assert_close(result%quadratic_width_decades, 1.0_dp, 2.0e-14_dp, &
    "Tek crossing quadratic width hatalı.")
  call assert_close(result%tail_width_decades, 1.0_dp, 2.0e-14_dp, &
    "Tek crossing tail width hatalı.")

  ! E) r=t-2, t=[0,4] iki crossing taşır. İki tail integralinin her biri 1,
  ! merkez quadratic integral 1/3 olduğundan total 7/3'tür.
  result = integrate_huber_linear_interval( &
    4.0_dp, -2.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
  call assert_close(result%integral, 7.0_dp/3.0_dp, 3.0e-14_dp, &
    "İki Huber threshold crossing integral'i hatalı.")
  call assert_close(result%quadratic_width_decades, 2.0_dp, 2.0e-14_dp, &
    "İki crossing quadratic width hatalı.")
  call assert_close(result%tail_width_decades, 2.0_dp, 2.0e-14_dp, &
    "İki crossing tail width hatalı.")

  ! F) Nearly-constant variance dalı d->0 quadratic limitine süreklidir.
  variance_perturbation = 0.5_dp*sqrt(sqrt(epsilon(1.0_dp)))
  result = integrate_huber_linear_interval(1.0_dp, 0.5_dp, 0.5_dp, &
    1.0_dp, 1.0_dp + variance_perturbation, 1.0_dp)
  call assert_true(result%valid, &
    "Nearly-constant Huber variance numerical failure üretti.")
  call assert_close(result%integral, 0.125_dp, 2.0e-5_dp, &
    "Nearly-constant Huber variance limit'i süreksiz.")

  ! G) Nonconstant v=1+t ve positive-tail r=2+0.5t için bağımsız u=v
  ! substitution: integral r/sqrt(v) dt hesaplanıp c|z|-c^2/2 uygulanır.
  expected = (1.0_dp/3.0_dp)*(2.0_dp**1.5_dp - 1.0_dp) + &
    3.0_dp*(sqrt(2.0_dp) - 1.0_dp) - 0.5_dp
  result = integrate_huber_linear_interval( &
    1.0_dp, 2.0_dp, 2.5_dp, 1.0_dp, 2.0_dp, 1.0_dp)
  call assert_close(result%integral, expected, 3.0e-14_dp, &
    "Nonconstant-variance Huber tail integral hatalı.")
  call assert_close(result%tail_width_decades, 1.0_dp, 2.0e-14_dp, &
    "Nonconstant-variance tail width hatalı.")

  ! Near-tangent threshold kökü duplicate crossing üretmemeli. Burada
  ! r=1+q, v=0.75+3q ve r^2-v=(q-0.5)^2; tek tangent dışındaki tüm width tail.
  result = integrate_huber_linear_interval( &
    1.0_dp, 1.0_dp, 2.0_dp, 0.75_dp, 3.75_dp, 1.0_dp)
  call assert_true(result%valid, "Near-tangent Huber root işlenemedi.")
  call assert_close(result%tail_width_decades, 1.0_dp, 2.0e-13_dp, &
    "Near-tangent root tail width'i bozdu.")

  result = integrate_huber_linear_interval( &
    1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp)
  call assert_true(.not. result%valid .and. &
    result%status == WEIGHTED_INTEGRAL_INVALID_INPUT, &
    "Nonpositive Huber c reddedilmedi.")

  print *, "V0.8.4 analytical standardized-Huber integralleri doğrulandı."
end program test_tts_huber_integrals
