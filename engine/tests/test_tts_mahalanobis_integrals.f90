program test_tts_mahalanobis_integrals
  use tms_kinds, only : dp
  use tms_bivariate_covariance_integrals, only : &
    bivariate_covariance_interval_integral_t, &
    integrate_bivariate_covariance_linear_interval, &
    BIVARIATE_INTEGRAL_ILL_CONDITIONED
  use tms_weighted_piecewise_integrals, only : &
    weighted_interval_integral_t, integrate_weighted_linear_interval
  use tms_tts_covariance_types, only : tts_covariance_matrix_2x2_t
  use tms_tts_covariance_validation, only : &
    interpolate_covariance_matrix_entries
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  type :: oracle_case_t
    real(dp) :: sr0
    real(dp) :: sr1
    real(dp) :: lr0
    real(dp) :: lr1
    real(dp) :: sv0
    real(dp) :: sv1
    real(dp) :: lv0
    real(dp) :: lv1
    real(dp) :: c0
    real(dp) :: c1
  end type oracle_case_t

  type(bivariate_covariance_interval_integral_t) :: integral
  type(weighted_interval_integral_t) :: storage_integral
  type(weighted_interval_integral_t) :: loss_integral
  type(tts_covariance_matrix_2x2_t) :: matrix0
  type(tts_covariance_matrix_2x2_t) :: matrix1
  type(tts_covariance_matrix_2x2_t) :: interpolated
  logical :: valid
  integer :: i

  ! A) Sabit covariance limitinde integrand residual lineerliğinden quadratic
  ! olur; production rational primitive independent adaptive oracle ile eşleşir.
  call check_case("constant covariance", 1.7_dp, &
    0.2_dp, 1.1_dp, -0.4_dp, 0.7_dp, &
    2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 0.4_dp, 0.4_dp, 2.0e-11_dp)

  ! B) c=0 diagonal limitinde yeni integral iki mevcut scalar weighted exact
  ! integralinin toplamına indirgenmelidir.
  integral = integrate_bivariate_covariance_linear_interval(1.3_dp, &
    0.1_dp, 0.8_dp, -0.2_dp, 0.5_dp, 1.5_dp, 1.5_dp, &
    2.5_dp, 2.5_dp, 0.0_dp, 0.0_dp)
  storage_integral = integrate_weighted_linear_interval( &
    1.3_dp, 0.1_dp, 0.8_dp, 1.5_dp, 1.5_dp)
  loss_integral = integrate_weighted_linear_interval( &
    1.3_dp, -0.2_dp, 0.5_dp, 2.5_dp, 2.5_dp)
  call assert_true(integral%valid, "Diagonal constant covariance çözülemedi.")
  call assert_close(integral%mahalanobis_integral, &
    storage_integral%integral + loss_integral%integral, 2.0e-13_dp, &
    "c=0 Mahalanobis scalar weighted toplamına indirgenmedi.")
  call assert_close(integral%mahalanobis_integral, &
    integral%matched_diagonal_integral, 2.0e-13_dp, &
    "c=0 matched-diagonal/Mahalanobis invariant bozuldu.")

  ! C-E) Nonconstant matrix ile pozitif ve negatif covariance yönleri,
  ! doğrudan adjugate integrand kullanan independent adaptive oracle'a karşıdır.
  call check_case("nonconstant covariance", 0.9_dp, &
    -0.3_dp, 0.9_dp, 0.8_dp, -0.1_dp, &
    1.2_dp, 2.0_dp, 2.8_dp, 2.1_dp, 0.2_dp, 0.5_dp, 3.0e-10_dp)
  call check_case("positive covariance", 1.1_dp, &
    0.4_dp, 1.0_dp, 0.2_dp, 1.3_dp, &
    2.2_dp, 2.8_dp, 1.9_dp, 2.5_dp, 0.6_dp, 0.7_dp, 3.0e-10_dp)
  call check_case("negative covariance", 1.1_dp, &
    0.4_dp, 1.0_dp, 0.2_dp, 1.3_dp, &
    2.2_dp, 2.8_dp, 1.9_dp, 2.5_dp, -0.6_dp, -0.3_dp, 3.0e-10_dp)

  ! F) q2≈0 nearly-linear determinant ve G) gerçek quadratic denominator
  ! discriminant dalları cancellation'a rağmen oracle ile uyumlu kalmalıdır.
  call check_case("nearly linear determinant", 0.8_dp, &
    -0.5_dp, 0.7_dp, 0.9_dp, -0.2_dp, &
    1.5_dp, 1.9_dp, 2.0_dp, 2.1_dp, 0.25_dp, &
    0.25_dp + 0.2_dp*(1.0_dp - 1.0e-10_dp), 2.0e-6_dp)
  call check_case("quadratic denominator", 1.4_dp, &
    -0.1_dp, 1.2_dp, 0.7_dp, -0.6_dp, &
    1.4_dp, 2.2_dp, 3.0_dp, 2.6_dp, 0.15_dp, 0.25_dp, 5.0e-10_dp)

  ! H) Matrix slope'u machine-small olduğunda constant-covariance limitine
  ! sürekli yaklaşmalı ve arbitrary resampling gerektirmemelidir.
  call check_case("near constant matrix slope", 2.0_dp, &
    0.3_dp, -0.8_dp, 0.5_dp, 1.1_dp, &
    2.0_dp, 2.0_dp*(1.0_dp + 1.0e-12_dp), &
    1.5_dp, 1.5_dp*(1.0_dp - 1.0e-12_dp), &
    0.2_dp, 0.2_dp*(1.0_dp + 1.0e-12_dp), 2.0e-9_dp)

  ! I) SPD cone convexdir: endpoint matris elemanlarının x=log10(f) üzerinde
  ! lineer interpolation'ı seçili alpha noktalarında SPD kalmalıdır.
  matrix0%storage_variance_pa2 = 2.0_dp
  matrix0%loss_variance_pa2 = 1.5_dp
  matrix0%storage_loss_covariance_pa2 = 0.4_dp
  matrix1%storage_variance_pa2 = 4.0_dp
  matrix1%loss_variance_pa2 = 3.0_dp
  matrix1%storage_loss_covariance_pa2 = -0.7_dp
  do i = 0, 10
    call interpolate_covariance_matrix_entries(matrix0, matrix1, &
      real(i, dp)/10.0_dp, interpolated, valid)
    call assert_true(valid, "Convex covariance interpolation SPD'yi bozdu.")
  end do

  ! J) rho≈1 nedeniyle machine-safe reciprocal condition sağlamayan interval
  ! explicit status ile reddedilir; jitter/pseudoinverse uygulanmaz.
  integral = integrate_bivariate_covariance_linear_interval(1.0_dp, &
    1.0_dp, 0.0_dp, 0.5_dp, -0.5_dp, 1.0_dp, 1.0_dp, &
    1.0_dp, 1.0_dp, 1.0_dp - 1.0e-10_dp, 1.0_dp - 1.0e-10_dp)
  call assert_true(.not. integral%valid .and. &
    integral%status == BIVARIATE_INTEGRAL_ILL_CONDITIONED, &
    "Near-singular interval explicit conditioning status vermedi.")

  print *, "V0.8.5 analytical Mahalanobis interval integralleri doğrulandı."

contains

  subroutine check_case(label, h, sr0, sr1, lr0, lr1, &
      sv0, sv1, lv0, lv1, c0, c1, tolerance)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: h
    real(dp), intent(in) :: sr0
    real(dp), intent(in) :: sr1
    real(dp), intent(in) :: lr0
    real(dp), intent(in) :: lr1
    real(dp), intent(in) :: sv0
    real(dp), intent(in) :: sv1
    real(dp), intent(in) :: lv0
    real(dp), intent(in) :: lv1
    real(dp), intent(in) :: c0
    real(dp), intent(in) :: c1
    real(dp), intent(in) :: tolerance
    type(oracle_case_t) :: data
    type(bivariate_covariance_interval_integral_t) :: actual
    real(dp) :: expected

    data = oracle_case_t(sr0, sr1, lr0, lr1, sv0, sv1, &
      lv0, lv1, c0, c1)
    actual = integrate_bivariate_covariance_linear_interval(h, &
      sr0, sr1, lr0, lr1, sv0, sv1, lv0, lv1, c0, c1)
    expected = h*adaptive_integral(0.0_dp, 1.0_dp, data, 1.0e-12_dp, 24)
    call assert_true(actual%valid, trim(label)// &
      ": production analytical integral geçersiz.")
    call assert_close(actual%mahalanobis_integral, expected, tolerance, &
      trim(label)//": independent adaptive oracle uyuşmazlığı.")
  end subroutine check_case

  function adaptive_integral(a, b, data, tolerance, depth) result(value)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    type(oracle_case_t), intent(in) :: data
    real(dp), intent(in) :: tolerance
    integer, intent(in) :: depth
    real(dp) :: value
    real(dp) :: midpoint
    real(dp) :: whole

    midpoint = 0.5_dp*(a + b)
    whole = (b - a) * (oracle_integrand(a, data) + &
      4.0_dp*oracle_integrand(midpoint, data) + &
      oracle_integrand(b, data))/6.0_dp
    value = adaptive_refine(a, b, oracle_integrand(a, data), &
      oracle_integrand(midpoint, data), oracle_integrand(b, data), &
      whole, tolerance, depth, data)
  end function adaptive_integral

  recursive function adaptive_refine( &
      a, b, fa, fm, fb, whole, tolerance, depth, data) result(value)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: fa
    real(dp), intent(in) :: fm
    real(dp), intent(in) :: fb
    real(dp), intent(in) :: whole
    real(dp), intent(in) :: tolerance
    integer, intent(in) :: depth
    type(oracle_case_t), intent(in) :: data
    real(dp) :: value
    real(dp) :: left_midpoint
    real(dp) :: right_midpoint
    real(dp) :: midpoint
    real(dp) :: flm
    real(dp) :: frm
    real(dp) :: left
    real(dp) :: right

    midpoint = 0.5_dp*(a + b)
    left_midpoint = 0.5_dp*(a + midpoint)
    right_midpoint = 0.5_dp*(midpoint + b)
    flm = oracle_integrand(left_midpoint, data)
    frm = oracle_integrand(right_midpoint, data)
    left = (midpoint - a)*(fa + 4.0_dp*flm + fm)/6.0_dp
    right = (b - midpoint)*(fm + 4.0_dp*frm + fb)/6.0_dp
    if (depth <= 0 .or. abs(left + right - whole) <= &
        15.0_dp*tolerance) then
      value = left + right + (left + right - whole)/15.0_dp
    else
      value = adaptive_refine(a, midpoint, fa, flm, fm, left, &
        0.5_dp*tolerance, depth - 1, data) + &
        adaptive_refine(midpoint, b, fm, frm, fb, right, &
        0.5_dp*tolerance, depth - 1, data)
    end if
  end function adaptive_refine

  function oracle_integrand(q, data) result(value)
    real(dp), intent(in) :: q
    type(oracle_case_t), intent(in) :: data
    real(dp) :: value
    real(dp) :: c
    real(dp) :: determinant
    real(dp) :: lr
    real(dp) :: lv
    real(dp) :: sr
    real(dp) :: sv

    sr = data%sr0 + q*(data%sr1 - data%sr0)
    lr = data%lr0 + q*(data%lr1 - data%lr0)
    sv = data%sv0 + q*(data%sv1 - data%sv0)
    lv = data%lv0 + q*(data%lv1 - data%lv0)
    c = data%c0 + q*(data%c1 - data%c0)
    determinant = sv*lv - c*c
    value = (lv*sr*sr - 2.0_dp*c*sr*lr + sv*lr*lr)/determinant
  end function oracle_integrand

end program test_tts_mahalanobis_integrals
