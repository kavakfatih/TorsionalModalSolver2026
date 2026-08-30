module tms_bivariate_covariance_integrals
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_weighted_piecewise_integrals, only : &
    weighted_interval_integral_t, integrate_weighted_linear_interval
  implicit none
  private

  integer, parameter, public :: BIVARIATE_INTEGRAL_SUCCESS = 0
  integer, parameter, public :: BIVARIATE_INTEGRAL_INVALID_INPUT = 1
  integer, parameter, public :: BIVARIATE_INTEGRAL_INVALID_COVARIANCE = 2
  integer, parameter, public :: BIVARIATE_INTEGRAL_ILL_CONDITIONED = 3
  integer, parameter, public :: BIVARIATE_INTEGRAL_NUMERICAL_FAILURE = 4

  !> Lineer residual ve entrywise-lineer 2x2 covariance üzerindeki exact
  !! interval integralidir. İki integral de standardized residual karesi ile
  !! log10-frequency decade genişliğinin çarpımıdır.
  type, public :: bivariate_covariance_interval_integral_t
    integer :: status = BIVARIATE_INTEGRAL_INVALID_INPUT
    logical :: valid = .false.
    real(dp) :: mahalanobis_integral = 0.0_dp
    real(dp) :: matched_diagonal_integral = 0.0_dp
  end type bivariate_covariance_interval_integral_t

  public :: integrate_bivariate_covariance_linear_interval

contains

  !> r_s(t), r_l(t) ve Sigma(t) elemanlarının t in [0,h] üzerinde lineer
  !! olduğu modelde Mahalanobis integrandını grid kullanmadan entegre eder.
  !! P3/Q2 polynomial division ve discriminant'a göre log/atan primitives
  !! kullanılır. h [decade], residual log10-modulus, covariance boyutsuzdur.
  !! Endpoint covariance'ları SPD ve machine-safe condition'da olmalıdır.
  pure function integrate_bivariate_covariance_linear_interval( &
      h, storage_r0, storage_r1, loss_r0, loss_r1, &
      storage_v0, storage_v1, loss_v0, loss_v1, &
      covariance0, covariance1) result(result)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: storage_r0
    real(dp), intent(in) :: storage_r1
    real(dp), intent(in) :: loss_r0
    real(dp), intent(in) :: loss_r1
    real(dp), intent(in) :: storage_v0
    real(dp), intent(in) :: storage_v1
    real(dp), intent(in) :: loss_v0
    real(dp), intent(in) :: loss_v1
    real(dp), intent(in) :: covariance0
    real(dp), intent(in) :: covariance1
    type(bivariate_covariance_interval_integral_t) :: result

    type(weighted_interval_integral_t) :: storage_diagonal
    type(weighted_interval_integral_t) :: loss_diagonal
    real(dp) :: a0
    real(dp) :: a1
    real(dp) :: b0
    real(dp) :: b1
    real(dp) :: c0
    real(dp) :: c1
    real(dp) :: da
    real(dp) :: db
    real(dp) :: dc
    real(dp) :: determinant_coefficients(0:2)
    real(dp) :: loss_delta
    real(dp) :: loss_start
    real(dp) :: numerator_coefficients(0:3)
    real(dp) :: rational_integral
    real(dp) :: residual_scale
    real(dp) :: storage_delta
    real(dp) :: storage_start
    real(dp) :: covariance_scale
    logical :: covariance_valid
    logical :: covariance_well_conditioned

    if (.not. inputs_are_finite(h, storage_r0, storage_r1, &
        loss_r0, loss_r1, storage_v0, storage_v1, loss_v0, loss_v1, &
        covariance0, covariance1) .or. h <= 0.0_dp) return
    call interval_covariance_is_safe(storage_v0, storage_v1, loss_v0, &
      loss_v1, covariance0, covariance1, covariance_valid, &
      covariance_well_conditioned)
    if (.not. covariance_valid) then
      result%status = BIVARIATE_INTEGRAL_INVALID_COVARIANCE
      return
    end if
    if (.not. covariance_well_conditioned) then
      result%status = BIVARIATE_INTEGRAL_ILL_CONDITIONED
      return
    end if

    storage_diagonal = integrate_weighted_linear_interval( &
      h, storage_r0, storage_r1, storage_v0, storage_v1)
    loss_diagonal = integrate_weighted_linear_interval( &
      h, loss_r0, loss_r1, loss_v0, loss_v1)
    if (.not. storage_diagonal%valid .or. .not. loss_diagonal%valid) then
      result%status = BIVARIATE_INTEGRAL_NUMERICAL_FAILURE
      return
    end if
    result%matched_diagonal_integral = storage_diagonal%integral + &
      loss_diagonal%integral

    covariance_scale = max(storage_v0, storage_v1, loss_v0, loss_v1, &
      abs(covariance0), abs(covariance1))
    residual_scale = max(abs(storage_r0), abs(storage_r1), &
      abs(loss_r0), abs(loss_r1))
    if (residual_scale <= 0.0_dp) then
      result%status = BIVARIATE_INTEGRAL_SUCCESS
      result%valid = .true.
      return
    end if

    a0 = storage_v0/covariance_scale
    a1 = storage_v1/covariance_scale
    b0 = loss_v0/covariance_scale
    b1 = loss_v1/covariance_scale
    c0 = covariance0/covariance_scale
    c1 = covariance1/covariance_scale
    da = a1 - a0
    db = b1 - b0
    dc = c1 - c0
    storage_start = storage_r0/residual_scale
    storage_delta = (storage_r1 - storage_r0)/residual_scale
    loss_start = loss_r0/residual_scale
    loss_delta = (loss_r1 - loss_r0)/residual_scale

    call build_rational_coefficients(a0, da, b0, db, c0, dc, &
      storage_start, storage_delta, loss_start, loss_delta, &
      numerator_coefficients, determinant_coefficients)
    rational_integral = integrate_cubic_over_quadratic( &
      numerator_coefficients, determinant_coefficients)
    if (.not. ieee_is_finite(rational_integral)) then
      result%status = BIVARIATE_INTEGRAL_NUMERICAL_FAILURE
      return
    end if
    result%mahalanobis_integral = h * &
      (residual_scale/sqrt(covariance_scale))**2*rational_integral
    if (.not. ieee_is_finite(result%mahalanobis_integral)) then
      result%status = BIVARIATE_INTEGRAL_NUMERICAL_FAILURE
      result%mahalanobis_integral = 0.0_dp
      return
    end if
    if (result%mahalanobis_integral < 0.0_dp .and. &
        abs(result%mahalanobis_integral) <= &
          512.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
            result%matched_diagonal_integral)) then
      result%mahalanobis_integral = 0.0_dp
    end if
    if (result%mahalanobis_integral < 0.0_dp) then
      result%status = BIVARIATE_INTEGRAL_NUMERICAL_FAILURE
      result%mahalanobis_integral = 0.0_dp
      return
    end if
    result%status = BIVARIATE_INTEGRAL_SUCCESS
    result%valid = .true.
  end function integrate_bivariate_covariance_linear_interval

  pure subroutine build_rational_coefficients( &
      a0, da, b0, db, c0, dc, s0, ds, l0, dl, p, q)
    real(dp), intent(in) :: a0
    real(dp), intent(in) :: da
    real(dp), intent(in) :: b0
    real(dp), intent(in) :: db
    real(dp), intent(in) :: c0
    real(dp), intent(in) :: dc
    real(dp), intent(in) :: s0
    real(dp), intent(in) :: ds
    real(dp), intent(in) :: l0
    real(dp), intent(in) :: dl
    real(dp), intent(out) :: p(0:3)
    real(dp), intent(out) :: q(0:2)

    real(dp) :: s2(0:2)
    real(dp) :: l2(0:2)
    real(dp) :: sl(0:2)

    s2 = [s0*s0, 2.0_dp*s0*ds, ds*ds]
    l2 = [l0*l0, 2.0_dp*l0*dl, dl*dl]
    sl = [s0*l0, s0*dl + ds*l0, ds*dl]
    p(0) = b0*s2(0) - 2.0_dp*c0*sl(0) + a0*l2(0)
    p(1) = b0*s2(1) + db*s2(0) - &
      2.0_dp*(c0*sl(1) + dc*sl(0)) + a0*l2(1) + da*l2(0)
    p(2) = b0*s2(2) + db*s2(1) - &
      2.0_dp*(c0*sl(2) + dc*sl(1)) + a0*l2(2) + da*l2(1)
    p(3) = db*s2(2) - 2.0_dp*dc*sl(2) + da*l2(2)
    q(0) = a0*b0 - c0*c0
    q(1) = a0*db + da*b0 - 2.0_dp*c0*dc
    q(2) = da*db - dc*dc
  end subroutine build_rational_coefficients

  !> P3(q)/Q2(q), q in [0,1] integralini polynomial division ile exact
  !! değerlendirir. Nearly-linear Q2 dalı yalnız quadratic coefficient
  !! makine hassasiyetinde ihmal edilebilir olduğunda kullanılır.
  pure function integrate_cubic_over_quadratic(p, q) result(value)
    real(dp), intent(in) :: p(0:3)
    real(dp), intent(in) :: q(0:2)
    real(dp) :: value

    real(dp) :: alpha
    real(dp) :: beta
    real(dp) :: denominator_integral
    real(dp) :: denominator_scale
    real(dp) :: linear_coefficient
    real(dp) :: logarithmic_change
    real(dp) :: remainder0
    real(dp) :: remainder1

    denominator_scale = max(abs(q(0)), abs(q(1)), abs(q(2)), tiny(1.0_dp))
    if (abs(q(2)) <= sqrt(epsilon(1.0_dp))*denominator_scale) then
      value = integrate_cubic_over_linear(p, q(0), q(1))
      return
    end if

    linear_coefficient = p(3)/q(2)
    beta = (p(2) - linear_coefficient*q(1))/q(2)
    remainder1 = p(1) - linear_coefficient*q(0) - beta*q(1)
    remainder0 = p(0) - beta*q(0)
    alpha = remainder1/(2.0_dp*q(2))
    beta = remainder0 - alpha*q(1)
    logarithmic_change = log(q(0) + q(1) + q(2)) - log(q(0))
    denominator_integral = integrate_inverse_quadratic(q)
    value = 0.5_dp*linear_coefficient + &
      (p(2) - linear_coefficient*q(1))/q(2) + &
      alpha*logarithmic_change + beta*denominator_integral
  end function integrate_cubic_over_quadratic

  !> ∫0^1 1/(a q^2+b q+c)dq integralini discriminant'a göre atan, log
  !! veya repeated-root primitive ile değerlendirir. Q interval boyunca
  !! pozitifliği çağıran yordamda doğrulanmıştır.
  pure function integrate_inverse_quadratic(q) result(value)
    real(dp), intent(in) :: q(0:2)
    real(dp) :: value

    real(dp) :: angle
    real(dp) :: discriminant
    real(dp) :: discriminant_scale
    real(dp) :: root_discriminant
    real(dp) :: z0
    real(dp) :: z1

    discriminant = 4.0_dp*q(2)*q(0) - q(1)*q(1)
    discriminant_scale = max(abs(4.0_dp*q(2)*q(0)), &
      q(1)*q(1), tiny(1.0_dp))
    if (discriminant > &
        128.0_dp*epsilon(1.0_dp)*discriminant_scale) then
      root_discriminant = sqrt(discriminant)
      angle = atan2(2.0_dp*q(2)*root_discriminant, &
        discriminant + q(1)*(q(1) + 2.0_dp*q(2)))
      value = 2.0_dp*angle/root_discriminant
    else if (discriminant < &
        -128.0_dp*epsilon(1.0_dp)*discriminant_scale) then
      root_discriminant = sqrt(-discriminant)
      z0 = q(1)
      z1 = 2.0_dp*q(2) + q(1)
      value = (log(abs((z1 - root_discriminant) / &
        (z1 + root_discriminant))) - &
        log(abs((z0 - root_discriminant) / &
        (z0 + root_discriminant))))/root_discriminant
    else
      value = 2.0_dp/q(1) - 2.0_dp/(2.0_dp*q(2) + q(1))
    end if
  end function integrate_inverse_quadratic

  !> Cubic polynomial / positive linear denominator integralini inverse
  !! linear moments ile analitik hesaplar. Küçük relative slope için
  !! convergent power series cancellation'ı önler; bu bir quadrature değildir.
  pure function integrate_cubic_over_linear(p, q0, q1) result(value)
    real(dp), intent(in) :: p(0:3)
    real(dp), intent(in) :: q0
    real(dp), intent(in) :: q1
    real(dp) :: value

    real(dp) :: moments(0:3)
    real(dp) :: relative_slope
    real(dp) :: power
    real(dp) :: term
    integer :: i
    integer :: k

    relative_slope = q1/q0
    if (abs(relative_slope) < 0.5_dp) then
      do i = 0, 3
        moments(i) = 0.0_dp
        power = 1.0_dp
        do k = 0, 512
          term = power/real(i + k + 1, dp)
          moments(i) = moments(i) + term
          if (abs(term) <= epsilon(1.0_dp) * &
              max(1.0_dp, abs(moments(i)))) exit
          power = -power*relative_slope
        end do
      end do
    else
      moments(0) = (log(q0 + q1) - log(q0))/relative_slope
      do i = 1, 3
        moments(i) = (1.0_dp/real(i, dp) - moments(i - 1)) / &
          relative_slope
      end do
    end if
    value = sum(p*moments)/q0
  end function integrate_cubic_over_linear

  pure subroutine interval_covariance_is_safe( &
      a0, a1, b0, b1, c0, c1, valid, well_conditioned)
    real(dp), intent(in) :: a0
    real(dp), intent(in) :: a1
    real(dp), intent(in) :: b0
    real(dp), intent(in) :: b1
    real(dp), intent(in) :: c0
    real(dp), intent(in) :: c1
    logical, intent(out) :: valid
    logical, intent(out) :: well_conditioned

    logical :: point_valid
    logical :: point_well_conditioned

    call covariance_point_is_safe(a0, b0, c0, valid, well_conditioned)
    if (.not. valid .or. .not. well_conditioned) return
    call covariance_point_is_safe(a1, b1, c1, point_valid, &
      point_well_conditioned)
    valid = valid .and. point_valid
    well_conditioned = well_conditioned .and. point_well_conditioned
    if (.not. valid .or. .not. well_conditioned) return
    call covariance_point_is_safe(0.5_dp*(a0 + a1), &
      0.5_dp*(b0 + b1), 0.5_dp*(c0 + c1), point_valid, &
      point_well_conditioned)
    valid = valid .and. point_valid
    well_conditioned = well_conditioned .and. point_well_conditioned
  end subroutine interval_covariance_is_safe

  pure subroutine covariance_point_is_safe(a, b, c, valid, well_conditioned)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: c
    logical, intent(out) :: valid
    logical, intent(out) :: well_conditioned

    real(dp) :: an
    real(dp) :: bn
    real(dp) :: cn
    real(dp) :: determinant
    real(dp) :: eigenvalue_maximum
    real(dp) :: eigenvalue_minimum
    real(dp) :: scale
    real(dp) :: tolerance

    valid = .false.
    well_conditioned = .false.
    if (a <= 0.0_dp .or. b <= 0.0_dp) return
    scale = max(a, b, abs(c))
    an = a/scale
    bn = b/scale
    cn = c/scale
    determinant = an*bn - cn*cn
    tolerance = 32.0_dp*epsilon(1.0_dp)* &
      max(an*bn, cn*cn, tiny(1.0_dp))
    if (determinant <= tolerance) return
    valid = .true.
    eigenvalue_maximum = 0.5_dp*(an + bn + &
      sqrt((an - bn)*(an - bn) + 4.0_dp*cn*cn))
    eigenvalue_minimum = determinant/eigenvalue_maximum
    well_conditioned = eigenvalue_minimum/eigenvalue_maximum > &
      sqrt(epsilon(1.0_dp))
  end subroutine covariance_point_is_safe

  pure function inputs_are_finite( &
      h, sr0, sr1, lr0, lr1, sv0, sv1, lv0, lv1, c0, c1) result(valid)
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
    logical :: valid

    valid = ieee_is_finite(h) .and. ieee_is_finite(sr0) .and. &
      ieee_is_finite(sr1) .and. ieee_is_finite(lr0) .and. &
      ieee_is_finite(lr1) .and. ieee_is_finite(sv0) .and. &
      ieee_is_finite(sv1) .and. ieee_is_finite(lv0) .and. &
      ieee_is_finite(lv1) .and. ieee_is_finite(c0) .and. &
      ieee_is_finite(c1)
  end function inputs_are_finite

end module tms_bivariate_covariance_integrals
