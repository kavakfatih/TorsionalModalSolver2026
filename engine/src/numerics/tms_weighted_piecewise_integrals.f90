module tms_weighted_piecewise_integrals
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: WEIGHTED_INTEGRAL_SUCCESS = 0
  integer, parameter, public :: WEIGHTED_INTEGRAL_INVALID_INPUT = 1
  integer, parameter, public :: WEIGHTED_INTEGRAL_NUMERICAL_FAILURE = 2

  !> Lineer residual ve lineer variance üzerindeki analitik integral sonucudur.
  !! Integral boyutsuz standardized-residual karesinin log10-frequency
  !! genişliğiyle çarpımıdır; width birimi decade'dir.
  type, public :: weighted_interval_integral_t
    integer :: status = WEIGHTED_INTEGRAL_INVALID_INPUT
    logical :: valid = .false.
    real(dp) :: integral = 0.0_dp
  end type weighted_interval_integral_t

  !> Standardized Huber interval integralini ve regime genişliklerini taşır.
  !! Genişlikler log10-frequency decade, integral rho_c(z)*decade birimindedir.
  type, public :: huber_interval_integral_t
    integer :: status = WEIGHTED_INTEGRAL_INVALID_INPUT
    logical :: valid = .false.
    real(dp) :: integral = 0.0_dp
    real(dp) :: total_width_decades = 0.0_dp
    real(dp) :: quadratic_width_decades = 0.0_dp
    real(dp) :: tail_width_decades = 0.0_dp
  end type huber_interval_integral_t

  public :: integrate_weighted_linear_interval
  public :: integrate_huber_linear_interval

contains

  !> r(t)=r0+(r1-r0)t/h ve v(t)=v0+(v1-v0)t/h modeli için
  !! integral_0^h r(t)^2/v(t) dt değerini analitik hesaplar. Residual
  !! log10-modulus, variance log10-modulus varyansı ve h decade'dir.
  !! Varsayım: h>0 ve iki variance endpoint'i sonlu/pozitiftir.
  pure function integrate_weighted_linear_interval( &
      h, r0, r1, v0, v1) result(result)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: r0
    real(dp), intent(in) :: r1
    real(dp), intent(in) :: v0
    real(dp), intent(in) :: v1
    type(weighted_interval_integral_t) :: result

    real(dp) :: a0
    real(dp) :: a1
    real(dp) :: a2
    real(dp) :: base_variance
    real(dp) :: end_variance
    real(dp) :: log_ratio
    real(dp) :: normalized_delta
    real(dp) :: normalized_r0
    real(dp) :: normalized_r_delta
    real(dp) :: residual_scale
    real(dp) :: start_residual
    real(dp) :: end_residual
    real(dp) :: scaled_integral

    if (.not. interval_inputs_are_valid(h, r0, r1, v0, v1)) return

    ! Büyük variance endpoint'inden başlanması p=(v_end/v_base)-1 değerini
    ! [-1,0] aralığında tutar; çok büyük variance oranında overflow önlenir.
    if (v0 >= v1) then
      base_variance = v0
      end_variance = v1
      start_residual = r0
      end_residual = r1
    else
      base_variance = v1
      end_variance = v0
      start_residual = r1
      end_residual = r0
    end if
    normalized_delta = end_variance/base_variance - 1.0_dp
    log_ratio = log(end_variance) - log(base_variance)
    a0 = inverse_linear_variance_moment(0, normalized_delta, log_ratio)
    a1 = inverse_linear_variance_moment(1, normalized_delta, log_ratio)
    a2 = inverse_linear_variance_moment(2, normalized_delta, log_ratio)

    residual_scale = max(abs(start_residual), abs(end_residual))
    if (residual_scale <= 0.0_dp) then
      result%status = WEIGHTED_INTEGRAL_SUCCESS
      result%valid = .true.
      return
    end if
    normalized_r0 = start_residual/residual_scale
    normalized_r_delta = (end_residual - start_residual)/residual_scale
    scaled_integral = normalized_r0*normalized_r0*a0 + &
      2.0_dp*normalized_r0*normalized_r_delta*a1 + &
      normalized_r_delta*normalized_r_delta*a2
    result%integral = h*(residual_scale/sqrt(base_variance))**2 * &
      scaled_integral
    if (.not. ieee_is_finite(result%integral)) then
      result%status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
      result%integral = 0.0_dp
      return
    end if
    if (result%integral < 0.0_dp) then
      if (abs(result%integral) <= 256.0_dp*epsilon(1.0_dp) * &
          max(1.0_dp, abs(h))) then
        result%integral = 0.0_dp
      else
        result%status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
        result%integral = 0.0_dp
        return
      end if
    end if
    result%status = WEIGHTED_INTEGRAL_SUCCESS
    result%valid = .true.
  end function integrate_weighted_linear_interval

  !> Standardized residual z=r/sqrt(v) için Huber rho_c(z) integralini
  !! grid kullanmadan hesaplar. r^2-c^2*v=0 crossing'leri kararlı quadratic
  !! köklerle bulunur; quadratic ve tail bölgeleri ayrı analitik entegre edilir.
  !! c boyutsuz tuning constant'tır, kabul limiti değildir.
  pure function integrate_huber_linear_interval( &
      h, r0, r1, v0, v1, huber_c) result(result)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: r0
    real(dp), intent(in) :: r1
    real(dp), intent(in) :: v0
    real(dp), intent(in) :: v1
    real(dp), intent(in) :: huber_c
    type(huber_interval_integral_t) :: result

    real(dp) :: knots(4)
    real(dp) :: q0
    real(dp) :: q1
    real(dp) :: qm
    real(dp) :: local_h
    real(dp) :: local_r0
    real(dp) :: local_r1
    real(dp) :: local_v0
    real(dp) :: local_v1
    real(dp) :: midpoint_r
    real(dp) :: midpoint_v
    real(dp) :: signed_integral
    real(dp) :: local_value
    integer :: i
    integer :: knot_count
    type(weighted_interval_integral_t) :: weighted

    if (.not. interval_inputs_are_valid(h, r0, r1, v0, v1)) return
    if (.not. ieee_is_finite(huber_c) .or. huber_c <= 0.0_dp) return
    call build_huber_regime_knots(r0, r1, v0, v1, huber_c, &
      knots, knot_count, result%status)
    if (result%status /= WEIGHTED_INTEGRAL_SUCCESS) return

    result%total_width_decades = h
    do i = 1, knot_count - 1
      q0 = knots(i)
      q1 = knots(i + 1)
      if (q1 <= q0) cycle
      qm = 0.5_dp*(q0 + q1)
      local_h = h*(q1 - q0)
      local_r0 = linear_value(r0, r1, q0)
      local_r1 = linear_value(r0, r1, q1)
      local_v0 = linear_value(v0, v1, q0)
      local_v1 = linear_value(v0, v1, q1)
      midpoint_r = linear_value(r0, r1, qm)
      midpoint_v = linear_value(v0, v1, qm)
      if (abs(midpoint_r) <= huber_c*sqrt(midpoint_v)) then
        weighted = integrate_weighted_linear_interval( &
          local_h, local_r0, local_r1, local_v0, local_v1)
        if (.not. weighted%valid) then
          result%status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
          return
        end if
        local_value = 0.5_dp*weighted%integral
        result%quadratic_width_decades = &
          result%quadratic_width_decades + local_h
      else
        signed_integral = integrate_linear_over_sqrt_variance( &
          local_h, local_r0, local_r1, local_v0, local_v1)
        if (.not. ieee_is_finite(signed_integral)) then
          result%status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
          return
        end if
        local_value = huber_c*abs(signed_integral) - &
          0.5_dp*huber_c*huber_c*local_h
        result%tail_width_decades = result%tail_width_decades + local_h
      end if
      if (.not. ieee_is_finite(local_value)) then
        result%status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
        return
      end if
      if (local_value < 0.0_dp .and. &
          abs(local_value) <= 256.0_dp*epsilon(1.0_dp) * &
            max(1.0_dp, local_h*huber_c*huber_c)) local_value = 0.0_dp
      if (local_value < 0.0_dp) then
        result%status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
        return
      end if
      result%integral = result%integral + local_value
    end do
    if (.not. ieee_is_finite(result%integral)) then
      result%status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
      return
    end if
    result%status = WEIGHTED_INTEGRAL_SUCCESS
    result%valid = .true.
  end function integrate_huber_linear_interval

  pure function interval_inputs_are_valid(h, r0, r1, v0, v1) result(valid)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: r0
    real(dp), intent(in) :: r1
    real(dp), intent(in) :: v0
    real(dp), intent(in) :: v1
    logical :: valid

    valid = ieee_is_finite(h) .and. h > 0.0_dp .and. &
      ieee_is_finite(r0) .and. ieee_is_finite(r1) .and. &
      ieee_is_finite(v0) .and. ieee_is_finite(v1) .and. &
      v0 > 0.0_dp .and. v1 > 0.0_dp
  end function interval_inputs_are_valid

  !> A_n=integral_0^1 q^n/(1+p*q)dq momentini hesaplar. Küçük p için
  !! alternating power series, diğer durumda log-ratio recurrence kullanılır;
  !! eşik yalnız makine hassasiyeti ölçeğinden türetilir.
  pure function inverse_linear_variance_moment( &
      order, p, log_ratio) result(moment)
    integer, intent(in) :: order
    real(dp), intent(in) :: p
    real(dp), intent(in) :: log_ratio
    real(dp) :: moment

    real(dp) :: power
    real(dp) :: term
    real(dp) :: a0
    real(dp) :: a1
    integer :: k

    if (abs(p) <= sqrt(sqrt(epsilon(1.0_dp)))) then
      moment = 0.0_dp
      power = 1.0_dp
      do k = 0, 64
        term = power/real(order + k + 1, dp)
        moment = moment + term
        if (abs(term) <= epsilon(1.0_dp)*max(1.0_dp, abs(moment))) exit
        power = -power*p
      end do
      return
    end if

    a0 = log_ratio/p
    if (order == 0) then
      moment = a0
      return
    end if
    a1 = (1.0_dp - a0)/p
    if (order == 1) then
      moment = a1
    else
      moment = (0.5_dp - a1)/p
    end if
  end function inverse_linear_variance_moment

  !> Lineer r'nin sqrt(lineer v) ile oranını exact entegre eder. Büyük
  !! variance endpoint'i base seçilir; B0 ve B1 momentlerinin rationalized
  !! biçimi d->0 limitinde de süreklidir.
  pure function integrate_linear_over_sqrt_variance( &
      h, r0, r1, v0, v1) result(value)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: r0
    real(dp), intent(in) :: r1
    real(dp), intent(in) :: v0
    real(dp), intent(in) :: v1
    real(dp) :: value

    real(dp) :: b0
    real(dp) :: b1
    real(dp) :: base_variance
    real(dp) :: end_variance
    real(dp) :: residual_delta
    real(dp) :: start_residual
    real(dp) :: end_residual
    real(dp) :: square_root_ratio

    if (v0 >= v1) then
      base_variance = v0
      end_variance = v1
      start_residual = r0
      end_residual = r1
    else
      base_variance = v1
      end_variance = v0
      start_residual = r1
      end_residual = r0
    end if
    square_root_ratio = sqrt(end_variance/base_variance)
    b0 = 2.0_dp/(square_root_ratio + 1.0_dp)
    b1 = 2.0_dp*(square_root_ratio + 2.0_dp) / &
      (3.0_dp*(square_root_ratio + 1.0_dp)**2)
    residual_delta = end_residual - start_residual
    value = h*(start_residual*b0 + residual_delta*b1) / &
      sqrt(base_variance)
  end function integrate_linear_over_sqrt_variance

  pure subroutine build_huber_regime_knots( &
      r0, r1, v0, v1, huber_c, knots, knot_count, status)
    real(dp), intent(in) :: r0
    real(dp), intent(in) :: r1
    real(dp), intent(in) :: v0
    real(dp), intent(in) :: v1
    real(dp), intent(in) :: huber_c
    real(dp), intent(out) :: knots(4)
    integer, intent(out) :: knot_count
    integer, intent(out) :: status

    real(dp) :: a
    real(dp) :: b
    real(dp) :: c0
    real(dp) :: coefficient_scale
    real(dp) :: discriminant
    real(dp) :: discriminant_scale
    real(dp) :: dr
    real(dp) :: q_first
    real(dp) :: q_second
    real(dp) :: root_scale
    real(dp) :: scale
    real(dp) :: square_root_discriminant
    real(dp) :: threshold0
    real(dp) :: threshold1

    knots = 0.0_dp
    knots(1) = 0.0_dp
    knots(2) = 1.0_dp
    knot_count = 2
    status = WEIGHTED_INTEGRAL_SUCCESS

    threshold0 = huber_c*sqrt(v0)
    threshold1 = huber_c*sqrt(v1)
    if (.not. ieee_is_finite(threshold0) .or. &
        .not. ieee_is_finite(threshold1)) then
      status = WEIGHTED_INTEGRAL_NUMERICAL_FAILURE
      return
    end if
    scale = max(abs(r0), abs(r1), threshold0, threshold1)
    if (scale <= 0.0_dp) return
    dr = (r1 - r0)/scale
    threshold0 = (threshold0/scale)**2
    threshold1 = (threshold1/scale)**2
    c0 = (r0/scale)**2 - threshold0
    b = 2.0_dp*(r0/scale)*dr - (threshold1 - threshold0)
    a = dr*dr
    coefficient_scale = max(abs(a), abs(b), abs(c0), tiny(1.0_dp))
    root_scale = 128.0_dp*epsilon(1.0_dp)*coefficient_scale

    if (abs(a) <= root_scale) then
      if (abs(b) > root_scale) then
        q_first = -c0/b
        call append_interior_root(knots, knot_count, q_first)
      end if
    else
      discriminant = b*b - 4.0_dp*a*c0
      discriminant_scale = max(b*b, abs(4.0_dp*a*c0), tiny(1.0_dp))
      if (abs(discriminant) <= 128.0_dp*epsilon(1.0_dp) * &
          discriminant_scale) discriminant = 0.0_dp
      if (discriminant >= 0.0_dp) then
        square_root_discriminant = sqrt(discriminant)
        q_first = -0.5_dp*(b + sign(square_root_discriminant, b))/a
        call append_interior_root(knots, knot_count, q_first)
        if (abs(q_first) > root_scale) then
          q_second = c0/(a*q_first)
        else
          q_second = (-b + square_root_discriminant)/(2.0_dp*a)
        end if
        call append_interior_root(knots, knot_count, q_second)
      end if
    end if
    call sort_small_knot_array(knots, knot_count)
  end subroutine build_huber_regime_knots

  pure subroutine append_interior_root(knots, knot_count, root)
    real(dp), intent(inout) :: knots(4)
    integer, intent(inout) :: knot_count
    real(dp), intent(in) :: root
    integer :: i

    if (.not. ieee_is_finite(root)) return
    if (root <= 0.0_dp .or. root >= 1.0_dp) return
    if (root_is_endpoint_equivalent(root, 0.0_dp) .or. &
        root_is_endpoint_equivalent(root, 1.0_dp)) return
    do i = 1, knot_count
      if (root_is_endpoint_equivalent(root, knots(i))) return
    end do
    if (knot_count < size(knots)) then
      knot_count = knot_count + 1
      knots(knot_count) = root
    end if
  end subroutine append_interior_root

  pure subroutine sort_small_knot_array(knots, knot_count)
    real(dp), intent(inout) :: knots(4)
    integer, intent(in) :: knot_count
    real(dp) :: key
    integer :: i
    integer :: j

    do i = 2, knot_count
      key = knots(i)
      j = i - 1
      do while (j >= 1)
        if (knots(j) <= key) exit
        knots(j + 1) = knots(j)
        j = j - 1
      end do
      knots(j + 1) = key
    end do
  end subroutine sort_small_knot_array

  pure function root_is_endpoint_equivalent(a, b) result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    equivalent = abs(a - b) <= 128.0_dp*epsilon(1.0_dp) * &
      max(1.0_dp, abs(a), abs(b))
  end function root_is_endpoint_equivalent

  pure function linear_value(value0, value1, coordinate) result(value)
    real(dp), intent(in) :: value0
    real(dp), intent(in) :: value1
    real(dp), intent(in) :: coordinate
    real(dp) :: value

    value = value0 + coordinate*(value1 - value0)
  end function linear_value

end module tms_weighted_piecewise_integrals
