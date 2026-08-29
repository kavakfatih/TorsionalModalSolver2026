module tms_scalar_minimizer
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: SCALAR_MINIMIZER_SUCCESS = 0
  integer, parameter, public :: SCALAR_MINIMIZER_INVALID_INPUT = 1
  integer, parameter, public :: SCALAR_MINIMIZER_INVALID_BRACKET = 2
  integer, parameter, public :: SCALAR_MINIMIZER_NONFINITE_OBJECTIVE = 3
  integer, parameter, public :: SCALAR_MINIMIZER_MAX_ITERATIONS = 4

  real(dp), parameter, public :: DEFAULT_SCALAR_MINIMIZER_TOLERANCE = &
    8.0_dp*sqrt(epsilon(1.0_dp))
  integer, parameter, public :: DEFAULT_SCALAR_MINIMIZER_MAX_ITERATIONS = 200

  !> Tek değişkenli minimizasyon sonucunu ve sayısal izini taşır. x boyutsuz
  !! veya caller tarafından tanımlanan birimdedir; f objective birimindedir.
  !! Bounds, Brent'e verilen gerçek üç noktalı interior bracket'ın dış uçlarıdır.
  type, public :: scalar_minimizer_result_t
    integer :: status = SCALAR_MINIMIZER_INVALID_INPUT
    logical :: converged = .false.
    real(dp) :: x_minimum = 0.0_dp
    real(dp) :: f_minimum = huge(1.0_dp)
    real(dp) :: lower_bound = 0.0_dp
    real(dp) :: upper_bound = 0.0_dp
    integer :: iteration_count = 0
    integer :: function_evaluation_count = 0
  end type scalar_minimizer_result_t

  abstract interface
    !> Yan etkisiz scalar objective arayüzüdür. Girdi x ve çıktı f sonlu
    !! olmalıdır; fiziksel anlam ve birimler caller sözleşmesine aittir.
    pure function scalar_objective_interface(x) result(value)
      import :: dp
      real(dp), intent(in) :: x
      real(dp) :: value
    end function scalar_objective_interface
  end interface

  public :: minimize_scalar_brent
  public :: is_valid_scalar_minimum_bracket

contains

  !> Üç sonlu noktada a<b<c ve f(a)>f(b)<f(c) koşulunu doğrular. Bu koşul
  !! interior minimum'un gerçek bir Brent bracket'ı içinde olduğunu gösterir;
  !! boundary minimum veya iki noktalı interval başarı sayılmaz.
  pure function is_valid_scalar_minimum_bracket( &
      lower_x, middle_x, upper_x, lower_f, middle_f, upper_f) result(valid)
    real(dp), intent(in) :: lower_x
    real(dp), intent(in) :: middle_x
    real(dp), intent(in) :: upper_x
    real(dp), intent(in) :: lower_f
    real(dp), intent(in) :: middle_f
    real(dp), intent(in) :: upper_f
    logical :: valid
    real(dp) :: objective_scale
    real(dp) :: strict_drop_tolerance

    objective_scale = max(1.0_dp, abs(lower_f), abs(middle_f), abs(upper_f))
    strict_drop_tolerance = 64.0_dp*epsilon(1.0_dp)*objective_scale
    valid = ieee_is_finite(lower_x) .and. ieee_is_finite(middle_x) .and. &
      ieee_is_finite(upper_x) .and. ieee_is_finite(lower_f) .and. &
      ieee_is_finite(middle_f) .and. ieee_is_finite(upper_f) .and. &
      lower_x < middle_x .and. middle_x < upper_x .and. &
      lower_f - middle_f > strict_drop_tolerance .and. &
      upper_f - middle_f > strict_drop_tolerance
  end function is_valid_scalar_minimum_bracket

  !> Brent'in safeguarded parabolic/golden-section algoritmasıyla sonlu bir
  !! scalar objective'i minimize eder. Model varsayımı, verilen üç noktanın
  !! gerçek interior minimum bracket'ı oluşturmasıdır. Absolute ve relative
  !! tolerance, x ölçeğinde numerical stopping ölçütüdür; deney belirsizliği
  !! veya mühendislik kabul eşiği değildir. Geçersiz bracket ve non-finite
  !! objective programı durdurmaz, açık status döndürür.
  pure function minimize_scalar_brent( &
      objective, lower_x, middle_x, upper_x, absolute_tolerance, &
      relative_tolerance, maximum_iterations) result(result)
    procedure(scalar_objective_interface) :: objective
    real(dp), intent(in) :: lower_x
    real(dp), intent(in) :: middle_x
    real(dp), intent(in) :: upper_x
    real(dp), intent(in), optional :: absolute_tolerance
    real(dp), intent(in), optional :: relative_tolerance
    integer, intent(in), optional :: maximum_iterations
    type(scalar_minimizer_result_t) :: result

    real(dp), parameter :: golden_ratio_complement = &
      0.3819660112501051518_dp
    real(dp) :: a
    real(dp) :: b
    real(dp) :: d
    real(dp) :: e
    real(dp) :: e_previous
    real(dp) :: f_lower
    real(dp) :: f_upper
    real(dp) :: fv
    real(dp) :: fw
    real(dp) :: fx
    real(dp) :: fu
    real(dp) :: midpoint
    real(dp) :: p
    real(dp) :: q
    real(dp) :: r
    real(dp) :: tolerance_absolute
    real(dp) :: tolerance_relative
    real(dp) :: tolerance_x
    real(dp) :: tolerance_x2
    real(dp) :: u
    real(dp) :: v
    real(dp) :: w
    real(dp) :: x
    integer :: iteration
    integer :: iteration_limit

    result%lower_bound = lower_x
    result%upper_bound = upper_x
    tolerance_absolute = DEFAULT_SCALAR_MINIMIZER_TOLERANCE
    tolerance_relative = DEFAULT_SCALAR_MINIMIZER_TOLERANCE
    iteration_limit = DEFAULT_SCALAR_MINIMIZER_MAX_ITERATIONS
    if (present(absolute_tolerance)) tolerance_absolute = absolute_tolerance
    if (present(relative_tolerance)) tolerance_relative = relative_tolerance
    if (present(maximum_iterations)) iteration_limit = maximum_iterations

    if (.not. ieee_is_finite(tolerance_absolute) .or. &
        .not. ieee_is_finite(tolerance_relative) .or. &
        tolerance_absolute <= 0.0_dp .or. tolerance_relative < 0.0_dp .or. &
        iteration_limit <= 0) then
      result%status = SCALAR_MINIMIZER_INVALID_INPUT
      return
    end if

    f_lower = objective(lower_x)
    fx = objective(middle_x)
    f_upper = objective(upper_x)
    result%function_evaluation_count = 3
    result%x_minimum = middle_x
    result%f_minimum = fx
    if (.not. ieee_is_finite(f_lower) .or. .not. ieee_is_finite(fx) .or. &
        .not. ieee_is_finite(f_upper)) then
      result%status = SCALAR_MINIMIZER_NONFINITE_OBJECTIVE
      return
    end if
    if (.not. is_valid_scalar_minimum_bracket( &
        lower_x, middle_x, upper_x, f_lower, fx, f_upper)) then
      result%status = SCALAR_MINIMIZER_INVALID_BRACKET
      return
    end if

    a = lower_x
    b = upper_x
    x = middle_x
    w = x
    v = x
    fw = fx
    fv = fx
    d = 0.0_dp
    e = 0.0_dp

    do iteration = 1, iteration_limit
      result%iteration_count = iteration
      midpoint = 0.5_dp*(a + b)
      tolerance_x = tolerance_absolute + tolerance_relative*abs(x)
      tolerance_x2 = 2.0_dp*tolerance_x
      if (abs(x - midpoint) <= tolerance_x2 - 0.5_dp*(b - a)) then
        result%status = SCALAR_MINIMIZER_SUCCESS
        result%converged = .true.
        result%x_minimum = x
        result%f_minimum = fx
        return
      end if

      if (abs(e) > tolerance_x) then
        r = (x - w)*(fx - fv)
        q = (x - v)*(fx - fw)
        p = (x - v)*q - (x - w)*r
        q = 2.0_dp*(q - r)
        if (q > 0.0_dp) p = -p
        q = abs(q)
        e_previous = e
        e = d
        if (abs(p) >= abs(0.5_dp*q*e_previous) .or. &
            p <= q*(a - x) .or. p >= q*(b - x)) then
          if (x >= midpoint) then
            e = a - x
          else
            e = b - x
          end if
          d = golden_ratio_complement*e
        else
          d = p/q
          u = x + d
          if (u - a < tolerance_x2 .or. b - u < tolerance_x2) then
            d = sign(tolerance_x, midpoint - x)
          end if
        end if
      else
        if (x >= midpoint) then
          e = a - x
        else
          e = b - x
        end if
        d = golden_ratio_complement*e
      end if

      if (abs(d) >= tolerance_x) then
        u = x + d
      else
        u = x + sign(tolerance_x, d)
      end if
      fu = objective(u)
      result%function_evaluation_count = &
        result%function_evaluation_count + 1
      if (.not. ieee_is_finite(fu)) then
        result%status = SCALAR_MINIMIZER_NONFINITE_OBJECTIVE
        result%x_minimum = x
        result%f_minimum = fx
        return
      end if

      if (fu <= fx) then
        if (u >= x) then
          a = x
        else
          b = x
        end if
        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu
      else
        if (u < x) then
          a = u
        else
          b = u
        end if
        if (fu <= fw .or. scalar_locations_are_equivalent(w, x)) then
          v = w
          fv = fw
          w = u
          fw = fu
        else if (fu <= fv .or. scalar_locations_are_equivalent(v, x) .or. &
            scalar_locations_are_equivalent(v, w)) then
          v = u
          fv = fu
        end if
      end if
    end do

    result%status = SCALAR_MINIMIZER_MAX_ITERATIONS
    result%x_minimum = x
    result%f_minimum = fx
  end function minimize_scalar_brent

  !> Brent state noktalarının yalnız floating-point temsil düzeyinde aynı
  !! konumu taşıyıp taşımadığını sınar. Bu numerical bookkeeping toleransıdır;
  !! physical veya experimental eşdeğerlik anlamına gelmez.
  pure elemental function scalar_locations_are_equivalent(a, b) &
      result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    equivalent = abs(a - b) <= 16.0_dp*epsilon(1.0_dp) * &
      max(1.0_dp, abs(a), abs(b))
  end function scalar_locations_are_equivalent

end module tms_scalar_minimizer
