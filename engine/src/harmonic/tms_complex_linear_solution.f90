module tms_complex_linear_solution
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> Benzersiz ve çalışma hassasiyetinde güvenilir çözüm bulundu durumudur.
  integer, parameter, public :: COMPLEX_SOLVE_SOLVED = 1

  !> Çözüm ve hata sınırları mevcut olmakla birlikte coefficient matrix'in
  !! çalışma hassasiyetinde kötü koşullu olduğu uyarı durumudur.
  integer, parameter, public :: COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED = 2

  !> Exact singular factorization nedeniyle benzersiz çözüm bulunmadı durumudur.
  integer, parameter, public :: COMPLEX_SOLVE_SINGULAR = 3

  !> Backend'den bağımsız complex doğrusal çözüm ve sayısal tanıları taşır.
  !!
  !! Matematiksel açıklama: SOLVED ve SOLVED_ILL_CONDITIONED durumlarında
  !! response X, A*X=B çözümüdür; reciprocal_condition_number A'nın boyutsuz
  !! RCOND tahminidir. FERR, BERR ve relative_residual her RHS için boyutsuz
  !! tanılardır. SINGULAR durumunda benzersiz X bulunmadığından response,
  !! FERR, BERR ve residual dizileri özellikle allocate edilmez.
  !! Fiziksel açıklama: Harmonic kullanımda X complex peak angular response
  !! [rad] olur. Storage private tutulur ve getter'lar bağımsız kopya döndürür.
  type, public :: complex_linear_solution_t
    private
    integer :: status = 0
    integer :: equation_count = 0
    integer :: rhs_count = 0
    logical :: response_available = .false.
    complex(dp), allocatable :: response(:, :)
    real(dp) :: reciprocal_condition_number = 0.0_dp
    real(dp), allocatable :: forward_error_bounds(:)
    real(dp), allocatable :: backward_errors(:)
    real(dp), allocatable :: relative_residuals(:)
    character(len=:), allocatable :: backend_identity
  end type complex_linear_solution_t

  public :: initialize_available_complex_linear_solution
  public :: initialize_singular_complex_linear_solution
  public :: validate_complex_linear_solution
  public :: get_complex_linear_solution_status
  public :: has_complex_linear_response
  public :: get_complex_linear_solution_equation_count
  public :: get_complex_linear_solution_rhs_count
  public :: get_complex_linear_response
  public :: get_complex_linear_reciprocal_condition_number
  public :: get_complex_linear_forward_error_bounds
  public :: get_complex_linear_backward_errors
  public :: get_complex_linear_relative_residuals
  public :: get_complex_linear_backend_identity

contains

  !> Çözümü bulunan backend sonucunu private solution contract'ına kopyalar.
  !!
  !! Matematiksel açıklama: response n x nrhs, FERR/BERR/residual uzunluğu
  !! nrhs olmalıdır. RCOND, hata tanıları ve residual boyutsuzdur. SOLVED için
  !! RCOND>=epsilon; SOLVED_ILL_CONDITIONED için RCOND<epsilon beklenir.
  !! Fiziksel açıklama: Harmonic response complex peak açıdır [rad].
  !! Varsayımlar: Bütün sayısal sonuçlar sonlu, error değerleri negatif olmayan
  !! ve status çözüm mevcut durumlarından biri olmalıdır.
  pure subroutine initialize_available_complex_linear_solution( &
      solution, status, response, reciprocal_condition_number, &
      forward_error_bounds, backward_errors, relative_residuals, &
      backend_identity)
    type(complex_linear_solution_t), intent(out) :: solution
    integer, intent(in) :: status
    complex(dp), intent(in) :: response(:, :)
    real(dp), intent(in) :: reciprocal_condition_number
    real(dp), intent(in) :: forward_error_bounds(:)
    real(dp), intent(in) :: backward_errors(:)
    real(dp), intent(in) :: relative_residuals(:)
    character(len=*), intent(in) :: backend_identity

    solution%status = status
    solution%equation_count = size(response, 1)
    solution%rhs_count = size(response, 2)
    solution%response_available = .true.
    solution%response = response
    solution%reciprocal_condition_number = reciprocal_condition_number
    solution%forward_error_bounds = forward_error_bounds
    solution%backward_errors = backward_errors
    solution%relative_residuals = relative_residuals
    solution%backend_identity = trim(backend_identity)
    call validate_complex_linear_solution(solution)
  end subroutine initialize_available_complex_linear_solution

  !> Exact singular backend sonucunu çözüm uydurmadan başlatır.
  !!
  !! Matematiksel açıklama: Exact singular A için benzersiz A*X=B çözümü yoktur.
  !! Bu nedenle yalnız n, nrhs, SINGULAR status, RCOND=0 ve backend kimliği
  !! saklanır; X/FERR/BERR/residual allocate edilmez. Sayılar boyutsuz n/nrhs
  !! ile RCOND, harmonic fizik bağlamında denklem sayılarıdır.
  pure subroutine initialize_singular_complex_linear_solution( &
      solution, equation_count, rhs_count, reciprocal_condition_number, &
      backend_identity)
    type(complex_linear_solution_t), intent(out) :: solution
    integer, intent(in) :: equation_count
    integer, intent(in) :: rhs_count
    real(dp), intent(in) :: reciprocal_condition_number
    character(len=*), intent(in) :: backend_identity

    solution%status = COMPLEX_SOLVE_SINGULAR
    solution%equation_count = equation_count
    solution%rhs_count = rhs_count
    solution%response_available = .false.
    solution%reciprocal_condition_number = reciprocal_condition_number
    solution%backend_identity = trim(backend_identity)
    call validate_complex_linear_solution(solution)
  end subroutine initialize_singular_complex_linear_solution

  !> Complex solution storage'ının status-bağımlı bütünlüğünü doğrular.
  !!
  !! SOLVED/ILL_CONDITIONED durumlarında response ve tüm per-RHS tanılar mevcut
  !! olmalıdır. SINGULAR durumda bunların hiçbiri mevcut olamaz; böylece stale
  !! LAPACK çalışma verisi fiziksel response gibi sunulamaz. Girdi solution,
  !! çıktı yoktur; sözleşme hatası error stop üretir.
  pure subroutine validate_complex_linear_solution(solution)
    type(complex_linear_solution_t), intent(in) :: solution

    logical :: available_status

    if (solution%equation_count <= 0 .or. solution%rhs_count <= 0) then
      error stop "Complex doğrusal sonuç pozitif denklem ve RHS sayısı taşımalıdır."
    end if
    if (.not. allocated(solution%backend_identity) .or. &
        len_trim(solution%backend_identity) == 0) then
      error stop "Complex doğrusal sonuç backend kimliğini taşımalıdır."
    end if
    if (.not. ieee_is_finite(solution%reciprocal_condition_number) .or. &
        solution%reciprocal_condition_number < 0.0_dp) then
      error stop "Complex doğrusal sonuç RCOND değeri sonlu ve negatif olmayan olmalıdır."
    end if

    available_status = solution%status == COMPLEX_SOLVE_SOLVED .or. &
      solution%status == COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED

    if (available_status) then
      call validate_available_storage(solution)
      if (solution%status == COMPLEX_SOLVE_SOLVED .and. &
          solution%reciprocal_condition_number < epsilon(1.0_dp)) then
        error stop "SOLVED sonucu çalışma hassasiyetinin altında RCOND taşıyamaz."
      end if
      if (solution%status == COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED .and. &
          solution%reciprocal_condition_number >= epsilon(1.0_dp)) then
        error stop "ILL_CONDITIONED sonucu machine epsilon altında RCOND gerektirir."
      end if
    else if (solution%status == COMPLEX_SOLVE_SINGULAR) then
      if (solution%response_available .or. allocated(solution%response) .or. &
          allocated(solution%forward_error_bounds) .or. &
          allocated(solution%backward_errors) .or. &
          allocated(solution%relative_residuals)) then
        error stop "SINGULAR sonuç response veya error bound uyduramaz."
      end if
      if (solution%reciprocal_condition_number > 0.0_dp) then
        error stop "Exact SINGULAR sonuç RCOND=0 taşımalıdır."
      end if
    else
      error stop "Complex doğrusal sonuç bilinmeyen solver status içeriyor."
    end if
  end subroutine validate_complex_linear_solution

  !> Çözüm mevcut durumlarının dizi boyutu, sonluluk ve negatif olmama
  !! koşullarını doğrular.
  pure subroutine validate_available_storage(solution)
    type(complex_linear_solution_t), intent(in) :: solution

    if (.not. solution%response_available .or. &
        .not. allocated(solution%response) .or. &
        .not. allocated(solution%forward_error_bounds) .or. &
        .not. allocated(solution%backward_errors) .or. &
        .not. allocated(solution%relative_residuals)) then
      error stop "Çözülmüş complex sonuç response ve tüm tanıları taşımalıdır."
    end if
    if (size(solution%response, 1) /= solution%equation_count .or. &
        size(solution%response, 2) /= solution%rhs_count .or. &
        size(solution%forward_error_bounds) /= solution%rhs_count .or. &
        size(solution%backward_errors) /= solution%rhs_count .or. &
        size(solution%relative_residuals) /= solution%rhs_count) then
      error stop "Complex response ve per-RHS tanı boyutları uyumsuz."
    end if
    if (.not. all_complex_finite(solution%response) .or. &
        .not. all(ieee_is_finite(solution%forward_error_bounds)) .or. &
        .not. all(ieee_is_finite(solution%backward_errors)) .or. &
        .not. all(ieee_is_finite(solution%relative_residuals))) then
      error stop "Complex çözüm ve tanılar yalnız sonlu değerler içermelidir."
    end if
    if (any(solution%forward_error_bounds < 0.0_dp) .or. &
        any(solution%backward_errors < 0.0_dp) .or. &
        any(solution%relative_residuals < 0.0_dp)) then
      error stop "Complex çözüm error ve residual tanıları negatif olamaz."
    end if
  end subroutine validate_available_storage

  !> Solver status kimliğini boyutsuz integer olarak döndürür.
  pure function get_complex_linear_solution_status(solution) result(status)
    type(complex_linear_solution_t), intent(in) :: solution
    integer :: status

    call validate_complex_linear_solution(solution)
    status = solution%status
  end function get_complex_linear_solution_status

  !> Benzersiz complex response'un kullanılabilir olup olmadığını döndürür.
  pure function has_complex_linear_response(solution) result(is_available)
    type(complex_linear_solution_t), intent(in) :: solution
    logical :: is_available

    call validate_complex_linear_solution(solution)
    is_available = solution%response_available
  end function has_complex_linear_response

  !> Sonuç denklem sayısını boyutsuz olarak döndürür.
  pure function get_complex_linear_solution_equation_count(solution) &
      result(equation_count)
    type(complex_linear_solution_t), intent(in) :: solution
    integer :: equation_count

    call validate_complex_linear_solution(solution)
    equation_count = solution%equation_count
  end function get_complex_linear_solution_equation_count

  !> Sonuç right-hand side sayısını boyutsuz olarak döndürür.
  pure function get_complex_linear_solution_rhs_count(solution) &
      result(rhs_count)
    type(complex_linear_solution_t), intent(in) :: solution
    integer :: rhs_count

    call validate_complex_linear_solution(solution)
    rhs_count = solution%rhs_count
  end function get_complex_linear_solution_rhs_count

  !> Complex response X'in bağımsız kopyasını döndürür.
  !! Harmonic kullanımda çıktı complex peak angular response [rad] matrisidir.
  !! SINGULAR durumda benzersiz X olmadığından clean diagnostic üretir.
  pure function get_complex_linear_response(solution) result(response)
    type(complex_linear_solution_t), intent(in) :: solution
    complex(dp), allocatable :: response(:, :)

    call validate_complex_linear_solution(solution)
    call require_available_response(solution)
    response = solution%response
  end function get_complex_linear_response

  !> A'nın boyutsuz reciprocal condition number tahminini döndürür.
  !! RCOND numerical conditioning tanısıdır; fiziksel resonance detector olarak
  !! yorumlanmamalıdır. SINGULAR sonuçta değer exact sıfırdır.
  pure function get_complex_linear_reciprocal_condition_number(solution) &
      result(rcond)
    type(complex_linear_solution_t), intent(in) :: solution
    real(dp) :: rcond

    call validate_complex_linear_solution(solution)
    rcond = solution%reciprocal_condition_number
  end function get_complex_linear_reciprocal_condition_number

  !> Her RHS için boyutsuz estimated forward error bound kopyasını döndürür.
  !! SINGULAR durumda error bound hesaplanmadığından clean diagnostic üretir.
  pure function get_complex_linear_forward_error_bounds(solution) result(ferr)
    type(complex_linear_solution_t), intent(in) :: solution
    real(dp), allocatable :: ferr(:)

    call validate_complex_linear_solution(solution)
    call require_available_response(solution)
    ferr = solution%forward_error_bounds
  end function get_complex_linear_forward_error_bounds

  !> Her RHS için boyutsuz componentwise backward error kopyasını döndürür.
  !! SINGULAR durumda error değeri hesaplanmadığından clean diagnostic üretir.
  pure function get_complex_linear_backward_errors(solution) result(berr)
    type(complex_linear_solution_t), intent(in) :: solution
    real(dp), allocatable :: berr(:)

    call validate_complex_linear_solution(solution)
    call require_available_response(solution)
    berr = solution%backward_errors
  end function get_complex_linear_backward_errors

  !> Her RHS için backend-independent boyutsuz relative residual döndürür.
  !! SINGULAR durumda çözüm bulunmadığından residual tanımlı değildir.
  pure function get_complex_linear_relative_residuals(solution) &
      result(residuals)
    type(complex_linear_solution_t), intent(in) :: solution
    real(dp), allocatable :: residuals(:)

    call validate_complex_linear_solution(solution)
    call require_available_response(solution)
    residuals = solution%relative_residuals
  end function get_complex_linear_relative_residuals

  !> Çözümü üreten numerical backend'in bağımsız kimlik kopyasını döndürür.
  pure function get_complex_linear_backend_identity(solution) result(identity)
    type(complex_linear_solution_t), intent(in) :: solution
    character(len=:), allocatable :: identity

    call validate_complex_linear_solution(solution)
    identity = solution%backend_identity
  end function get_complex_linear_backend_identity

  !> Status'un gerçek response taşıdığını doğrular; SINGULAR response erişimini
  !! programlama hatası olarak clean diagnostic ile reddeder.
  pure subroutine require_available_response(solution)
    type(complex_linear_solution_t), intent(in) :: solution

    if (.not. solution%response_available) then
      error stop "SINGULAR complex doğrusal sonuç benzersiz response taşımıyor."
    end if
  end subroutine require_available_response

  !> Complex dizinin reel ve sanal bileşen sonluluğunu değerlendirir.
  pure function all_complex_finite(values) result(is_finite)
    complex(dp), intent(in) :: values(:, :)
    logical :: is_finite

    is_finite = all(ieee_is_finite(real(values, dp))) .and. &
      all(ieee_is_finite(aimag(values)))
  end function all_complex_finite

end module tms_complex_linear_solution
