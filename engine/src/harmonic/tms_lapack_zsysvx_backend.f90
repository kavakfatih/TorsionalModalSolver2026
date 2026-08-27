module tms_lapack_zsysvx_backend
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_complex_linear_problem, only : complex_linear_problem_t, &
    validate_complex_linear_problem, get_complex_linear_problem_order, &
    get_complex_linear_rhs_count, get_complex_linear_coefficient_matrix, &
    get_complex_linear_right_hand_sides
  use tms_complex_linear_solution, only : COMPLEX_SOLVE_SOLVED, &
    COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED, COMPLEX_SOLVE_SINGULAR
  implicit none
  private

  character(len=1), parameter :: compute_factorization = 'N'
  character(len=1), parameter :: symmetric_triangle = 'U'
  character(len=*), parameter :: backend_name = 'LAPACK ZSYSVX (LP64)'

  public :: solve_with_lapack_zsysvx

  ! LAPACK Fortran ABI için explicit interface kullanılır. bind(C) kullanılmaz;
  ! çünkü ZSYSVX C routine'i değil, Fortran LAPACK routine'idir. complex(dp)
  ! DOUBLE COMPLEX, real(dp) DOUBLE PRECISION ve default INTEGER LP64 32-bit
  ! INTEGER ile eşleşmelidir.
  interface
    subroutine zsysvx( &
        fact, uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &
        rcond, ferr, berr, work, lwork, rwork, info)
      import :: dp
      implicit none
      character(len=1), intent(in) :: fact
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      integer, intent(in) :: lda
      integer, intent(in) :: ldaf
      integer, intent(in) :: ldb
      integer, intent(in) :: ldx
      integer, intent(in) :: lwork
      complex(dp), intent(in) :: a(lda, *)
      complex(dp), intent(inout) :: af(ldaf, *)
      integer, intent(inout) :: ipiv(*)
      complex(dp), intent(in) :: b(ldb, *)
      complex(dp), intent(out) :: x(ldx, *)
      real(dp), intent(out) :: rcond
      real(dp), intent(out) :: ferr(*)
      real(dp), intent(out) :: berr(*)
      complex(dp), intent(inout) :: work(*)
      real(dp), intent(out) :: rwork(*)
      integer, intent(out) :: info
    end subroutine zsysvx
  end interface

contains

  !> Complex symmetric A*X=B problemini LAPACK ZSYSVX ile çözer.
  !!
  !! Matematiksel model: FACT='N' ile A=U*D*U^T diagonal-pivoting
  !! factorization'ı backend içinde oluşturulur; UPLO='U' yalnız üst üçgeni
  !! factorization için seçer. ZSYSVX condition estimate, iterative refinement,
  !! FERR ve BERR üretir. Matris complex symmetric A^T=A'dır; Hermitian olması
  !! gerekmez.
  !! Fiziksel açıklama: Harmonic kullanımında A=Z [N*m/rad], B peak torque
  !! [N*m], X peak angular response [rad] olur. Bu alt katman fizik modelini
  !! bilmez ve multiple RHS'yi native biçimde korur.
  !! Çıktılar: Backend-neutral status, yalnız çözüm mevcut durumlarında X ile
  !! per-RHS FERR/BERR, her durumda anlamlı RCOND ve backend kimliğidir.
  !! Exact singular durumda X/FERR/BERR allocate edilmez. Workspace LWORK=-1
  !! sorgusuyla belirlenir; authoritative problem A/B değerleri değiştirilmez.
  subroutine solve_with_lapack_zsysvx( &
      problem, status, response, reciprocal_condition_number, &
      forward_error_bounds, backward_errors, backend_identity)
    type(complex_linear_problem_t), intent(in) :: problem
    integer, intent(out) :: status
    complex(dp), allocatable, intent(out) :: response(:, :)
    real(dp), intent(out) :: reciprocal_condition_number
    real(dp), allocatable, intent(out) :: forward_error_bounds(:)
    real(dp), allocatable, intent(out) :: backward_errors(:)
    character(len=:), allocatable, intent(out) :: backend_identity

    complex(dp), allocatable :: coefficient_work(:, :)
    complex(dp), allocatable :: factor_work(:, :)
    complex(dp), allocatable :: response_work(:, :)
    complex(dp), allocatable :: right_hand_side_work(:, :)
    complex(dp), allocatable :: work(:)
    complex(dp) :: workspace_query(1)
    integer, allocatable :: pivot_indices(:)
    real(dp), allocatable :: backward_error_work(:)
    real(dp), allocatable :: forward_error_work(:)
    real(dp), allocatable :: real_work(:)
    real(dp) :: rcond_work
    integer :: info
    integer :: lwork
    integer :: problem_order
    integer :: rhs_count

    call validate_lapack_complex_abi()
    call validate_complex_linear_problem(problem)
    problem_order = get_complex_linear_problem_order(problem)
    rhs_count = get_complex_linear_rhs_count(problem)
    backend_identity = backend_name

    allocate(factor_work(problem_order, problem_order))
    allocate(response_work(problem_order, rhs_count))
    allocate(pivot_indices(problem_order))
    allocate(forward_error_work(rhs_count))
    allocate(backward_error_work(rhs_count))
    allocate(real_work(max(1, problem_order)))

    call reset_backend_work_arrays( &
      problem, coefficient_work, right_hand_side_work, factor_work, &
      pivot_indices, response_work, rcond_work, forward_error_work, &
      backward_error_work, real_work)

    lwork = -1
    workspace_query = cmplx(0.0_dp, 0.0_dp, kind=dp)
    call zsysvx( &
      compute_factorization, symmetric_triangle, problem_order, rhs_count, &
      coefficient_work, problem_order, factor_work, problem_order, &
      pivot_indices, right_hand_side_work, problem_order, response_work, &
      problem_order, rcond_work, forward_error_work, backward_error_work, &
      workspace_query, lwork, real_work, info)
    call require_successful_workspace_query(info)
    lwork = decode_workspace_size(workspace_query(1), problem_order)
    allocate(work(lwork))

    ! Workspace query'nin hiçbir çalışma dizisini koruduğu varsayılmaz.
    ! Actual solve authoritative problem içindeki bağımsız A/B kopyalarından
    ! ve yeniden sıfırlanmış factor/diagnostic storage'dan başlar.
    call reset_backend_work_arrays( &
      problem, coefficient_work, right_hand_side_work, factor_work, &
      pivot_indices, response_work, rcond_work, forward_error_work, &
      backward_error_work, real_work)
    work = cmplx(0.0_dp, 0.0_dp, kind=dp)

    call zsysvx( &
      compute_factorization, symmetric_triangle, problem_order, rhs_count, &
      coefficient_work, problem_order, factor_work, problem_order, &
      pivot_indices, right_hand_side_work, problem_order, response_work, &
      problem_order, rcond_work, forward_error_work, backward_error_work, &
      work, lwork, real_work, info)

    call translate_zsysvx_result( &
      info, problem_order, response_work, rcond_work, forward_error_work, &
      backward_error_work, status, response, reciprocal_condition_number, &
      forward_error_bounds, backward_errors)
  end subroutine solve_with_lapack_zsysvx

  !> LAPACK çalışma dizilerini actual problem girdilerinden yeniden başlatır.
  !! Girdi A/B private problem kopyalarıdır; çıktı backend'e ait mutable A/B,
  !! AF, IPIV, X ve diagnostic çalışma alanlarıdır. Fiziksel birimler A/B
  !! tarafından taşınır; sıfırlanan storage henüz fiziksel sonuç değildir.
  subroutine reset_backend_work_arrays( &
      problem, coefficient_work, right_hand_side_work, factor_work, &
      pivot_indices, response_work, rcond_work, forward_error_work, &
      backward_error_work, real_work)
    type(complex_linear_problem_t), intent(in) :: problem
    complex(dp), allocatable, intent(out) :: coefficient_work(:, :)
    complex(dp), allocatable, intent(out) :: right_hand_side_work(:, :)
    complex(dp), intent(out) :: factor_work(:, :)
    integer, intent(out) :: pivot_indices(:)
    complex(dp), intent(out) :: response_work(:, :)
    real(dp), intent(out) :: rcond_work
    real(dp), intent(out) :: forward_error_work(:)
    real(dp), intent(out) :: backward_error_work(:)
    real(dp), intent(out) :: real_work(:)

    coefficient_work = get_complex_linear_coefficient_matrix(problem)
    right_hand_side_work = get_complex_linear_right_hand_sides(problem)
    factor_work = cmplx(0.0_dp, 0.0_dp, kind=dp)
    pivot_indices = 0
    response_work = cmplx(0.0_dp, 0.0_dp, kind=dp)
    rcond_work = 0.0_dp
    forward_error_work = 0.0_dp
    backward_error_work = 0.0_dp
    real_work = 0.0_dp
  end subroutine reset_backend_work_arrays

  !> LAPACK ABI tür varsayımlarını çalıştırma anında doğrular.
  !! dp Fortran DOUBLE PRECISION, complex(dp) DOUBLE COMPLEX ve default INTEGER
  !! CMake ile seçilen LP64 32-bit INTEGER sözleşmesiyle eşleşmelidir.
  subroutine validate_lapack_complex_abi()
    if (dp /= kind(0.0d0)) then
      error stop "TMS26 dp türü LAPACK DOUBLE PRECISION ABI ile uyumsuz."
    end if
    if (storage_size(cmplx(0.0_dp, 0.0_dp, kind=dp)) /= 128) then
      error stop "TMS26 complex(dp) türü LAPACK DOUBLE COMPLEX ABI ile uyumsuz."
    end if
    if (storage_size(0) /= 32) then
      error stop "TMS26 ZSYSVX backend 32-bit LP64 LAPACK INTEGER gerektirir."
    end if
  end subroutine validate_lapack_complex_abi

  !> ZSYSVX workspace sorgusunun yalnız başarılı INFO=0 ile dönmesini ister.
  !! Workspace query factorization yapmadığından pozitif INFO de internal
  !! backend contract hatasıdır; public analysis status olarak yorumlanmaz.
  subroutine require_successful_workspace_query(info)
    integer, intent(in) :: info

    if (info < 0) then
      error stop "LAPACK ZSYSVX workspace sorgusunda geçersiz argüman oluştu."
    end if
    if (info > 0) then
      error stop "LAPACK ZSYSVX workspace sorgusu beklenmeyen status döndürdü."
    end if
  end subroutine require_successful_workspace_query

  !> Complex workspace query değerinin reel kısmını güvenli LP64 LWORK
  !! uzunluğuna çevirir.
  !!
  !! Matematiksel açıklama: ZSYSVX için minimum LWORK=max(1,2*N)'dir; optimal
  !! değer WORK(1)'in reel kısmındadır. Reel ve sanal bileşen sonlu olmalı,
  !! optimal değer minimumu sağlamalı ve default INTEGER aralığında kalmalıdır.
  function decode_workspace_size(workspace_value, problem_order) &
      result(workspace_size)
    complex(dp), intent(in) :: workspace_value
    integer, intent(in) :: problem_order
    integer :: workspace_size

    integer :: minimum_workspace
    real(dp) :: workspace_real

    if (real(problem_order, dp) > &
        real(huge(minimum_workspace), dp) / 2.0_dp) then
      error stop "LAPACK ZSYSVX problem boyutu LP64 workspace aralığını aşıyor."
    end if
    minimum_workspace = max(1, 2 * problem_order)
    workspace_real = real(workspace_value, dp)
    if (.not. ieee_is_finite(workspace_real) .or. &
        .not. ieee_is_finite(aimag(workspace_value)) .or. &
        workspace_real < real(minimum_workspace, dp)) then
      error stop "LAPACK ZSYSVX geçerli optimal workspace boyutu döndürmedi."
    end if
    if (workspace_real > real(huge(workspace_size), dp)) then
      error stop "LAPACK ZSYSVX workspace boyutu LP64 INTEGER aralığını aşıyor."
    end if

    workspace_size = ceiling(workspace_real)
  end function decode_workspace_size

  !> ZSYSVX INFO kodunu backend-neutral status ve status'a uygun sonuçlara
  !! dönüştürür.
  !!
  !! INFO=0 SOLVED; 1..N exact factor singularity nedeniyle SINGULAR;
  !! INFO=N+1 çözüm ve hata sınırları mevcut SOLVED_ILL_CONDITIONED durumudur.
  !! INFO<0 veya INFO>N+1 backend/programming contract hatasıdır. Singular
  !! durumda LAPACK X/FERR/BERR tanımlamadığından bu diziler kopyalanmaz.
  subroutine translate_zsysvx_result( &
      info, problem_order, response_work, rcond_work, forward_error_work, &
      backward_error_work, status, response, reciprocal_condition_number, &
      forward_error_bounds, backward_errors)
    integer, intent(in) :: info
    integer, intent(in) :: problem_order
    complex(dp), intent(in) :: response_work(:, :)
    real(dp), intent(in) :: rcond_work
    real(dp), intent(in) :: forward_error_work(:)
    real(dp), intent(in) :: backward_error_work(:)
    integer, intent(out) :: status
    complex(dp), allocatable, intent(out) :: response(:, :)
    real(dp), intent(out) :: reciprocal_condition_number
    real(dp), allocatable, intent(out) :: forward_error_bounds(:)
    real(dp), allocatable, intent(out) :: backward_errors(:)

    if (info < 0) then
      error stop "LAPACK ZSYSVX çözüm çağrısında geçersiz argüman oluştu."
    end if
    if (info > problem_order + 1) then
      error stop "LAPACK ZSYSVX belgelenmeyen bir INFO status döndürdü."
    end if

    if (info >= 1 .and. info <= problem_order) then
      status = COMPLEX_SOLVE_SINGULAR
      reciprocal_condition_number = 0.0_dp
      return
    end if

    if (.not. ieee_is_finite(rcond_work) .or. rcond_work < 0.0_dp) then
      error stop "LAPACK ZSYSVX sonlu ve negatif olmayan RCOND döndürmelidir."
    end if
    if (info == 0) then
      status = COMPLEX_SOLVE_SOLVED
    else if (info == problem_order + 1) then
      status = COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED
    else
      error stop "LAPACK ZSYSVX çözülemeyen internal status üretti."
    end if

    reciprocal_condition_number = rcond_work
    response = response_work
    forward_error_bounds = forward_error_work
    backward_errors = backward_error_work
  end subroutine translate_zsysvx_result

end module tms_lapack_zsysvx_backend
