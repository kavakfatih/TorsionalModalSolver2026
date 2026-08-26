module tms_lapack_dsygv_backend
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_generalized_eigen_problem, only : generalized_eigen_problem_t, &
    validate_generalized_eigen_problem, get_generalized_eigen_problem_order, &
    get_generalized_eigen_stiffness, get_generalized_eigen_mass
  use tms_eigen_solution, only : eigen_solution_t, initialize_eigen_solution
  implicit none
  private

  integer, parameter :: dsygv_problem_type = 1
  character(len=1), parameter :: calculate_eigenvectors = 'V'
  character(len=1), parameter :: symmetric_triangle = 'U'
  character(len=*), parameter :: backend_identity = 'LAPACK DSYGV (LP64)'

  public :: solve_with_lapack_dsygv

  ! LAPACK Fortran ABI için explicit interface tanımlanır. bind(C) kullanılmaz;
  ! çünkü DSYGV C arayüzü değil, Fortran LAPACK routine'idir. real(dp) proje
  ! yapılandırmasında DOUBLE PRECISION, default INTEGER ise LP64 32-bit INTEGER
  ! ile eşleşmelidir.
  interface
    subroutine dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
      import :: dp
      implicit none
      integer, intent(in) :: itype
      character(len=1), intent(in) :: jobz
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n
      integer, intent(in) :: lda
      integer, intent(in) :: ldb
      real(dp), intent(inout) :: a(lda, *)
      real(dp), intent(inout) :: b(ldb, *)
      real(dp), intent(out) :: w(*)
      real(dp), intent(inout) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine dsygv
  end interface

contains

  !> Reel simetrik-definite generalized eigenproblem'i LAPACK DSYGV ile çözer.
  !!
  !! Fiziksel açıklama: Torsional reduced K [N*m/rad] ile reduced M [kg*m^2]
  !! kullanılarak doğal mode shape'ler ve lambda=omega^2 [1/s^2] elde edilir.
  !! Matematiksel model: ITYPE=1 için K*phi=lambda*M*phi; JOBZ='V' özvektörleri
  !! ister ve UPLO='U' bütün backend boyunca üst simetrik üçgeni seçer. K
  !! singular/PSD olabilir ve rigid-body mode oluşturabilir; M SPD olmalıdır.
  !! Girdi doğrulanmış generalized_eigen_problem_t, çıktı artan lambda sırasıyla
  !! eşleşmiş eigenpair'ler taşıyan eigen_solution_t değeridir.
  !! Varsayımlar ve geçerlilik: DSYGV input A/B matrislerini overwrite ettiği
  !! için yalnız bağımsız çalışma kopyaları kullanılır. LWORK=-1 sorgusundan
  !! sonra original problemden yeni kopyalar yüklenir. Bu external-library
  !! çağrısı yan etkisiz Fortran yordamı olmadığından procedure PURE değildir.
  function solve_with_lapack_dsygv(problem) result(solution)
    type(generalized_eigen_problem_t), intent(in) :: problem
    type(eigen_solution_t) :: solution

    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: mass_work(:, :)
    real(dp), allocatable :: stiffness_work(:, :)
    real(dp), allocatable :: work(:)
    real(dp) :: workspace_query(1)
    integer :: info
    integer :: lwork
    integer :: problem_order

    call validate_lapack_abi()
    call validate_generalized_eigen_problem(problem)
    problem_order = get_generalized_eigen_problem_order(problem)

    if (problem_order <= 0) then
      error stop "LAPACK DSYGV en az bir aktif serbestlik derecesi gerektirir."
    end if

    allocate(eigenvalues(problem_order))
    stiffness_work = get_generalized_eigen_stiffness(problem)
    mass_work = get_generalized_eigen_mass(problem)

    lwork = -1
    call dsygv( &
      dsygv_problem_type, calculate_eigenvectors, symmetric_triangle, &
      problem_order, stiffness_work, problem_order, mass_work, problem_order, &
      eigenvalues, workspace_query, lwork, info)
    call require_successful_dsygv(info, problem_order, .true.)
    lwork = decode_workspace_size(workspace_query(1))
    allocate(work(lwork))

    ! Workspace query'nin A/B üzerinde değişiklik yapmayacağı varsayılmaz.
    ! Actual solve her zaman problem içindeki korunmuş original K/M'den başlar.
    stiffness_work = get_generalized_eigen_stiffness(problem)
    mass_work = get_generalized_eigen_mass(problem)
    call dsygv( &
      dsygv_problem_type, calculate_eigenvectors, symmetric_triangle, &
      problem_order, stiffness_work, problem_order, mass_work, problem_order, &
      eigenvalues, work, lwork, info)
    call require_successful_dsygv(info, problem_order, .false.)

    ! DSYGV, JOBZ='V' ve ITYPE=1 için A sütunlarında M-normalized
    ! eigenvectors ve W içinde artan eigenvalues döndürür.
    call initialize_eigen_solution( &
      solution, eigenvalues, stiffness_work, backend_identity)
  end function solve_with_lapack_dsygv

  !> LAPACK ABI tür varsayımlarını çalıştırma anında doğrular.
  !! dp, Fortran DOUBLE PRECISION türüyle; default INTEGER ise CMake'te seçilen
  !! LP64 LAPACK INTEGER=32-bit sözleşmesiyle eşleşmelidir. Çıktı üretmez.
  subroutine validate_lapack_abi()
    if (dp /= kind(0.0d0)) then
      error stop "TMS26 dp türü LAPACK DOUBLE PRECISION ABI ile uyumsuz."
    end if
    if (storage_size(0) /= 32) then
      error stop "TMS26 DSYGV backend 32-bit LP64 LAPACK INTEGER gerektirir."
    end if
  end subroutine validate_lapack_abi

  !> DSYGV workspace sorgusunun döndürdüğü reel boyutu güvenli tamsayıya çevirir.
  !! Girdi boyutsuz optimal LWORK değeridir. Çıktı pozitif default INTEGER
  !! workspace uzunluğudur; NaN, Inf, sıfır, negatif veya LP64 aralığını aşan
  !! değerler backend diagnostic'i ile reddedilir.
  function decode_workspace_size(workspace_value) result(workspace_size)
    real(dp), intent(in) :: workspace_value
    integer :: workspace_size

    if (.not. ieee_is_finite(workspace_value) .or. workspace_value < 1.0_dp) then
      error stop "LAPACK DSYGV geçerli bir optimal workspace boyutu döndürmedi."
    end if
    if (workspace_value > real(huge(workspace_size), dp)) then
      error stop "LAPACK DSYGV workspace boyutu LP64 INTEGER aralığını aşıyor."
    end if

    workspace_size = ceiling(workspace_value)
  end function decode_workspace_size

  !> LAPACK DSYGV INFO kodunu anlamlı TMS26 diagnostic'ine dönüştürür.
  !! INFO=0 başarıdır. INFO<0 backend çağrı sözleşmesi hatası, 1..N eigen
  !! convergence failure, INFO>N ise M'nin leading principal minor düzeyinde
  !! positive definite olmadığını bildirir. Raw INFO public API'ye sızdırılmaz.
  subroutine require_successful_dsygv(info, problem_order, workspace_phase)
    integer, intent(in) :: info
    integer, intent(in) :: problem_order
    logical, intent(in) :: workspace_phase

    if (info == 0) return

    if (info < 0) then
      if (workspace_phase) then
        error stop &
          "LAPACK DSYGV workspace sorgusunda geçersiz backend argümanı oluştu."
      else
        error stop "LAPACK DSYGV çözüm çağrısında geçersiz backend argümanı oluştu."
      end if
    end if

    if (info <= problem_order) then
      error stop "LAPACK DSYGV genelleştirilmiş özdeğer çözümünde yakınsamadı."
    end if

    error stop "Genelleştirilmiş özdeğer probleminde M positive definite değil."
  end subroutine require_successful_dsygv

end module tms_lapack_dsygv_backend
