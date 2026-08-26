module tms_eigen_solution
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> Backend-neutral reel genelleştirilmiş özdeğer çözümünü taşır.
  !!
  !! Matematiksel açıklama: Her eigenvalue lambda_i ile eigenvectors(:,i)
  !! aynı eigenpair'i oluşturur ve değerler artan sıradadır. Özvektörlerin satır
  !! sayısı problem DOF sayısı, sütun sayısı mode sayısıdır.
  !! Fiziksel açıklama: Torsional modal problemde lambda [1/s^2] ve özvektörler
  !! modal genlik taşır. Frequency conversion, rigid-mode classification,
  !! residual ve physical recovery bu düşük seviye sonucun sorumluluğu değildir.
  !! Tüm depolama private tutulur ve getters bağımsız kopya döndürür.
  type, public :: eigen_solution_t
    private
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: eigenvectors(:, :)
    character(len=:), allocatable :: backend_identity
  end type eigen_solution_t

  public :: initialize_eigen_solution
  public :: validate_eigen_solution
  public :: get_eigen_mode_count
  public :: get_eigenvalues
  public :: get_eigenvectors
  public :: get_eigen_backend_identity

contains

  !> Backend sonucunu eşleşmiş eigenpair'ler halinde başlatır.
  !!
  !! Matematiksel açıklama: eigenvectors(:,i), eigenvalues(i) değerine aittir;
  !! eigenvalue dizisi azalmayan sırada olmalıdır. Girdiler lambda [1/s^2],
  !! karşılık gelen modal genlik matrisi ve boyutsuz backend kimliğidir. Çıktı
  !! private ve bağımsız kopyalar taşıyan eigen_solution_t değeridir.
  !! Varsayımlar: En az bir DOF ve bir mode bulunur; bütün sayılar sonludur.
  pure subroutine initialize_eigen_solution( &
      solution, eigenvalues, eigenvectors, backend_identity)
    type(eigen_solution_t), intent(out) :: solution
    real(dp), intent(in) :: eigenvalues(:)
    real(dp), intent(in) :: eigenvectors(:, :)
    character(len=*), intent(in) :: backend_identity

    solution%eigenvalues = eigenvalues
    solution%eigenvectors = eigenvectors
    solution%backend_identity = trim(backend_identity)
    call validate_eigen_solution(solution)
  end subroutine initialize_eigen_solution

  !> Backend-neutral eigen solution veri bütünlüğünü doğrular.
  !!
  !! Matematiksel açıklama: Mode sayısı eigenvalue uzunluğudur; eigenvector
  !! sütun sayısı buna eşit, satır sayısı pozitif olmalıdır. Lambda değerleri
  !! [1/s^2] ve bütün modal genlikler sonlu, lambda sırası azalmayan olmalıdır.
  !! Individual eigenvector sign veya repeated-eigenvalue baz seçimi burada
  !! doğrulanmaz; ikisi de fiziksel olarak tekil değildir.
  pure subroutine validate_eigen_solution(solution)
    type(eigen_solution_t), intent(in) :: solution

    integer :: mode_index

    if (.not. allocated(solution%eigenvalues) .or. &
        .not. allocated(solution%eigenvectors) .or. &
        .not. allocated(solution%backend_identity)) then
      error stop "Özdeğer çözümü kullanılmadan önce başlatılmalıdır."
    end if

    if (size(solution%eigenvalues) == 0 .or. &
        size(solution%eigenvectors, 1) == 0) then
      error stop "Özdeğer çözümü en az bir DOF ve bir mode içermelidir."
    end if
    if (size(solution%eigenvectors, 2) /= size(solution%eigenvalues)) then
      error stop "Özdeğer ve özvektör mode sayıları uyumsuz."
    end if
    if (size(solution%eigenvalues) > size(solution%eigenvectors, 1)) then
      error stop "Özdeğer çözümü fiziksel DOF sayısından fazla mode içeremez."
    end if
    if (len_trim(solution%backend_identity) == 0) then
      error stop "Özdeğer çözümü anlamlı bir backend kimliği taşımalıdır."
    end if

    if (.not. all(ieee_is_finite(solution%eigenvalues))) then
      error stop "Özdeğer çözümündeki lambda değerleri sonlu olmalıdır."
    end if
    if (.not. all(ieee_is_finite(solution%eigenvectors))) then
      error stop "Özdeğer çözümündeki modal genlikler sonlu olmalıdır."
    end if

    do mode_index = 2, size(solution%eigenvalues)
      if (solution%eigenvalues(mode_index) < &
          solution%eigenvalues(mode_index - 1)) then
        error stop "Özdeğerler azalmayan sırada saklanmalıdır."
      end if
    end do
  end subroutine validate_eigen_solution

  !> Backend'in döndürdüğü eigenpair sayısını boyutsuz olarak verir.
  pure function get_eigen_mode_count(solution) result(mode_count)
    type(eigen_solution_t), intent(in) :: solution
    integer :: mode_count

    call validate_eigen_solution(solution)
    mode_count = size(solution%eigenvalues)
  end function get_eigen_mode_count

  !> Artan sıradaki eigenvalue değerlerinin bağımsız kopyasını döndürür.
  !! Çıktı lambda(:) [1/s^2] olup eigenvectors getter'ının aynı sütun sırasıyla
  !! eşleşir.
  pure function get_eigenvalues(solution) result(eigenvalues)
    type(eigen_solution_t), intent(in) :: solution
    real(dp), allocatable :: eigenvalues(:)

    call validate_eigen_solution(solution)
    eigenvalues = solution%eigenvalues
  end function get_eigenvalues

  !> Modal eigenvector matrisinin bağımsız kopyasını döndürür.
  !! Çıktıda her sütun aynı indeksli lambda değerine ait modal genlik vektörüdür.
  !! İşaret arbitrary olabilir; phi ve -phi aynı fiziksel mode shape'i temsil eder.
  pure function get_eigenvectors(solution) result(eigenvectors)
    type(eigen_solution_t), intent(in) :: solution
    real(dp), allocatable :: eigenvectors(:, :)

    call validate_eigen_solution(solution)
    eigenvectors = solution%eigenvectors
  end function get_eigenvectors

  !> Çözümü üreten numerical backend'in bağımsız metin kopyasını döndürür.
  !! Kimlik tanısal metadata'dır; modal fiziği veya storage biçimini belirlemez.
  pure function get_eigen_backend_identity(solution) result(backend_identity)
    type(eigen_solution_t), intent(in) :: solution
    character(len=:), allocatable :: backend_identity

    call validate_eigen_solution(solution)
    backend_identity = solution%backend_identity
  end function get_eigen_backend_identity

end module tms_eigen_solution
