module tms_complex_linear_problem
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> Complex simetri kontrolünde makine epsilonu ile çarpılan boyutsuz
  !! yuvarlama hata katsayısıdır.
  real(dp), parameter :: SYMMETRY_TOLERANCE_MULTIPLIER = 128.0_dp

  !> Backend'den bağımsız complex symmetric doğrusal denklem problemini taşır.
  !!
  !! Fiziksel açıklama: V0.6 harmonic torsional analizinde coefficient matrix
  !! dynamic stiffness Z [N*m/rad], right-hand sides ise bir veya daha fazla
  !! complex peak torque genliği B [N*m] taşır. Çözüm X complex açısal genlik
  !! [rad] olacaktır.
  !! Matematiksel açıklama: A*X=B problemi tanımlanır. A kare, pozitif
  !! mertebeli, sonlu ve complex symmetric olmalıdır: A^T=A. Hermitian olma
  !! koşulu aranmaz; genel durumda A^H/=A geçerlidir. B'nin satır sayısı A'nın
  !! mertebesine eşit ve sütun sayısı en az bir olmalıdır.
  !! Depolama private ve bağımsız kopyalara dayalıdır. Böylece factorization
  !! backend'i çalışma dizilerini değiştirse bile authoritative A/B korunur.
  type, public :: complex_linear_problem_t
    private
    complex(dp), allocatable :: coefficient_matrix(:, :)
    complex(dp), allocatable :: right_hand_sides(:, :)
  end type complex_linear_problem_t

  public :: create_complex_linear_problem
  public :: validate_complex_linear_problem
  public :: get_complex_linear_problem_order
  public :: get_complex_linear_rhs_count
  public :: get_complex_linear_coefficient_matrix
  public :: get_complex_linear_right_hand_sides

contains

  !> Complex symmetric doğrusal denklem problemini bağımsız kopyalarla kurar.
  !!
  !! Matematiksel model A*X=B'dir. Girdiler aynı fiziksel problemde A [çıktı
  !! birimi/girdi birimi] ve B [çıktı birimi] katsayılarıdır; V0.6 harmonic
  !! kullanımında bunlar sırasıyla [N*m/rad] ve [N*m] olur. Çıktı private
  !! depolamalı complex_linear_problem_t değeridir.
  !! Varsayımlar: n>0, nrhs>=1, A n x n, B n x nrhs, bütün reel/sanal
  !! bileşenler sonlu ve A ölçeğe duyarlı tolerans içinde A^T=A olmalıdır.
  pure function create_complex_linear_problem( &
      coefficient_matrix, right_hand_sides) result(problem)
    complex(dp), intent(in) :: coefficient_matrix(:, :)
    complex(dp), intent(in) :: right_hand_sides(:, :)
    type(complex_linear_problem_t) :: problem

    problem%coefficient_matrix = coefficient_matrix
    problem%right_hand_sides = right_hand_sides
    call validate_complex_linear_problem(problem)
  end function create_complex_linear_problem

  !> Complex doğrusal problem boyut, sonluluk ve complex simetri koşullarını
  !! doğrular.
  !!
  !! Matematiksel açıklama: n>0, nrhs>=1, A n x n, B n x nrhs ve
  !! max|A-A^T|/max|A| <= 128*epsilon koşulları uygulanır. Transpose conjugate
  !! değildir; bu nedenle routine Hermitian simetriyi değil complex simetriyi
  !! sınar. Girdi problem, çıktı yoktur; geçersiz sözleşme error stop üretir.
  pure subroutine validate_complex_linear_problem(problem)
    type(complex_linear_problem_t), intent(in) :: problem

    integer :: problem_order

    if (.not. allocated(problem%coefficient_matrix) .or. &
        .not. allocated(problem%right_hand_sides)) then
      error stop "Complex doğrusal problem kullanılmadan önce oluşturulmalıdır."
    end if

    problem_order = size(problem%coefficient_matrix, 1)
    if (problem_order <= 0) then
      error stop "Complex doğrusal problem en az bir denklem içermelidir."
    end if
    if (size(problem%coefficient_matrix, 2) /= problem_order) then
      error stop "Complex doğrusal problem coefficient matrix'i kare olmalıdır."
    end if
    if (size(problem%right_hand_sides, 1) /= problem_order) then
      error stop "Complex doğrusal problem A ve B satır boyutları uyumsuz."
    end if
    if (size(problem%right_hand_sides, 2) <= 0) then
      error stop "Complex doğrusal problem en az bir right-hand side içermelidir."
    end if

    if (.not. all_complex_finite(problem%coefficient_matrix)) then
      error stop "Complex doğrusal problem A matrisi yalnız sonlu değerler içermelidir."
    end if
    if (.not. all_complex_finite(problem%right_hand_sides)) then
      error stop "Complex doğrusal problem B matrisi yalnız sonlu değerler içermelidir."
    end if
    if (.not. is_complex_symmetric(problem%coefficient_matrix)) then
      error stop "Complex doğrusal problem A^T=A complex simetri koşulunu sağlamalıdır."
    end if
  end subroutine validate_complex_linear_problem

  !> Doğrusal problemdeki denklem sayısını boyutsuz olarak döndürür.
  pure function get_complex_linear_problem_order(problem) result(order)
    type(complex_linear_problem_t), intent(in) :: problem
    integer :: order

    call validate_complex_linear_problem(problem)
    order = size(problem%coefficient_matrix, 1)
  end function get_complex_linear_problem_order

  !> Doğrusal problemdeki bağımsız right-hand side sayısını döndürür.
  pure function get_complex_linear_rhs_count(problem) result(rhs_count)
    type(complex_linear_problem_t), intent(in) :: problem
    integer :: rhs_count

    call validate_complex_linear_problem(problem)
    rhs_count = size(problem%right_hand_sides, 2)
  end function get_complex_linear_rhs_count

  !> Coefficient matrix A'nın bağımsız complex kopyasını döndürür.
  !! V0.6 harmonic kullanımında çıktı tam mantıksal Z [N*m/rad] matrisidir;
  !! yalnız upper triangle kullanan backend'den bağımsız olarak iki üçgen de
  !! residual ve tanı hesapları için korunur.
  pure function get_complex_linear_coefficient_matrix(problem) result(matrix)
    type(complex_linear_problem_t), intent(in) :: problem
    complex(dp), allocatable :: matrix(:, :)

    call validate_complex_linear_problem(problem)
    matrix = problem%coefficient_matrix
  end function get_complex_linear_coefficient_matrix

  !> Right-hand side B'nin bağımsız complex kopyasını döndürür.
  !! V0.6 harmonic kullanımında her sütun bir peak torque load case'i [N*m]
  !! temsil eder.
  pure function get_complex_linear_right_hand_sides(problem) result(rhs)
    type(complex_linear_problem_t), intent(in) :: problem
    complex(dp), allocatable :: rhs(:, :)

    call validate_complex_linear_problem(problem)
    rhs = problem%right_hand_sides
  end function get_complex_linear_right_hand_sides

  !> Complex dizinin bütün reel ve sanal bileşenlerinin IEEE-sonlu olduğunu
  !! değerlendirir. Çıktı boyutsuz logical değerdir.
  pure function all_complex_finite(values) result(is_finite)
    complex(dp), intent(in) :: values(:, :)
    logical :: is_finite

    is_finite = all(ieee_is_finite(real(values, dp))) .and. &
      all(ieee_is_finite(aimag(values)))
  end function all_complex_finite

  !> Kare complex matrisin transpose simetrisini ölçeğe duyarlı doğrular.
  !! Matematiksel ölçüt max|A-A^T|/max|A| <= 128*epsilon biçimindedir.
  !! Conjugate transpose kullanılmaz. Exact sıfır matris simetrik kabul edilir.
  pure function is_complex_symmetric(matrix) result(is_symmetric)
    complex(dp), intent(in) :: matrix(:, :)
    logical :: is_symmetric

    real(dp) :: asymmetry_scale
    real(dp) :: matrix_scale

    matrix_scale = maxval(abs(matrix))
    if (matrix_scale <= 0.0_dp) then
      is_symmetric = .true.
      return
    end if

    asymmetry_scale = maxval(abs(matrix - transpose(matrix)))
    is_symmetric = asymmetry_scale / matrix_scale <= &
      SYMMETRY_TOLERANCE_MULTIPLIER * epsilon(1.0_dp)
  end function is_complex_symmetric

end module tms_complex_linear_problem
