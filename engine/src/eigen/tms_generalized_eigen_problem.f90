module tms_generalized_eigen_problem
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> Simetri doğrulamasında izin verilen boyutsuz yuvarlama hata çarpanıdır.
  !! Matris ölçeğine göre kullanılan göreli sınır
  !! 128*epsilon(real64) biçimindedir; sabit bir fiziksel tolerans değildir.
  real(dp), parameter :: symmetry_tolerance_multiplier = 128.0_dp

  !> Reel simetrik genelleştirilmiş özdeğer probleminin bağımsız veri kopyasını
  !! taşır.
  !!
  !! Fiziksel açıklama: Torsional modal analizde stiffness katsayıları K
  !! [N*m/rad], mass katsayıları ise polar kütle ataleti M [kg*m^2]
  !! birimindedir. Bu alt katman geometri, malzeme, düğüm ve constraint bilmez.
  !! Matematiksel açıklama: K*phi=lambda*M*phi problemi tanımlanır. K ve M aynı
  !! boyutta reel simetrik kare matrislerdir. K singular veya positive
  !! semidefinite olabilir; bu durum rigid-body mode üretebilir. M'nin positive
  !! definite olması gerekir, ancak bu koşul LAPACK DSYGV INFO semantiğiyle
  !! backend içinde doğrulanır.
  !! Depolama private ve kopyaya dayalıdır; böylece LAPACK'in input matrislerini
  !! overwrite etmesi çağıranın özgün K/M değerlerini değiştiremez.
  type, public :: generalized_eigen_problem_t
    private
    real(dp), allocatable :: stiffness(:, :)
    real(dp), allocatable :: mass(:, :)
  end type generalized_eigen_problem_t

  public :: create_generalized_eigen_problem
  public :: validate_generalized_eigen_problem
  public :: get_generalized_eigen_problem_order
  public :: get_generalized_eigen_stiffness
  public :: get_generalized_eigen_mass

contains

  !> Reel simetrik genelleştirilmiş özdeğer problemi oluşturur.
  !!
  !! Fiziksel açıklama: Girdi K katsayıları torsional rijitlik [N*m/rad], M
  !! katsayıları polar kütle ataleti [kg*m^2] taşır. Çıktı aynı fiziksel sistemi
  !! temsil eden, girdilerden bağımsız generalized_eigen_problem_t değeridir.
  !! Matematiksel açıklama: Problem K*phi=lambda*M*phi biçimindedir. Matrisler
  !! kare, eş boyutlu, sonlu ve ölçeğe bağlı tolerans içinde simetrik olmalıdır.
  !! Varsayımlar ve geçerlilik: K singular olabilir. M-SPD kontrolü burada
  !! yinelenmez; DSYGV faktorizasyon bilgisi authoritative doğrulamadır. Sıfır
  !! boyutlu problem yapısal olarak saklanabilir, fakat eigensolver'a verilmez.
  pure function create_generalized_eigen_problem(stiffness, mass) &
      result(problem)
    real(dp), intent(in) :: stiffness(:, :)
    real(dp), intent(in) :: mass(:, :)
    type(generalized_eigen_problem_t) :: problem

    problem%stiffness = stiffness
    problem%mass = mass
    call validate_generalized_eigen_problem(problem)
  end function create_generalized_eigen_problem

  !> Genelleştirilmiş özdeğer probleminin yapısal önkoşullarını doğrular.
  !!
  !! Matematiksel açıklama: K ve M kare, aynı n mertebesinde, IEEE-sonlu ve
  !! göreli simetri hatası 128*epsilon sınırını aşmayacak biçimde simetrik
  !! olmalıdır. Simetri ölçütü max|A-A^T|/max|A| değeridir; sıfır matris için
  !! exact symmetry aranır. Girdi problemidir; çıktı üretmez. K'nın singular
  !! olması hata değildir ve M-SPD testi LAPACK backend'e bırakılır.
  pure subroutine validate_generalized_eigen_problem(problem)
    type(generalized_eigen_problem_t), intent(in) :: problem

    if (.not. allocated(problem%stiffness) .or. &
        .not. allocated(problem%mass)) then
      error stop "Genelleştirilmiş özdeğer problemi önce oluşturulmalıdır."
    end if

    if (size(problem%stiffness, 1) /= size(problem%stiffness, 2)) then
      error stop "Genelleştirilmiş özdeğer probleminde K kare olmalıdır."
    end if
    if (size(problem%mass, 1) /= size(problem%mass, 2)) then
      error stop "Genelleştirilmiş özdeğer probleminde M kare olmalıdır."
    end if
    if (any(shape(problem%stiffness) /= shape(problem%mass))) then
      error stop "Genelleştirilmiş özdeğer probleminde K ve M eş boyutlu olmalıdır."
    end if

    if (.not. all(ieee_is_finite(problem%stiffness))) then
      error stop "Genelleştirilmiş özdeğer probleminde K sonlu olmalıdır."
    end if
    if (.not. all(ieee_is_finite(problem%mass))) then
      error stop "Genelleştirilmiş özdeğer probleminde M sonlu olmalıdır."
    end if

    if (.not. is_symmetric_matrix(problem%stiffness)) then
      error stop "Genelleştirilmiş özdeğer probleminde K simetrik olmalıdır."
    end if
    if (.not. is_symmetric_matrix(problem%mass)) then
      error stop "Genelleştirilmiş özdeğer probleminde M simetrik olmalıdır."
    end if
  end subroutine validate_generalized_eigen_problem

  !> Genelleştirilmiş özdeğer probleminin matris mertebesini döndürür.
  !! Girdi doğrulanmış problemdir; çıktı K ve M'nin boyutsuz ortak satır/sütun
  !! sayısı n'dir. Tamamen kısıtlanmış yapısal problem için sıfır olabilir.
  pure function get_generalized_eigen_problem_order(problem) result(order)
    type(generalized_eigen_problem_t), intent(in) :: problem
    integer :: order

    call validate_generalized_eigen_problem(problem)
    order = size(problem%stiffness, 1)
  end function get_generalized_eigen_problem_order

  !> Rijitlik matrisi K'nın bağımsız katsayı kopyasını döndürür.
  !! Çıktı [N*m/rad] birimindeki n x n reel matristir. Dönen değer private
  !! problem depolamasını veya çağıranın özgün reduced sistemini değiştiremez.
  pure function get_generalized_eigen_stiffness(problem) result(stiffness)
    type(generalized_eigen_problem_t), intent(in) :: problem
    real(dp), allocatable :: stiffness(:, :)

    call validate_generalized_eigen_problem(problem)
    stiffness = problem%stiffness
  end function get_generalized_eigen_stiffness

  !> Polar atalet matrisi M'nin bağımsız katsayı kopyasını döndürür.
  !! Çıktı [kg*m^2] birimindeki n x n reel matristir. Dönen değer LAPACK çalışma
  !! matrisi olarak değiştirilebilir; problem içindeki özgün M korunur.
  pure function get_generalized_eigen_mass(problem) result(mass)
    type(generalized_eigen_problem_t), intent(in) :: problem
    real(dp), allocatable :: mass(:, :)

    call validate_generalized_eigen_problem(problem)
    mass = problem%mass
  end function get_generalized_eigen_mass

  !> Reel kare matrisin ölçek-bağımsız simetri koşulunu değerlendirir.
  !! Matematiksel ölçüt max|A-A^T|/max|A| <= 128*epsilon biçimindedir.
  !! Girdi üst katmanın fiziksel birimindeki sonlu matristir; çıktı boyutsuz
  !! logical sonuçtur. Sıfır boyutlu ve exact sıfır matris simetrik kabul edilir.
  pure function is_symmetric_matrix(matrix) result(is_symmetric)
    real(dp), intent(in) :: matrix(:, :)
    logical :: is_symmetric

    real(dp) :: asymmetry_scale
    real(dp) :: matrix_scale
    real(dp) :: relative_asymmetry

    if (size(matrix, 1) == 0) then
      is_symmetric = .true.
      return
    end if

    asymmetry_scale = maxval(abs(matrix - transpose(matrix)))
    matrix_scale = maxval(abs(matrix))

    if (.not. matrix_scale > 0.0_dp) then
      ! max|A| sıfırsa sonlu matrisin bütün katsayıları ve asimetrisi sıfırdır.
      is_symmetric = .true.
      return
    end if

    relative_asymmetry = asymmetry_scale / matrix_scale
    is_symmetric = relative_asymmetry <= &
      symmetry_tolerance_multiplier * epsilon(1.0_dp)
  end function is_symmetric_matrix

end module tms_generalized_eigen_problem
