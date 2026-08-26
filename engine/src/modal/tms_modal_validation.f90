module tms_modal_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> Otomatik rijit-cisim özdeğer toleransında makine epsilonu ile çarpılan
  !! boyutsuz varsayılan güvenlik katsayısıdır. Sabit bir Hz eşiği değildir.
  real(dp), parameter, public :: &
    AUTO_RIGID_TOLERANCE_MULTIPLIER = 100.0_dp

  !> DSYGV'nin M-normalizasyonunu yeniden ölçeklemeden kabul etmek için
  !! kullanılan makine-duyarlı boyutsuz eşik katsayısıdır.
  real(dp), parameter :: MASS_NORMALIZATION_EPSILON_MULTIPLIER = 100.0_dp

  !> Modal bazın M-ortogonalliği için platformlar arası boyutsuz üst sınırdır.
  real(dp), parameter :: MASS_ORTHOGONALITY_BASE_TOLERANCE = 1.0e-10_dp

  public :: calculate_auto_rigid_eigenvalue_tolerance
  public :: is_rigid_eigenvalue
  public :: normalize_and_validate_mass_modes
  public :: calculate_relative_eigenpair_residuals
  public :: calculate_mass_orthogonality_error

contains

  !> Tek bir özdeğerin ölçeğe duyarlı tolerans altında rijit-cisim modu olup
  !! olmadığını belirler ve anlamlı negatif özdeğeri reddeder.
  !!
  !! Fiziksel açıklama: Serbest torsional sistemin exact sıfır özdeğeri sonlu
  !! hassasiyette küçük negatif veya pozitif çıkabilir; bu durum instability
  !! değildir. Toleransın belirgin biçimde altındaki negatif lambda ise V0.5
  !! lineer modelinde negatif rijitlik enerjisi/kararsızlık tanısıdır.
  !! Matematiksel model: lambda < -tau reddedilir; -tau <= lambda <= tau için
  !! sonuç true, lambda > tau için false olur. Girdiler lambda ve tau [1/s^2],
  !! çıktı boyutsuz logical sınıftır. Her iki girdi sonlu, tau negatif olmayan
  !! olmalıdır; sabit Hz eşiği kullanılmaz.
  pure function is_rigid_eigenvalue(eigenvalue, tolerance) result(is_rigid)
    real(dp), intent(in) :: eigenvalue
    real(dp), intent(in) :: tolerance
    logical :: is_rigid

    if (.not. ieee_is_finite(eigenvalue) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) then
      error stop "Modal özdeğer ve rijit tolerans sonlu, tolerans negatif olmayan olmalıdır."
    end if
    if (eigenvalue < -tolerance) then
      error stop &
        "Modal sistem AUTO toleransın altında negatif, kararsız bir özdeğer içeriyor."
    end if

    is_rigid = eigenvalue <= tolerance
  end function is_rigid_eigenvalue

  !> Rijit-cisim modlarını ölçekten bağımsız sınıflandırmak için özdeğer
  !! toleransını hesaplar.
  !!
  !! Fiziksel açıklama: Serbest sistemlerde tam sıfır olması gereken rijit-cisim
  !! özdeğerleri sonlu hassasiyet nedeniyle küçük pozitif veya negatif çıkabilir.
  !! Matematiksel model:
  !!   tau_lambda = c*epsilon*max(max_i|lambda_i|, ||K||_inf/||M||_inf).
  !! Girdiler: K [N*m/rad], M [kg*m^2], lambda [1/s^2] ve boyutsuz pozitif c.
  !! Çıktı: lambda ile aynı [1/s^2] birimindeki negatif olmayan toleranstır.
  !! Varsayımlar: K ve M eş boyutlu, kare ve sonlu; M normu pozitiftir. Bu
  !! yordam sabit bir frekans eşiği kullanmaz ve fiziksel kararsızlığı gizlemez.
  pure function calculate_auto_rigid_eigenvalue_tolerance( &
      stiffness, mass, eigenvalues, multiplier) result(tolerance)
    real(dp), intent(in) :: stiffness(:, :)
    real(dp), intent(in) :: mass(:, :)
    real(dp), intent(in) :: eigenvalues(:)
    real(dp), intent(in) :: multiplier
    real(dp) :: tolerance

    real(dp) :: eigenvalue_scale
    real(dp) :: mass_norm
    real(dp) :: matrix_scale
    real(dp) :: scaled_epsilon
    real(dp) :: stiffness_norm

    call validate_modal_matrix_pair(stiffness, mass)
    if (size(eigenvalues) <= 0 .or. &
        size(eigenvalues) > size(stiffness, 1) .or. &
        .not. all(ieee_is_finite(eigenvalues))) then
      error stop "Özdeğer dizisi 1..N mode aralığında ve sonlu olmalıdır."
    end if
    if (.not. ieee_is_finite(multiplier) .or. multiplier <= 0.0_dp) then
      error stop "Rijit mod tolerans katsayısı sonlu ve pozitif olmalıdır."
    end if

    stiffness_norm = matrix_infinity_norm(stiffness)
    mass_norm = matrix_infinity_norm(mass)
    if (mass_norm <= 0.0_dp) then
      error stop "Otomatik rijit mod toleransı pozitif M normu gerektirir."
    end if

    eigenvalue_scale = maxval(abs(eigenvalues))
    matrix_scale = safe_nonnegative_ratio(stiffness_norm, mass_norm)
    scaled_epsilon = multiplier * epsilon(1.0_dp)
    if (.not. ieee_is_finite(scaled_epsilon)) then
      error stop "Rijit mod tolerans katsayısı temsil aralığını aşamaz."
    end if

    tolerance = safe_nonnegative_product( &
      max(eigenvalue_scale, matrix_scale), scaled_epsilon)
    if (.not. ieee_is_finite(tolerance)) then
      error stop "Otomatik rijit mod özdeğer toleransı sonlu olmalıdır."
    end if
  end function calculate_auto_rigid_eigenvalue_tolerance

  !> Genelleştirilmiş özvektörleri önce M-normalizasyonu açısından doğrular,
  !! yalnız gerekli sütunlarda sayısal yeniden ölçekleme uygular ve son olarak
  !! tüm modal bazın M-ortogonalliğini sınar.
  !!
  !! Fiziksel açıklama: Kinetik enerji metriği M olduğundan her mod için
  !! phi_i^T*M*phi_i=1 ve farklı modlar için phi_i^T*M*phi_j=0 beklenir.
  !! Matematiksel model: Norm sapması O(epsilon) içindeyse DSYGV sonucu aynen
  !! korunur; daha büyük sonlu sapmada phi<-phi/sqrt(phi^T*M*phi) temizliği
  !! uygulanır. Ardından ||Phi^T*M*Phi-I||_F kontrol edilir.
  !! Girdiler: M [kg*m^2] ve sütunları reduced mode shape olan Phi. Çıktı:
  !! Gerektiğinde kütle-normalize edilmiş aynı boyutlu Phi matrisidir.
  !! Varsayımlar: M SPD ve Phi sonludur. İşaret keyfîdir ve değiştirilmez.
  pure subroutine normalize_and_validate_mass_modes(mass, mode_shapes)
    real(dp), intent(in) :: mass(:, :)
    real(dp), intent(inout) :: mode_shapes(:, :)

    real(dp) :: allowed_orthogonality_error
    real(dp) :: cleanup_trigger
    real(dp) :: mass_norm_squared
    real(dp) :: orthogonality_error
    integer :: mode_index

    call validate_mass_and_modes(mass, mode_shapes)

    cleanup_trigger = MASS_NORMALIZATION_EPSILON_MULTIPLIER * &
      epsilon(1.0_dp) * real(max(1, size(mass, 1)), dp)

    do mode_index = 1, size(mode_shapes, 2)
      mass_norm_squared = dot_product( &
        mode_shapes(:, mode_index), &
        matmul(mass, mode_shapes(:, mode_index)))
      if (.not. ieee_is_finite(mass_norm_squared) .or. &
          mass_norm_squared <= 0.0_dp) then
        error stop "Modal özvektörün M-normu sonlu ve pozitif olmalıdır."
      end if

      if (abs(mass_norm_squared - 1.0_dp) > cleanup_trigger) then
        mode_shapes(:, mode_index) = &
          mode_shapes(:, mode_index) / sqrt(mass_norm_squared)
      end if
    end do

    orthogonality_error = &
      calculate_mass_orthogonality_error(mass, mode_shapes)
    allowed_orthogonality_error = MASS_ORTHOGONALITY_BASE_TOLERANCE * &
      max(1.0_dp, sqrt(real(size(mode_shapes, 2), dp)))
    if (orthogonality_error > allowed_orthogonality_error) then
      error stop "Modal baz Phi^T*M*Phi=I kütle ortogonalliğini sağlamıyor."
    end if
  end subroutine normalize_and_validate_mass_modes

  !> Her özçift için boyutsuz göreli genelleştirilmiş özdeğer kalıntısını
  !! hesaplar.
  !!
  !! Fiziksel açıklama: Kalıntı, hesaplanan modun K*phi=lambda*M*phi denge
  !! eşitliğini ne ölçüde sağladığını gösterir.
  !! Matematiksel model:
  !! rho=||K*phi-lambda*M*phi||_2 /
  !!     ((||K||_inf+|lambda|*||M||_inf)*||phi||_2).
  !! Girdiler: K [N*m/rad], M [kg*m^2], lambda [1/s^2] ve reduced Phi.
  !! Çıktı: Her mode için boyutsuz, sonlu ve negatif olmayan rho dizisidir.
  !! Varsayımlar: Tüm girdiler sonludur. Payda ve pay birlikte sıfırsa tam sıfır
  !! operatörünün doğru özçifti kabul edilerek rho=0 döndürülür.
  pure function calculate_relative_eigenpair_residuals( &
      stiffness, mass, eigenvalues, mode_shapes) result(residuals)
    real(dp), intent(in) :: stiffness(:, :)
    real(dp), intent(in) :: mass(:, :)
    real(dp), intent(in) :: eigenvalues(:)
    real(dp), intent(in) :: mode_shapes(:, :)
    real(dp), allocatable :: residuals(:)

    real(dp) :: denominator
    real(dp) :: mass_norm
    real(dp) :: mode_norm
    real(dp) :: operator_norm
    real(dp) :: residual_norm
    real(dp) :: stiffness_norm
    real(dp), allocatable :: residual_vector(:)
    integer :: mode_index

    call validate_modal_matrix_pair(stiffness, mass)
    call validate_mass_and_modes(mass, mode_shapes)
    if (size(eigenvalues) /= size(mode_shapes, 2) .or. &
        .not. all(ieee_is_finite(eigenvalues))) then
      error stop "Kalıntı hesabında özdeğer ve mode sayıları uyumlu olmalıdır."
    end if

    stiffness_norm = matrix_infinity_norm(stiffness)
    mass_norm = matrix_infinity_norm(mass)
    allocate(residuals(size(eigenvalues)))

    do mode_index = 1, size(eigenvalues)
      residual_vector = matmul(stiffness, mode_shapes(:, mode_index)) - &
        eigenvalues(mode_index) * &
        matmul(mass, mode_shapes(:, mode_index))
      if (.not. all(ieee_is_finite(residual_vector))) then
        error stop "Özçift kalıntısı sonlu sayı aralığında kalmalıdır."
      end if

      residual_norm = stable_vector_two_norm(residual_vector)
      mode_norm = stable_vector_two_norm(mode_shapes(:, mode_index))
      operator_norm = safe_nonnegative_sum( &
        stiffness_norm, safe_nonnegative_product( &
        abs(eigenvalues(mode_index)), mass_norm))
      denominator = safe_nonnegative_product(operator_norm, mode_norm)

      if (denominator <= 0.0_dp) then
        if (residual_norm <= 0.0_dp) then
          residuals(mode_index) = 0.0_dp
        else
          error stop "Özçift göreli kalıntısının paydası sıfırken payı pozitiftir."
        end if
      else
        residuals(mode_index) = residual_norm / denominator
      end if

      if (.not. ieee_is_finite(residuals(mode_index)) .or. &
          residuals(mode_index) < 0.0_dp) then
        error stop "Özçift göreli kalıntısı sonlu ve negatif olmayan olmalıdır."
      end if
    end do
  end function calculate_relative_eigenpair_residuals

  !> Modal bazın global kütle ortogonalliği hatasını hesaplar.
  !!
  !! Fiziksel açıklama: M-normalize modların kinetik enerji metriğinde birim ve
  !! birbirine dik olması gerekir.
  !! Matematiksel model: E_M=||Phi^T*M*Phi-I||_F.
  !! Girdiler: M [kg*m^2] ve reduced modal genlik matrisi Phi. Çıktı: Boyutsuz
  !! ve negatif olmayan Frobenius hata normudur.
  !! Varsayımlar: M kare, Phi satır sayısı M boyutuna eşit ve girdiler sonludur.
  pure function calculate_mass_orthogonality_error(mass, mode_shapes) &
      result(error_norm)
    real(dp), intent(in) :: mass(:, :)
    real(dp), intent(in) :: mode_shapes(:, :)
    real(dp) :: error_norm

    real(dp), allocatable :: gram_error(:, :)
    integer :: mode_index

    call validate_mass_and_modes(mass, mode_shapes)
    allocate(gram_error(size(mode_shapes, 2), size(mode_shapes, 2)))
    gram_error = matmul(transpose(mode_shapes), matmul(mass, mode_shapes))
    do mode_index = 1, size(mode_shapes, 2)
      gram_error(mode_index, mode_index) = &
        gram_error(mode_index, mode_index) - 1.0_dp
    end do

    if (.not. all(ieee_is_finite(gram_error))) then
      error stop "Kütle ortogonalliği Gram matrisi sonlu olmalıdır."
    end if
    error_norm = stable_matrix_frobenius_norm(gram_error)
  end function calculate_mass_orthogonality_error

  !> K ve M matrislerinin kare, eş boyutlu, pozitif boyutlu ve sonlu olduğunu
  !! doğrular. Fiziksel birimler üst seviye modal yordam sözleşmesine aittir.
  pure subroutine validate_modal_matrix_pair(stiffness, mass)
    real(dp), intent(in) :: stiffness(:, :)
    real(dp), intent(in) :: mass(:, :)

    if (size(stiffness, 1) <= 0 .or. &
        size(stiffness, 1) /= size(stiffness, 2) .or. &
        size(mass, 1) /= size(mass, 2) .or. &
        any(shape(stiffness) /= shape(mass))) then
      error stop "Modal doğrulama eş boyutlu, pozitif boyutlu kare K/M gerektirir."
    end if
    if (.not. all(ieee_is_finite(stiffness)) .or. &
        .not. all(ieee_is_finite(mass))) then
      error stop "Modal doğrulama yalnız sonlu K/M katsayıları kabul eder."
    end if
  end subroutine validate_modal_matrix_pair

  !> M matrisi ile modal şekil matrisinin boyut ve sonluluk bütünlüğünü sınar.
  pure subroutine validate_mass_and_modes(mass, mode_shapes)
    real(dp), intent(in) :: mass(:, :)
    real(dp), intent(in) :: mode_shapes(:, :)

    if (size(mass, 1) <= 0 .or. size(mass, 1) /= size(mass, 2) .or. &
        size(mode_shapes, 1) /= size(mass, 1)) then
      error stop "Modal şekil satır sayısı kare M matrisi boyutuyla uyuşmalıdır."
    end if
    if (.not. all(ieee_is_finite(mass)) .or. &
        .not. all(ieee_is_finite(mode_shapes))) then
      error stop "M matrisi ve modal şekiller yalnız sonlu değerler içermelidir."
    end if
  end subroutine validate_mass_and_modes

  !> Sonlu bir matrisin sonsuz normunu ölçekli satır toplamlarıyla hesaplar.
  pure function matrix_infinity_norm(matrix) result(norm_value)
    real(dp), intent(in) :: matrix(:, :)
    real(dp) :: norm_value

    real(dp) :: row_norm
    real(dp) :: row_scale
    integer :: row

    if (.not. all(ieee_is_finite(matrix))) then
      error stop "Matris normu yalnız sonlu katsayılarla hesaplanabilir."
    end if
    norm_value = 0.0_dp
    do row = 1, size(matrix, 1)
      if (size(matrix, 2) == 0) cycle
      row_scale = maxval(abs(matrix(row, :)))
      if (row_scale <= 0.0_dp) cycle
      row_norm = safe_nonnegative_product( &
        row_scale, sum(abs(matrix(row, :)) / row_scale))
      norm_value = max(norm_value, row_norm)
    end do
  end function matrix_infinity_norm

  !> Sonlu vektörün iki normunu taşma riskini azaltan ölçekli biçimde hesaplar.
  pure function stable_vector_two_norm(vector) result(norm_value)
    real(dp), intent(in) :: vector(:)
    real(dp) :: norm_value

    real(dp) :: scale

    if (.not. all(ieee_is_finite(vector))) then
      error stop "Vektör normu yalnız sonlu bileşenlerle hesaplanabilir."
    end if
    if (size(vector) == 0) then
      norm_value = 0.0_dp
      return
    end if
    scale = maxval(abs(vector))
    if (scale <= 0.0_dp) then
      norm_value = 0.0_dp
    else
      norm_value = safe_nonnegative_product( &
        scale, sqrt(sum((vector / scale)**2)))
    end if
  end function stable_vector_two_norm

  !> Sonlu matrisin Frobenius normunu ölçekli biçimde hesaplar.
  pure function stable_matrix_frobenius_norm(matrix) result(norm_value)
    real(dp), intent(in) :: matrix(:, :)
    real(dp) :: norm_value

    real(dp) :: scale

    if (.not. all(ieee_is_finite(matrix))) then
      error stop "Frobenius normu yalnız sonlu katsayılarla hesaplanabilir."
    end if
    if (size(matrix) == 0) then
      norm_value = 0.0_dp
      return
    end if
    scale = maxval(abs(matrix))
    if (scale <= 0.0_dp) then
      norm_value = 0.0_dp
    else
      norm_value = safe_nonnegative_product( &
        scale, sqrt(sum((matrix / scale)**2)))
    end if
  end function stable_matrix_frobenius_norm

  !> Negatif olmayan iki sonlu sayının taşmaya dirençli çarpımını döndürür.
  pure function safe_nonnegative_product(first, second) result(product_value)
    real(dp), intent(in) :: first
    real(dp), intent(in) :: second
    real(dp) :: product_value

    if (.not. ieee_is_finite(first) .or. .not. ieee_is_finite(second) .or. &
        first < 0.0_dp .or. second < 0.0_dp) then
      error stop "Güvenli çarpım sonlu ve negatif olmayan girdiler gerektirir."
    end if
    if (first <= 0.0_dp .or. second <= 0.0_dp) then
      product_value = 0.0_dp
    else if (first > huge(1.0_dp) / second) then
      product_value = huge(1.0_dp)
    else
      product_value = first * second
    end if
  end function safe_nonnegative_product

  !> Negatif olmayan iki sonlu sayının taşmaya dirençli toplamını döndürür.
  pure function safe_nonnegative_sum(first, second) result(sum_value)
    real(dp), intent(in) :: first
    real(dp), intent(in) :: second
    real(dp) :: sum_value

    if (.not. ieee_is_finite(first) .or. .not. ieee_is_finite(second) .or. &
        first < 0.0_dp .or. second < 0.0_dp) then
      error stop "Güvenli toplam sonlu ve negatif olmayan girdiler gerektirir."
    end if
    if (first > huge(1.0_dp) - second) then
      sum_value = huge(1.0_dp)
    else
      sum_value = first + second
    end if
  end function safe_nonnegative_sum

  !> Negatif olmayan iki sonlu sayının taşmaya dirençli oranını döndürür.
  pure function safe_nonnegative_ratio(numerator, denominator) &
      result(ratio_value)
    real(dp), intent(in) :: numerator
    real(dp), intent(in) :: denominator
    real(dp) :: ratio_value

    if (.not. ieee_is_finite(numerator) .or. &
        .not. ieee_is_finite(denominator) .or. numerator < 0.0_dp .or. &
        denominator <= 0.0_dp) then
      error stop "Güvenli oran sonlu, negatif olmayan pay ve pozitif payda ister."
    end if
    if (numerator <= 0.0_dp) then
      ratio_value = 0.0_dp
    else if (denominator < 1.0_dp .and. &
        numerator > huge(1.0_dp) * denominator) then
      ratio_value = huge(1.0_dp)
    else
      ratio_value = numerator / denominator
    end if
  end function safe_nonnegative_ratio

end module tms_modal_validation
