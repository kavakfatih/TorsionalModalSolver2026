module tms_modal_result
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : pi
  implicit none
  private

  !> Modal sonuç alanları arasındaki cebirsel bağıntıların doğrulama katsayısıdır.
  real(dp), parameter :: RESULT_RELATION_EPSILON_MULTIPLIER = 256.0_dp

  !> Rijit-cisim torsional modunun boyutsuz sınıflandırma kimliğidir.
  integer, parameter, public :: RIGID_MODE = 1

  !> Pozitif geri çağırıcı rijitlik taşıyan elastik modun boyutsuz
  !! sınıflandırma kimliğidir.
  integer, parameter, public :: ELASTIC_MODE = 2

  !> Backend'den bağımsız doğrulanmış modal analiz sonucudur.
  !!
  !! Fiziksel sözleşme: Sonuç lineer, sönümsüz ve frozen-property torsional
  !! analiz içindir. Frekanslar Hz, açısal frekanslar rad/s, özdeğerler 1/s^2
  !! cinsindedir. Reduced ve physical mode shape sütunları aynı mode sırasını
  !! kullanır; physical şekiller constrained bileşenlerde sıfırdır.
  !! Depolama private tutulur. Tüm diziler getter yordamlarıyla bağımsız kopya
  !! olarak döndürülür; böylece istemci sonucu veya backend iç durumunu mutasyona
  !! uğratamaz.
  type, public :: modal_result_t
    private
    integer :: number_of_modes = 0
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: angular_frequencies_rad_s(:)
    real(dp), allocatable :: frequencies_hz(:)
    integer, allocatable :: mode_classifications(:)
    real(dp), allocatable :: relative_residuals(:)
    real(dp), allocatable :: reduced_mode_shapes(:, :)
    real(dp), allocatable :: physical_mode_shapes(:, :)
    real(dp) :: mass_orthogonality_error = 0.0_dp
    real(dp) :: rigid_eigenvalue_tolerance = 0.0_dp
    character(len=:), allocatable :: solver_backend_identity
    logical :: linear_model = .true.
    logical :: undamped_model = .true.
    logical :: frozen_property_model = .true.
  end type modal_result_t

  public :: create_modal_result
  public :: validate_modal_result
  public :: get_modal_mode_count
  public :: get_modal_eigenvalues
  public :: get_modal_angular_frequencies_rad_s
  public :: get_modal_frequencies_hz
  public :: get_modal_mode_classifications
  public :: get_modal_relative_residuals
  public :: get_modal_reduced_mode_shapes
  public :: get_modal_physical_mode_shapes
  public :: get_modal_mass_orthogonality_error
  public :: get_modal_rigid_eigenvalue_tolerance
  public :: get_modal_solver_backend_identity
  public :: is_modal_result_linear
  public :: is_modal_result_undamped
  public :: is_modal_result_frozen_property

contains

  !> Birbirine karşılık gelen backend-neutral modal sonuç dizilerini tek ve
  !! değiştirilemez bir sonuç sözleşmesinde birleştirir.
  !!
  !! Fiziksel açıklama: Her mode için lambda, omega, f, sınıf, residual ve mode
  !! shape aynı sütun/indis altında tutulur. Physical mode shape, reduced şeklin
  !! constraint-aware recovery sonucudur.
  !! Matematiksel açıklama: lambda [1/s^2], omega=sqrt(max(lambda,0)) [rad/s]
  !! ve f=omega/(2*pi) [Hz] eşleşmesi çağıran analiz katmanında kurulmuştur.
  !! Girdiler: SI frekans dizileri, boyutsuz sınıf/residual değerleri,
  !! mass-normalized modal genlikler, boyutsuz ortogonallik hatası ve [1/s^2]
  !! rijit toleransıdır. Modal genlik ölçeği genel olarak boyutsuz değildir.
  !! Çıktı: Private depolamalı modal_result_t kopyasıdır.
  !! Varsayımlar: En az bir mode vardır; tüm diziler sonlu ve boyutça tutarlıdır.
  pure function create_modal_result( &
      eigenvalues, angular_frequencies_rad_s, frequencies_hz, &
      mode_classifications, relative_residuals, reduced_mode_shapes, &
      physical_mode_shapes, mass_orthogonality_error, &
      rigid_eigenvalue_tolerance, solver_backend_identity) result(modal_result)
    real(dp), intent(in) :: eigenvalues(:)
    real(dp), intent(in) :: angular_frequencies_rad_s(:)
    real(dp), intent(in) :: frequencies_hz(:)
    integer, intent(in) :: mode_classifications(:)
    real(dp), intent(in) :: relative_residuals(:)
    real(dp), intent(in) :: reduced_mode_shapes(:, :)
    real(dp), intent(in) :: physical_mode_shapes(:, :)
    real(dp), intent(in) :: mass_orthogonality_error
    real(dp), intent(in) :: rigid_eigenvalue_tolerance
    character(len=*), intent(in) :: solver_backend_identity
    type(modal_result_t) :: modal_result

    modal_result%number_of_modes = size(eigenvalues)
    modal_result%eigenvalues = eigenvalues
    modal_result%angular_frequencies_rad_s = angular_frequencies_rad_s
    modal_result%frequencies_hz = frequencies_hz
    modal_result%mode_classifications = mode_classifications
    modal_result%relative_residuals = relative_residuals
    modal_result%reduced_mode_shapes = reduced_mode_shapes
    modal_result%physical_mode_shapes = physical_mode_shapes
    modal_result%mass_orthogonality_error = mass_orthogonality_error
    modal_result%rigid_eigenvalue_tolerance = rigid_eigenvalue_tolerance
    modal_result%solver_backend_identity = trim(solver_backend_identity)
    modal_result%linear_model = .true.
    modal_result%undamped_model = .true.
    modal_result%frozen_property_model = .true.

    call validate_modal_result(modal_result)
  end function create_modal_result

  !> Modal sonuç depolamasının boyut, sonluluk, sınıflandırma ve metadata
  !! bütünlüğünü doğrular.
  !!
  !! Matematiksel açıklama: Tüm mode-bağımlı dizilerin uzunluğu ve her mode
  !! shape matrisinin sütun sayısı number_of_modes olmalıdır. Frekans, residual,
  !! ortogonallik hatası ve tolerans negatif olamaz.
  !! Girdi: modal_result_t. Çıktı üretmez; geçersiz sonuç error stop ile
  !! reddedilir. Fizik denklemini yeniden çözmez.
  pure subroutine validate_modal_result(modal_result)
    type(modal_result_t), intent(in) :: modal_result

    integer :: mode_count
    integer :: mode_index
    real(dp) :: expected_angular_frequency
    real(dp) :: expected_frequency
    real(dp) :: relation_tolerance

    mode_count = modal_result%number_of_modes
    if (mode_count <= 0) then
      error stop "Modal sonuç en az bir mode içermelidir."
    end if
    if (.not. allocated(modal_result%eigenvalues) .or. &
        .not. allocated(modal_result%angular_frequencies_rad_s) .or. &
        .not. allocated(modal_result%frequencies_hz) .or. &
        .not. allocated(modal_result%mode_classifications) .or. &
        .not. allocated(modal_result%relative_residuals) .or. &
        .not. allocated(modal_result%reduced_mode_shapes) .or. &
        .not. allocated(modal_result%physical_mode_shapes) .or. &
        .not. allocated(modal_result%solver_backend_identity)) then
      error stop "Modal sonuç kullanılmadan önce bütün alanları başlatılmalıdır."
    end if

    if (size(modal_result%eigenvalues) /= mode_count .or. &
        size(modal_result%angular_frequencies_rad_s) /= mode_count .or. &
        size(modal_result%frequencies_hz) /= mode_count .or. &
        size(modal_result%mode_classifications) /= mode_count .or. &
        size(modal_result%relative_residuals) /= mode_count .or. &
        size(modal_result%reduced_mode_shapes, 2) /= mode_count .or. &
        size(modal_result%physical_mode_shapes, 2) /= mode_count) then
      error stop "Modal sonuç dizileri aynı mode sayısını taşımalıdır."
    end if

    if (size(modal_result%reduced_mode_shapes, 1) <= 0 .or. &
        size(modal_result%physical_mode_shapes, 1) < &
        size(modal_result%reduced_mode_shapes, 1)) then
      error stop "Physical mode shape boyutu reduced aktif DOF boyutundan az olamaz."
    end if
    if (.not. all(ieee_is_finite(modal_result%eigenvalues)) .or. &
        .not. all(ieee_is_finite( &
        modal_result%angular_frequencies_rad_s)) .or. &
        .not. all(ieee_is_finite(modal_result%frequencies_hz)) .or. &
        .not. all(ieee_is_finite(modal_result%relative_residuals)) .or. &
        .not. all(ieee_is_finite(modal_result%reduced_mode_shapes)) .or. &
        .not. all(ieee_is_finite(modal_result%physical_mode_shapes))) then
      error stop "Modal sonuç yalnız sonlu sayısal değerler içermelidir."
    end if

    if (any(modal_result%angular_frequencies_rad_s < 0.0_dp) .or. &
        any(modal_result%frequencies_hz < 0.0_dp) .or. &
        any(modal_result%relative_residuals < 0.0_dp)) then
      error stop "Modal frekanslar ve göreli kalıntılar negatif olamaz."
    end if
    if (any(modal_result%mode_classifications /= RIGID_MODE .and. &
        modal_result%mode_classifications /= ELASTIC_MODE)) then
      error stop "Modal sonuç bilinmeyen bir mode sınıflandırması içeriyor."
    end if
    if (.not. ieee_is_finite(modal_result%mass_orthogonality_error) .or. &
        modal_result%mass_orthogonality_error < 0.0_dp .or. &
        .not. ieee_is_finite(modal_result%rigid_eigenvalue_tolerance) .or. &
        modal_result%rigid_eigenvalue_tolerance < 0.0_dp) then
      error stop "Modal diagnostic değerleri sonlu ve negatif olmayan olmalıdır."
    end if
    if (len_trim(modal_result%solver_backend_identity) == 0) then
      error stop "Modal sonuç solver backend kimliğini taşımalıdır."
    end if
    if (.not. modal_result%linear_model .or. &
        .not. modal_result%undamped_model .or. &
        .not. modal_result%frozen_property_model) then
      error stop "V0.5 modal sonuç lineer, sönümsüz ve frozen-property olmalıdır."
    end if

    relation_tolerance = RESULT_RELATION_EPSILON_MULTIPLIER * epsilon(1.0_dp)
    do mode_index = 1, mode_count
      if (mode_index > 1) then
        if (modal_result%eigenvalues(mode_index) < &
            modal_result%eigenvalues(mode_index - 1)) then
          error stop "Modal sonuç özdeğerleri azalmayan sırada olmalıdır."
        end if
      end if
      if (modal_result%eigenvalues(mode_index) < &
          -modal_result%rigid_eigenvalue_tolerance) then
        error stop "Modal sonuç anlamlı negatif özdeğer içeremez."
      end if

      if (modal_result%eigenvalues(mode_index) <= &
          modal_result%rigid_eigenvalue_tolerance) then
        if (modal_result%mode_classifications(mode_index) /= RIGID_MODE) then
          error stop "Rijit tolerans içindeki özdeğer RIGID_MODE olmalıdır."
        end if
        expected_angular_frequency = 0.0_dp
      else
        if (modal_result%mode_classifications(mode_index) /= ELASTIC_MODE) then
          error stop "Pozitif elastik özdeğer ELASTIC_MODE olmalıdır."
        end if
        expected_angular_frequency = &
          sqrt(modal_result%eigenvalues(mode_index))
      end if
      expected_frequency = expected_angular_frequency / (2.0_dp * pi)

      if (abs(modal_result%angular_frequencies_rad_s(mode_index) - &
          expected_angular_frequency) > relation_tolerance * &
          max(1.0_dp, expected_angular_frequency) .or. &
          abs(modal_result%frequencies_hz(mode_index) - expected_frequency) > &
          relation_tolerance * max(1.0_dp, expected_frequency)) then
        error stop "Modal sonuç lambda, omega ve frekans bağıntıları tutarsız."
      end if
    end do
  end subroutine validate_modal_result

  !> Modal sonuçtaki boyutsuz mode sayısını döndürür.
  pure function get_modal_mode_count(modal_result) result(mode_count)
    type(modal_result_t), intent(in) :: modal_result
    integer :: mode_count

    call validate_modal_result(modal_result)
    mode_count = modal_result%number_of_modes
  end function get_modal_mode_count

  !> Artan mode sırasındaki özdeğerlerin [1/s^2] bağımsız kopyasını döndürür.
  pure function get_modal_eigenvalues(modal_result) result(eigenvalues)
    type(modal_result_t), intent(in) :: modal_result
    real(dp), allocatable :: eigenvalues(:)

    call validate_modal_result(modal_result)
    eigenvalues = modal_result%eigenvalues
  end function get_modal_eigenvalues

  !> Artan mode sırasındaki omega [rad/s] değerlerinin bağımsız kopyasını verir.
  pure function get_modal_angular_frequencies_rad_s(modal_result) &
      result(angular_frequencies_rad_s)
    type(modal_result_t), intent(in) :: modal_result
    real(dp), allocatable :: angular_frequencies_rad_s(:)

    call validate_modal_result(modal_result)
    angular_frequencies_rad_s = modal_result%angular_frequencies_rad_s
  end function get_modal_angular_frequencies_rad_s

  !> Artan mode sırasındaki doğal frekansların [Hz] bağımsız kopyasını verir.
  pure function get_modal_frequencies_hz(modal_result) result(frequencies_hz)
    type(modal_result_t), intent(in) :: modal_result
    real(dp), allocatable :: frequencies_hz(:)

    call validate_modal_result(modal_result)
    frequencies_hz = modal_result%frequencies_hz
  end function get_modal_frequencies_hz

  !> Her mode için RIGID_MODE veya ELASTIC_MODE kimliklerinin kopyasını verir.
  pure function get_modal_mode_classifications(modal_result) &
      result(mode_classifications)
    type(modal_result_t), intent(in) :: modal_result
    integer, allocatable :: mode_classifications(:)

    call validate_modal_result(modal_result)
    mode_classifications = modal_result%mode_classifications
  end function get_modal_mode_classifications

  !> Her mode için boyutsuz göreli özçift kalıntılarının kopyasını döndürür.
  pure function get_modal_relative_residuals(modal_result) result(residuals)
    type(modal_result_t), intent(in) :: modal_result
    real(dp), allocatable :: residuals(:)

    call validate_modal_result(modal_result)
    residuals = modal_result%relative_residuals
  end function get_modal_relative_residuals

  !> Active Equation ID satırlı, mode sütunlu reduced şekillerin kopyasını verir.
  pure function get_modal_reduced_mode_shapes(modal_result) result(mode_shapes)
    type(modal_result_t), intent(in) :: modal_result
    real(dp), allocatable :: mode_shapes(:, :)

    call validate_modal_result(modal_result)
    mode_shapes = modal_result%reduced_mode_shapes
  end function get_modal_reduced_mode_shapes

  !> Physical DOF satırlı, mode sütunlu full şekillerin bağımsız kopyasını verir.
  pure function get_modal_physical_mode_shapes(modal_result) result(mode_shapes)
    type(modal_result_t), intent(in) :: modal_result
    real(dp), allocatable :: mode_shapes(:, :)

    call validate_modal_result(modal_result)
    mode_shapes = modal_result%physical_mode_shapes
  end function get_modal_physical_mode_shapes

  !> `||Phi^T*M*Phi-I||_F` boyutsuz global diagnostic değerini döndürür.
  pure function get_modal_mass_orthogonality_error(modal_result) &
      result(error_norm)
    type(modal_result_t), intent(in) :: modal_result
    real(dp) :: error_norm

    call validate_modal_result(modal_result)
    error_norm = modal_result%mass_orthogonality_error
  end function get_modal_mass_orthogonality_error

  !> Analizde kullanılan ölçek-duyarlı lambda toleransını [1/s^2] döndürür.
  pure function get_modal_rigid_eigenvalue_tolerance(modal_result) &
      result(tolerance)
    type(modal_result_t), intent(in) :: modal_result
    real(dp) :: tolerance

    call validate_modal_result(modal_result)
    tolerance = modal_result%rigid_eigenvalue_tolerance
  end function get_modal_rigid_eigenvalue_tolerance

  !> Backend kimliğinin çağıran tarafından değiştirilemeyen karakter kopyasını
  !! döndürür; DSYGV ayrıntısı modal hesap API'sinin türlerine sızmaz.
  pure function get_modal_solver_backend_identity(modal_result) &
      result(identity)
    type(modal_result_t), intent(in) :: modal_result
    character(len=:), allocatable :: identity

    call validate_modal_result(modal_result)
    identity = modal_result%solver_backend_identity
  end function get_modal_solver_backend_identity

  !> Sonucun lineer modal model kullandığını doğrulayan metadata bayrağını verir.
  pure function is_modal_result_linear(modal_result) result(is_linear)
    type(modal_result_t), intent(in) :: modal_result
    logical :: is_linear

    call validate_modal_result(modal_result)
    is_linear = modal_result%linear_model
  end function is_modal_result_linear

  !> Sonucun sönümsüz modal model kullandığını belirten metadata bayrağını verir.
  pure function is_modal_result_undamped(modal_result) result(is_undamped)
    type(modal_result_t), intent(in) :: modal_result
    logical :: is_undamped

    call validate_modal_result(modal_result)
    is_undamped = modal_result%undamped_model
  end function is_modal_result_undamped

  !> Sonucun sabitlenmiş malzeme özelliği kullandığını belirten metadata'yı verir.
  pure function is_modal_result_frozen_property(modal_result) &
      result(is_frozen_property)
    type(modal_result_t), intent(in) :: modal_result
    logical :: is_frozen_property

    call validate_modal_result(modal_result)
    is_frozen_property = modal_result%frozen_property_model
  end function is_modal_result_frozen_property

end module tms_modal_result
