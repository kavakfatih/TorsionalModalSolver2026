module tms_tabulated_dynamic_modulus_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_dynamic_modulus, only : dynamic_shear_modulus, &
    calculate_loss_factor
  use tms_material_frequency, only : material_frequency_point
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    validate_dynamic_material_metadata
  use tms_dynamic_modulus_provider, only : dynamic_modulus_provider_t, &
    dynamic_modulus_evaluation_t, LINEAR_FREQUENCY, &
    LINEAR_LOG_FREQUENCY, are_machine_equivalent
  implicit none
  private

  !> Tek measured isotherm için doğrulanmış G'(f), G''(f) tablosunu taşır.
  !! Authoritative input değerleri private dizilere kopyalanır; caller dizisi
  !! daha sonra değişse bile provider constitutive durumu değişmez.
  type, extends(dynamic_modulus_provider_t), public :: &
      tabulated_dynamic_modulus_provider_t
    private
    real(dp), allocatable :: frequencies_hz(:)
    real(dp), allocatable :: storage_moduli_pa(:)
    real(dp), allocatable :: loss_moduli_pa(:)
    integer :: interpolation_policy = LINEAR_FREQUENCY
    type(dynamic_material_metadata_t) :: metadata
  contains
    procedure, public :: evaluate => evaluate_tabulated_provider
    procedure, public :: get_metadata => get_tabulated_metadata
    procedure, public :: validate => validate_tabulated_provider
  end type tabulated_dynamic_modulus_provider_t

  public :: create_tabulated_dynamic_modulus_provider
  public :: get_tabulated_frequency_points
  public :: get_tabulated_interpolation_policy

contains

  !> Canonical SI birimli deney noktalarından bağımsız-kopyalı provider
  !! oluşturur.
  !! Girdiler: En az iki strictly increasing f [Hz], G'>0 ve G''>=0 [Pa]
  !! noktası; tek operating-state metadata'sı ve seçimlik interpolation
  !! policy. Çıktı: Doğrulanmış tabulated provider.
  !! Varsayımlar ve sınırlar: Her noktanın sıcaklığı metadata isotherm'iyle
  !! machine tolerance içinde aynı olmalıdır. Default LINEAR_FREQUENCY'dir;
  !! extrapolation, temperature interpolation, smoothing ve curve fitting yoktur.
  pure function create_tabulated_dynamic_modulus_provider( &
      points, metadata, interpolation_policy) result(provider)
    type(material_frequency_point), intent(in) :: points(:)
    type(dynamic_material_metadata_t), intent(in) :: metadata
    integer, intent(in), optional :: interpolation_policy
    type(tabulated_dynamic_modulus_provider_t) :: provider

    integer :: point_index

    allocate(provider%frequencies_hz(size(points)))
    allocate(provider%storage_moduli_pa(size(points)))
    allocate(provider%loss_moduli_pa(size(points)))
    do point_index = 1, size(points)
      provider%frequencies_hz(point_index) = points(point_index)%frequency
      provider%storage_moduli_pa(point_index) = &
        points(point_index)%storage_modulus
      provider%loss_moduli_pa(point_index) = points(point_index)%loss_modulus
    end do
    provider%metadata = metadata
    if (present(interpolation_policy)) then
      provider%interpolation_policy = interpolation_policy
    end if

    call validate_points_against_metadata(points, metadata)
    call provider%validate()
  end function create_tabulated_dynamic_modulus_provider

  !> f [Hz] ve aynı measured-isotherm T [K] için G' ve G'' [Pa] üretir.
  !! Matematiksel model LINEAR_FREQUENCY için
  !! alpha=(f-f1)/(f2-f1), LINEAR_LOG_FREQUENCY için
  !! alpha=(log10(f)-log10(f1))/(log10(f2)-log10(f1)); iki durumda da
  !! G=(1-alpha)G1+alpha*G2 ayrı ayrı G' ve G'' üzerine uygulanır.
  !! Girdi/çıktı birimleri: f [Hz], T [K], G'/G'' [Pa], alpha [-].
  !! Varsayımlar ve sınırlar: G modülünün logaritması alınmaz; exact point
  !! machine tolerance ile bulunur. Extrapolation ve temperature interpolation
  !! error stop ile reddedilir; negatif G'' clamp edilmez. Ayrıntılar:
  !! docs/mathematics/dynamic_modulus_interpolation.md.
  pure function evaluate_tabulated_provider( &
      self, frequency_hz, temperature_k) result(evaluation)
    class(tabulated_dynamic_modulus_provider_t), intent(in) :: self
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: temperature_k
    type(dynamic_modulus_evaluation_t) :: evaluation

    integer :: lower_index
    integer :: point_index
    real(dp) :: alpha
    real(dp) :: interpolation_coordinate
    real(dp) :: lower_coordinate
    real(dp) :: upper_coordinate

    call self%validate()
    if (.not. ieee_is_finite(frequency_hz) .or. frequency_hz <= 0.0_dp) then
      error stop "Material provider sorgu frekansı sonlu ve pozitif olmalıdır."
    end if
    if (.not. ieee_is_finite(temperature_k) .or. temperature_k <= 0.0_dp) then
      error stop "Material provider sorgu sıcaklığı sonlu ve pozitif olmalıdır."
    end if
    if (.not. are_machine_equivalent( &
        temperature_k, self%metadata%dataset_temperature_k)) then
      error stop "Operating sıcaklığı provider measured isotherm'iyle eşleşmiyor."
    end if

    evaluation%interpolation_policy = self%interpolation_policy
    do point_index = 1, size(self%frequencies_hz)
      if (are_machine_equivalent( &
          frequency_hz, self%frequencies_hz(point_index))) then
        evaluation%modulus = dynamic_shear_modulus( &
          storage_modulus=self%storage_moduli_pa(point_index), &
          loss_modulus=self%loss_moduli_pa(point_index), &
          frequency=frequency_hz, temperature=temperature_k)
        evaluation%exact_table_point = .true.
        evaluation%lower_frequency_hz = self%frequencies_hz(point_index)
        evaluation%upper_frequency_hz = self%frequencies_hz(point_index)
        evaluation%interpolation_alpha = 0.0_dp
        return
      end if
    end do

    if (frequency_hz < self%frequencies_hz(1) .or. &
        frequency_hz > self%frequencies_hz(size(self%frequencies_hz))) then
      error stop "Material provider frekans extrapolation'ına izin vermez."
    end if

    lower_index = 0
    do point_index = 1, size(self%frequencies_hz) - 1
      if (frequency_hz > self%frequencies_hz(point_index) .and. &
          frequency_hz < self%frequencies_hz(point_index + 1)) then
        lower_index = point_index
        exit
      end if
    end do
    if (lower_index == 0) then
      error stop "Material provider interpolation bracket'i bulunamadı."
    end if

    select case (self%interpolation_policy)
    case (LINEAR_FREQUENCY)
      lower_coordinate = self%frequencies_hz(lower_index)
      upper_coordinate = self%frequencies_hz(lower_index + 1)
      interpolation_coordinate = frequency_hz
    case (LINEAR_LOG_FREQUENCY)
      lower_coordinate = log10(self%frequencies_hz(lower_index))
      upper_coordinate = log10(self%frequencies_hz(lower_index + 1))
      interpolation_coordinate = log10(frequency_hz)
    case default
      error stop "Material provider interpolation policy desteklenmiyor."
    end select

    alpha = (interpolation_coordinate - lower_coordinate) / &
      (upper_coordinate - lower_coordinate)
    evaluation%modulus = dynamic_shear_modulus( &
      storage_modulus=(1.0_dp - alpha) * &
        self%storage_moduli_pa(lower_index) + alpha * &
        self%storage_moduli_pa(lower_index + 1), &
      loss_modulus=(1.0_dp - alpha) * &
        self%loss_moduli_pa(lower_index) + alpha * &
        self%loss_moduli_pa(lower_index + 1), &
      frequency=frequency_hz, temperature=temperature_k)
    evaluation%exact_table_point = .false.
    evaluation%lower_frequency_hz = self%frequencies_hz(lower_index)
    evaluation%upper_frequency_hz = self%frequencies_hz(lower_index + 1)
    evaluation%interpolation_alpha = alpha

    ! calculate_loss_factor çağrısı output passivity ve sonluluk invariantlarını
    ! derived quantity üzerinden de doğrular; tan(delta) burada saklanmaz.
    alpha = calculate_loss_factor(evaluation%modulus)
  end function evaluate_tabulated_provider

  !> Provider dataset metadata'sının bağımsız kopyasını döndürür.
  pure function get_tabulated_metadata(self) result(metadata)
    class(tabulated_dynamic_modulus_provider_t), intent(in) :: self
    type(dynamic_material_metadata_t) :: metadata

    call self%validate()
    metadata = self%metadata
  end function get_tabulated_metadata

  !> Private tablonun dataset, monotonluk, passivity ve policy invariantlarını
  !! doğrular. En az iki nokta, sonlu ve strictly increasing f>0 [Hz],
  !! G'>0 [Pa], G''>=0 [Pa] ve geçerli tek-isotherm metadata gerektirir.
  pure subroutine validate_tabulated_provider(self)
    class(tabulated_dynamic_modulus_provider_t), intent(in) :: self

    integer :: point_index

    if (.not. allocated(self%frequencies_hz) .or. &
        .not. allocated(self%storage_moduli_pa) .or. &
        .not. allocated(self%loss_moduli_pa)) then
      error stop "Tabulated dynamic modulus provider başlatılmamış."
    end if
    if (size(self%frequencies_hz) < 2 .or. &
        size(self%storage_moduli_pa) /= size(self%frequencies_hz) .or. &
        size(self%loss_moduli_pa) /= size(self%frequencies_hz)) then
      error stop "Tabulated provider en az iki eş boyutlu data point gerektirir."
    end if
    if (self%interpolation_policy /= LINEAR_FREQUENCY .and. &
        self%interpolation_policy /= LINEAR_LOG_FREQUENCY) then
      error stop "Tabulated provider interpolation policy geçersiz."
    end if

    call validate_dynamic_material_metadata(self%metadata)
    do point_index = 1, size(self%frequencies_hz)
      if (.not. ieee_is_finite(self%frequencies_hz(point_index)) .or. &
          self%frequencies_hz(point_index) <= 0.0_dp) then
        error stop "Tabulated provider frekansları sonlu ve pozitif olmalıdır."
      end if
      if (.not. ieee_is_finite(self%storage_moduli_pa(point_index)) .or. &
          self%storage_moduli_pa(point_index) <= 0.0_dp) then
        error stop "Tabulated provider G' değerleri sonlu ve pozitif olmalıdır."
      end if
      if (.not. ieee_is_finite(self%loss_moduli_pa(point_index)) .or. &
          self%loss_moduli_pa(point_index) < 0.0_dp) then
        error stop "Tabulated provider G'' negatif veya sonlu olmayan olamaz."
      end if
      if (point_index > 1) then
        if (self%frequencies_hz(point_index) <= &
            self%frequencies_hz(point_index - 1) .or. &
            are_machine_equivalent( &
              self%frequencies_hz(point_index), &
              self%frequencies_hz(point_index - 1))) then
          error stop "Provider frekansları strictly increasing ve distinct olmalıdır."
        end if
      end if
    end do
  end subroutine validate_tabulated_provider

  !> Private tablonun material_frequency_point dizisi olarak bağımsız
  !! kopyasını döndürür. Caller bu kopyayı değiştirse provider etkilenmez.
  pure function get_tabulated_frequency_points(provider) result(points)
    type(tabulated_dynamic_modulus_provider_t), intent(in) :: provider
    type(material_frequency_point), allocatable :: points(:)

    integer :: point_index

    call provider%validate()
    allocate(points(size(provider%frequencies_hz)))
    do point_index = 1, size(points)
      points(point_index)%frequency = provider%frequencies_hz(point_index)
      points(point_index)%temperature = &
        provider%metadata%dataset_temperature_k
      points(point_index)%storage_modulus = &
        provider%storage_moduli_pa(point_index)
      points(point_index)%loss_modulus = provider%loss_moduli_pa(point_index)
    end do
  end function get_tabulated_frequency_points

  !> Provider interpolation policy kimliğini [-] döndürür.
  pure function get_tabulated_interpolation_policy(provider) result(policy)
    type(tabulated_dynamic_modulus_provider_t), intent(in) :: provider
    integer :: policy

    call provider%validate()
    policy = provider%interpolation_policy
  end function get_tabulated_interpolation_policy

  !> Constructor input noktalarının dataset temperature ve constitutive
  !! alanlarını birlikte doğrular. Tüm noktalar tek measured isotherm'e aittir.
  pure subroutine validate_points_against_metadata(points, metadata)
    type(material_frequency_point), intent(in) :: points(:)
    type(dynamic_material_metadata_t), intent(in) :: metadata

    integer :: point_index

    call validate_dynamic_material_metadata(metadata)
    do point_index = 1, size(points)
      if (.not. ieee_is_finite(points(point_index)%temperature) .or. &
          .not. are_machine_equivalent( &
            points(point_index)%temperature, &
            metadata%dataset_temperature_k)) then
        error stop "Tüm frequency point'ler dataset isotherm sıcaklığında olmalıdır."
      end if
    end do
  end subroutine validate_points_against_metadata

end module tms_tabulated_dynamic_modulus_provider
