module tms_tabulated_temperature_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_temperature_shift_types, only : temperature_shift_domain_t, &
    temperature_shift_evaluation_t, TABULATED_LOG10_SHIFT
  use tms_temperature_shift_provider, only : temperature_shift_provider_t, &
    validate_temperature_shift_domain, temperature_is_in_shift_domain, &
    are_temperature_shift_values_machine_equivalent, &
    calculate_shift_factor_from_log10
  implicit none
  private

  !> Deneysel sıcaklık noktalarında s=log10(a_T) verisini saklar. Input
  !! dizileri private alana kopyalandığından caller mutation provider'ı
  !! değiştirmez. Shift değerlerinde monotonicity fizik kuralı olarak zorlanmaz.
  type, extends(temperature_shift_provider_t), public :: &
      tabulated_log10_shift_provider_t
    private
    real(dp), allocatable :: temperature_points_k(:)
    real(dp), allocatable :: log10_a_t_points(:)
    type(temperature_shift_domain_t) :: domain
  contains
    procedure, public :: evaluate => evaluate_tabulated_temperature_shift
    procedure, public :: validate => validate_tabulated_temperature_shift
    procedure, public :: get_domain => get_tabulated_temperature_domain
  end type tabulated_log10_shift_provider_t

  public :: create_tabulated_temperature_shift_provider
  public :: get_tabulated_shift_temperature_points
  public :: get_tabulated_log10_shift_points

contains

  !> En az iki temperature [K] ve s=log10(a_T) [-] noktasından doğrulanmış
  !! provider oluşturur. T noktaları sonlu, pozitif, distinct ve strictly
  !! increasing olmalıdır. T_ref [K] tabloda explicit bulunmalı ve o noktada
  !! s machine tolerance içinde sıfır olmalıdır. Extrapolation uygulanmaz.
  pure function create_tabulated_temperature_shift_provider( &
      temperature_points_k, log10_a_t_points, reference_temperature_k) &
      result(provider)
    real(dp), intent(in) :: temperature_points_k(:)
    real(dp), intent(in) :: log10_a_t_points(:)
    real(dp), intent(in) :: reference_temperature_k
    type(tabulated_log10_shift_provider_t) :: provider

    provider%temperature_points_k = temperature_points_k
    provider%log10_a_t_points = log10_a_t_points
    if (size(temperature_points_k) >= 1) then
      provider%domain%minimum_temperature_k = temperature_points_k(1)
      provider%domain%maximum_temperature_k = &
        temperature_points_k(size(temperature_points_k))
    end if
    provider%domain%reference_temperature_k = reference_temperature_k

    call provider%validate()
  end function create_tabulated_temperature_shift_provider

  !> T [K] ekseninde s=log10(a_T) değerini lineer interpolate eder:
  !! alpha=(T-T1)/(T2-T1), s=(1-alpha)s1+alpha*s2; ardından a_T=10^s.
  !! Doğrudan a_T interpolation yapılmaz. Çıktı boyutsuz shift ve temperature
  !! bracket izidir. Validated kapalı domain dışında extrapolation reddedilir.
  pure function evaluate_tabulated_temperature_shift(self, temperature_k) &
      result(evaluation)
    class(tabulated_log10_shift_provider_t), intent(in) :: self
    real(dp), intent(in) :: temperature_k
    type(temperature_shift_evaluation_t) :: evaluation

    integer :: lower_index
    integer :: point_index
    real(dp) :: alpha

    call self%validate()
    if (.not. temperature_is_in_shift_domain(temperature_k, self%domain)) then
      error stop "Tabulated temperature-shift query domain dışında."
    end if

    evaluation%shift_model_kind = TABULATED_LOG10_SHIFT
    evaluation%operating_temperature_k = temperature_k
    evaluation%reference_temperature_k = &
      self%domain%reference_temperature_k
    evaluation%has_temperature_bracket = .true.

    do point_index = 1, size(self%temperature_points_k)
      if (are_temperature_shift_values_machine_equivalent( &
          temperature_k, self%temperature_points_k(point_index))) then
        if (are_temperature_shift_values_machine_equivalent( &
            self%temperature_points_k(point_index), &
            self%domain%reference_temperature_k)) then
          ! Reference shift machine-zero girdisi canonical identity'ye çevrilir.
          evaluation%log10_a_t = 0.0_dp
        else
          evaluation%log10_a_t = self%log10_a_t_points(point_index)
        end if
        evaluation%a_t = calculate_shift_factor_from_log10( &
          evaluation%log10_a_t)
        evaluation%exact_temperature_point = .true.
        evaluation%lower_temperature_k = &
          self%temperature_points_k(point_index)
        evaluation%upper_temperature_k = &
          self%temperature_points_k(point_index)
        evaluation%interpolation_alpha = 0.0_dp
        return
      end if
    end do

    lower_index = 0
    do point_index = 1, size(self%temperature_points_k) - 1
      if (temperature_k > self%temperature_points_k(point_index) .and. &
          temperature_k < self%temperature_points_k(point_index + 1)) then
        lower_index = point_index
        exit
      end if
    end do
    if (lower_index == 0) then
      error stop "Tabulated temperature-shift interpolation bracket'i bulunamadı."
    end if

    alpha = (temperature_k - self%temperature_points_k(lower_index)) / &
      (self%temperature_points_k(lower_index + 1) - &
        self%temperature_points_k(lower_index))
    evaluation%log10_a_t = (1.0_dp - alpha) * &
      self%log10_a_t_points(lower_index) + alpha * &
      self%log10_a_t_points(lower_index + 1)
    evaluation%a_t = calculate_shift_factor_from_log10( &
      evaluation%log10_a_t)
    evaluation%exact_temperature_point = .false.
    evaluation%lower_temperature_k = &
      self%temperature_points_k(lower_index)
    evaluation%upper_temperature_k = &
      self%temperature_points_k(lower_index + 1)
    evaluation%interpolation_alpha = alpha
  end function evaluate_tabulated_temperature_shift

  !> Private tablonun boyut, sonluluk, reference identity ve temperature
  !! ordering invariantlarını doğrular. s değerleri boyutsuz ve sonlu olmalı;
  !! measured shift verisine monotonicity zorlanmaz.
  pure subroutine validate_tabulated_temperature_shift(self)
    class(tabulated_log10_shift_provider_t), intent(in) :: self

    integer :: point_index
    integer :: reference_index

    if (.not. allocated(self%temperature_points_k) .or. &
        .not. allocated(self%log10_a_t_points)) then
      error stop "Tabulated temperature-shift provider başlatılmamış."
    end if
    if (size(self%temperature_points_k) < 2 .or. &
        size(self%temperature_points_k) /= size(self%log10_a_t_points)) then
      error stop "Tabulated temperature shift en az iki eş boyutlu nokta ister."
    end if

    call validate_temperature_shift_domain(self%domain)
    reference_index = 0
    do point_index = 1, size(self%temperature_points_k)
      if (.not. ieee_is_finite(self%temperature_points_k(point_index)) .or. &
          self%temperature_points_k(point_index) <= 0.0_dp) then
        error stop "Tabulated shift sıcaklıkları sonlu ve pozitif olmalıdır."
      end if
      if (.not. ieee_is_finite(self%log10_a_t_points(point_index))) then
        error stop "Tabulated log10(a_T) değerleri sonlu olmalıdır."
      end if
      if (point_index > 1) then
        if (self%temperature_points_k(point_index) <= &
            self%temperature_points_k(point_index - 1) .or. &
            are_temperature_shift_values_machine_equivalent( &
              self%temperature_points_k(point_index), &
              self%temperature_points_k(point_index - 1))) then
          error stop "Tabulated shift sıcaklıkları strictly increasing olmalıdır."
        end if
      end if
      if (are_temperature_shift_values_machine_equivalent( &
          self%temperature_points_k(point_index), &
          self%domain%reference_temperature_k)) then
        reference_index = point_index
      end if
    end do

    if (reference_index == 0) then
      error stop "Reference sıcaklığı tabulated shift dataset'inde bulunmalıdır."
    end if
    if (.not. are_temperature_shift_values_machine_equivalent( &
        self%log10_a_t_points(reference_index), 0.0_dp)) then
      error stop "Reference noktada log10(a_T) machine tolerance içinde sıfır olmalıdır."
    end if
  end subroutine validate_tabulated_temperature_shift

  !> Tabulated provider'ın kapalı T_min/T_max ve T_ref [K] domain'ini döndürür.
  pure function get_tabulated_temperature_domain(self) result(domain)
    class(tabulated_log10_shift_provider_t), intent(in) :: self
    type(temperature_shift_domain_t) :: domain

    call self%validate()
    domain = self%domain
  end function get_tabulated_temperature_domain

  !> Private temperature [K] noktalarının bağımsız kopyasını döndürür.
  pure function get_tabulated_shift_temperature_points(provider) result(points)
    type(tabulated_log10_shift_provider_t), intent(in) :: provider
    real(dp), allocatable :: points(:)

    call provider%validate()
    points = provider%temperature_points_k
  end function get_tabulated_shift_temperature_points

  !> Private s=log10(a_T) [-] noktalarının bağımsız kopyasını döndürür.
  pure function get_tabulated_log10_shift_points(provider) result(points)
    type(tabulated_log10_shift_provider_t), intent(in) :: provider
    real(dp), allocatable :: points(:)

    call provider%validate()
    points = provider%log10_a_t_points
  end function get_tabulated_log10_shift_points

end module tms_tabulated_temperature_shift
