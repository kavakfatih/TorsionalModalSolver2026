module tms_wlf_temperature_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_temperature_shift_types, only : temperature_shift_domain_t, &
    temperature_shift_evaluation_t, WLF_TEMPERATURE_SHIFT
  use tms_temperature_shift_provider, only : temperature_shift_provider_t, &
    validate_temperature_shift_domain, temperature_is_in_shift_domain, &
    are_temperature_shift_values_machine_equivalent, &
    calculate_shift_factor_from_log10
  implicit none
  private

  !> Williams-Landel-Ferry sıcaklık kaydırma modelidir. TMS26 convention'ı
  !! s=log10(a_T)=-C1(T-T_ref)/(C2+T-T_ref) biçimindedir. C1 boyutsuz,
  !! C2 ve tüm sıcaklıklar K'dir; validated domain WLF pole'ünü geçemez.
  type, extends(temperature_shift_provider_t), public :: &
      wlf_temperature_shift_provider_t
    private
    real(dp) :: c1 = 0.0_dp
    real(dp) :: c2_k = 0.0_dp
    type(temperature_shift_domain_t) :: domain
  contains
    procedure, public :: evaluate => evaluate_wlf_temperature_shift
    procedure, public :: validate => validate_wlf_temperature_shift
    procedure, public :: get_domain => get_wlf_temperature_domain
  end type wlf_temperature_shift_provider_t

  public :: create_wlf_temperature_shift_provider

contains

  !> C1>0 [-], C2>0 [K], T_ref [K] ve explicit [T_min,T_max] [K] ile WLF
  !! provider oluşturur. Domain boyunca C2+(T-T_ref)>0 zorunludur; böylece
  !! WLF pole veya pole'ün karşı tarafındaki fiziksel olmayan kullanım reddedilir.
  pure function create_wlf_temperature_shift_provider( &
      c1, c2_k, reference_temperature_k, minimum_temperature_k, &
      maximum_temperature_k) result(provider)
    real(dp), intent(in) :: c1
    real(dp), intent(in) :: c2_k
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: minimum_temperature_k
    real(dp), intent(in) :: maximum_temperature_k
    type(wlf_temperature_shift_provider_t) :: provider

    provider%c1 = c1
    provider%c2_k = c2_k
    provider%domain = temperature_shift_domain_t( &
      minimum_temperature_k=minimum_temperature_k, &
      maximum_temperature_k=maximum_temperature_k, &
      reference_temperature_k=reference_temperature_k)
    call provider%validate()
  end function create_wlf_temperature_shift_provider

  !> Operating T [K] için WLF denklemiyle s=log10(a_T) hesaplar:
  !! s=-C1*dT/(C2+dT), a_T=10^s. Çıktılar boyutsuzdur. Model conventional
  !! positive C1/C2, thermorheologically-simple yatay kaydırma ve externally
  !! prescribed sıcaklık varsayar; validated domain dışında extrapolation yoktur.
  pure function evaluate_wlf_temperature_shift(self, temperature_k) &
      result(evaluation)
    class(wlf_temperature_shift_provider_t), intent(in) :: self
    real(dp), intent(in) :: temperature_k
    type(temperature_shift_evaluation_t) :: evaluation

    real(dp) :: denominator_k
    real(dp) :: temperature_difference_k

    call self%validate()
    if (.not. temperature_is_in_shift_domain(temperature_k, self%domain)) then
      error stop "WLF operating sıcaklığı validated domain dışında."
    end if

    temperature_difference_k = temperature_k - &
      self%domain%reference_temperature_k
    denominator_k = self%c2_k + temperature_difference_k
    if (.not. ieee_is_finite(denominator_k) .or. &
        denominator_k <= 0.0_dp) then
      error stop "WLF denominator pozitif ve sonlu olmalıdır."
    end if

    evaluation%shift_model_kind = WLF_TEMPERATURE_SHIFT
    evaluation%operating_temperature_k = temperature_k
    evaluation%reference_temperature_k = &
      self%domain%reference_temperature_k
    if (are_temperature_shift_values_machine_equivalent( &
        temperature_k, self%domain%reference_temperature_k)) then
      evaluation%log10_a_t = 0.0_dp
      evaluation%a_t = 1.0_dp
      return
    end if
    evaluation%log10_a_t = &
      (-temperature_difference_k / denominator_k) * self%c1
    if (.not. ieee_is_finite(evaluation%log10_a_t)) then
      error stop "WLF log10(a_T) sonucu sonlu olmalıdır."
    end if
    evaluation%a_t = calculate_shift_factor_from_log10( &
      evaluation%log10_a_t)
  end function evaluate_wlf_temperature_shift

  !> WLF parametre ve domain invariantlarını doğrular. C1>0 [-], C2>0 [K],
  !! geçerli mutlak sıcaklık domain'i ve domain boyunca D(T)=C2+T-T_ref>0
  !! gerekir. Endpoint sonuçlarının sonlu olması da kontrol edilir.
  pure subroutine validate_wlf_temperature_shift(self)
    class(wlf_temperature_shift_provider_t), intent(in) :: self

    real(dp) :: maximum_denominator_k
    real(dp) :: minimum_denominator_k
    real(dp) :: pole_separation_scale_k
    real(dp) :: shift_at_maximum
    real(dp) :: shift_at_minimum

    if (.not. ieee_is_finite(self%c1) .or. self%c1 <= 0.0_dp) then
      error stop "WLF C1 katsayısı sonlu ve pozitif olmalıdır."
    end if
    if (.not. ieee_is_finite(self%c2_k) .or. self%c2_k <= 0.0_dp) then
      error stop "WLF C2 katsayısı sonlu ve pozitif [K] olmalıdır."
    end if
    call validate_temperature_shift_domain(self%domain)

    minimum_denominator_k = self%c2_k + &
      (self%domain%minimum_temperature_k - &
        self%domain%reference_temperature_k)
    maximum_denominator_k = self%c2_k + &
      (self%domain%maximum_temperature_k - &
        self%domain%reference_temperature_k)
    pole_separation_scale_k = max( &
      abs(self%c2_k), &
      abs(self%domain%minimum_temperature_k- &
        self%domain%reference_temperature_k), tiny(1.0_dp))
    if (.not. ieee_is_finite(minimum_denominator_k) .or. &
        .not. ieee_is_finite(maximum_denominator_k) .or. &
        minimum_denominator_k <= 0.0_dp .or. &
        minimum_denominator_k <= &
          64.0_dp*epsilon(1.0_dp)*pole_separation_scale_k .or. &
        maximum_denominator_k <= 0.0_dp) then
      error stop "Validated WLF domain pole'ü içeremez veya geçemez."
    end if

    shift_at_minimum = (-(&
      self%domain%minimum_temperature_k - &
      self%domain%reference_temperature_k) / minimum_denominator_k) * self%c1
    shift_at_maximum = (-(&
      self%domain%maximum_temperature_k - &
      self%domain%reference_temperature_k) / maximum_denominator_k) * self%c1
    if (.not. ieee_is_finite(shift_at_minimum) .or. &
        .not. ieee_is_finite(shift_at_maximum)) then
      error stop "WLF domain endpoint log10(a_T) değerleri sonlu olmalıdır."
    end if
  end subroutine validate_wlf_temperature_shift

  !> WLF provider'ın explicit T_min/T_max ve T_ref [K] domain'ini döndürür.
  pure function get_wlf_temperature_domain(self) result(domain)
    class(wlf_temperature_shift_provider_t), intent(in) :: self
    type(temperature_shift_domain_t) :: domain

    call self%validate()
    domain = self%domain
  end function get_wlf_temperature_domain

end module tms_wlf_temperature_shift
