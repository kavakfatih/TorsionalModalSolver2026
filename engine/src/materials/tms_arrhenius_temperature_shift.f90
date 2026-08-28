module tms_arrhenius_temperature_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : universal_gas_constant_j_per_mol_k
  use tms_temperature_shift_types, only : temperature_shift_domain_t, &
    temperature_shift_evaluation_t, ARRHENIUS_TEMPERATURE_SHIFT
  use tms_temperature_shift_provider, only : temperature_shift_provider_t, &
    validate_temperature_shift_domain, temperature_is_in_shift_domain, &
    are_temperature_shift_values_machine_equivalent, &
    calculate_shift_factor_from_log10
  implicit none
  private

  !> Thermally activated Arrhenius temperature-shift modelidir. TMS26
  !! convention'ında s=log10(a_T)=Ea/[ln(10)R]*(1/T-1/T_ref).
  !! Ea [J/mol], R [J/(mol K)] ve sıcaklıklar [K] cinsindedir.
  type, extends(temperature_shift_provider_t), public :: &
      arrhenius_temperature_shift_provider_t
    private
    real(dp) :: activation_energy_j_per_mol = 0.0_dp
    type(temperature_shift_domain_t) :: domain
  contains
    procedure, public :: evaluate => evaluate_arrhenius_temperature_shift
    procedure, public :: validate => validate_arrhenius_temperature_shift
    procedure, public :: get_domain => get_arrhenius_temperature_domain
  end type arrhenius_temperature_shift_provider_t

  public :: create_arrhenius_temperature_shift_provider

contains

  !> Ea>0 [J/mol], T_ref [K] ve explicit [T_min,T_max] [K] ile Arrhenius
  !! provider oluşturur. Matematiksel formül domain dışında finite olsa bile
  !! extrapolation yapılmaz; kullanım aralığı constructor'da sabitlenir.
  pure function create_arrhenius_temperature_shift_provider( &
      activation_energy_j_per_mol, reference_temperature_k, &
      minimum_temperature_k, maximum_temperature_k) result(provider)
    real(dp), intent(in) :: activation_energy_j_per_mol
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: minimum_temperature_k
    real(dp), intent(in) :: maximum_temperature_k
    type(arrhenius_temperature_shift_provider_t) :: provider

    provider%activation_energy_j_per_mol = activation_energy_j_per_mol
    provider%domain = temperature_shift_domain_t( &
      minimum_temperature_k=minimum_temperature_k, &
      maximum_temperature_k=maximum_temperature_k, &
      reference_temperature_k=reference_temperature_k)
    call provider%validate()
  end function create_arrhenius_temperature_shift_provider

  !> Operating T [K] için
  !! s=Ea/[ln(10)R]*(1/T-1/T_ref) ve a_T=10^s hesaplar. Ea [J/mol],
  !! R [J/(mol K)], s ve a_T boyutsuzdur. Model tek activation energy,
  !! horizontal shifting ve externally prescribed sıcaklık varsayar.
  pure function evaluate_arrhenius_temperature_shift(self, temperature_k) &
      result(evaluation)
    class(arrhenius_temperature_shift_provider_t), intent(in) :: self
    real(dp), intent(in) :: temperature_k
    type(temperature_shift_evaluation_t) :: evaluation

    real(dp) :: inverse_temperature_difference

    call self%validate()
    if (.not. temperature_is_in_shift_domain(temperature_k, self%domain)) then
      error stop "Arrhenius operating sıcaklığı validated domain dışında."
    end if

    ! (T_ref-T)/(T*T_ref) cebirsel biçimi yerine ardışık bölme, ara
    ! T*T_ref taşmasını önler ve T=T_ref kimliğinde tam sıfır üretir.
    inverse_temperature_difference = &
      ((self%domain%reference_temperature_k - temperature_k) / &
        temperature_k) / self%domain%reference_temperature_k
    evaluation%shift_model_kind = ARRHENIUS_TEMPERATURE_SHIFT
    evaluation%operating_temperature_k = temperature_k
    evaluation%reference_temperature_k = &
      self%domain%reference_temperature_k
    if (are_temperature_shift_values_machine_equivalent( &
        temperature_k, self%domain%reference_temperature_k)) then
      evaluation%log10_a_t = 0.0_dp
      evaluation%a_t = 1.0_dp
      return
    end if
    evaluation%log10_a_t = inverse_temperature_difference * &
      (self%activation_energy_j_per_mol / &
        universal_gas_constant_j_per_mol_k) / log(10.0_dp)
    if (.not. ieee_is_finite(evaluation%log10_a_t)) then
      error stop "Arrhenius log10(a_T) sonucu sonlu olmalıdır."
    end if
    evaluation%a_t = calculate_shift_factor_from_log10( &
      evaluation%log10_a_t)
  end function evaluate_arrhenius_temperature_shift

  !> Arrhenius parameter ve domain invariantlarını doğrular. Ea sonlu ve
  !! pozitif [J/mol], bütün mutlak sıcaklıklar sonlu ve pozitif [K] olmalı;
  !! endpoint log10(a_T) değerleri finite kalmalıdır.
  pure subroutine validate_arrhenius_temperature_shift(self)
    class(arrhenius_temperature_shift_provider_t), intent(in) :: self

    real(dp) :: shift_at_maximum
    real(dp) :: shift_at_minimum

    if (.not. ieee_is_finite(self%activation_energy_j_per_mol) .or. &
        self%activation_energy_j_per_mol <= 0.0_dp) then
      error stop "Arrhenius activation energy sonlu ve pozitif [J/mol] olmalıdır."
    end if
    call validate_temperature_shift_domain(self%domain)

    shift_at_minimum = &
      (((self%domain%reference_temperature_k - &
        self%domain%minimum_temperature_k) / &
        self%domain%minimum_temperature_k) / &
        self%domain%reference_temperature_k) * &
      (self%activation_energy_j_per_mol / &
        universal_gas_constant_j_per_mol_k) / log(10.0_dp)
    shift_at_maximum = &
      (((self%domain%reference_temperature_k - &
        self%domain%maximum_temperature_k) / &
        self%domain%maximum_temperature_k) / &
        self%domain%reference_temperature_k) * &
      (self%activation_energy_j_per_mol / &
        universal_gas_constant_j_per_mol_k) / log(10.0_dp)
    if (.not. ieee_is_finite(shift_at_minimum) .or. &
        .not. ieee_is_finite(shift_at_maximum)) then
      error stop "Arrhenius domain endpoint log10(a_T) değerleri sonlu olmalıdır."
    end if
  end subroutine validate_arrhenius_temperature_shift

  !> Arrhenius provider'ın explicit T_min/T_max ve T_ref [K] domain'ini döndürür.
  pure function get_arrhenius_temperature_domain(self) result(domain)
    class(arrhenius_temperature_shift_provider_t), intent(in) :: self
    type(temperature_shift_domain_t) :: domain

    call self%validate()
    domain = self%domain
  end function get_arrhenius_temperature_domain

end module tms_arrhenius_temperature_shift
