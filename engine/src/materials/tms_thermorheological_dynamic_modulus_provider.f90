module tms_thermorheological_dynamic_modulus_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t
  use tms_dynamic_modulus_provider, only : dynamic_modulus_provider_t, &
    dynamic_modulus_evaluation_t, are_machine_equivalent
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t, get_tabulated_frequency_domain
  use tms_temperature_shift_types, only : temperature_shift_domain_t, &
    temperature_shift_evaluation_t
  use tms_temperature_shift_provider, only : temperature_shift_provider_t, &
    evaluate_temperature_shift
  implicit none
  private

  !> Doğrulanmış V0.7 reference master curve ile bağımsız temperature-shift
  !! modelini bileştirir. Bu concrete provider mevcut dynamic modulus provider
  !! sınırını uygular; harmonic solver veya torsional binding özel davranışı
  !! içermez.
  type, extends(dynamic_modulus_provider_t), public :: &
      thermorheological_dynamic_modulus_provider_t
    private
    type(tabulated_dynamic_modulus_provider_t) :: master_curve
    class(temperature_shift_provider_t), allocatable :: shift_provider
  contains
    procedure, public :: evaluate => evaluate_thermorheological_provider
    procedure, public :: get_metadata => get_thermorheological_metadata
    procedure, public :: validate => validate_thermorheological_provider
  end type thermorheological_dynamic_modulus_provider_t

  public :: create_thermorheological_dynamic_modulus_provider

contains

  !> Validated reference master curve ile shift provider'ın bağımsız
  !! kopyalarını birleştirir.
  !! Girdi master curve G'/G'' [Pa] ve f_r [Hz] tablosu; shift provider ise
  !! T [K] -> log10(a_T) dönüşümüdür. Master metadata sıcaklığı ile shift
  !! reference sıcaklığı machine-equivalent olmalıdır.
  !! Çıktı horizontal-only, linear-viscoelastic runtime provider'dır.
  function create_thermorheological_dynamic_modulus_provider( &
      master_curve, shift_provider) result(provider)
    type(tabulated_dynamic_modulus_provider_t), intent(in) :: master_curve
    class(temperature_shift_provider_t), intent(in) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: provider

    provider%master_curve = master_curve
    allocate(provider%shift_provider, source=shift_provider)
    call provider%validate()
  end function create_thermorheological_dynamic_modulus_provider

  !> Fiziksel f [Hz] ve externally prescribed T [K] için horizontal
  !! temperature-shifted G*(f,T)=G'+iG'' [Pa] durumunu hesaplar.
  !! Matematiksel model: s=log10(a_T), log10(f_r)=log10(f)+s ve master curve
  !! lookup'u f_r'dedir. Önce log-domain doğrulanır; a_T*f doğrudan
  !! hesaplanmaz. Returned modulus physical f,T'yi, interpolation trace ise
  !! reduced lookup f_r'yi taşır.
  !! Varsayımlar ve sınırlar: Thermorheological simplicity önceden
  !! doğrulanmıştır; yalnız horizontal shift vardır. Temperature/frequency
  !! extrapolation, vertical shift, self-heating ve state interpolation yoktur.
  pure function evaluate_thermorheological_provider( &
      self, frequency_hz, temperature_k) result(evaluation)
    class(thermorheological_dynamic_modulus_provider_t), intent(in) :: self
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: temperature_k
    type(dynamic_modulus_evaluation_t) :: evaluation

    type(dynamic_modulus_evaluation_t) :: master_evaluation
    type(temperature_shift_evaluation_t) :: shift_evaluation
    real(dp) :: log_lookup_frequency
    real(dp) :: log_maximum_frequency
    real(dp) :: log_minimum_frequency
    real(dp) :: lookup_frequency_hz
    real(dp) :: maximum_frequency_hz
    real(dp) :: minimum_frequency_hz

    call self%validate()
    if (.not. ieee_is_finite(frequency_hz) .or. frequency_hz <= 0.0_dp) then
      error stop "Thermorheological physical frequency sonlu ve pozitif olmalıdır."
    end if
    if (.not. ieee_is_finite(temperature_k) .or. temperature_k <= 0.0_dp) then
      error stop "Thermorheological operating temperature pozitif K olmalıdır."
    end if

    shift_evaluation = evaluate_temperature_shift( &
      self%shift_provider, temperature_k)
    call get_tabulated_frequency_domain( &
      self%master_curve, minimum_frequency_hz, maximum_frequency_hz)
    log_minimum_frequency = log10(minimum_frequency_hz)
    log_maximum_frequency = log10(maximum_frequency_hz)
    log_lookup_frequency = log10(frequency_hz)+shift_evaluation%log10_a_t
    if (.not. ieee_is_finite(log_lookup_frequency)) then
      error stop "Reduced frequency logarithması sonlu değildir."
    end if

    if (log_lookup_frequency < log_minimum_frequency .and. &
        .not. are_log_coordinates_equivalent( &
          log_lookup_frequency, log_minimum_frequency)) then
      error stop "Reduced frequency master-curve alt domain'inin dışındadır."
    end if
    if (log_lookup_frequency > log_maximum_frequency .and. &
        .not. are_log_coordinates_equivalent( &
          log_lookup_frequency, log_maximum_frequency)) then
      error stop "Reduced frequency master-curve üst domain'inin dışındadır."
    end if

    ! Endpoint equality yalnız floating-point representation seviyesinde
    ! canonicalize edilir; fiziksel endpoint clamp veya extrapolation değildir.
    if (are_log_coordinates_equivalent( &
        log_lookup_frequency, log_minimum_frequency)) then
      lookup_frequency_hz = minimum_frequency_hz
    else if (are_log_coordinates_equivalent( &
        log_lookup_frequency, log_maximum_frequency)) then
      lookup_frequency_hz = maximum_frequency_hz
    else if (are_machine_equivalent( &
        shift_evaluation%operating_temperature_k, &
        shift_evaluation%reference_temperature_k) .and. &
        abs(shift_evaluation%log10_a_t) <= 64.0_dp*epsilon(1.0_dp)) then
      lookup_frequency_hz = frequency_hz
    else
      lookup_frequency_hz = 10.0_dp**log_lookup_frequency
    end if
    if (.not. ieee_is_finite(lookup_frequency_hz) .or. &
        lookup_frequency_hz <= 0.0_dp) then
      error stop "Reduced frequency sonlu ve pozitif olarak oluşturulamadı."
    end if

    master_evaluation = self%master_curve%evaluate( &
      lookup_frequency_hz, shift_evaluation%reference_temperature_k)
    evaluation = master_evaluation
    evaluation%modulus%frequency = frequency_hz
    evaluation%modulus%temperature = temperature_k
    evaluation%physical_frequency_hz = frequency_hz
    evaluation%lookup_frequency_hz = lookup_frequency_hz
    evaluation%temperature_shift_applied = .true.
    evaluation%shift_model_kind = shift_evaluation%shift_model_kind
    evaluation%reference_temperature_k = &
      shift_evaluation%reference_temperature_k
    evaluation%log10_a_t = shift_evaluation%log10_a_t
    evaluation%a_t = shift_evaluation%a_t
    evaluation%has_temperature_bracket = &
      shift_evaluation%has_temperature_bracket
    evaluation%shift_exact_temperature_point = &
      shift_evaluation%exact_temperature_point
    evaluation%lower_temperature_k = &
      shift_evaluation%lower_temperature_k
    evaluation%upper_temperature_k = &
      shift_evaluation%upper_temperature_k
    evaluation%temperature_interpolation_alpha = &
      shift_evaluation%interpolation_alpha
  end function evaluate_thermorheological_provider

  !> Underlying master-curve/reference-state metadata'sının bağımsız
  !! kopyasını döndürür. dataset_temperature_k operating temperature değil,
  !! reference/master-curve sıcaklığı [K] anlamındadır.
  pure function get_thermorheological_metadata(self) result(metadata)
    class(thermorheological_dynamic_modulus_provider_t), intent(in) :: self
    type(dynamic_material_metadata_t) :: metadata

    call self%validate()
    metadata = self%master_curve%get_metadata()
  end function get_thermorheological_metadata

  !> Master curve, shift provider ve ortak reference-temperature
  !! invariantlarını doğrular. Fiziksel query veya interpolation yapmaz.
  pure subroutine validate_thermorheological_provider(self)
    class(thermorheological_dynamic_modulus_provider_t), intent(in) :: self

    type(dynamic_material_metadata_t) :: metadata
    type(temperature_shift_domain_t) :: shift_domain

    call self%master_curve%validate()
    if (.not. allocated(self%shift_provider)) then
      error stop "Thermorheological provider temperature-shift modeli içermiyor."
    end if
    call self%shift_provider%validate()
    metadata = self%master_curve%get_metadata()
    shift_domain = self%shift_provider%get_domain()
    if (.not. are_machine_equivalent( &
        metadata%dataset_temperature_k, &
        shift_domain%reference_temperature_k)) then
      error stop "Master curve ve shift reference sıcaklıkları eşleşmiyor."
    end if
  end subroutine validate_thermorheological_provider

  !> İki log10-frequency koordinatının yalnız machine representation
  !! seviyesinde eşdeğerliğini sınar. Girdiler log10(Hz), çıktı logical'dır;
  !! deneysel bandwidth toleransı değildir.
  pure elemental function are_log_coordinates_equivalent(a, b) &
      result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    equivalent = abs(a-b) <= &
      64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(a), abs(b))
  end function are_log_coordinates_equivalent

end module tms_thermorheological_dynamic_modulus_provider
