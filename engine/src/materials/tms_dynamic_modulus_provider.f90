module tms_dynamic_modulus_provider
  use tms_kinds, only : dp
  use tms_dynamic_modulus, only : dynamic_shear_modulus
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t
  use tms_temperature_shift_types, only : TEMPERATURE_SHIFT_NONE
  implicit none
  private

  integer, parameter, public :: LINEAR_FREQUENCY = 1
  integer, parameter, public :: LINEAR_LOG_FREQUENCY = 2

  !> Tek bir provider sorgusunun constitutive sonucu ile interpolation
  !! izini birlikte taşır. G' ve G'' primary değerlerdir; kayıp faktörü
  !! bunlardan ayrıca türetilir.
  type, public :: dynamic_modulus_evaluation_t
    type(dynamic_shear_modulus) :: modulus
    !> Harmonic excitation/query frekansı [Hz]. Returned modulus içindeki
    !! frequency alanıyla aynı fiziksel koordinattır.
    real(dp) :: physical_frequency_hz = 0.0_dp
    !> Constitutive tablo interpolation'ında kullanılan frekans [Hz].
    !! V0.7 provider'da physical frequency'ye, shifted provider'da f_r'ye eşittir.
    real(dp) :: lookup_frequency_hz = 0.0_dp
    integer :: interpolation_policy = 0
    logical :: exact_table_point = .false.
    !> Aşağıdaki bracket master-curve lookup frequency eksenindedir.
    real(dp) :: lower_frequency_hz = 0.0_dp
    real(dp) :: upper_frequency_hz = 0.0_dp
    real(dp) :: interpolation_alpha = 0.0_dp

    !> Temperature-shift context'i. Unshifted V0.7 için false, model NONE,
    !! log10(a_T)=0, a_T=1 ve lookup=physical anlamlı default'ları kullanılır.
    logical :: temperature_shift_applied = .false.
    integer :: shift_model_kind = TEMPERATURE_SHIFT_NONE
    real(dp) :: reference_temperature_k = 0.0_dp
    real(dp) :: log10_a_t = 0.0_dp
    real(dp) :: a_t = 1.0_dp
    logical :: has_temperature_bracket = .false.
    logical :: shift_exact_temperature_point = .false.
    real(dp) :: lower_temperature_k = 0.0_dp
    real(dp) :: upper_temperature_k = 0.0_dp
    real(dp) :: temperature_interpolation_alpha = 0.0_dp
  end type dynamic_modulus_evaluation_t

  !> Harmonic solver ile constitutive material model arasındaki küçük ve
  !! genişletilebilir provider sınırıdır. Solver tablo, master curve veya
  !! gelecekteki shift-model ayrıntısını bilmeden aynı sorgu sözleşmesini
  !! kullanır.
  type, abstract, public :: dynamic_modulus_provider_t
  contains
    procedure(provider_evaluate_interface), deferred, public :: evaluate
    procedure(provider_metadata_interface), deferred, public :: get_metadata
    procedure(provider_validate_interface), deferred, public :: validate
  end type dynamic_modulus_provider_t

  abstract interface
    !> f [Hz] ve T [K] çalışma noktasında G' ve G'' [Pa] ile interpolation
    !! izini döndüren saf provider arayüzüdür.
    pure function provider_evaluate_interface( &
        self, frequency_hz, temperature_k) result(evaluation)
      import :: dp, dynamic_modulus_evaluation_t, dynamic_modulus_provider_t
      class(dynamic_modulus_provider_t), intent(in) :: self
      real(dp), intent(in) :: frequency_hz
      real(dp), intent(in) :: temperature_k
      type(dynamic_modulus_evaluation_t) :: evaluation
    end function provider_evaluate_interface

    !> Dataset-level metadata'nın bağımsız kopyasını döndüren arayüzdür.
    pure function provider_metadata_interface(self) result(metadata)
      import :: dynamic_material_metadata_t, dynamic_modulus_provider_t
      class(dynamic_modulus_provider_t), intent(in) :: self
      type(dynamic_material_metadata_t) :: metadata
    end function provider_metadata_interface

    !> Provider iç durumunu doğrulayan saf arayüzdür.
    pure subroutine provider_validate_interface(self)
      import :: dynamic_modulus_provider_t
      class(dynamic_modulus_provider_t), intent(in) :: self
    end subroutine provider_validate_interface
  end interface

  public :: query_dynamic_shear_modulus
  public :: evaluate_dynamic_shear_modulus
  public :: get_dynamic_modulus_provider_metadata
  public :: validate_dynamic_modulus_provider
  public :: are_machine_equivalent

contains

  !> Provider'dan doğrudan constitutive G*(f,T)=G'(f,T)+iG''(f,T) durumunu alır.
  !! Girdiler physical frequency [Hz] ve operating temperature [K]; çıktı
  !! G', G'' [Pa] ile aynı physical f [Hz] ve T [K] alanlarıdır. Isotherm
  !! eşlemesi, temperature shift ve extrapolation politikası concrete provider
  !! sözleşmesinin sorumluluğundadır.
  pure function query_dynamic_shear_modulus( &
      provider, frequency_hz, temperature_k) result(modulus)
    class(dynamic_modulus_provider_t), intent(in) :: provider
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: temperature_k
    type(dynamic_shear_modulus) :: modulus

    type(dynamic_modulus_evaluation_t) :: evaluation

    evaluation = provider%evaluate(frequency_hz, temperature_k)
    modulus = evaluation%modulus
  end function query_dynamic_shear_modulus

  !> Provider sorgusunun G'/G'' sonucu ile exact/interpolated bracket izini
  !! birlikte döndürür. Fizik ve birimler query_dynamic_shear_modulus ile
  !! aynıdır; alpha boyutsuzdur.
  pure function evaluate_dynamic_shear_modulus( &
      provider, frequency_hz, temperature_k) result(evaluation)
    class(dynamic_modulus_provider_t), intent(in) :: provider
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: temperature_k
    type(dynamic_modulus_evaluation_t) :: evaluation

    evaluation = provider%evaluate(frequency_hz, temperature_k)
  end function evaluate_dynamic_shear_modulus

  !> Provider'ın dataset-level deney metadata'sının bağımsız kopyasını verir.
  pure function get_dynamic_modulus_provider_metadata(provider) &
      result(metadata)
    class(dynamic_modulus_provider_t), intent(in) :: provider
    type(dynamic_material_metadata_t) :: metadata

    metadata = provider%get_metadata()
  end function get_dynamic_modulus_provider_metadata

  !> Concrete provider'ın kendi veri ve operating-state invariantlarını
  !! doğrulatır. Fizik hesabı veya interpolation yapmaz.
  pure subroutine validate_dynamic_modulus_provider(provider)
    class(dynamic_modulus_provider_t), intent(in) :: provider

    call provider%validate()
  end subroutine validate_dynamic_modulus_provider

  !> İki sonlu fiziksel büyüklüğün yalnız floating-point representation
  !! düzeyinde eşdeğer olup olmadığını sınar.
  !! Matematiksel model: |a-b| <= 64*epsilon*max(|a|,|b|,tiny).
  !! Girdiler aynı SI birimindeki a ve b; çıktı boyutsuz logical değerdir.
  !! Bu tolerans deney belirsizliği, fiziksel bandwidth veya kullanıcı
  !! toleransı değildir; yalnız exact table-point ve isotherm eşlemesi içindir.
  pure elemental function are_machine_equivalent(a, b) result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    real(dp) :: scale

    scale = max(abs(a), abs(b), tiny(1.0_dp))
    equivalent = abs(a - b) <= 64.0_dp*epsilon(1.0_dp)*scale
  end function are_machine_equivalent

end module tms_dynamic_modulus_provider
