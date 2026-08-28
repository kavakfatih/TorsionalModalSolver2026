module tms_temperature_shift_provider
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_temperature_shift_types, only : temperature_shift_domain_t, &
    temperature_shift_evaluation_t, TABULATED_LOG10_SHIFT, &
    WLF_TEMPERATURE_SHIFT, ARRHENIUS_TEMPERATURE_SHIFT
  implicit none
  private

  !> Temperature-shift constitutive ayrıntısını thermorheological material
  !! provider'dan ayıran küçük arayüzdür. Concrete model yalnız T [K] girdisini
  !! a_T ve log10(a_T) boyutsuz sonuçlarına dönüştürür; frequency shifting veya
  !! master-curve lookup bu katmanın sorumluluğunda değildir.
  type, abstract, public :: temperature_shift_provider_t
  contains
    procedure(shift_evaluate_interface), deferred, public :: evaluate
    procedure(shift_validate_interface), deferred, public :: validate
    procedure(shift_domain_interface), deferred, public :: get_domain
  end type temperature_shift_provider_t

  abstract interface
    !> Operating T [K] için authoritative log10(a_T) ve türetilmiş a_T [-]
    !! değerlerini döndüren saf arayüzdür. Query, provider'ın doğrulanmış
    !! sıcaklık domain'i dışında extrapolation yapamaz.
    pure function shift_evaluate_interface(self, temperature_k) &
        result(evaluation)
      import :: dp, temperature_shift_evaluation_t, &
        temperature_shift_provider_t
      class(temperature_shift_provider_t), intent(in) :: self
      real(dp), intent(in) :: temperature_k
      type(temperature_shift_evaluation_t) :: evaluation
    end function shift_evaluate_interface

    !> Concrete provider parametre ve sıcaklık domain'ini doğrulayan saf
    !! arayüzdür. Fiziksel hesap sonucu üretmez.
    pure subroutine shift_validate_interface(self)
      import :: temperature_shift_provider_t
      class(temperature_shift_provider_t), intent(in) :: self
    end subroutine shift_validate_interface

    !> Minimum, maksimum ve reference sıcaklıklarını [K] bağımsız kopya
    !! halinde döndüren saf arayüzdür.
    pure function shift_domain_interface(self) result(domain)
      import :: temperature_shift_domain_t, temperature_shift_provider_t
      class(temperature_shift_provider_t), intent(in) :: self
      type(temperature_shift_domain_t) :: domain
    end function shift_domain_interface
  end interface

  public :: evaluate_temperature_shift
  public :: validate_temperature_shift_provider
  public :: get_temperature_shift_domain
  public :: validate_temperature_shift_evaluation
  public :: validate_temperature_shift_domain
  public :: temperature_is_in_shift_domain
  public :: are_temperature_shift_values_machine_equivalent
  public :: calculate_shift_factor_from_log10

contains

  !> Concrete provider'ı operating T [K] değerinde sorgular. Çıktı
  !! authoritative log10(a_T), derived a_T [-] ve model izidir. TMS26 yatay
  !! kaydırma convention'ı f_r=a_T*f bu sonucun tüketildiği katmanda uygulanır.
  pure function evaluate_temperature_shift(provider, temperature_k) &
      result(evaluation)
    class(temperature_shift_provider_t), intent(in) :: provider
    real(dp), intent(in) :: temperature_k
    type(temperature_shift_evaluation_t) :: evaluation

    evaluation = provider%evaluate(temperature_k)
    call validate_temperature_shift_evaluation( &
      provider, temperature_k, evaluation)
  end function evaluate_temperature_shift

  !> Concrete provider'ın model parametreleri ve explicit T domain'ini
  !! doğrulatır. Girdi provider, çıktı yoktur; geçersiz durum error stop üretir.
  pure subroutine validate_temperature_shift_provider(provider)
    class(temperature_shift_provider_t), intent(in) :: provider

    call provider%validate()
  end subroutine validate_temperature_shift_provider

  !> Provider'ın T_min, T_max ve T_ref [K] değerlerini bağımsız kopya olarak
  !! döndürür. Fonksiyon constitutive kaydırma hesabı yapmaz.
  pure function get_temperature_shift_domain(provider) result(domain)
    class(temperature_shift_provider_t), intent(in) :: provider
    type(temperature_shift_domain_t) :: domain

    domain = provider%get_domain()
  end function get_temperature_shift_domain

  !> Concrete provider evaluation'ının ortak runtime sözleşmesini doğrular.
  !! Girdi query/operating/reference sıcaklıkları [K], s=log10(a_T) ve a_T
  !! boyutsuzdur; tabulated bracket sıcaklıkları [K], alpha boyutsuzdur.
  !! Model: a_T=10^s ve T_ref identity'sinde s=0, a_T=1. Bu yordam shift
  !! denklemini yeniden hesaplamaz; provider sınırından çıkan iz tutarlılığını
  !! lookup yapılmadan önce doğrular.
  pure subroutine validate_temperature_shift_evaluation( &
      provider, temperature_k, evaluation)
    class(temperature_shift_provider_t), intent(in) :: provider
    real(dp), intent(in) :: temperature_k
    type(temperature_shift_evaluation_t), intent(in) :: evaluation

    type(temperature_shift_domain_t) :: domain
    real(dp) :: expected_alpha

    call provider%validate()
    domain = provider%get_domain()
    call validate_temperature_shift_domain(domain)
    if (.not. temperature_is_in_shift_domain(temperature_k, domain)) then
      error stop "Temperature-shift evaluation sorgusu validated domain dışında."
    end if
    if (.not. ieee_is_finite(evaluation%operating_temperature_k) .or. &
        .not. are_temperature_shift_values_machine_equivalent( &
          evaluation%operating_temperature_k, temperature_k)) then
      error stop "Temperature-shift evaluation operating sıcaklığı hatalı."
    end if
    if (.not. ieee_is_finite(evaluation%reference_temperature_k) .or. &
        .not. are_temperature_shift_values_machine_equivalent( &
          evaluation%reference_temperature_k, &
          domain%reference_temperature_k)) then
      error stop "Temperature-shift evaluation reference sıcaklığı hatalı."
    end if
    if (.not. ieee_is_finite(evaluation%log10_a_t) .or. &
        .not. ieee_is_finite(evaluation%a_t) .or. &
        evaluation%a_t <= 0.0_dp .or. &
        .not. are_temperature_shift_values_machine_equivalent( &
          log10(evaluation%a_t), evaluation%log10_a_t)) then
      error stop "Temperature-shift evaluation s ve a_T değerleri tutarsız."
    end if
    if (are_temperature_shift_values_machine_equivalent( &
        temperature_k, domain%reference_temperature_k)) then
      if (.not. are_temperature_shift_values_machine_equivalent( &
          evaluation%log10_a_t, 0.0_dp) .or. &
          .not. are_temperature_shift_values_machine_equivalent( &
            evaluation%a_t, 1.0_dp)) then
        error stop "Temperature-shift evaluation T_ref identity'sini sağlamıyor."
      end if
    end if

    select case (evaluation%shift_model_kind)
    case (TABULATED_LOG10_SHIFT)
      if (.not. evaluation%has_temperature_bracket .or. &
          .not. ieee_is_finite(evaluation%lower_temperature_k) .or. &
          .not. ieee_is_finite(evaluation%upper_temperature_k) .or. &
          .not. ieee_is_finite(evaluation%interpolation_alpha) .or. &
          evaluation%lower_temperature_k <= 0.0_dp .or. &
          evaluation%upper_temperature_k < &
            evaluation%lower_temperature_k .or. &
          evaluation%interpolation_alpha < 0.0_dp .or. &
          evaluation%interpolation_alpha > 1.0_dp) then
        error stop "Tabulated temperature-shift evaluation bracket'i geçersiz."
      end if
      if (evaluation%exact_temperature_point) then
        if (.not. are_temperature_shift_values_machine_equivalent( &
            evaluation%lower_temperature_k, &
            evaluation%upper_temperature_k) .or. &
            .not. are_temperature_shift_values_machine_equivalent( &
              evaluation%operating_temperature_k, &
              evaluation%lower_temperature_k) .or. &
            .not. are_temperature_shift_values_machine_equivalent( &
              evaluation%interpolation_alpha, 0.0_dp)) then
          error stop "Exact tabulated temperature-shift evaluation hatalı."
        end if
      else
        if (evaluation%lower_temperature_k >= &
            evaluation%upper_temperature_k .or. &
            evaluation%operating_temperature_k <= &
              evaluation%lower_temperature_k .or. &
            evaluation%operating_temperature_k >= &
              evaluation%upper_temperature_k) then
          error stop "Interpolated temperature-shift evaluation bracket dışında."
        end if
        expected_alpha = (evaluation%operating_temperature_k - &
          evaluation%lower_temperature_k) / &
          (evaluation%upper_temperature_k - &
            evaluation%lower_temperature_k)
        if (.not. are_temperature_shift_values_machine_equivalent( &
            evaluation%interpolation_alpha, expected_alpha)) then
          error stop "Temperature-shift interpolation alpha değeri tutarsız."
        end if
      end if
    case (WLF_TEMPERATURE_SHIFT, ARRHENIUS_TEMPERATURE_SHIFT)
      if (evaluation%has_temperature_bracket .or. &
          evaluation%exact_temperature_point .or. &
          .not. are_temperature_shift_values_machine_equivalent( &
            evaluation%lower_temperature_k, 0.0_dp) .or. &
          .not. are_temperature_shift_values_machine_equivalent( &
            evaluation%upper_temperature_k, 0.0_dp) .or. &
          .not. are_temperature_shift_values_machine_equivalent( &
            evaluation%interpolation_alpha, 0.0_dp)) then
        error stop "Analytical temperature-shift evaluation table bracket taşıyor."
      end if
    case default
      error stop "Temperature-shift evaluation model kimliği geçersiz."
    end select
  end subroutine validate_temperature_shift_evaluation

  !> Temperature-shift kullanım domain'ini doğrular. T_min, T_max ve T_ref
  !! sonlu ve pozitif [K], T_max>T_min ve T_ref kapalı domain içinde olmalıdır.
  !! Geçersiz girdi extrapolation'a veya fiziksel olmayan mutlak sıcaklığa yol
  !! açacağından error stop ile reddedilir.
  pure subroutine validate_temperature_shift_domain(domain)
    type(temperature_shift_domain_t), intent(in) :: domain

    if (.not. ieee_is_finite(domain%minimum_temperature_k) .or. &
        domain%minimum_temperature_k <= 0.0_dp) then
      error stop "Temperature-shift minimum sıcaklığı sonlu ve pozitif olmalıdır."
    end if
    if (.not. ieee_is_finite(domain%maximum_temperature_k) .or. &
        domain%maximum_temperature_k <= domain%minimum_temperature_k) then
      error stop "Temperature-shift maksimum sıcaklığı minimumdan büyük olmalıdır."
    end if
    if (.not. ieee_is_finite(domain%reference_temperature_k) .or. &
        domain%reference_temperature_k <= 0.0_dp) then
      error stop "Temperature-shift reference sıcaklığı sonlu ve pozitif olmalıdır."
    end if
    if (.not. temperature_is_in_shift_domain( &
        domain%reference_temperature_k, domain)) then
      error stop "Temperature-shift reference sıcaklığı validated domain dışında."
    end if
  end subroutine validate_temperature_shift_domain

  !> T [K] değerinin kapalı validated domain içinde olup olmadığını machine
  !! representation toleransıyla sınar. Bu tolerans deney belirsizliği veya
  !! extrapolation izni değildir; yalnız endpoint round-off'unu ayırır.
  pure function temperature_is_in_shift_domain(temperature_k, domain) &
      result(is_in_domain)
    real(dp), intent(in) :: temperature_k
    type(temperature_shift_domain_t), intent(in) :: domain
    logical :: is_in_domain

    is_in_domain = ieee_is_finite(temperature_k) .and. &
      temperature_k > 0.0_dp .and. &
      (temperature_k > domain%minimum_temperature_k .or. &
        are_temperature_shift_values_machine_equivalent( &
          temperature_k, domain%minimum_temperature_k)) .and. &
      (temperature_k < domain%maximum_temperature_k .or. &
        are_temperature_shift_values_machine_equivalent( &
          temperature_k, domain%maximum_temperature_k))
  end function temperature_is_in_shift_domain

  !> Aynı fiziksel büyüklüğü temsil eden iki sonlu sayının yalnız floating-
  !! point düzeyinde eşdeğerliğini sınar. Model:
  !! |a-b|<=64 epsilon max(1,|a|,|b|). Girdiler aynı birimde, çıktı logical'dır.
  pure elemental function are_temperature_shift_values_machine_equivalent( &
      a, b) result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    real(dp) :: scale

    scale = max(1.0_dp, abs(a), abs(b))
    equivalent = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
      abs(a - b) <= 64.0_dp*epsilon(1.0_dp)*scale
  end function are_temperature_shift_values_machine_equivalent

  !> Authoritative s=log10(a_T) değerinden boyutsuz a_T=10^s üretir.
  !! Girdi s boyutsuzdur; çıktı pozitif ve sonlu a_T'dir. Representable dp
  !! aralığını aşan derived değerler sessiz overflow/underflow yerine clean
  !! error üretir. Reduced-frequency domain kontrolü s üzerinden yapılmalıdır.
  pure function calculate_shift_factor_from_log10(log10_a_t) result(a_t)
    real(dp), intent(in) :: log10_a_t
    real(dp) :: a_t

    real(dp) :: maximum_exponent
    real(dp) :: minimum_exponent

    if (.not. ieee_is_finite(log10_a_t)) then
      error stop "log10(a_T) sonlu olmalıdır."
    end if
    maximum_exponent = log10(huge(1.0_dp))
    minimum_exponent = log10(tiny(1.0_dp))
    if (log10_a_t > maximum_exponent .or. &
        log10_a_t < minimum_exponent) then
      error stop "Derived a_T çalışma hassasiyetinin sonlu aralığını aşıyor."
    end if

    a_t = 10.0_dp**log10_a_t
    if (.not. ieee_is_finite(a_t) .or. a_t <= 0.0_dp) then
      error stop "Derived a_T sonlu ve pozitif olmalıdır."
    end if
  end function calculate_shift_factor_from_log10

end module tms_temperature_shift_provider
