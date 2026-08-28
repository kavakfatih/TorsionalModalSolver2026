module tms_dynamic_torsional_property_binding
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_geometry, only : rubber_geometry_t, &
    calculate_annular_bush_torsion_geometry_factor
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    DYNAMIC_DEFORMATION_MODE_SHEAR
  use tms_dynamic_modulus_provider, only : dynamic_modulus_provider_t, &
    dynamic_modulus_evaluation_t, validate_dynamic_modulus_provider, &
    get_dynamic_modulus_provider_metadata
  use tms_dynamic_torsional_stiffness, only : &
    complex_torsional_stiffness_t, &
    calculate_dynamic_torsional_stiffness_from_geometry_factor
  implicit none
  private

  !> Bir provider sorgusunun eleman-level G* -> K* mapping sonucunu ve
  !! interpolation izini birlikte taşır.
  type, public :: dynamic_torsional_property_state_t
    integer :: element_id = 0
    !> Geriye uyumlu physical query frequency [Hz] alanı.
    real(dp) :: frequency_hz = 0.0_dp
    !> Explicit physical excitation frequency [Hz].
    real(dp) :: physical_frequency_hz = 0.0_dp
    !> Master curve üzerinde kullanılan reduced/lookup frequency [Hz].
    real(dp) :: lookup_frequency_hz = 0.0_dp
    !> Geriye uyumlu operating temperature [K] alanı.
    real(dp) :: temperature_k = 0.0_dp
    !> Explicit externally prescribed operating temperature [K].
    real(dp) :: operating_temperature_k = 0.0_dp
    real(dp) :: storage_modulus_pa = 0.0_dp
    real(dp) :: loss_modulus_pa = 0.0_dp
    real(dp) :: loss_factor = 0.0_dp
    real(dp) :: storage_stiffness_nm_per_rad = 0.0_dp
    real(dp) :: loss_stiffness_nm_per_rad = 0.0_dp
    integer :: interpolation_policy = 0
    logical :: exact_table_point = .false.
    real(dp) :: lower_frequency_hz = 0.0_dp
    real(dp) :: upper_frequency_hz = 0.0_dp
    real(dp) :: interpolation_alpha = 0.0_dp
    logical :: temperature_shift_applied = .false.
    integer :: shift_model_kind = 0
    real(dp) :: reference_temperature_k = 0.0_dp
    real(dp) :: log10_a_t = 0.0_dp
    real(dp) :: a_t = 1.0_dp
    logical :: has_temperature_bracket = .false.
    logical :: shift_exact_temperature_point = .false.
    real(dp) :: lower_temperature_k = 0.0_dp
    real(dp) :: upper_temperature_k = 0.0_dp
    real(dp) :: temperature_interpolation_alpha = 0.0_dp
  end type dynamic_torsional_property_state_t

  !> Bir torsional elemanı dinamik shear-modulus provider'ına bağlar.
  !! Provider ve önceden hesaplanmış bonded-annular C_theta private kopya
  !! olarak saklanır; material tablosu torsional_element_t içine gömülmez.
  type, public :: dynamic_torsional_property_binding_t
    private
    integer :: element_id = 0
    real(dp) :: geometry_factor_m3 = 0.0_dp
    class(dynamic_modulus_provider_t), allocatable :: provider
  end type dynamic_torsional_property_binding_t

  public :: create_dynamic_torsional_property_binding
  public :: validate_dynamic_torsional_property_binding
  public :: evaluate_dynamic_torsional_property
  public :: get_dynamic_binding_element_id
  public :: get_dynamic_binding_geometry_factor
  public :: get_dynamic_binding_metadata

contains

  !> Eleman-provider bağını oluşturur ve frequency-independent annular
  !! geometri katsayısını yalnız bir kez hesaplar.
  !! Fiziksel model: C_theta=4*pi*L*ri^2*ro^2/(ro^2-ri^2) [m^3], rijit hub
  !! ve ring arasındaki homojen bonded annular elastomer içindir.
  !! Girdiler: Pozitif element ID [-], direct dynamic SHEAR provider ve
  !! ri/ro/L [m] geometri. Çıktı: Bağımsız provider kopyası taşıyan binding.
  !! Varsayımlar ve sınırlar: TENSILE mode torsional shear mapping'e sessizce
  !! dönüştürülmez. E->G dönüşümü, distributed inertia ve bond compliance yoktur.
  function create_dynamic_torsional_property_binding( &
      element_id, provider, rubber) result(binding)
    integer, intent(in) :: element_id
    class(dynamic_modulus_provider_t), intent(in) :: provider
    type(rubber_geometry_t), intent(in) :: rubber
    type(dynamic_torsional_property_binding_t) :: binding

    type(dynamic_material_metadata_t) :: metadata

    if (element_id <= 0) then
      error stop "Dynamic torsional binding eleman kimliği pozitif olmalıdır."
    end if
    call validate_dynamic_modulus_provider(provider)
    metadata = get_dynamic_modulus_provider_metadata(provider)
    if (metadata%deformation_mode /= DYNAMIC_DEFORMATION_MODE_SHEAR) then
      error stop "Direct torsional material binding yalnız SHEAR dataset kabul eder."
    end if

    binding%element_id = element_id
    binding%geometry_factor_m3 = &
      calculate_annular_bush_torsion_geometry_factor( &
        inner_radius=rubber%inner_radius_m, &
        outer_radius=rubber%outer_radius_m, &
        axial_length=rubber%axial_length_m)
    allocate(binding%provider, source=provider)
    call validate_dynamic_torsional_property_binding(binding)
  end function create_dynamic_torsional_property_binding

  !> Binding'in kimlik, geometri, provider ve direct-SHEAR invariantlarını
  !! doğrular. C_theta [m^3] sonlu ve pozitif olmalıdır.
  pure subroutine validate_dynamic_torsional_property_binding(binding)
    type(dynamic_torsional_property_binding_t), intent(in) :: binding

    type(dynamic_material_metadata_t) :: metadata

    if (binding%element_id <= 0) then
      error stop "Dynamic torsional binding eleman kimliği pozitif olmalıdır."
    end if
    if (.not. ieee_is_finite(binding%geometry_factor_m3) .or. &
        binding%geometry_factor_m3 <= 0.0_dp) then
      error stop "Dynamic torsional binding C_theta değeri geçersiz."
    end if
    if (.not. allocated(binding%provider)) then
      error stop "Dynamic torsional binding provider içermelidir."
    end if
    call validate_dynamic_modulus_provider(binding%provider)
    metadata = get_dynamic_modulus_provider_metadata(binding%provider)
    if (metadata%deformation_mode /= DYNAMIC_DEFORMATION_MODE_SHEAR) then
      error stop "Dynamic torsional binding SHEAR dataset gerektirir."
    end if
  end subroutine validate_dynamic_torsional_property_binding

  !> Binding'in provider durumunu sorgular ve bonded-annular K* katsayılarına
  !! dönüştürür.
  !! Matematiksel model: G*=G'+iG'', K'=C_theta*G', K''=C_theta*G'' ve
  !! tan(delta)=G''/G'=K''/K'. Girdiler physical f [Hz] ile externally
  !! prescribed operating T [K]; çıktı G'/G'' [Pa], K'/K'' [N*m/rad] ve
  !! provider interpolation/shift trace'idir.
  !! Varsayımlar ve sınırlar: Provider kendi validated domain'leri dışında
  !! extrapolation yapmaz; pasif modelde G'>0, G''>=0, K'>0 ve K''>=0 korunur.
  !! Ayrıntılar: docs/architecture/V0.8_thermorheological_runtime.md.
  pure function evaluate_dynamic_torsional_property( &
      binding, frequency_hz, temperature_k) result(state)
    type(dynamic_torsional_property_binding_t), intent(in) :: binding
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: temperature_k
    type(dynamic_torsional_property_state_t) :: state

    type(dynamic_modulus_evaluation_t) :: modulus_evaluation
    type(complex_torsional_stiffness_t) :: stiffness

    call validate_dynamic_torsional_property_binding(binding)
    modulus_evaluation = binding%provider%evaluate( &
      frequency_hz, temperature_k)
    stiffness = calculate_dynamic_torsional_stiffness_from_geometry_factor( &
      modulus_evaluation%modulus, binding%geometry_factor_m3)

    state%element_id = binding%element_id
    state%frequency_hz = modulus_evaluation%modulus%frequency
    state%physical_frequency_hz = &
      modulus_evaluation%physical_frequency_hz
    state%lookup_frequency_hz = modulus_evaluation%lookup_frequency_hz
    state%temperature_k = modulus_evaluation%modulus%temperature
    state%operating_temperature_k = &
      modulus_evaluation%modulus%temperature
    state%storage_modulus_pa = &
      modulus_evaluation%modulus%storage_modulus
    state%loss_modulus_pa = modulus_evaluation%modulus%loss_modulus
    state%loss_factor = stiffness%loss_factor
    state%storage_stiffness_nm_per_rad = stiffness%storage_stiffness
    state%loss_stiffness_nm_per_rad = stiffness%loss_stiffness
    state%interpolation_policy = modulus_evaluation%interpolation_policy
    state%exact_table_point = modulus_evaluation%exact_table_point
    state%lower_frequency_hz = modulus_evaluation%lower_frequency_hz
    state%upper_frequency_hz = modulus_evaluation%upper_frequency_hz
    state%interpolation_alpha = modulus_evaluation%interpolation_alpha
    state%temperature_shift_applied = &
      modulus_evaluation%temperature_shift_applied
    state%shift_model_kind = modulus_evaluation%shift_model_kind
    state%reference_temperature_k = &
      modulus_evaluation%reference_temperature_k
    state%log10_a_t = modulus_evaluation%log10_a_t
    state%a_t = modulus_evaluation%a_t
    state%has_temperature_bracket = &
      modulus_evaluation%has_temperature_bracket
    state%shift_exact_temperature_point = &
      modulus_evaluation%shift_exact_temperature_point
    state%lower_temperature_k = modulus_evaluation%lower_temperature_k
    state%upper_temperature_k = modulus_evaluation%upper_temperature_k
    state%temperature_interpolation_alpha = &
      modulus_evaluation%temperature_interpolation_alpha
  end function evaluate_dynamic_torsional_property

  !> Binding'in authoritative torsional eleman kimliğini [-] döndürür.
  pure function get_dynamic_binding_element_id(binding) result(element_id)
    type(dynamic_torsional_property_binding_t), intent(in) :: binding
    integer :: element_id

    call validate_dynamic_torsional_property_binding(binding)
    element_id = binding%element_id
  end function get_dynamic_binding_element_id

  !> Bir kez hazırlanmış bonded-annular C_theta [m^3] değerini döndürür.
  pure function get_dynamic_binding_geometry_factor(binding) &
      result(geometry_factor_m3)
    type(dynamic_torsional_property_binding_t), intent(in) :: binding
    real(dp) :: geometry_factor_m3

    call validate_dynamic_torsional_property_binding(binding)
    geometry_factor_m3 = binding%geometry_factor_m3
  end function get_dynamic_binding_geometry_factor

  !> Binding provider'ının dataset-level metadata kopyasını döndürür.
  pure function get_dynamic_binding_metadata(binding) result(metadata)
    type(dynamic_torsional_property_binding_t), intent(in) :: binding
    type(dynamic_material_metadata_t) :: metadata

    call validate_dynamic_torsional_property_binding(binding)
    metadata = get_dynamic_modulus_provider_metadata(binding%provider)
  end function get_dynamic_binding_metadata

end module tms_dynamic_torsional_property_binding
