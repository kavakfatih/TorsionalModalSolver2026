module tms_material_state_trace
  use tms_kinds, only : dp
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t
  use tms_dynamic_torsional_property_binding, only : &
    dynamic_torsional_property_binding_t, &
    dynamic_torsional_property_state_t, &
    get_dynamic_binding_element_id, &
    get_dynamic_binding_geometry_factor, get_dynamic_binding_metadata
  implicit none
  private

  !> Sonuçta bir kez saklanan dataset/binding izidir. Büyük metadata her
  !! frequency point'te tekrarlanmaz.
  type, public :: material_binding_trace_t
    integer :: element_id = 0
    real(dp) :: geometry_factor_m3 = 0.0_dp
    type(dynamic_material_metadata_t) :: metadata
  end type material_binding_trace_t

  !> Tek dynamic binding ve requested frequency için constitutive/mapped
  !! malzeme durumunu taşır. Complex çözüm daha sonra singular olsa bile bu
  !! iz korunur.
  type, public :: material_state_trace_t
    integer :: element_id = 0
    !> Geriye uyumlu physical query frequency [Hz].
    real(dp) :: frequency_hz = 0.0_dp
    !> Explicit physical excitation frequency [Hz].
    real(dp) :: physical_frequency_hz = 0.0_dp
    !> Master-curve interpolation coordinate'i olan reduced frequency [Hz].
    real(dp) :: lookup_frequency_hz = 0.0_dp
    !> Geriye uyumlu externally prescribed operating temperature [K].
    real(dp) :: temperature_k = 0.0_dp
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
  end type material_state_trace_t

  public :: create_material_binding_trace
  public :: create_material_state_trace

contains

  !> Binding kimliği, frequency-independent C_theta [m^3] ve dataset
  !! metadata'sını sonuç-level bağımsız bir kayda kopyalar. Fizik hesabı yapmaz.
  pure function create_material_binding_trace(binding) result(trace)
    type(dynamic_torsional_property_binding_t), intent(in) :: binding
    type(material_binding_trace_t) :: trace

    trace%element_id = get_dynamic_binding_element_id(binding)
    trace%geometry_factor_m3 = get_dynamic_binding_geometry_factor(binding)
    trace%metadata = get_dynamic_binding_metadata(binding)
  end function create_material_binding_trace

  !> Önceden prevalidate edilmiş provider/binding durumunu immutable sonuç
  !! izine dönüştürür. G'/G'' [Pa], K'/K'' [N*m/rad], f [Hz], T [K] ve
  !! boyutsuz interpolation bilgisini aynen korur; yeni interpolation yapmaz.
  pure function create_material_state_trace(state) result(trace)
    type(dynamic_torsional_property_state_t), intent(in) :: state
    type(material_state_trace_t) :: trace

    trace%element_id = state%element_id
    trace%frequency_hz = state%frequency_hz
    trace%physical_frequency_hz = state%physical_frequency_hz
    trace%lookup_frequency_hz = state%lookup_frequency_hz
    trace%temperature_k = state%temperature_k
    trace%operating_temperature_k = state%operating_temperature_k
    trace%storage_modulus_pa = state%storage_modulus_pa
    trace%loss_modulus_pa = state%loss_modulus_pa
    trace%loss_factor = state%loss_factor
    trace%storage_stiffness_nm_per_rad = &
      state%storage_stiffness_nm_per_rad
    trace%loss_stiffness_nm_per_rad = state%loss_stiffness_nm_per_rad
    trace%interpolation_policy = state%interpolation_policy
    trace%exact_table_point = state%exact_table_point
    trace%lower_frequency_hz = state%lower_frequency_hz
    trace%upper_frequency_hz = state%upper_frequency_hz
    trace%interpolation_alpha = state%interpolation_alpha
    trace%temperature_shift_applied = state%temperature_shift_applied
    trace%shift_model_kind = state%shift_model_kind
    trace%reference_temperature_k = state%reference_temperature_k
    trace%log10_a_t = state%log10_a_t
    trace%a_t = state%a_t
    trace%has_temperature_bracket = state%has_temperature_bracket
    trace%shift_exact_temperature_point = &
      state%shift_exact_temperature_point
    trace%lower_temperature_k = state%lower_temperature_k
    trace%upper_temperature_k = state%upper_temperature_k
    trace%temperature_interpolation_alpha = &
      state%temperature_interpolation_alpha
  end function create_material_state_trace

end module tms_material_state_trace
