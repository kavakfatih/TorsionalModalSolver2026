module tms_tts_runtime_export
  use tms_kinds, only : dp
  use tms_material_frequency, only : material_frequency_point
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    DYNAMIC_DEFORMATION_MODE_SHEAR
  use tms_dynamic_modulus_provider, only : LINEAR_LOG_FREQUENCY
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t, &
    create_tabulated_dynamic_modulus_provider
  use tms_tabulated_temperature_shift, only : &
    tabulated_log10_shift_provider_t, &
    create_tabulated_temperature_shift_provider
  use tms_thermorheological_dynamic_modulus_provider, only : &
    thermorheological_dynamic_modulus_provider_t, &
    create_thermorheological_dynamic_modulus_provider
  use tms_tts_types, only : tts_identification_result_t, tts_text_length
  implicit none
  private

  !> V0.8.1 distilled table ve empirical shift table'ını mevcut V0.8.0
  !! concrete provider nesneleri olarak taşır. Parallel provider mimarisi
  !! oluşturmaz; bütün nesneler mevcut runtime API'lerini doğrudan kullanır.
  type, public :: tts_runtime_export_t
    logical :: available = .false.
    character(len=tts_text_length) :: message = ""
    type(tabulated_dynamic_modulus_provider_t) :: master_curve_provider
    type(tabulated_log10_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: &
      thermorheological_provider
  end type tts_runtime_export_t

  public :: create_tts_runtime_export

contains

  !> Başarılı identification sonucundan V0.8.0 tabulated modulus, tabulated
  !! log10(a_T) ve thermorheological provider'larını kurar. Master frekansları
  !! [Hz], modüller [Pa], sıcaklıklar [K] ve strain'ler boyutsuzdur. Export
  !! readiness yoksa error-stop yerine available=false döner.
  function create_tts_runtime_export(result) result(runtime_export)
    type(tts_identification_result_t), intent(in) :: result
    type(tts_runtime_export_t) :: runtime_export

    type(dynamic_material_metadata_t) :: metadata
    type(material_frequency_point), allocatable :: points(:)
    real(dp), allocatable :: shifts(:)
    real(dp), allocatable :: temperatures(:)
    integer :: i

    if (.not. result%runtime_export_ready .or. &
        .not. allocated(result%runtime_master_table) .or. &
        .not. allocated(result%empirical_shifts)) then
      runtime_export%message = "Identification sonucu runtime export'a hazır değil."
      return
    end if
    if (size(result%runtime_master_table) < 2 .or. &
        size(result%empirical_shifts) < 2) then
      runtime_export%message = "Runtime provider tabloları en az iki nokta ister."
      return
    end if

    allocate(points(size(result%runtime_master_table)))
    do i = 1, size(points)
      points(i)%frequency = &
        result%runtime_master_table(i)%reduced_frequency_hz
      points(i)%storage_modulus = &
        result%runtime_master_table(i)%storage_modulus_pa
      points(i)%loss_modulus = &
        result%runtime_master_table(i)%loss_modulus_pa
      points(i)%temperature = result%reference_temperature_k
    end do
    metadata = make_runtime_metadata(result)
    runtime_export%master_curve_provider = &
      create_tabulated_dynamic_modulus_provider(points, metadata, &
        LINEAR_LOG_FREQUENCY)

    allocate(temperatures(size(result%empirical_shifts)))
    allocate(shifts(size(result%empirical_shifts)))
    do i = 1, size(result%empirical_shifts)
      temperatures(i) = result%empirical_shifts(i)%temperature_k
      shifts(i) = result%empirical_shifts(i)%log10_a_t
    end do
    runtime_export%shift_provider = &
      create_tabulated_temperature_shift_provider(temperatures, shifts, &
        result%reference_temperature_k)
    runtime_export%thermorheological_provider = &
      create_thermorheological_dynamic_modulus_provider( &
        runtime_export%master_curve_provider, runtime_export%shift_provider)
    runtime_export%available = .true.
    runtime_export%message = "V0.8.0 runtime provider export hazır."
  end function create_tts_runtime_export

  pure function make_runtime_metadata(result) result(metadata)
    type(tts_identification_result_t), intent(in) :: result
    type(dynamic_material_metadata_t) :: metadata

    metadata%dataset_identifier = &
      trim(result%source_family%family_identifier)//"-V0.8.1-MASTER"
    metadata%material_identifier = &
      result%source_family%common_state%material_identifier
    metadata%dataset_temperature_k = result%reference_temperature_k
    metadata%has_dynamic_shear_strain_amplitude = .true.
    metadata%dynamic_shear_strain_amplitude = &
      result%source_family%common_state%dynamic_strain_amplitude_ratio
    metadata%has_static_shear_prestrain = .true.
    metadata%static_shear_prestrain = &
      result%source_family%common_state%static_prestrain_ratio
    metadata%deformation_mode = DYNAMIC_DEFORMATION_MODE_SHEAR
    metadata%has_conditioning_state = .true.
    metadata%conditioning_state = &
      result%source_family%common_state%conditioning_description
    metadata%has_material_state = .true.
    metadata%material_state = &
      result%source_family%common_state%batch_state_identifier
    metadata%has_test_method_source = .true.
    metadata%test_method_source = &
      trim(result%source_family%common_state%test_method)//" | "// &
      trim(result%source_family%common_state%source_metadata)
    metadata%has_notes = .true.
    metadata%notes = &
      "V0.8.1 measured-isotherm identification; no smoothing/averaging."
  end function make_runtime_metadata

end module tms_tts_runtime_export
