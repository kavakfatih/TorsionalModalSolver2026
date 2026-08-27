module tms_dynamic_material_metadata
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: DYNAMIC_DEFORMATION_MODE_SHEAR = 1
  integer, parameter, public :: DYNAMIC_DEFORMATION_MODE_TENSILE = 2
  integer, parameter, public :: dynamic_metadata_text_length = 256

  !> Tek bir dinamik malzeme veri kümesinin deney ve operating-state
  !! izlenebilirliğini taşır. Bir provider yalnız bu tek durumu temsil eder;
  !! farklı sıcaklık, genlik, prestrain veya conditioning durumları ayrı veri
  !! kümeleri olarak oluşturulmalıdır. Bilinmeyen seçimlik alanlar magic sayı
  !! yerine açık has_* bayraklarıyla gösterilir.
  type, public :: dynamic_material_metadata_t
    character(len=dynamic_metadata_text_length) :: dataset_identifier = ""
    character(len=dynamic_metadata_text_length) :: material_identifier = ""

    logical :: has_specimen_identifier = .false.
    character(len=dynamic_metadata_text_length) :: specimen_identifier = ""

    !> Veri kümesinin ölçülmüş isotherm sıcaklığı [K]. V0.7 bu sıcaklıklar
    !! arasında interpolation veya shift işlemi yapmaz.
    real(dp) :: dataset_temperature_k = 0.0_dp

    !> Deneydeki dinamik shear strain amplitude gamma_a [-]. Alan yalnız
    !! has_dynamic_shear_strain_amplitude doğruysa kullanılabilir.
    logical :: has_dynamic_shear_strain_amplitude = .false.
    real(dp) :: dynamic_shear_strain_amplitude = 0.0_dp

    !> Deneyin statik shear prestrain değeri gamma_0 [-]. Alan yalnız
    !! has_static_shear_prestrain doğruysa kullanılabilir.
    logical :: has_static_shear_prestrain = .false.
    real(dp) :: static_shear_prestrain = 0.0_dp

    !> Ölçümün deformation mode kimliği [-]. Direct torsional binding yalnız
    !! DYNAMIC_DEFORMATION_MODE_SHEAR kabul eder.
    integer :: deformation_mode = 0

    logical :: has_conditioning_state = .false.
    character(len=dynamic_metadata_text_length) :: conditioning_state = ""

    logical :: has_material_state = .false.
    character(len=dynamic_metadata_text_length) :: material_state = ""

    logical :: has_test_method_source = .false.
    character(len=dynamic_metadata_text_length) :: test_method_source = ""

    logical :: has_standard_reference = .false.
    character(len=dynamic_metadata_text_length) :: standard_reference = ""

    logical :: has_notes = .false.
    character(len=dynamic_metadata_text_length) :: notes = ""
  end type dynamic_material_metadata_t

  public :: validate_dynamic_material_metadata

contains

  !> Dataset-level deney metadata'sının tek operating-state sözleşmesini
  !! doğrular.
  !! Girdi: Kimlikler, mutlak sıcaklık [K], seçimlik boyutsuz strain değerleri
  !! ve deformation mode taşıyan metadata.
  !! Çıktı: Yoktur; geçersiz veya availability bayrağıyla tutarsız alan error
  !! stop ile reddedilir.
  !! Varsayım ve sınır: Bu doğrulama ölçüm belirsizliği yorumlamaz; genlik,
  !! prestrain veya sıcaklık interpolation'ı yapmaz. Her provider tek deney
  !! durumunu temsil eder.
  pure subroutine validate_dynamic_material_metadata(metadata)
    type(dynamic_material_metadata_t), intent(in) :: metadata

    if (len_trim(metadata%dataset_identifier) == 0) then
      error stop "Dinamik malzeme dataset kimliği boş olamaz."
    end if
    if (len_trim(metadata%material_identifier) == 0) then
      error stop "Dinamik malzeme/compound kimliği boş olamaz."
    end if
    if (.not. ieee_is_finite(metadata%dataset_temperature_k) .or. &
        metadata%dataset_temperature_k <= 0.0_dp) then
      error stop "Dataset sıcaklığı sonlu ve pozitif Kelvin olmalıdır."
    end if

    select case (metadata%deformation_mode)
    case (DYNAMIC_DEFORMATION_MODE_SHEAR, &
          DYNAMIC_DEFORMATION_MODE_TENSILE)
      continue
    case default
      error stop "Dinamik malzeme deformation mode tanımsızdır."
    end select

    if (metadata%has_dynamic_shear_strain_amplitude) then
      if (.not. ieee_is_finite( &
          metadata%dynamic_shear_strain_amplitude) .or. &
          metadata%dynamic_shear_strain_amplitude <= 0.0_dp) then
        error stop "Dinamik shear strain amplitude sonlu ve pozitif olmalıdır."
      end if
    else if (.not. ieee_is_finite( &
        metadata%dynamic_shear_strain_amplitude) .or. &
        abs(metadata%dynamic_shear_strain_amplitude) > 0.0_dp) then
      error stop "Unavailable strain amplitude magic değer taşıyamaz."
    end if
    if (metadata%has_static_shear_prestrain) then
      if (.not. ieee_is_finite(metadata%static_shear_prestrain) .or. &
          metadata%static_shear_prestrain < 0.0_dp) then
        error stop "Statik shear prestrain sonlu ve negatif olmamalıdır."
      end if
    else if (.not. ieee_is_finite(metadata%static_shear_prestrain) .or. &
        abs(metadata%static_shear_prestrain) > 0.0_dp) then
      error stop "Unavailable static prestrain magic değer taşıyamaz."
    end if

    call validate_optional_text(metadata%has_specimen_identifier, &
      metadata%specimen_identifier, "Specimen kimliği")
    call validate_optional_text(metadata%has_conditioning_state, &
      metadata%conditioning_state, "Conditioning state")
    call validate_optional_text(metadata%has_material_state, &
      metadata%material_state, "Material state")
    call validate_optional_text(metadata%has_test_method_source, &
      metadata%test_method_source, "Test method/source")
    call validate_optional_text(metadata%has_standard_reference, &
      metadata%standard_reference, "Standard/reference")
    call validate_optional_text(metadata%has_notes, metadata%notes, "Not")
  end subroutine validate_dynamic_material_metadata

  !> Seçimlik metnin availability bayrağıyla tutarlı olduğunu doğrular.
  !! Bu yardımcı fizik veya matematik hesabı yapmaz.
  pure subroutine validate_optional_text(is_available, value, field_name)
    logical, intent(in) :: is_available
    character(len=*), intent(in) :: value
    character(len=*), intent(in) :: field_name

    if (is_available .and. len_trim(value) == 0) then
      error stop trim(field_name)//" available ise boş olamaz."
    end if
    if (.not. is_available .and. len_trim(value) /= 0) then
      error stop trim(field_name)//" availability bayrağı ile tutarsız."
    end if
  end subroutine validate_optional_text

end module tms_dynamic_material_metadata
