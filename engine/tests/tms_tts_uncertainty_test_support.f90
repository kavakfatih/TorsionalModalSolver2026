module tms_tts_uncertainty_test_support
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, &
    UNCERTAINTY_SOURCE_COMBINED_STANDARD
  implicit none
  private

  public :: make_relative_uncertainty_family
  public :: scale_uncertainty_family

contains

  !> Synthetic test fixture'ında u_G/G oranını sabit tutan explicit pointwise
  !! standard uncertainty overlay'i üretir. Bu helper production uncertainty
  !! tahmini yapmaz; yalnız bağımsız deterministic test girdisi kurar.
  function make_relative_uncertainty_family( &
      family, storage_relative_uncertainty, loss_relative_uncertainty) &
      result(uncertainty_family)
    type(tts_material_family_t), intent(in) :: family
    real(dp), intent(in) :: storage_relative_uncertainty
    real(dp), intent(in) :: loss_relative_uncertainty
    type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty_family

    integer :: i
    integer :: j

    uncertainty_family%family_identifier = "SYNTHETIC-POINTWISE-UNCERTAINTY"
    uncertainty_family%provenance = &
      "Deterministic combined-standard uncertainty test fixture"
    allocate(uncertainty_family%isotherms(size(family%isotherms)))
    do i = 1, size(family%isotherms)
      uncertainty_family%isotherms(i)%isotherm_identifier = &
        family%isotherms(i)%isotherm_identifier
      uncertainty_family%isotherms(i)%temperature_k = &
        family%isotherms(i)%temperature_k
      allocate(uncertainty_family%isotherms(i)%points( &
        size(family%isotherms(i)%points)))
      do j = 1, size(family%isotherms(i)%points)
        uncertainty_family%isotherms(i)%points(j)%temperature_k = &
          family%isotherms(i)%temperature_k
        uncertainty_family%isotherms(i)%points(j)%frequency_hz = &
          family%isotherms(i)%points(j)%frequency_hz
        uncertainty_family%isotherms(i)%points(j) &
          %storage_standard_uncertainty_pa = &
          storage_relative_uncertainty * &
          family%isotherms(i)%points(j)%storage_modulus_pa
        uncertainty_family%isotherms(i)%points(j) &
          %loss_standard_uncertainty_pa = &
          loss_relative_uncertainty * &
          family%isotherms(i)%points(j)%loss_modulus_pa
        uncertainty_family%isotherms(i)%points(j) &
          %storage_uncertainty_available = .true.
        uncertainty_family%isotherms(i)%points(j) &
          %loss_uncertainty_available = .true.
        uncertainty_family%isotherms(i)%points(j)%uncertainty_source = &
          UNCERTAINTY_SOURCE_COMBINED_STANDARD
        uncertainty_family%isotherms(i)%points(j)%source_metadata = &
          "Synthetic standard uncertainty"
      end do
    end do
  end function make_relative_uncertainty_family

  !> Bütün available standard uncertainty değerlerini aynı pozitif lambda ile
  !! ölçekler. Weighted minimizer scale-invariance test fixture'ı içindir.
  subroutine scale_uncertainty_family(uncertainty_family, scale_factor)
    type(tts_dynamic_modulus_uncertainty_family_t), intent(inout) :: &
      uncertainty_family
    real(dp), intent(in) :: scale_factor
    integer :: i
    integer :: j

    do i = 1, size(uncertainty_family%isotherms)
      do j = 1, size(uncertainty_family%isotherms(i)%points)
        uncertainty_family%isotherms(i)%points(j) &
          %storage_standard_uncertainty_pa = scale_factor * &
          uncertainty_family%isotherms(i)%points(j) &
            %storage_standard_uncertainty_pa
        uncertainty_family%isotherms(i)%points(j) &
          %loss_standard_uncertainty_pa = scale_factor * &
          uncertainty_family%isotherms(i)%points(j) &
            %loss_standard_uncertainty_pa
      end do
    end do
  end subroutine scale_uncertainty_family

end module tms_tts_uncertainty_test_support
