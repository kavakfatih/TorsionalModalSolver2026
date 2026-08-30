module tms_tts_covariance_test_support
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t
  use tms_tts_covariance_types, only : &
    tts_dynamic_modulus_covariance_family_t, COVARIANCE_SOURCE_DIRECT
  implicit none
  private

  public :: make_covariance_family
  public :: scale_covariance_inputs

contains

  !> V0.8.5 testleri için V0.8.4 u_G diyagonalleriyle tam uyumlu explicit
  !! point-local covariance overlay'i üretir. rho boyutsuz ve caller-defined
  !! test girdisidir; production covariance tahmini veya universal rubber
  !! correlation varsayımı değildir.
  function make_covariance_family( &
      family, uncertainty_family, correlation) result(covariance_family)
    type(tts_material_family_t), intent(in) :: family
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    real(dp), intent(in) :: correlation
    type(tts_dynamic_modulus_covariance_family_t) :: covariance_family

    real(dp) :: storage_uncertainty
    real(dp) :: loss_uncertainty
    integer :: i
    integer :: j

    covariance_family%family_identifier = "SYNTHETIC-POINT-COVARIANCE"
    covariance_family%provenance = &
      "Deterministic direct covariance validation fixture"
    allocate(covariance_family%isotherms(size(family%isotherms)))
    do i = 1, size(family%isotherms)
      covariance_family%isotherms(i)%isotherm_identifier = &
        family%isotherms(i)%isotherm_identifier
      covariance_family%isotherms(i)%temperature_k = &
        family%isotherms(i)%temperature_k
      allocate(covariance_family%isotherms(i)%points( &
        size(family%isotherms(i)%points)))
      do j = 1, size(family%isotherms(i)%points)
        storage_uncertainty = uncertainty_family%isotherms(i)%points(j) &
          %storage_standard_uncertainty_pa
        loss_uncertainty = uncertainty_family%isotherms(i)%points(j) &
          %loss_standard_uncertainty_pa
        covariance_family%isotherms(i)%points(j)%temperature_k = &
          family%isotherms(i)%temperature_k
        covariance_family%isotherms(i)%points(j)%frequency_hz = &
          family%isotherms(i)%points(j)%frequency_hz
        covariance_family%isotherms(i)%points(j)%covariance &
          %storage_variance_pa2 = storage_uncertainty**2
        covariance_family%isotherms(i)%points(j)%covariance &
          %loss_variance_pa2 = loss_uncertainty**2
        covariance_family%isotherms(i)%points(j)%covariance &
          %storage_loss_covariance_pa2 = correlation * &
          storage_uncertainty*loss_uncertainty
        covariance_family%isotherms(i)%points(j)%covariance_available = &
          uncertainty_family%isotherms(i)%points(j) &
            %storage_uncertainty_available .and. &
          uncertainty_family%isotherms(i)%points(j) &
            %loss_uncertainty_available
        covariance_family%isotherms(i)%points(j) &
          %covariance_source_kind = COVARIANCE_SOURCE_DIRECT
        covariance_family%isotherms(i)%points(j)%measurement_method = &
          "Synthetic same-point bivariate DMA covariance"
        covariance_family%isotherms(i)%points(j)%instrument_identifier = &
          "DMA-SYNTHETIC-01"
        covariance_family%isotherms(i)%points(j)%calibration_reference = &
          "Synthetic calibration provenance"
        covariance_family%isotherms(i)%points(j)%source_identifier = &
          "COV-SOURCE-01"
        covariance_family%isotherms(i)%points(j)%source_description = &
          "Test-only explicit covariance evidence"
      end do
    end do
  end function make_covariance_family

  !> Standard uncertainty'leri lambda, covariance matrislerini lambda^2 ile
  !! tutarlı ölçekler. Mahalanobis ve matched-diagonal matrix-scale invariant
  !! testleri içindir; physical measurement değerleri değişmez.
  subroutine scale_covariance_inputs( &
      uncertainty_family, covariance_family, scale_factor)
    type(tts_dynamic_modulus_uncertainty_family_t), intent(inout) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(inout) :: &
      covariance_family
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
        covariance_family%isotherms(i)%points(j)%covariance &
          %storage_variance_pa2 = scale_factor**2 * &
          covariance_family%isotherms(i)%points(j)%covariance &
            %storage_variance_pa2
        covariance_family%isotherms(i)%points(j)%covariance &
          %loss_variance_pa2 = scale_factor**2 * &
          covariance_family%isotherms(i)%points(j)%covariance &
            %loss_variance_pa2
        covariance_family%isotherms(i)%points(j)%covariance &
          %storage_loss_covariance_pa2 = scale_factor**2 * &
          covariance_family%isotherms(i)%points(j)%covariance &
            %storage_loss_covariance_pa2
      end do
    end do
  end subroutine scale_covariance_inputs

end module tms_tts_covariance_test_support
