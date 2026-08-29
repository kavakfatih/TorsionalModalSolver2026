module tms_tts_test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_isotherm_t, &
    tts_measurement_point_t, tts_empirical_shift_t, MEASUREMENT_VALID, &
    TTS_DEFORMATION_MODE_SHEAR
  implicit none
  private

  real(dp), parameter, public :: test_storage_intercept = 6.2_dp
  real(dp), parameter, public :: test_loss_intercept = 5.4_dp
  real(dp), parameter, public :: test_storage_slope = 0.25_dp
  real(dp), parameter, public :: test_loss_slope = 0.15_dp

  public :: make_exact_trs_family
  public :: populate_isotherm_from_log_grid
  public :: replace_loss_truth_shift
  public :: truth_storage_modulus
  public :: truth_loss_modulus
  public :: find_empirical_shift
  public :: assert_true
  public :: assert_close

contains

  function make_exact_trs_family( &
      temperatures_k, log10_shifts, storage_slope, loss_slope) result(family)
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: log10_shifts(:)
    real(dp), intent(in), optional :: storage_slope
    real(dp), intent(in), optional :: loss_slope
    type(tts_material_family_t) :: family

    real(dp), parameter :: default_grid(7) = &
      [-2.0_dp, -1.35_dp, -0.65_dp, 0.0_dp, 0.55_dp, 1.25_dp, 2.0_dp]
    real(dp) :: selected_loss_slope
    real(dp) :: selected_storage_slope
    integer :: i
    character(len=32) :: identifier

    if (size(temperatures_k) /= size(log10_shifts)) then
      error stop "Test family temperature/shift boyutları eşit olmalıdır."
    end if
    selected_storage_slope = test_storage_slope
    selected_loss_slope = test_loss_slope
    if (present(storage_slope)) selected_storage_slope = storage_slope
    if (present(loss_slope)) selected_loss_slope = loss_slope

    family%family_identifier = "SYNTHETIC-TTS-FAMILY"
    family%common_state%material_identifier = "SYNTHETIC-ELASTOMER"
    family%common_state%batch_state_identifier = "BATCH-EXACT-01"
    family%common_state%dynamic_strain_amplitude_ratio = 0.005_dp
    family%common_state%static_prestrain_ratio = 0.02_dp
    family%common_state%deformation_mode = TTS_DEFORMATION_MODE_SHEAR
    family%common_state%conditioning_description = "10-cycle preconditioning"
    family%common_state%test_method = "Synthetic DMA shear sweep"
    family%common_state%source_metadata = "TMS26 deterministic test data"
    allocate(family%isotherms(size(temperatures_k)))
    do i = 1, size(temperatures_k)
      write(identifier, '("ISO-",I0)') i
      family%isotherms(i)%isotherm_identifier = trim(identifier)
      family%isotherms(i)%temperature_k = temperatures_k(i)
      write(identifier, '("SPECIMEN-",I0)') i
      family%isotherms(i)%specimen_identifier = trim(identifier)
      write(identifier, '("SOURCE-",I0)') i
      family%isotherms(i)%source_identifier = trim(identifier)
      call populate_isotherm_from_log_grid(family%isotherms(i), &
        default_grid, log10_shifts(i), selected_storage_slope, &
        selected_loss_slope)
    end do
  end function make_exact_trs_family

  subroutine populate_isotherm_from_log_grid( &
      isotherm, log10_frequency_grid, truth_shift, storage_slope, loss_slope)
    type(tts_isotherm_t), intent(inout) :: isotherm
    real(dp), intent(in) :: log10_frequency_grid(:)
    real(dp), intent(in) :: truth_shift
    real(dp), intent(in) :: storage_slope
    real(dp), intent(in) :: loss_slope
    integer :: i

    if (allocated(isotherm%points)) deallocate(isotherm%points)
    allocate(isotherm%points(size(log10_frequency_grid)))
    do i = 1, size(log10_frequency_grid)
      isotherm%points(i)%frequency_hz = 10.0_dp**log10_frequency_grid(i)
      isotherm%points(i)%storage_modulus_pa = 10.0_dp**( &
        test_storage_intercept + storage_slope * &
        (log10_frequency_grid(i) + truth_shift))
      isotherm%points(i)%loss_modulus_pa = 10.0_dp**( &
        test_loss_intercept + loss_slope * &
        (log10_frequency_grid(i) + truth_shift))
      isotherm%points(i)%storage_quality = MEASUREMENT_VALID
      isotherm%points(i)%loss_quality = MEASUREMENT_VALID
    end do
  end subroutine populate_isotherm_from_log_grid

  subroutine replace_loss_truth_shift(isotherm, truth_shift, loss_slope)
    type(tts_isotherm_t), intent(inout) :: isotherm
    real(dp), intent(in) :: truth_shift
    real(dp), intent(in), optional :: loss_slope
    real(dp) :: selected_slope
    integer :: i

    selected_slope = test_loss_slope
    if (present(loss_slope)) selected_slope = loss_slope
    do i = 1, size(isotherm%points)
      isotherm%points(i)%loss_modulus_pa = 10.0_dp**( &
        test_loss_intercept + selected_slope * &
        (log10(isotherm%points(i)%frequency_hz) + truth_shift))
    end do
  end subroutine replace_loss_truth_shift

  pure function truth_storage_modulus( &
      reduced_frequency_hz, storage_slope) result(value)
    real(dp), intent(in) :: reduced_frequency_hz
    real(dp), intent(in), optional :: storage_slope
    real(dp) :: value
    real(dp) :: selected_slope

    selected_slope = test_storage_slope
    if (present(storage_slope)) selected_slope = storage_slope
    value = 10.0_dp**(test_storage_intercept + &
      selected_slope*log10(reduced_frequency_hz))
  end function truth_storage_modulus

  pure function truth_loss_modulus( &
      reduced_frequency_hz, loss_slope) result(value)
    real(dp), intent(in) :: reduced_frequency_hz
    real(dp), intent(in), optional :: loss_slope
    real(dp) :: value
    real(dp) :: selected_slope

    selected_slope = test_loss_slope
    if (present(loss_slope)) selected_slope = loss_slope
    value = 10.0_dp**(test_loss_intercept + &
      selected_slope*log10(reduced_frequency_hz))
  end function truth_loss_modulus

  pure function find_empirical_shift(shifts, source_index) result(value)
    type(tts_empirical_shift_t), intent(in) :: shifts(:)
    integer, intent(in) :: source_index
    real(dp) :: value
    integer :: i

    value = huge(1.0_dp)
    do i = 1, size(shifts)
      if (shifts(i)%source_isotherm_index == source_index) then
        value = shifts(i)%log10_a_t
        return
      end if
    end do
  end function find_empirical_shift

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message

    if (.not. condition) error stop message
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message
    real(dp) :: scale

    scale = max(1.0_dp, abs(expected))
    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp .or. &
        abs(actual - expected) > tolerance*scale) then
      error stop message
    end if
  end subroutine assert_close

end module tms_tts_test_support
