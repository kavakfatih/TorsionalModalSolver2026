module tms_tts_test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : pi
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
  public :: make_generalized_maxwell_trs_family
  public :: populate_isotherm_from_log_grid
  public :: replace_loss_truth_shift
  public :: truth_storage_modulus
  public :: truth_loss_modulus
  public :: generalized_maxwell_storage_modulus
  public :: generalized_maxwell_loss_modulus
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

  !> Deterministic generalized-Maxwell analytical truth ile curved exact-TRS
  !! family üretir. Model üretim material implementation'ı değildir; yalnız
  !! bağımsız test denklemidir. f [Hz], T [K], G [Pa], tau [s] ve
  !! f_r=10^s*f kullanılır. Grid -3...+3 log10(Hz) aralığında 0.1 decade'dir.
  function make_generalized_maxwell_trs_family( &
      temperatures_k, log10_shifts) result(family)
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: log10_shifts(:)
    type(tts_material_family_t) :: family

    real(dp) :: log10_frequency_grid(61)
    integer :: i
    integer :: isotherm_index

    if (size(temperatures_k) /= size(log10_shifts)) then
      error stop "Generalized-Maxwell temperature/shift boyutları eşit olmalı."
    end if
    family = make_exact_trs_family(temperatures_k, log10_shifts)
    family%family_identifier = "SYNTHETIC-GENERALIZED-MAXWELL-TTS"
    family%common_state%batch_state_identifier = "BATCH-MAXWELL-01"
    family%common_state%test_method = &
      "Independent generalized-Maxwell analytical benchmark"
    do i = 1, size(log10_frequency_grid)
      log10_frequency_grid(i) = -3.0_dp + 0.1_dp*real(i - 1, dp)
    end do
    do isotherm_index = 1, size(family%isotherms)
      call populate_generalized_maxwell_isotherm( &
        family%isotherms(isotherm_index), log10_frequency_grid, &
        log10_shifts(isotherm_index))
    end do
  end function make_generalized_maxwell_trs_family

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

  !> Measured physical f [Hz] noktasını f_r=10^s*f ile bağımsız analytical
  !! generalized-Maxwell master'a taşır ve G'/G'' [Pa] üretir. Horizontal-only
  !! exact TRS varsayılır; production identification yordamı çağrılmaz.
  subroutine populate_generalized_maxwell_isotherm( &
      isotherm, log10_frequency_grid, truth_shift)
    type(tts_isotherm_t), intent(inout) :: isotherm
    real(dp), intent(in) :: log10_frequency_grid(:)
    real(dp), intent(in) :: truth_shift

    real(dp) :: physical_frequency_hz
    real(dp) :: reduced_frequency_hz
    integer :: i

    if (allocated(isotherm%points)) deallocate(isotherm%points)
    allocate(isotherm%points(size(log10_frequency_grid)))
    do i = 1, size(log10_frequency_grid)
      physical_frequency_hz = 10.0_dp**log10_frequency_grid(i)
      reduced_frequency_hz = &
        10.0_dp**(log10_frequency_grid(i) + truth_shift)
      isotherm%points(i)%frequency_hz = physical_frequency_hz
      isotherm%points(i)%storage_modulus_pa = &
        generalized_maxwell_storage_modulus(reduced_frequency_hz)
      isotherm%points(i)%loss_modulus_pa = &
        generalized_maxwell_loss_modulus(reduced_frequency_hz)
      isotherm%points(i)%storage_quality = MEASUREMENT_VALID
      isotherm%points(i)%loss_quality = MEASUREMENT_VALID
    end do
  end subroutine populate_generalized_maxwell_isotherm

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

  !> Üç relaxation branch'li generalized-Maxwell modelinin storage shear
  !! modulus'unu hesaplar. Model:
  !! G'=G_inf+sum[G_k*(omega*tau_k)^2/(1+(omega*tau_k)^2)].
  !! Girdi f [Hz], omega=2*pi*f [rad/s]; G_inf/G_k ve çıktı [Pa], tau [s].
  !! Parametreler yalnız deterministic validation truth'udur.
  pure function generalized_maxwell_storage_modulus(frequency_hz) &
      result(storage_modulus_pa)
    real(dp), intent(in) :: frequency_hz
    real(dp) :: storage_modulus_pa

    real(dp), parameter :: equilibrium_modulus_pa = 0.8e6_dp
    real(dp), parameter :: branch_moduli_pa(3) = &
      [0.5e6_dp, 1.0e6_dp, 2.0e6_dp]
    real(dp), parameter :: relaxation_times_s(3) = &
      [1.0e-3_dp, 1.0e-1_dp, 1.0e1_dp]
    real(dp) :: omega_tau
    integer :: branch_index

    storage_modulus_pa = equilibrium_modulus_pa
    do branch_index = 1, size(branch_moduli_pa)
      omega_tau = 2.0_dp*pi*frequency_hz* &
        relaxation_times_s(branch_index)
      storage_modulus_pa = storage_modulus_pa + &
        branch_moduli_pa(branch_index)*omega_tau**2 / &
        (1.0_dp + omega_tau**2)
    end do
  end function generalized_maxwell_storage_modulus

  !> Üç relaxation branch'li generalized-Maxwell modelinin loss shear
  !! modulus'unu hesaplar. Model:
  !! G''=sum[G_k*(omega*tau_k)/(1+(omega*tau_k)^2)].
  !! Girdi f [Hz], omega=2*pi*f [rad/s]; G_k ve çıktı [Pa], tau [s].
  !! Parametreler yalnız deterministic validation truth'udur.
  pure function generalized_maxwell_loss_modulus(frequency_hz) &
      result(loss_modulus_pa)
    real(dp), intent(in) :: frequency_hz
    real(dp) :: loss_modulus_pa

    real(dp), parameter :: branch_moduli_pa(3) = &
      [0.5e6_dp, 1.0e6_dp, 2.0e6_dp]
    real(dp), parameter :: relaxation_times_s(3) = &
      [1.0e-3_dp, 1.0e-1_dp, 1.0e1_dp]
    real(dp) :: omega_tau
    integer :: branch_index

    loss_modulus_pa = 0.0_dp
    do branch_index = 1, size(branch_moduli_pa)
      omega_tau = 2.0_dp*pi*frequency_hz* &
        relaxation_times_s(branch_index)
      loss_modulus_pa = loss_modulus_pa + &
        branch_moduli_pa(branch_index)*omega_tau / &
        (1.0_dp + omega_tau**2)
    end do
  end function generalized_maxwell_loss_modulus

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
