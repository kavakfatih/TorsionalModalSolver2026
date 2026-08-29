module tms_tts_shift_law_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_identification_result_t
  use tms_tts_shift_law_types, only : tts_pair_shift_observation_t, &
    tts_wlf_fit_configuration_t, tts_arrhenius_fit_result_t, &
    tts_wlf_fit_result_t, tts_shift_law_identification_result_t, &
    SHIFT_LAW_FIT_SUCCESS, SHIFT_LAW_FIT_INVALID_INPUT, &
    SHIFT_LAW_FIT_INSUFFICIENT_DATA
  use tms_tts_shift_law_fit, only : fit_tts_arrhenius_shift_law, &
    fit_tts_wlf_shift_law, predict_tts_arrhenius_pair_shift, &
    predict_tts_wlf_pair_shift, validate_tts_pair_shift_observations
  implicit none
  private

  integer, parameter :: minimum_loto_temperature_count = 5

  public :: fit_tts_shift_laws
  public :: extract_tts_pair_shift_observations

contains

  !> V0.8.1 identification sonucunun adjacent pair delta_s gözlemlerini
  !! kullanarak Arrhenius ve WLF approximation'larını üretir. Empirical
  !! shifts, master cloud ve runtime table yalnız intent(in) kaynaktır ve
  !! mutate edilmez; cumulative absolute s değerleri fit observation değildir.
  !! Equal-weight pair fit sonrası en az beş measured temperature varsa LOTO
  !! predictive diagnostics hesaplanır. Automatic best-model seçilmez.
  pure function fit_tts_shift_laws( &
      identification, wlf_configuration) result(result)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_wlf_fit_configuration_t), intent(in), optional :: &
      wlf_configuration
    type(tts_shift_law_identification_result_t) :: result

    type(tts_pair_shift_observation_t), allocatable :: observations(:)
    type(tts_wlf_fit_configuration_t) :: settings
    integer :: extraction_status
    character(len=256) :: extraction_message

    if (present(wlf_configuration)) settings = wlf_configuration
    call extract_tts_pair_shift_observations(identification, observations, &
      extraction_status, extraction_message)
    if (extraction_status /= SHIFT_LAW_FIT_SUCCESS) then
      result%status = extraction_status
      result%message = extraction_message
      return
    end if

    result%pair_observations_available = .true.
    result%pair_observation_count = size(observations)
    result%arrhenius = fit_tts_arrhenius_shift_law(observations, &
      identification%reference_temperature_k)
    result%wlf = fit_tts_wlf_shift_law(observations, &
      identification%reference_temperature_k, settings)

    call attach_arrhenius_loto(observations, &
      identification%reference_temperature_k, result%arrhenius)
    call attach_wlf_loto(observations, &
      identification%reference_temperature_k, settings, result%wlf)

    result%status = SHIFT_LAW_FIT_SUCCESS
    result%message = "Empirical adjacent-pair shift-law evidence üretildi."
  end function fit_tts_shift_laws

  !> V0.8.1 pair result provenance'ını temperature endpoint'leriyle
  !! birleştirerek delta_s=s_j-s_i observation dizisi üretir. Yalnız
  !! pair_shift_results okunur; empirical cumulative shift table bilinçli
  !! olarak bu extraction'ın girdisi değildir. Failure status ile döner.
  pure subroutine extract_tts_pair_shift_observations( &
      identification, observations, status, message)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_pair_shift_observation_t), allocatable, intent(out) :: &
      observations(:)
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    real(dp) :: maximum_temperature_k
    real(dp) :: minimum_temperature_k
    integer :: i
    integer :: index_i
    integer :: index_j
    integer :: validation_status
    character(len=256) :: validation_message

    status = SHIFT_LAW_FIT_INVALID_INPUT
    message = ""
    if (.not. allocated(identification%source_family%isotherms) .or. &
        .not. allocated(identification%pair_shift_results)) then
      status = SHIFT_LAW_FIT_INSUFFICIENT_DATA
      message = "Identification adjacent pair results içermiyor."
      return
    end if
    if (size(identification%pair_shift_results) == 0) then
      status = SHIFT_LAW_FIT_INSUFFICIENT_DATA
      message = "Identification en az bir adjacent pair result içermeli."
      return
    end if
    if (identification%reference_isotherm_index < 1 .or. &
        identification%reference_isotherm_index > &
          size(identification%source_family%isotherms)) then
      message = "Identification measured reference indeksi geçersiz."
      return
    end if
    if (.not. temperatures_are_machine_equivalent( &
        identification%reference_temperature_k, &
        identification%source_family%isotherms( &
          identification%reference_isotherm_index)%temperature_k)) then
      message = "Identification reference temperature provenance tutarsız."
      return
    end if

    allocate(observations(size(identification%pair_shift_results)))
    do i = 1, size(observations)
      index_i = identification%pair_shift_results(i) &
        %reference_isotherm_index
      index_j = identification%pair_shift_results(i)%moving_isotherm_index
      if (.not. identification%pair_shift_results(i)%shift_available .or. &
          index_i < 1 .or. index_i > &
            size(identification%source_family%isotherms) .or. &
          index_j < 1 .or. index_j > &
            size(identification%source_family%isotherms) .or. &
          index_i == index_j) then
        deallocate(observations)
        message = "Adjacent pair shift/provenance fit için geçersiz."
        return
      end if
      observations(i)%isotherm_i_index = index_i
      observations(i)%isotherm_j_index = index_j
      observations(i)%temperature_i_k = identification%source_family &
        %isotherms(index_i)%temperature_k
      observations(i)%temperature_j_k = identification%source_family &
        %isotherms(index_j)%temperature_k
      observations(i)%delta_s_j_minus_i = &
        identification%pair_shift_results(i)%delta_s
    end do

    call validate_tts_pair_shift_observations(observations, &
      identification%reference_temperature_k, 1, minimum_temperature_k, &
      maximum_temperature_k, validation_status, validation_message)
    if (validation_status /= SHIFT_LAW_FIT_SUCCESS) then
      deallocate(observations)
      status = validation_status
      message = validation_message
      return
    end if

    status = SHIFT_LAW_FIT_SUCCESS
    message = "Adjacent empirical pair observations hazır."
  end subroutine extract_tts_pair_shift_observations

  pure subroutine attach_arrhenius_loto( &
      observations, reference_temperature_k, fit)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: reference_temperature_k
    type(tts_arrhenius_fit_result_t), intent(inout) :: fit

    type(tts_arrhenius_fit_result_t) :: training_fit
    type(tts_pair_shift_observation_t), allocatable :: held_out(:)
    type(tts_pair_shift_observation_t), allocatable :: training(:)
    real(dp), allocatable :: residual_pool(:)
    real(dp), allocatable :: temperatures(:)
    real(dp), allocatable :: fold_residuals(:)
    integer :: fold_index
    integer :: i
    integer :: residual_count
    logical :: all_folds_available

    if (.not. fit%fit_available) return
    temperatures = collect_unique_temperatures(observations)
    allocate(fit%loto_diagnostics(size(temperatures)))
    do fold_index = 1, size(temperatures)
      fit%loto_diagnostics(fold_index)%held_out_temperature_k = &
        temperatures(fold_index)
    end do
    if (size(temperatures) < minimum_loto_temperature_count) return

    allocate(residual_pool(2*size(observations)), source=0.0_dp)
    residual_count = 0
    all_folds_available = .true.
    do fold_index = 1, size(temperatures)
      call partition_loto_observations(observations, temperatures(fold_index), &
        training, held_out)
      fit%loto_diagnostics(fold_index)%training_observation_count = &
        size(training)
      fit%loto_diagnostics(fold_index)%validation_observation_count = &
        size(held_out)
      training_fit = fit_tts_arrhenius_shift_law(training, &
        reference_temperature_k)
      if (.not. training_fit%fit_available .or. size(held_out) == 0) then
        all_folds_available = .false.
        cycle
      end if
      allocate(fold_residuals(size(held_out)))
      do i = 1, size(held_out)
        fold_residuals(i) = held_out(i)%delta_s_j_minus_i - &
          predict_tts_arrhenius_pair_shift(training_fit%beta_k, held_out(i))
      end do
      call set_fold_metrics(fold_residuals, &
        fit%loto_diagnostics(fold_index))
      residual_pool(residual_count+1:residual_count+size(fold_residuals)) = &
        fold_residuals
      residual_count = residual_count + size(fold_residuals)
      deallocate(fold_residuals)
    end do
    if (.not. all_folds_available .or. residual_count == 0) return
    fit%predictive_validation_available = .true.
    fit%loto_rmse = sqrt(sum(residual_pool(1:residual_count)**2) / &
      real(residual_count, dp))
    fit%loto_max_abs_residual = &
      maxval(abs(residual_pool(1:residual_count)))
  end subroutine attach_arrhenius_loto

  pure subroutine attach_wlf_loto( &
      observations, reference_temperature_k, settings, fit)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: reference_temperature_k
    type(tts_wlf_fit_configuration_t), intent(in) :: settings
    type(tts_wlf_fit_result_t), intent(inout) :: fit

    type(tts_wlf_fit_result_t) :: training_fit
    type(tts_pair_shift_observation_t), allocatable :: held_out(:)
    type(tts_pair_shift_observation_t), allocatable :: training(:)
    real(dp), allocatable :: residual_pool(:)
    real(dp), allocatable :: temperatures(:)
    real(dp), allocatable :: fold_residuals(:)
    integer :: fold_index
    integer :: i
    integer :: residual_count
    logical :: all_folds_available

    if (.not. fit%fit_available .or. .not. fit%parameter_identifiable) return
    temperatures = collect_unique_temperatures(observations)
    allocate(fit%loto_diagnostics(size(temperatures)))
    do fold_index = 1, size(temperatures)
      fit%loto_diagnostics(fold_index)%held_out_temperature_k = &
        temperatures(fold_index)
    end do
    if (size(temperatures) < minimum_loto_temperature_count) return

    allocate(residual_pool(2*size(observations)), source=0.0_dp)
    residual_count = 0
    all_folds_available = .true.
    do fold_index = 1, size(temperatures)
      call partition_loto_observations(observations, temperatures(fold_index), &
        training, held_out)
      fit%loto_diagnostics(fold_index)%training_observation_count = &
        size(training)
      fit%loto_diagnostics(fold_index)%validation_observation_count = &
        size(held_out)
      training_fit = fit_tts_wlf_shift_law(training, &
        reference_temperature_k, settings)
      if (.not. training_fit%fit_available .or. &
          .not. training_fit%parameter_identifiable .or. &
          size(held_out) == 0) then
        all_folds_available = .false.
        cycle
      end if
      if (.not. held_out_domain_is_pole_safe(held_out, &
          reference_temperature_k, training_fit%c2_k)) then
        all_folds_available = .false.
        cycle
      end if
      allocate(fold_residuals(size(held_out)))
      do i = 1, size(held_out)
        fold_residuals(i) = held_out(i)%delta_s_j_minus_i - &
          predict_tts_wlf_pair_shift(training_fit%c1, training_fit%c2_k, &
            reference_temperature_k, held_out(i))
      end do
      call set_fold_metrics(fold_residuals, &
        fit%loto_diagnostics(fold_index))
      residual_pool(residual_count+1:residual_count+size(fold_residuals)) = &
        fold_residuals
      residual_count = residual_count + size(fold_residuals)
      deallocate(fold_residuals)
    end do
    if (.not. all_folds_available .or. residual_count == 0) return
    fit%predictive_validation_available = .true.
    fit%loto_rmse = sqrt(sum(residual_pool(1:residual_count)**2) / &
      real(residual_count, dp))
    fit%loto_max_abs_residual = &
      maxval(abs(residual_pool(1:residual_count)))
  end subroutine attach_wlf_loto

  pure function collect_unique_temperatures(observations) result(temperatures)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), allocatable :: temperatures(:)

    real(dp), allocatable :: candidates(:)
    real(dp) :: key
    integer :: count
    integer :: i
    integer :: j

    allocate(candidates(2*size(observations)))
    do i = 1, size(observations)
      candidates(2*i-1) = observations(i)%temperature_i_k
      candidates(2*i) = observations(i)%temperature_j_k
    end do
    do i = 2, size(candidates)
      key = candidates(i)
      j = i-1
      do while (j >= 1)
        if (candidates(j) <= key) exit
        candidates(j+1) = candidates(j)
        j = j-1
      end do
      candidates(j+1) = key
    end do

    count = 0
    do i = 1, size(candidates)
      if (count == 0) then
        count = 1
        candidates(count) = candidates(i)
      else if (.not. temperatures_are_machine_equivalent( &
          candidates(i), candidates(count))) then
        count = count+1
        candidates(count) = candidates(i)
      end if
    end do
    allocate(temperatures(count))
    temperatures = candidates(1:count)
  end function collect_unique_temperatures

  pure subroutine partition_loto_observations( &
      observations, held_temperature_k, training, held_out)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: held_temperature_k
    type(tts_pair_shift_observation_t), allocatable, intent(out) :: training(:)
    type(tts_pair_shift_observation_t), allocatable, intent(out) :: held_out(:)

    logical, allocatable :: held_mask(:)
    integer :: held_count
    integer :: i
    integer :: held_position
    integer :: training_position

    allocate(held_mask(size(observations)))
    do i = 1, size(observations)
      held_mask(i) = temperatures_are_machine_equivalent( &
        observations(i)%temperature_i_k, held_temperature_k) .or. &
        temperatures_are_machine_equivalent( &
          observations(i)%temperature_j_k, held_temperature_k)
    end do
    held_count = count(held_mask)
    allocate(training(size(observations)-held_count))
    allocate(held_out(held_count))
    held_position = 0
    training_position = 0
    do i = 1, size(observations)
      if (held_mask(i)) then
        held_position = held_position+1
        held_out(held_position) = observations(i)
      else
        training_position = training_position+1
        training(training_position) = observations(i)
      end if
    end do
  end subroutine partition_loto_observations

  pure subroutine set_fold_metrics(residuals, diagnostic)
    use tms_tts_shift_law_types, only : tts_loto_fold_diagnostic_t
    real(dp), intent(in) :: residuals(:)
    type(tts_loto_fold_diagnostic_t), intent(inout) :: diagnostic

    diagnostic%available = size(residuals) > 0
    if (.not. diagnostic%available) return
    diagnostic%pair_rmse = sqrt(sum(residuals**2)/real(size(residuals), dp))
    diagnostic%pair_max_abs_residual = maxval(abs(residuals))
    diagnostic%pair_mean_residual = sum(residuals)/real(size(residuals), dp)
  end subroutine set_fold_metrics

  pure function held_out_domain_is_pole_safe( &
      held_out, reference_temperature_k, c2_k) result(safe)
    type(tts_pair_shift_observation_t), intent(in) :: held_out(:)
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: c2_k
    logical :: safe
    integer :: i

    safe = ieee_is_finite(c2_k) .and. c2_k > 0.0_dp
    if (.not. safe) return
    do i = 1, size(held_out)
      if (c2_k+held_out(i)%temperature_i_k-reference_temperature_k <= &
          0.0_dp .or. &
          c2_k+held_out(i)%temperature_j_k-reference_temperature_k <= &
          0.0_dp) then
        safe = .false.
        return
      end if
    end do
  end function held_out_domain_is_pole_safe

  pure elemental function temperatures_are_machine_equivalent(a, b) &
      result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    equivalent = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
      abs(a-b) <= 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(a), abs(b))
  end function temperatures_are_machine_equivalent

end module tms_tts_shift_law_validation
