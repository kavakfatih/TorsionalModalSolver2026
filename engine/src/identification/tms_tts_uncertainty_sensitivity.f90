module tms_tts_uncertainty_sensitivity
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_identification_result_t, &
    TTS_IDENTIFICATION_SUCCESS
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, &
    tts_uncertainty_sensitivity_configuration_t, &
    tts_uncertainty_sensitivity_result_t, &
    tts_uncertainty_validation_result_t, TTS_UNCERTAINTY_SUCCESS, &
    TTS_UNCERTAINTY_INVALID_INPUT, TTS_UNCERTAINTY_NO_SUPPORT, &
    TTS_WEIGHTED_L2_OBJECTIVE, TTS_STANDARDIZED_HUBER_OBJECTIVE
  use tms_tts_uncertainty_validation, only : &
    validate_tts_uncertainty_family
  use tms_tts_weighted_pair_shift, only : &
    identify_tts_uncertainty_pair_shift
  implicit none
  private

  public :: analyze_tts_uncertainty_sensitivity

contains

  !> Existing V0.8.1 identification'ın adjacent pair baseline shift'lerini
  !! değiştirmeden weighted-L2 ve standardized-Huber sensitivity evidence'ı
  !! üretir. Input intent(in)'dir; empirical shifts, master cloud ve runtime
  !! table yeniden kurulmaz. Başarı engineering acceptance/TRS kanıtı değildir.
  pure function analyze_tts_uncertainty_sensitivity( &
      identification, uncertainty_family, configuration) result(sensitivity)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_uncertainty_sensitivity_configuration_t), intent(in), &
      optional :: configuration
    type(tts_uncertainty_sensitivity_result_t) :: sensitivity

    type(tts_uncertainty_sensitivity_configuration_t) :: settings
    type(tts_uncertainty_validation_result_t) :: validation
    integer :: i
    integer :: reference_index
    integer :: moving_index

    if (present(configuration)) settings = configuration
    sensitivity%configuration = settings
    sensitivity%uncertainty_family_identifier = &
      uncertainty_family%family_identifier
    sensitivity%uncertainty_provenance = uncertainty_family%provenance
    if (.not. sensitivity_configuration_is_valid(settings)) then
      sensitivity%message = "Uncertainty sensitivity ayarları geçersiz."
      return
    end if
    if (identification%status /= TTS_IDENTIFICATION_SUCCESS .or. &
        .not. identification%shift_chain_available .or. &
        .not. allocated(identification%source_family%isotherms) .or. &
        .not. allocated(identification%pair_shift_results)) then
      sensitivity%message = "Başarılı V0.8.1 identification sonucu gerekli."
      return
    end if
    validation = validate_tts_uncertainty_family( &
      identification%source_family, uncertainty_family)
    if (.not. validation%valid) then
      sensitivity%status = validation%status
      sensitivity%message = validation%message
      return
    end if

    allocate(sensitivity%pair_results( &
      size(identification%pair_shift_results)))
    do i = 1, size(identification%pair_shift_results)
      reference_index = identification%pair_shift_results(i) &
        %reference_isotherm_index
      moving_index = identification%pair_shift_results(i) &
        %moving_isotherm_index
      if (reference_index < 1 .or. moving_index < 1 .or. &
          reference_index > size(identification%source_family%isotherms) .or. &
          moving_index > size(identification%source_family%isotherms)) then
        sensitivity%status = TTS_UNCERTAINTY_INVALID_INPUT
        sensitivity%message = "Baseline pair source indeksleri geçersiz."
        return
      end if
      call populate_pair_identity(sensitivity, identification, i, &
        reference_index, moving_index)
      sensitivity%pair_results(i)%weighted = &
        identify_tts_uncertainty_pair_shift( &
          identification%source_family%isotherms(reference_index), &
          identification%source_family%isotherms(moving_index), &
          uncertainty_family, TTS_WEIGHTED_L2_OBJECTIVE, &
          settings%pair_shift, settings%huber_c)
      sensitivity%pair_results(i)%huber = &
        identify_tts_uncertainty_pair_shift( &
          identification%source_family%isotherms(reference_index), &
          identification%source_family%isotherms(moving_index), &
          uncertainty_family, TTS_STANDARDIZED_HUBER_OBJECTIVE, &
          settings%pair_shift, settings%huber_c)
      call populate_pair_deltas(sensitivity, i)
    end do

    if (sensitivity%weighted_available_count == 0 .and. &
        sensitivity%huber_available_count == 0) then
      sensitivity%status = TTS_UNCERTAINTY_NO_SUPPORT
      sensitivity%message = "Pair-level weighted/Huber support bulunamadı."
      return
    end if
    sensitivity%status = TTS_UNCERTAINTY_SUCCESS
    sensitivity%message = &
      "Baseline korunarak uncertainty-weighted sensitivity üretildi."
  end function analyze_tts_uncertainty_sensitivity

  pure subroutine populate_pair_identity( &
      sensitivity, identification, pair_index, reference_index, moving_index)
    type(tts_uncertainty_sensitivity_result_t), intent(inout) :: sensitivity
    type(tts_identification_result_t), intent(in) :: identification
    integer, intent(in) :: pair_index
    integer, intent(in) :: reference_index
    integer, intent(in) :: moving_index

    sensitivity%pair_results(pair_index)%reference_isotherm_index = &
      reference_index
    sensitivity%pair_results(pair_index)%moving_isotherm_index = moving_index
    sensitivity%pair_results(pair_index)%reference_isotherm_identifier = &
      identification%source_family%isotherms(reference_index) &
        %isotherm_identifier
    sensitivity%pair_results(pair_index)%moving_isotherm_identifier = &
      identification%source_family%isotherms(moving_index)%isotherm_identifier
    sensitivity%pair_results(pair_index)%reference_temperature_k = &
      identification%source_family%isotherms(reference_index)%temperature_k
    sensitivity%pair_results(pair_index)%moving_temperature_k = &
      identification%source_family%isotherms(moving_index)%temperature_k
    sensitivity%pair_results(pair_index)%baseline_status = &
      identification%pair_shift_results(pair_index)%status
    sensitivity%pair_results(pair_index)%baseline_shift_available = &
      identification%pair_shift_results(pair_index)%shift_available
    sensitivity%pair_results(pair_index)%baseline_shift = &
      identification%pair_shift_results(pair_index)%delta_s
    sensitivity%pair_results(pair_index)%huber_c = &
      sensitivity%configuration%huber_c
    sensitivity%pair_results(pair_index)%uncertainty_provenance = &
      sensitivity%uncertainty_provenance
    if (sensitivity%pair_results(pair_index)%baseline_shift_available) &
      sensitivity%baseline_available_count = &
        sensitivity%baseline_available_count + 1
  end subroutine populate_pair_identity

  pure subroutine populate_pair_deltas(sensitivity, pair_index)
    type(tts_uncertainty_sensitivity_result_t), intent(inout) :: sensitivity
    integer, intent(in) :: pair_index

    if (sensitivity%pair_results(pair_index)%weighted%shift_available) then
      sensitivity%weighted_available_count = &
        sensitivity%weighted_available_count + 1
      if (sensitivity%pair_results(pair_index)%baseline_shift_available) then
        sensitivity%pair_results(pair_index)%weighted_delta_available = .true.
        sensitivity%pair_results(pair_index)%delta_weighted_vs_baseline = &
          sensitivity%pair_results(pair_index)%weighted%shift - &
          sensitivity%pair_results(pair_index)%baseline_shift
      end if
    end if
    if (sensitivity%pair_results(pair_index)%huber%shift_available) then
      sensitivity%huber_available_count = sensitivity%huber_available_count + 1
      if (sensitivity%pair_results(pair_index)%baseline_shift_available) then
        sensitivity%pair_results(pair_index) &
          %huber_baseline_delta_available = .true.
        sensitivity%pair_results(pair_index)%delta_huber_vs_baseline = &
          sensitivity%pair_results(pair_index)%huber%shift - &
          sensitivity%pair_results(pair_index)%baseline_shift
      end if
    end if
    if (sensitivity%pair_results(pair_index)%weighted%shift_available .and. &
        sensitivity%pair_results(pair_index)%huber%shift_available) then
      sensitivity%pair_results(pair_index) &
        %huber_weighted_delta_available = .true.
      sensitivity%pair_results(pair_index)%delta_huber_vs_weighted = &
        sensitivity%pair_results(pair_index)%huber%shift - &
        sensitivity%pair_results(pair_index)%weighted%shift
    end if
  end subroutine populate_pair_deltas

  pure function sensitivity_configuration_is_valid(settings) result(valid)
    type(tts_uncertainty_sensitivity_configuration_t), intent(in) :: settings
    logical :: valid

    valid = ieee_is_finite(settings%huber_c) .and. &
      settings%huber_c > 0.0_dp .and. &
      settings%pair_shift%coarse_scan_point_count >= 3 .and. &
      ieee_is_finite(settings%pair_shift%absolute_tolerance) .and. &
      ieee_is_finite(settings%pair_shift%relative_tolerance) .and. &
      settings%pair_shift%absolute_tolerance > 0.0_dp .and. &
      settings%pair_shift%relative_tolerance >= 0.0_dp .and. &
      settings%pair_shift%maximum_iterations > 0
  end function sensitivity_configuration_is_valid

end module tms_tts_uncertainty_sensitivity
