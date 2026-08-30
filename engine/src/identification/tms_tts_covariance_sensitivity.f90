module tms_tts_covariance_sensitivity
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_identification_result_t, &
    TTS_IDENTIFICATION_SUCCESS, are_tts_values_machine_equivalent
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, TTS_WEIGHTED_L2_OBJECTIVE, &
    DEFAULT_TTS_HUBER_C
  use tms_tts_weighted_pair_shift, only : &
    identify_tts_uncertainty_pair_shift
  use tms_tts_covariance_types, only : &
    tts_dynamic_modulus_covariance_family_t, &
    tts_covariance_sensitivity_configuration_t, &
    tts_covariance_sensitivity_result_t, &
    tts_covariance_validation_result_t, tts_covariance_pair_solution_t, &
    TTS_COVARIANCE_SUCCESS, TTS_COVARIANCE_INVALID_INPUT, &
    TTS_COVARIANCE_NO_BIVARIATE_SUPPORT
  use tms_tts_covariance_validation, only : &
    validate_tts_covariance_family
  use tms_tts_mahalanobis_pair_shift, only : &
    identify_tts_covariance_pair_shift
  implicit none
  private

  public :: analyze_tts_covariance_sensitivity

contains

  !> Existing V0.8.1 identification ve V0.8.4 uncertainty girdilerini
  !! değiştirmeden dört ayrı shift evidence'ı üretir: authoritative baseline,
  !! original weighted, O_B matched-diagonal ve Mahalanobis. Delta değerleri
  !! boyutsuz s=log10(a_T)'dir. Başarı preferred shift/TRS kabulü değildir.
  pure function analyze_tts_covariance_sensitivity( &
      identification, uncertainty_family, covariance_family, configuration) &
      result(sensitivity)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    type(tts_covariance_sensitivity_configuration_t), intent(in), &
      optional :: configuration
    type(tts_covariance_sensitivity_result_t) :: sensitivity

    type(tts_covariance_sensitivity_configuration_t) :: settings
    type(tts_covariance_validation_result_t) :: validation
    type(tts_covariance_pair_solution_t) :: covariance_solution
    integer :: i
    integer :: reference_index
    integer :: moving_index

    if (present(configuration)) settings = configuration
    sensitivity%configuration = settings
    sensitivity%covariance_family_identifier = &
      covariance_family%family_identifier
    sensitivity%covariance_provenance = covariance_family%provenance
    if (.not. configuration_is_valid(settings)) then
      sensitivity%message = "Covariance sensitivity ayarları geçersiz."
      return
    end if
    if (identification%status /= TTS_IDENTIFICATION_SUCCESS .or. &
        .not. identification%shift_chain_available .or. &
        .not. allocated(identification%source_family%isotherms) .or. &
        .not. allocated(identification%pair_shift_results)) then
      sensitivity%message = "Başarılı V0.8.1 identification sonucu gerekli."
      return
    end if
    validation = validate_tts_covariance_family( &
      identification%source_family, uncertainty_family, covariance_family)
    if (.not. validation%valid) then
      sensitivity%status = validation%status
      sensitivity%message = validation%message
      return
    end if

    sensitivity%pair_count = size(identification%pair_shift_results)
    allocate(sensitivity%pair_results(sensitivity%pair_count))
    do i = 1, sensitivity%pair_count
      reference_index = identification%pair_shift_results(i) &
        %reference_isotherm_index
      moving_index = identification%pair_shift_results(i) &
        %moving_isotherm_index
      if (reference_index < 1 .or. moving_index < 1 .or. &
          reference_index > size(identification%source_family%isotherms) .or. &
          moving_index > size(identification%source_family%isotherms)) then
        sensitivity%status = TTS_COVARIANCE_INVALID_INPUT
        sensitivity%message = "Baseline pair source indeksleri geçersiz."
        return
      end if
      call populate_pair_identity(sensitivity, identification, i, &
        reference_index, moving_index)
      sensitivity%pair_results(i)%weighted_original = &
        identify_tts_uncertainty_pair_shift( &
          identification%source_family%isotherms(reference_index), &
          identification%source_family%isotherms(moving_index), &
          uncertainty_family, TTS_WEIGHTED_L2_OBJECTIVE, &
          settings%pair_shift, DEFAULT_TTS_HUBER_C)
      covariance_solution = identify_tts_covariance_pair_shift( &
        identification%source_family%isotherms(reference_index), &
        identification%source_family%isotherms(moving_index), &
        uncertainty_family, covariance_family, settings%pair_shift)
      sensitivity%pair_results(i)%diagonal_matched = &
        covariance_solution%diagonal_matched
      sensitivity%pair_results(i)%mahalanobis = covariance_solution%mahalanobis
      call populate_pair_support_counts(sensitivity, covariance_family, i)
      call populate_pair_deltas(sensitivity, i)
    end do

    if (sensitivity%mahalanobis_available_count == 0 .and. &
        sensitivity%matched_diagonal_available_count == 0) then
      sensitivity%status = TTS_COVARIANCE_NO_BIVARIATE_SUPPORT
      sensitivity%message = "Pair-level bivariate common support bulunamadı."
      return
    end if
    sensitivity%status = TTS_COVARIANCE_SUCCESS
    sensitivity%message = &
      "Baseline korunarak point-local covariance sensitivity üretildi."
  end function analyze_tts_covariance_sensitivity

  pure subroutine populate_pair_identity( &
      sensitivity, identification, pair_index, reference_index, moving_index)
    type(tts_covariance_sensitivity_result_t), intent(inout) :: sensitivity
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
  end subroutine populate_pair_identity

  pure subroutine populate_pair_support_counts( &
      sensitivity, covariance_family, pair_index)
    type(tts_covariance_sensitivity_result_t), intent(inout) :: sensitivity
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    integer, intent(in) :: pair_index
    integer :: i
    integer :: j
    logical :: pair_temperature

    do i = 1, size(covariance_family%isotherms)
      pair_temperature = are_tts_values_machine_equivalent( &
        covariance_family%isotherms(i)%temperature_k, &
        sensitivity%pair_results(pair_index)%reference_temperature_k) .or. &
        are_tts_values_machine_equivalent( &
        covariance_family%isotherms(i)%temperature_k, &
        sensitivity%pair_results(pair_index)%moving_temperature_k)
      if (.not. pair_temperature .or. &
          .not. allocated(covariance_family%isotherms(i)%points)) cycle
      do j = 1, size(covariance_family%isotherms(i)%points)
        if (covariance_family%isotherms(i)%points(j)%covariance_available) then
          sensitivity%pair_results(pair_index)%covariance_point_count = &
            sensitivity%pair_results(pair_index)%covariance_point_count + 1
        else
          sensitivity%pair_results(pair_index)%covariance_gap_count = &
            sensitivity%pair_results(pair_index)%covariance_gap_count + 1
        end if
      end do
    end do
  end subroutine populate_pair_support_counts

  pure subroutine populate_pair_deltas(sensitivity, pair_index)
    type(tts_covariance_sensitivity_result_t), intent(inout) :: sensitivity
    integer, intent(in) :: pair_index

    if (sensitivity%pair_results(pair_index) &
        %diagonal_matched%shift_available) then
      sensitivity%matched_diagonal_available_count = &
        sensitivity%matched_diagonal_available_count + 1
      sensitivity%pair_results(pair_index)%bivariate_overlap_width_decades = &
        sensitivity%pair_results(pair_index) &
          %diagonal_matched%diagnostics%overlap_width_decades
    end if
    if (sensitivity%pair_results(pair_index)%mahalanobis%shift_available) then
      sensitivity%mahalanobis_available_count = &
        sensitivity%mahalanobis_available_count + 1
      sensitivity%pair_results(pair_index)%bivariate_overlap_width_decades = &
        sensitivity%pair_results(pair_index) &
          %mahalanobis%diagnostics%overlap_width_decades
    end if
    if (.not. sensitivity%pair_results(pair_index) &
        %diagonal_matched%shift_available .and. &
        .not. sensitivity%pair_results(pair_index) &
          %mahalanobis%shift_available) then
      sensitivity%unsupported_pair_count = &
        sensitivity%unsupported_pair_count + 1
    end if

    if (sensitivity%pair_results(pair_index) &
        %weighted_original%shift_available .and. &
        sensitivity%pair_results(pair_index) &
          %diagonal_matched%shift_available) then
      sensitivity%pair_results(pair_index)%delta_support_available = .true.
      sensitivity%pair_results(pair_index)%delta_support = &
        sensitivity%pair_results(pair_index)%diagonal_matched%shift - &
        sensitivity%pair_results(pair_index)%weighted_original%shift
    end if
    if (sensitivity%pair_results(pair_index) &
        %diagonal_matched%shift_available .and. &
        sensitivity%pair_results(pair_index)%mahalanobis%shift_available) then
      sensitivity%pair_results(pair_index)%delta_covariance_available = .true.
      sensitivity%pair_results(pair_index)%delta_covariance = &
        sensitivity%pair_results(pair_index)%mahalanobis%shift - &
        sensitivity%pair_results(pair_index)%diagonal_matched%shift
    end if
    if (sensitivity%pair_results(pair_index) &
        %weighted_original%shift_available .and. &
        sensitivity%pair_results(pair_index)%mahalanobis%shift_available) then
      sensitivity%pair_results(pair_index)%delta_total_available = .true.
      sensitivity%pair_results(pair_index)%delta_total = &
        sensitivity%pair_results(pair_index)%mahalanobis%shift - &
        sensitivity%pair_results(pair_index)%weighted_original%shift
    end if
    if (sensitivity%pair_results(pair_index)%baseline_shift_available .and. &
        sensitivity%pair_results(pair_index)%mahalanobis%shift_available) then
      sensitivity%pair_results(pair_index) &
        %delta_mahalanobis_vs_baseline_available = .true.
      sensitivity%pair_results(pair_index)%delta_mahalanobis_vs_baseline = &
        sensitivity%pair_results(pair_index)%mahalanobis%shift - &
        sensitivity%pair_results(pair_index)%baseline_shift
    end if
  end subroutine populate_pair_deltas

  pure function configuration_is_valid(settings) result(valid)
    type(tts_covariance_sensitivity_configuration_t), intent(in) :: settings
    logical :: valid

    valid = settings%pair_shift%coarse_scan_point_count >= 3 .and. &
      ieee_is_finite(settings%pair_shift%absolute_tolerance) .and. &
      ieee_is_finite(settings%pair_shift%relative_tolerance) .and. &
      settings%pair_shift%absolute_tolerance > 0.0_dp .and. &
      settings%pair_shift%relative_tolerance >= 0.0_dp .and. &
      settings%pair_shift%maximum_iterations > 0
  end function configuration_is_valid

end module tms_tts_covariance_sensitivity
