program test_tts_robust_shift_sensitivity
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_identification_result_t, TTS_IDENTIFICATION_SUCCESS
  use tms_tts_identification, only : identify_tts_master_curve
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, &
    tts_uncertainty_sensitivity_configuration_t, &
    tts_uncertainty_sensitivity_result_t, TTS_UNCERTAINTY_SUCCESS
  use tms_tts_uncertainty_sensitivity, only : &
    analyze_tts_uncertainty_sensitivity
  use tms_tts_test_support, only : make_generalized_maxwell_trs_family, &
    assert_true, assert_close
  use tms_tts_uncertainty_test_support, only : &
    make_relative_uncertainty_family
  implicit none

  type(tts_material_family_t) :: family
  type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty
  type(tts_identification_result_t) :: identification
  type(tts_identification_result_t) :: snapshot
  type(tts_uncertainty_sensitivity_configuration_t) :: configuration
  type(tts_uncertainty_sensitivity_result_t) :: sensitivity
  real(dp), parameter :: temperatures(2) = [293.15_dp, 313.15_dp]
  real(dp), parameter :: truth_shift = 0.65_dp
  integer, parameter :: anomaly_index = 31

  ! Clean exact-TRS fixture'ında optimum çevresindeki standardized residual
  ! quadratic regime'dedir; Huber ve weighted-L2 aynı argmin'i vermelidir.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  identification = identify_tts_master_curve(family, "ISO-1")
  call assert_true(identification%status == TTS_IDENTIFICATION_SUCCESS, &
    "Clean sensitivity fixture baseline identification başarısız.")
  uncertainty = make_relative_uncertainty_family(family, 0.05_dp, 0.05_dp)
  sensitivity = analyze_tts_uncertainty_sensitivity( &
    identification, uncertainty, configuration)
  call assert_true(sensitivity%status == TTS_UNCERTAINTY_SUCCESS .and. &
    sensitivity%weighted_available_count == 1 .and. &
    sensitivity%huber_available_count == 1, &
    "Clean weighted/Huber top-level sensitivity üretilemedi.")
  call assert_close(sensitivity%pair_results(1)%huber%shift, &
    sensitivity%pair_results(1)%weighted%shift, 2.0e-6_dp, &
    "Clean-data Huber ve weighted minimizer'ları ayrıldı.")
  call assert_true( &
    sensitivity%pair_results(1)%huber%production_diagnostics &
      %tail_fraction < 1.0e-8_dp, &
    "Clean exact-TRS optimum Huber tail üretti.")

  ! Localized gross residual anomaly ölçümden silinmez. Comparable pointwise
  ! uncertainty altında weighted L2 anomaly'den etkilenirken standardized
  ! Huber tail influence'ı sınırlar ve bu fixture'da known truth'e yaklaşır.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  family%isotherms(2)%points(anomaly_index)%storage_modulus_pa = &
    family%isotherms(2)%points(anomaly_index)%storage_modulus_pa*100.0_dp
  family%isotherms(2)%points(anomaly_index)%loss_modulus_pa = &
    family%isotherms(2)%points(anomaly_index)%loss_modulus_pa*100.0_dp
  identification = identify_tts_master_curve(family, "ISO-1")
  call assert_true(identification%status == TTS_IDENTIFICATION_SUCCESS, &
    "Anomaly fixture authoritative baseline identification başarısız.")
  snapshot = identification
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  sensitivity = analyze_tts_uncertainty_sensitivity( &
    identification, uncertainty, configuration)
  call assert_true(sensitivity%status == TTS_UNCERTAINTY_SUCCESS, &
    "Anomaly fixture sensitivity sonucu üretilemedi.")
  call assert_true(abs(sensitivity%pair_results(1)%huber%shift - truth_shift) < &
    abs(sensitivity%pair_results(1)%weighted%shift - truth_shift), &
    "Standardized Huber anomaly influence'ını bu fixture'da azaltmadı.")
  call assert_true( &
    sensitivity%pair_results(1)%huber%production_diagnostics &
      %tail_fraction > 0.0_dp, &
    "Localized anomaly Huber tail diagnostic'i üretmedi.")
  call assert_true(sensitivity%pair_results(1)%weighted_delta_available .and. &
    sensitivity%pair_results(1)%huber_baseline_delta_available .and. &
    sensitivity%pair_results(1)%huber_weighted_delta_available, &
    "Baseline/weighted/Huber sensitivity delta availability eksik.")
  call assert_close(sensitivity%pair_results(1)%delta_weighted_vs_baseline, &
    sensitivity%pair_results(1)%weighted%shift - &
      sensitivity%pair_results(1)%baseline_shift, 1.0e-15_dp, &
    "Weighted-baseline delta tanımı hatalı.")
  call assert_close(sensitivity%pair_results(1)%delta_huber_vs_baseline, &
    sensitivity%pair_results(1)%huber%shift - &
      sensitivity%pair_results(1)%baseline_shift, 1.0e-15_dp, &
    "Huber-baseline delta tanımı hatalı.")
  call assert_close(sensitivity%pair_results(1)%delta_huber_vs_weighted, &
    sensitivity%pair_results(1)%huber%shift - &
      sensitivity%pair_results(1)%weighted%shift, 1.0e-15_dp, &
    "Huber-weighted delta tanımı hatalı.")
  call assert_close(sensitivity%pair_results(1)%huber_c, 1.345_dp, &
    0.0_dp, "Exact Huber c provenance sonucu içinde korunmadı.")
  call assert_true(.not. sensitivity%cross_channel_covariance_modeled .and. &
    .not. sensitivity%cross_isotherm_covariance_modeled .and. &
    .not. sensitivity%pair_results(1)%cross_channel_covariance_modeled .and. &
    .not. sensitivity%pair_results(1)%cross_isotherm_covariance_modeled, &
    "V0.8.4 diagonal covariance assumptions açık değil.")

  ! Additive analiz authoritative V0.8.1 pair/empirical/master/runtime
  ! dizilerinin tek bitini değiştirmemeli ve anomaly point'i silmemelidir.
  call assert_identification_unchanged(identification, snapshot)

  print *, "V0.8.4 robust shift sensitivity ve immutability doğrulandı."

contains

  subroutine assert_identification_unchanged(actual, expected)
    type(tts_identification_result_t), intent(in) :: actual
    type(tts_identification_result_t), intent(in) :: expected
    integer :: i
    integer :: j

    call assert_true(size(actual%pair_shift_results) == &
      size(expected%pair_shift_results), "Pair result array boyutu değişti.")
    do i = 1, size(actual%pair_shift_results)
      call assert_true(actual%pair_shift_results(i)%status == &
        expected%pair_shift_results(i)%status .and. &
        actual%pair_shift_results(i)%reference_isotherm_identifier == &
        expected%pair_shift_results(i)%reference_isotherm_identifier .and. &
        actual%pair_shift_results(i)%moving_isotherm_identifier == &
        expected%pair_shift_results(i)%moving_isotherm_identifier .and. &
        actual%pair_shift_results(i)%reference_isotherm_index == &
        expected%pair_shift_results(i)%reference_isotherm_index .and. &
        actual%pair_shift_results(i)%moving_isotherm_index == &
        expected%pair_shift_results(i)%moving_isotherm_index .and. &
        actual%pair_shift_results(i)%production_channel == &
        expected%pair_shift_results(i)%production_channel .and. &
        (actual%pair_shift_results(i)%shift_available .eqv. &
        expected%pair_shift_results(i)%shift_available) .and. &
        (actual%pair_shift_results(i)%joint_shift_available .eqv. &
        expected%pair_shift_results(i)%joint_shift_available) .and. &
        (actual%pair_shift_results(i)%storage_shift_available .eqv. &
        expected%pair_shift_results(i)%storage_shift_available) .and. &
        (actual%pair_shift_results(i)%loss_shift_available .eqv. &
        expected%pair_shift_results(i)%loss_shift_available) .and. &
        actual%pair_shift_results(i)%iteration_count == &
        expected%pair_shift_results(i)%iteration_count .and. &
        actual%pair_shift_results(i)%evaluation_count == &
        expected%pair_shift_results(i)%evaluation_count, &
        "Pair result discrete alanları mutate edildi.")
      call assert_true(same_real_bits(actual%pair_shift_results(i)%delta_s, &
        expected%pair_shift_results(i)%delta_s) .and. &
        same_real_bits(actual%pair_shift_results(i)%delta_s_joint, &
        expected%pair_shift_results(i)%delta_s_joint) .and. &
        same_real_bits(actual%pair_shift_results(i)%delta_s_storage, &
        expected%pair_shift_results(i)%delta_s_storage) .and. &
        same_real_bits(actual%pair_shift_results(i)%delta_s_loss, &
        expected%pair_shift_results(i)%delta_s_loss) .and. &
        same_real_bits(actual%pair_shift_results(i) &
          %storage_loss_shift_discrepancy, expected%pair_shift_results(i) &
          %storage_loss_shift_discrepancy) .and. &
        same_real_bits(actual%pair_shift_results(i)%objective_minimum, &
        expected%pair_shift_results(i)%objective_minimum) .and. &
        same_real_bits(actual%pair_shift_results(i)%overlap_width_decades, &
        expected%pair_shift_results(i)%overlap_width_decades) .and. &
        same_real_bits(actual%pair_shift_results(i)%overlap_fraction, &
        expected%pair_shift_results(i)%overlap_fraction) .and. &
        same_real_bits(actual%pair_shift_results(i) &
          %storage_overlap_width_decades, expected%pair_shift_results(i) &
          %storage_overlap_width_decades) .and. &
        same_real_bits(actual%pair_shift_results(i) &
          %loss_overlap_width_decades, expected%pair_shift_results(i) &
          %loss_overlap_width_decades) .and. &
        same_real_bits(actual%pair_shift_results(i)%objective_curvature, &
        expected%pair_shift_results(i)%objective_curvature), &
        "Pair result real alanları mutate edildi.")
    end do

    call assert_true(size(actual%empirical_shifts) == &
      size(expected%empirical_shifts), "Empirical shift array boyutu değişti.")
    do i = 1, size(actual%empirical_shifts)
      call assert_true(actual%empirical_shifts(i)%source_isotherm_index == &
        expected%empirical_shifts(i)%source_isotherm_index .and. &
        actual%empirical_shifts(i)%source_isotherm_identifier == &
        expected%empirical_shifts(i)%source_isotherm_identifier .and. &
        same_real_bits(actual%empirical_shifts(i)%temperature_k, &
        expected%empirical_shifts(i)%temperature_k) .and. &
        same_real_bits(actual%empirical_shifts(i)%log10_a_t, &
        expected%empirical_shifts(i)%log10_a_t) .and. &
        same_real_bits(actual%empirical_shifts(i)%a_t, &
        expected%empirical_shifts(i)%a_t), &
        "Empirical shift array mutate edildi.")
    end do

    call assert_true(size(actual%master_cloud) == size(expected%master_cloud), &
      "Master cloud array boyutu değişti.")
    do i = 1, size(actual%master_cloud)
      call assert_true(actual%master_cloud(i)%source_isotherm_index == &
        expected%master_cloud(i)%source_isotherm_index .and. &
        actual%master_cloud(i)%source_point_index == &
        expected%master_cloud(i)%source_point_index .and. &
        actual%master_cloud(i)%source_isotherm_identifier == &
        expected%master_cloud(i)%source_isotherm_identifier .and. &
        actual%master_cloud(i)%specimen_identifier == &
        expected%master_cloud(i)%specimen_identifier .and. &
        actual%master_cloud(i)%source_identifier == &
        expected%master_cloud(i)%source_identifier .and. &
        same_real_bits(actual%master_cloud(i)%source_temperature_k, &
        expected%master_cloud(i)%source_temperature_k) .and. &
        same_real_bits(actual%master_cloud(i)%source_frequency_hz, &
        expected%master_cloud(i)%source_frequency_hz) .and. &
        same_real_bits(actual%master_cloud(i)%log10_a_t, &
        expected%master_cloud(i)%log10_a_t) .and. &
        same_real_bits(actual%master_cloud(i)%reduced_frequency_hz, &
        expected%master_cloud(i)%reduced_frequency_hz) .and. &
        same_real_bits(actual%master_cloud(i)%storage_modulus_pa, &
        expected%master_cloud(i)%storage_modulus_pa) .and. &
        same_real_bits(actual%master_cloud(i)%loss_modulus_pa, &
        expected%master_cloud(i)%loss_modulus_pa) .and. &
        actual%master_cloud(i)%storage_quality == &
        expected%master_cloud(i)%storage_quality .and. &
        actual%master_cloud(i)%loss_quality == &
        expected%master_cloud(i)%loss_quality .and. &
        (actual%master_cloud(i)%contributes_to_validation .eqv. &
        expected%master_cloud(i)%contributes_to_validation) .and. &
        (actual%master_cloud(i)%contributes_to_runtime_extension .eqv. &
        expected%master_cloud(i)%contributes_to_runtime_extension), &
        "Master cloud array mutate edildi.")
    end do

    call assert_true(size(actual%runtime_master_table) == &
      size(expected%runtime_master_table), &
      "Runtime master table boyutu değişti.")
    do i = 1, size(actual%runtime_master_table)
      call assert_true(same_real_bits( &
        actual%runtime_master_table(i)%reduced_frequency_hz, &
        expected%runtime_master_table(i)%reduced_frequency_hz) .and. &
        same_real_bits(actual%runtime_master_table(i)%storage_modulus_pa, &
        expected%runtime_master_table(i)%storage_modulus_pa) .and. &
        same_real_bits(actual%runtime_master_table(i)%loss_modulus_pa, &
        expected%runtime_master_table(i)%loss_modulus_pa) .and. &
        actual%runtime_master_table(i)%source_isotherm_index == &
        expected%runtime_master_table(i)%source_isotherm_index .and. &
        actual%runtime_master_table(i)%source_point_index == &
        expected%runtime_master_table(i)%source_point_index .and. &
        actual%runtime_master_table(i)%source_isotherm_identifier == &
        expected%runtime_master_table(i)%source_isotherm_identifier, &
        "Runtime master table mutate edildi.")
    end do

    call assert_true(size(actual%source_family%isotherms) == &
      size(expected%source_family%isotherms), &
      "Source family isotherm sayısı değişti.")
    call assert_true(actual%source_family%family_identifier == &
      expected%source_family%family_identifier .and. &
      actual%source_family%common_state%material_identifier == &
      expected%source_family%common_state%material_identifier .and. &
      actual%source_family%common_state%batch_state_identifier == &
      expected%source_family%common_state%batch_state_identifier .and. &
      actual%source_family%common_state%deformation_mode == &
      expected%source_family%common_state%deformation_mode .and. &
      actual%source_family%common_state%conditioning_description == &
      expected%source_family%common_state%conditioning_description .and. &
      actual%source_family%common_state%test_method == &
      expected%source_family%common_state%test_method .and. &
      actual%source_family%common_state%source_metadata == &
      expected%source_family%common_state%source_metadata .and. &
      same_real_bits(actual%source_family%common_state &
        %dynamic_strain_amplitude_ratio, expected%source_family%common_state &
        %dynamic_strain_amplitude_ratio) .and. &
      same_real_bits(actual%source_family%common_state%static_prestrain_ratio, &
        expected%source_family%common_state%static_prestrain_ratio), &
      "Authoritative source family/common state mutate edildi.")
    do i = 1, size(actual%source_family%isotherms)
      call assert_true(actual%source_family%isotherms(i) &
        %isotherm_identifier == expected%source_family%isotherms(i) &
        %isotherm_identifier .and. actual%source_family%isotherms(i) &
        %specimen_identifier == expected%source_family%isotherms(i) &
        %specimen_identifier .and. actual%source_family%isotherms(i) &
        %source_identifier == expected%source_family%isotherms(i) &
        %source_identifier .and. same_real_bits( &
        actual%source_family%isotherms(i)%temperature_k, &
        expected%source_family%isotherms(i)%temperature_k), &
        "Authoritative source isotherm metadata mutate edildi.")
      call assert_true(size(actual%source_family%isotherms(i)%points) == &
        size(expected%source_family%isotherms(i)%points), &
        "Source measurement point sayısı değişti.")
      do j = 1, size(actual%source_family%isotherms(i)%points)
        call assert_true(same_real_bits( &
          actual%source_family%isotherms(i)%points(j)%frequency_hz, &
          expected%source_family%isotherms(i)%points(j)%frequency_hz) .and. &
          same_real_bits(actual%source_family%isotherms(i)%points(j) &
            %storage_modulus_pa, expected%source_family%isotherms(i)%points(j) &
            %storage_modulus_pa) .and. same_real_bits( &
          actual%source_family%isotherms(i)%points(j)%loss_modulus_pa, &
          expected%source_family%isotherms(i)%points(j)%loss_modulus_pa) .and. &
          actual%source_family%isotherms(i)%points(j)%storage_quality == &
          expected%source_family%isotherms(i)%points(j)%storage_quality .and. &
          actual%source_family%isotherms(i)%points(j)%loss_quality == &
          expected%source_family%isotherms(i)%points(j)%loss_quality, &
          "Authoritative source measurement mutate/silindi.")
      end do
    end do
  end subroutine assert_identification_unchanged

  pure function same_real_bits(a, b) result(same)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: same

    same = transfer(a, 0_int64) == transfer(b, 0_int64)
  end function same_real_bits

end program test_tts_robust_shift_sensitivity
