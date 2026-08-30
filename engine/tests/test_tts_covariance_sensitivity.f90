program test_tts_covariance_sensitivity
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_identification_result_t, TTS_IDENTIFICATION_SUCCESS
  use tms_tts_identification, only : identify_tts_master_curve
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t
  use tms_tts_covariance_types, only : &
    tts_dynamic_modulus_covariance_family_t, &
    tts_covariance_sensitivity_configuration_t, &
    tts_covariance_sensitivity_result_t, TTS_COVARIANCE_SUCCESS
  use tms_tts_covariance_sensitivity, only : &
    analyze_tts_covariance_sensitivity
  use tms_tts_test_support, only : make_exact_trs_family, &
    make_generalized_maxwell_trs_family, replace_loss_truth_shift, &
    assert_true, assert_close
  use tms_tts_uncertainty_test_support, only : &
    make_relative_uncertainty_family
  use tms_tts_covariance_test_support, only : make_covariance_family
  implicit none

  type(tts_material_family_t) :: family
  type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty
  type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty_snapshot
  type(tts_dynamic_modulus_covariance_family_t) :: covariance
  type(tts_dynamic_modulus_covariance_family_t) :: covariance_snapshot
  type(tts_dynamic_modulus_covariance_family_t) :: correlated_covariance
  type(tts_identification_result_t) :: identification
  type(tts_identification_result_t) :: identification_snapshot
  type(tts_covariance_sensitivity_configuration_t) :: configuration
  type(tts_covariance_sensitivity_result_t) :: support_sensitivity
  type(tts_covariance_sensitivity_result_t) :: diagonal_sensitivity
  type(tts_covariance_sensitivity_result_t) :: correlated_sensitivity
  real(dp), parameter :: temperatures(2) = [293.15_dp, 313.15_dp]
  real(dp), parameter :: truth_shift = 0.65_dp
  integer :: i

  ! Support-decomposition fixture: original V0.8.4 uncertainty support geniştir
  ! ve high-frequency biased bölgeyi içerir. V0.8.5 covariance support'u bu
  ! bölgeyi gap yapar. Off-diagonal sıfırken matched ile Mahalanobis aynı kalır;
  ! gözlenen fark yalnız delta_support içinde olmalıdır.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  do i = 49, size(family%isotherms(2)%points)
    family%isotherms(2)%points(i)%storage_modulus_pa = &
      family%isotherms(2)%points(i)%storage_modulus_pa*10.0_dp**0.35_dp
    family%isotherms(2)%points(i)%loss_modulus_pa = &
      family%isotherms(2)%points(i)%loss_modulus_pa*10.0_dp**0.35_dp
  end do
  identification = identify_tts_master_curve(family, "ISO-1")
  call assert_true(identification%status == TTS_IDENTIFICATION_SUCCESS, &
    "Support fixture baseline identification başarısız.")
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  covariance = make_covariance_family(family, uncertainty, 0.0_dp)
  do i = 46, size(covariance%isotherms(1)%points)
    covariance%isotherms(1)%points(i)%covariance_available = .false.
    covariance%isotherms(2)%points(i)%covariance_available = .false.
  end do
  identification_snapshot = identification
  uncertainty_snapshot = uncertainty
  covariance_snapshot = covariance
  support_sensitivity = analyze_tts_covariance_sensitivity( &
    identification, uncertainty, covariance, configuration)
  call assert_true(support_sensitivity%status == TTS_COVARIANCE_SUCCESS .and. &
    support_sensitivity%mahalanobis_available_count == 1 .and. &
    support_sensitivity%matched_diagonal_available_count == 1, &
    "Support-decomposition top-level sensitivity üretilemedi.")
  call assert_true(support_sensitivity%pair_results(1) &
      %baseline_shift_available .and. &
    support_sensitivity%pair_results(1) &
      %weighted_original%shift_available .and. &
    support_sensitivity%pair_results(1) &
      %diagonal_matched%shift_available .and. &
    support_sensitivity%pair_results(1)%mahalanobis%shift_available, &
    "Dört ayrı shift hierarchy alanı birlikte üretilmedi.")
  call assert_close(support_sensitivity%pair_results(1)%baseline_shift, &
    identification%pair_shift_results(1)%delta_s, 0.0_dp, &
    "V0.8.1 authoritative baseline sonuçta korunmadı.")
  call assert_true(abs(support_sensitivity%pair_results(1)%delta_support) > &
    1.0e-4_dp, "Daralan covariance support delta_support üretmedi.")
  call assert_close(support_sensitivity%pair_results(1)%mahalanobis%shift, &
    support_sensitivity%pair_results(1)%diagonal_matched%shift, &
    5.0e-8_dp, "Zero covariance matched/Mahalanobis minimizer'ları ayrıldı.")
  call assert_close(support_sensitivity%pair_results(1)%delta_covariance, &
    0.0_dp, 5.0e-8_dp, &
    "Zero covariance fixture saf covariance delta üretti.")
  call assert_close(support_sensitivity%pair_results(1)%delta_total, &
    support_sensitivity%pair_results(1)%delta_support + &
      support_sensitivity%pair_results(1)%delta_covariance, &
    2.0e-14_dp, "Support+covariance total decomposition identity bozuldu.")
  call assert_true(support_sensitivity%pair_results(1) &
      %first_order_covariance_propagation .and. &
    support_sensitivity%pair_results(1)%cross_channel_covariance_modeled .and. &
    .not. support_sensitivity%pair_results(1) &
      %cross_isotherm_covariance_modeled .and. &
    .not. support_sensitivity%pair_results(1) &
      %cross_frequency_covariance_modeled .and. &
    .not. support_sensitivity%pair_results(1) &
      %automatic_regularization_used, &
    "V0.8.5 scope/provenance flags açık değil.")
  call assert_identification_unchanged(identification, &
    identification_snapshot)
  call assert_uncertainty_unchanged(uncertainty, uncertainty_snapshot)
  call assert_covariance_unchanged(covariance, covariance_snapshot)

  ! Pure-covariance fixture: original ve bivariate support birebir aynıdır.
  ! Storage 0.60, loss 0.20 shift tercih eder. Matched-diagonal V0.8.4 original
  ! ile çakışırken positive rho yalnız covariance geometry'siyle shift'i artırır.
  family = make_exact_trs_family(temperatures, [0.0_dp, 0.60_dp])
  call replace_loss_truth_shift(family%isotherms(2), 0.20_dp)
  identification = identify_tts_master_curve(family, "ISO-1")
  call assert_true(identification%status == TTS_IDENTIFICATION_SUCCESS, &
    "Pure-covariance fixture baseline identification başarısız.")
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  covariance = make_covariance_family(family, uncertainty, 0.0_dp)
  correlated_covariance = make_covariance_family(family, uncertainty, 0.5_dp)
  diagonal_sensitivity = analyze_tts_covariance_sensitivity( &
    identification, uncertainty, covariance, configuration)
  correlated_sensitivity = analyze_tts_covariance_sensitivity( &
    identification, uncertainty, correlated_covariance, configuration)
  call assert_true(diagonal_sensitivity%status == TTS_COVARIANCE_SUCCESS .and. &
    correlated_sensitivity%status == TTS_COVARIANCE_SUCCESS, &
    "Pure-covariance sensitivity çözümleri üretilemedi.")
  call assert_close(diagonal_sensitivity%pair_results(1)%delta_support, &
    0.0_dp, 7.0e-8_dp, &
    "Identical support fixture V0.8.4/matched support farkı üretti.")
  call assert_close(correlated_sensitivity%pair_results(1)%delta_support, &
    0.0_dp, 7.0e-8_dp, &
    "Correlation değişimi matched support delta'sını değiştirdi.")
  call assert_close(diagonal_sensitivity%pair_results(1)%delta_covariance, &
    0.0_dp, 7.0e-8_dp, &
    "Diagonal pure-covariance kontrolü sıfır delta vermedi.")
  call assert_true(correlated_sensitivity%pair_results(1)%delta_covariance > &
    1.0e-3_dp, &
    "Positive correlation predicted positive covariance shift delta üretmedi.")
  call assert_true(correlated_sensitivity%pair_results(1) &
      %mahalanobis%shift > &
    diagonal_sensitivity%pair_results(1)%mahalanobis%shift, &
    "Pure covariance effect analytical yönde değil.")
  call assert_close(correlated_sensitivity%pair_results(1)%delta_total, &
    correlated_sensitivity%pair_results(1)%delta_support + &
      correlated_sensitivity%pair_results(1)%delta_covariance, &
    2.0e-14_dp, "Pure covariance total delta identity bozuldu.")
  call assert_close(correlated_sensitivity%pair_results(1) &
      %delta_mahalanobis_vs_baseline, &
    correlated_sensitivity%pair_results(1)%mahalanobis%shift - &
      correlated_sensitivity%pair_results(1)%baseline_shift, &
    2.0e-14_dp, "Mahalanobis-baseline delta tanımı hatalı.")

  print *, "V0.8.5 covariance sensitivity decomposition doğrulandı."

contains

  subroutine assert_identification_unchanged(actual, expected)
    type(tts_identification_result_t), intent(in) :: actual
    type(tts_identification_result_t), intent(in) :: expected
    integer :: isotherm_index
    integer :: point_index
    integer :: result_index

    call assert_true(size(actual%pair_shift_results) == &
      size(expected%pair_shift_results), "Pair result boyutu mutate edildi.")
    do result_index = 1, size(actual%pair_shift_results)
      call assert_true(actual%pair_shift_results(result_index)%status == &
        expected%pair_shift_results(result_index)%status .and. &
        same_real_bits(actual%pair_shift_results(result_index)%delta_s, &
        expected%pair_shift_results(result_index)%delta_s) .and. &
        same_real_bits(actual%pair_shift_results(result_index) &
          %objective_minimum, expected%pair_shift_results(result_index) &
          %objective_minimum), "Authoritative pair result mutate edildi.")
    end do
    call assert_true(size(actual%empirical_shifts) == &
      size(expected%empirical_shifts), "Empirical shift boyutu mutate edildi.")
    do result_index = 1, size(actual%empirical_shifts)
      call assert_true(same_real_bits( &
        actual%empirical_shifts(result_index)%log10_a_t, &
        expected%empirical_shifts(result_index)%log10_a_t) .and. &
        same_real_bits(actual%empirical_shifts(result_index)%a_t, &
        expected%empirical_shifts(result_index)%a_t), &
        "Empirical shift mutate edildi.")
    end do
    call assert_true(size(actual%master_cloud) == size(expected%master_cloud), &
      "Master cloud boyutu mutate edildi.")
    do result_index = 1, size(actual%master_cloud)
      call assert_true(same_real_bits( &
        actual%master_cloud(result_index)%reduced_frequency_hz, &
        expected%master_cloud(result_index)%reduced_frequency_hz) .and. &
        same_real_bits(actual%master_cloud(result_index)%storage_modulus_pa, &
        expected%master_cloud(result_index)%storage_modulus_pa) .and. &
        same_real_bits(actual%master_cloud(result_index)%loss_modulus_pa, &
        expected%master_cloud(result_index)%loss_modulus_pa), &
        "Master cloud mutate edildi.")
    end do
    call assert_true(size(actual%runtime_master_table) == &
      size(expected%runtime_master_table), &
      "Runtime master table boyutu mutate edildi.")
    do result_index = 1, size(actual%runtime_master_table)
      call assert_true(same_real_bits(actual%runtime_master_table(result_index) &
        %reduced_frequency_hz, expected%runtime_master_table(result_index) &
        %reduced_frequency_hz) .and. same_real_bits( &
        actual%runtime_master_table(result_index)%storage_modulus_pa, &
        expected%runtime_master_table(result_index)%storage_modulus_pa) .and. &
        same_real_bits(actual%runtime_master_table(result_index) &
        %loss_modulus_pa, expected%runtime_master_table(result_index) &
        %loss_modulus_pa), "Runtime master table mutate edildi.")
    end do
    do isotherm_index = 1, size(actual%source_family%isotherms)
      call assert_true(same_real_bits(actual%source_family &
        %isotherms(isotherm_index)%temperature_k, expected%source_family &
        %isotherms(isotherm_index)%temperature_k), &
        "Source isotherm sıcaklığı mutate edildi.")
      do point_index = 1, size(actual%source_family &
          %isotherms(isotherm_index)%points)
        call assert_true(same_real_bits(actual%source_family &
          %isotherms(isotherm_index)%points(point_index)%frequency_hz, &
          expected%source_family%isotherms(isotherm_index) &
          %points(point_index)%frequency_hz) .and. same_real_bits( &
          actual%source_family%isotherms(isotherm_index) &
          %points(point_index)%storage_modulus_pa, expected%source_family &
          %isotherms(isotherm_index)%points(point_index) &
          %storage_modulus_pa) .and. same_real_bits(actual%source_family &
          %isotherms(isotherm_index)%points(point_index)%loss_modulus_pa, &
          expected%source_family%isotherms(isotherm_index) &
          %points(point_index)%loss_modulus_pa), &
          "Authoritative source family measurement mutate edildi.")
      end do
    end do
  end subroutine assert_identification_unchanged

  subroutine assert_uncertainty_unchanged(actual, expected)
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: actual
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: expected
    integer :: i
    integer :: j

    do i = 1, size(actual%isotherms)
      do j = 1, size(actual%isotherms(i)%points)
        call assert_true(same_real_bits(actual%isotherms(i)%points(j) &
          %storage_standard_uncertainty_pa, expected%isotherms(i)%points(j) &
          %storage_standard_uncertainty_pa) .and. same_real_bits( &
          actual%isotherms(i)%points(j)%loss_standard_uncertainty_pa, &
          expected%isotherms(i)%points(j)%loss_standard_uncertainty_pa) .and. &
          (actual%isotherms(i)%points(j)%storage_uncertainty_available .eqv. &
          expected%isotherms(i)%points(j)%storage_uncertainty_available) .and. &
          (actual%isotherms(i)%points(j)%loss_uncertainty_available .eqv. &
          expected%isotherms(i)%points(j)%loss_uncertainty_available), &
          "V0.8.4 uncertainty input mutate edildi.")
      end do
    end do
  end subroutine assert_uncertainty_unchanged

  subroutine assert_covariance_unchanged(actual, expected)
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: actual
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: expected
    integer :: i
    integer :: j

    do i = 1, size(actual%isotherms)
      do j = 1, size(actual%isotherms(i)%points)
        call assert_true(same_real_bits(actual%isotherms(i)%points(j) &
          %covariance%storage_variance_pa2, expected%isotherms(i)%points(j) &
          %covariance%storage_variance_pa2) .and. same_real_bits( &
          actual%isotherms(i)%points(j)%covariance%loss_variance_pa2, &
          expected%isotherms(i)%points(j)%covariance &
          %loss_variance_pa2) .and. same_real_bits(actual%isotherms(i) &
          %points(j)%covariance%storage_loss_covariance_pa2, &
          expected%isotherms(i)%points(j)%covariance &
          %storage_loss_covariance_pa2) .and. &
          (actual%isotherms(i)%points(j)%covariance_available .eqv. &
          expected%isotherms(i)%points(j)%covariance_available), &
          "V0.8.5 covariance input mutate edildi.")
      end do
    end do
  end subroutine assert_covariance_unchanged

  pure function same_real_bits(a, b) result(same)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: same

    same = transfer(a, 0_int64) == transfer(b, 0_int64)
  end function same_real_bits

end program test_tts_covariance_sensitivity
