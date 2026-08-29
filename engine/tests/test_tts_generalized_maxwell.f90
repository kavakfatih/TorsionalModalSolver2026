program test_tts_generalized_maxwell
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_identification_result_t, TTS_IDENTIFICATION_SUCCESS, &
    PAIR_SHIFT_SUCCESS
  use tms_tts_identification, only : identify_tts_master_curve
  use tms_tts_test_support, only : &
    make_generalized_maxwell_trs_family, &
    generalized_maxwell_storage_modulus, &
    generalized_maxwell_loss_modulus, find_empirical_shift, &
    assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures_k(3) = &
    [273.15_dp, 293.15_dp, 313.15_dp]
  real(dp), parameter :: truth_shifts(3) = [1.0_dp, 0.0_dp, -1.0_dp]
  type(tts_material_family_t) :: family
  type(tts_identification_result_t) :: result
  real(dp) :: high_region_slope
  real(dp) :: low_region_slope
  integer :: i

  family = make_generalized_maxwell_trs_family( &
    temperatures_k, truth_shifts)

  ! Synthetic input equations production identification'dan bağımsızdır.
  ! ISO-1'de f=0.1 Hz ve s=+1, f_r=1 Hz analytical master noktasına gider.
  call assert_close(family%isotherms(1)%points(21)%storage_modulus_pa, &
    generalized_maxwell_storage_modulus(1.0_dp), 2.0e-13_dp, &
    "Generalized-Maxwell storage test truth yanlış üretildi.")
  call assert_close(family%isotherms(1)%points(21)%loss_modulus_pa, &
    generalized_maxwell_loss_modulus(1.0_dp), 2.0e-13_dp, &
    "Generalized-Maxwell loss test truth yanlış üretildi.")
  do i = 1, size(family%isotherms)
    call assert_true(all(family%isotherms(i)%points%storage_modulus_pa > &
      0.0_dp), "Generalized-Maxwell G' pozitif değil.")
    call assert_true(all(family%isotherms(i)%points%loss_modulus_pa >= &
      0.0_dp), "Generalized-Maxwell G'' passive değil.")
  end do

  ! Farklı log-frequency bölgelerindeki analitik secant slope'ların ayrışması
  ! benchmark'ın straight power-law değil gerçekten curved olduğunu gösterir.
  low_region_slope = log10( &
    generalized_maxwell_storage_modulus(1.0e-2_dp) / &
    generalized_maxwell_storage_modulus(1.0e-3_dp))
  high_region_slope = log10( &
    generalized_maxwell_storage_modulus(1.0_dp) / &
    generalized_maxwell_storage_modulus(1.0e-1_dp))
  call assert_true(abs(low_region_slope - high_region_slope) > 0.1_dp, &
    "Generalized-Maxwell benchmark curvature göstermiyor.")

  result = identify_tts_master_curve(family, "ISO-2")
  call assert_true(result%status == TTS_IDENTIFICATION_SUCCESS .and. &
    result%runtime_export_ready, &
    "Curved exact-TRS identification/runtime table üretilemedi.")
  call assert_close(find_empirical_shift(result%empirical_shifts, 1), &
    truth_shifts(1), 2.0e-6_dp, &
    "Curved colder absolute shift doğru yönde recovery edilmedi.")
  call assert_close(find_empirical_shift(result%empirical_shifts, 2), &
    0.0_dp, 0.0_dp, "Curved reference shift sıfır değil.")
  call assert_close(find_empirical_shift(result%empirical_shifts, 3), &
    truth_shifts(3), 2.0e-6_dp, &
    "Curved hotter absolute shift doğru yönde recovery edilmedi.")
  do i = 1, size(result%pair_shift_results)
    call assert_true( &
      result%pair_shift_results(i)%status == PAIR_SHIFT_SUCCESS .and. &
      result%pair_shift_results(i)%joint_shift_available, &
      "Curved benchmark joint G'/G'' shift üretmedi.")
    call assert_close(result%pair_shift_results(i) &
      %storage_loss_shift_discrepancy, 0.0_dp, 2.0e-6_dp, &
      "Curved benchmark G'/G'' shift'leri uyuşmuyor.")
    call assert_true(ieee_is_finite(result%pair_shift_results(i) &
      %objective_curvature) .and. &
      result%pair_shift_results(i)%objective_curvature > 0.0_dp, &
      "Curved benchmark objective curvature pozitif/sonlu değil.")
    call assert_true(result%pair_shift_results(i)%objective_minimum < &
      1.0e-12_dp, "Curved exact-TRS pair residual yeterince küçük değil.")
  end do
  call assert_true(result%master_cloud_available .and. &
    size(result%master_cloud) == 183, &
    "Curved master cloud bütün source points'i korumadı.")
  call assert_true(result%master_cloud(1)%source_isotherm_index == 1 .and. &
    result%master_cloud(1)%source_point_index == 1 .and. &
    trim(result%master_cloud(1)%specimen_identifier) == "SPECIMEN-1", &
    "Curved master cloud provenance kayboldu.")
  call assert_true(size(result%runtime_master_table) > 61, &
    "Curved runtime stitching reference domain'ini genişletmedi.")

  print *, "V0.8.1 curved generalized-Maxwell exact-TRS doğrulandı."
end program test_tts_generalized_maxwell
