program test_tts_diagnostics
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_vgp_point_t, &
    tts_cole_cole_point_t, MEASUREMENT_VALID
  use tms_tts_diagnostics, only : build_tts_vgp_cloud, &
    build_tts_cole_cole_cloud
  use tms_tts_test_support, only : make_exact_trs_family, &
    assert_true, assert_close
  implicit none

  type(tts_material_family_t) :: family
  type(tts_vgp_point_t), allocatable :: vgp(:)
  type(tts_cole_cole_point_t), allocatable :: cole(:)
  integer :: zero_loss_index

  family = make_exact_trs_family([293.15_dp, 313.15_dp], [0.0_dp, -1.0_dp])
  family%isotherms(1)%points(1)%storage_modulus_pa = 3.0_dp
  family%isotherms(1)%points(1)%loss_modulus_pa = 4.0_dp
  family%isotherms(1)%points(2)%loss_modulus_pa = 0.0_dp
  family%isotherms(1)%points(2)%loss_quality = MEASUREMENT_VALID

  vgp = build_tts_vgp_cloud(family)
  cole = build_tts_cole_cole_cloud(family)
  call assert_true(size(vgp) == 14 .and. size(cole) == 14, &
    "VGP/Cole-Cole VALID complex point sayısı hatalı.")
  call assert_close(vgp(1)%complex_modulus_magnitude_pa, 5.0_dp, &
    1.0e-14_dp, "VGP |G*|=sqrt(G'^2+G''^2) sonucu hatalı.")
  call assert_close(vgp(1)%phase_angle_rad, atan2(4.0_dp, 3.0_dp), &
    1.0e-14_dp, "VGP delta=atan2(G'',G') sonucu hatalı.")
  call assert_close(cole(1)%storage_modulus_pa, 3.0_dp, 0.0_dp, &
    "Cole-Cole G' ordinate/abscissa eşlemesi hatalı.")
  call assert_close(cole(1)%loss_modulus_pa, 4.0_dp, 0.0_dp, &
    "Cole-Cole G'' ordinate/abscissa eşlemesi hatalı.")
  call assert_true(vgp(1)%source_isotherm_index == 1 .and. &
    vgp(1)%source_point_index == 1 .and. &
    abs(vgp(1)%source_temperature_k - 293.15_dp) < 1.0e-12_dp, &
    "VGP temperature/source provenance kayboldu.")

  zero_loss_index = 2
  call assert_close(vgp(zero_loss_index)%phase_angle_rad, 0.0_dp, &
    0.0_dp, "VALID zero-loss VGP cloud'dan çıkarıldı veya phase hatalı.")
  call assert_close(vgp(zero_loss_index)%complex_modulus_magnitude_pa, &
    family%isotherms(1)%points(2)%storage_modulus_pa, 1.0e-14_dp, &
    "Zero-loss VGP magnitude storage modulus'a eşit değil.")

  print *, "V0.8.1 VGP ve Cole-Cole point-cloud diagnostics doğrulandı."
end program test_tts_diagnostics
