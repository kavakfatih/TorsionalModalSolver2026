program test_hub_inertia
  use tms_kinds, only : dp
  use tms_geometry, only : hub_geometry_t
  use tms_inertia, only : calculate_annular_hub_properties
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-10_dp
  real(dp), parameter :: expected_volume_m3 = &
    0.0015079644737231012_dp
  real(dp), parameter :: expected_mass_kg = &
    11.762122895040189_dp
  real(dp), parameter :: expected_inertia_kg_m2 = &
    0.061163039054208994_dp

  type(hub_geometry_t) :: hub
  real(dp) :: volume_m3
  real(dp) :: mass_kg
  real(dp) :: polar_inertia_kg_m2

  ! Analitik referans; 20 mm delik yarıçaplı, 100 mm dış yarıçaplı,
  ! 50 mm genişliğinde ve 7800 kg/m^3 yoğunluklu homojen annüler göbektir.
  ! Üretim yordamının V, m ve J_h çıktıları bağımsız sabitlerle karşılaştırılır.
  hub = hub_geometry_t( &
    bore_radius_m=0.02_dp, &
    outer_radius_m=0.10_dp, &
    axial_length_m=0.05_dp)

  call calculate_annular_hub_properties( &
    hub, 7800.0_dp, volume_m3, mass_kg, polar_inertia_kg_m2)

  call assert_relative_close( &
    volume_m3, expected_volume_m3, &
    "Göbek hacmi bağımsız analitik referansla uyuşmuyor.")
  call assert_relative_close( &
    mass_kg, expected_mass_kg, &
    "Göbek kütlesi bağımsız analitik referansla uyuşmuyor.")
  call assert_relative_close( &
    polar_inertia_kg_m2, expected_inertia_kg_m2, &
    "Göbek polar kütle ataleti bağımsız analitik referansla uyuşmuyor.")

  print *, "Annüler göbek hacmi, kütlesi ve polar ataleti doğrulandı."

contains

  !> Hesaplanan pozitif fiziksel büyüklüğü analitik referansla karşılaştırır.
  !! Matematiksel açıklama: |actual-expected|/|expected| < 1e-10 koşulu
  !! uygulanır. Girdiler aynı SI biriminde olmalıdır; çıktı üretmez ve hata
  !! sınırı aşılırsa test error stop ile sonlanır.
  subroutine assert_relative_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (abs(actual - expected) / abs(expected) >= &
        maximum_relative_error) then
      error stop message
    end if
  end subroutine assert_relative_close

end program test_hub_inertia
