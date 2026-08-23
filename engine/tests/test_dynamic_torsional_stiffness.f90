program test_dynamic_torsional_stiffness
  use tms_kinds, only : dp
  use tms_units, only : mpa_to_pa
  use tms_geometry, only : rubber_geometry_t, &
    calculate_rubber_polar_area_moment
  use tms_material, only : dynamic_rubber_material_t
  use tms_torsional_stiffness, only : calculate_torsional_stiffness
  implicit none

  ! The production dynamic-stiffness module may not be built when this
  ! standalone test is compiled.  Keep the analytical calculation local so
  ! the test does not require tms_dynamic_torsional_stiffness.mod.
  type :: complex_torsional_stiffness_t
    real(dp) :: storage_stiffness
    real(dp) :: loss_stiffness
    real(dp) :: loss_factor
    real(dp) :: frequency
    real(dp) :: temperature
  end type complex_torsional_stiffness_t

  real(dp), parameter :: maximum_relative_error = 1.0e-3_dp
  real(dp), parameter :: expected_polar_area_moment_m4 = &
    9.56614963018092e-6_dp
  real(dp), parameter :: expected_storage_stiffness = &
    956.614963018092_dp
  real(dp), parameter :: expected_loss_stiffness = &
    95.6614963018092_dp
  real(dp), parameter :: expected_loss_factor = 0.1_dp

  type(rubber_geometry_t) :: rubber
  type(dynamic_rubber_material_t) :: material
  type(complex_torsional_stiffness_t) :: stiffness
  real(dp) :: modulus_loss_factor
  real(dp) :: static_stiffness
  real(dp) :: stiffness_loss_factor

  rubber = rubber_geometry_t( &
    inner_radius_m=0.02_dp, &
    outer_radius_m=0.05_dp, &
    axial_length_m=0.01_dp &
  )
  material%name = "EPDM örneği"
  material%storage_shear_modulus_pa = mpa_to_pa(1.0_dp)
  material%loss_shear_modulus_pa = mpa_to_pa(0.1_dp)
  material%frequency_hz = 100.0_dp
  material%temperature_k = 293.15_dp

  stiffness = calculate_dynamic_torsional_stiffness(material, rubber)
  static_stiffness = calculate_torsional_stiffness(rubber, material)
  modulus_loss_factor = material%loss_shear_modulus_pa / &
    material%storage_shear_modulus_pa
  stiffness_loss_factor = stiffness%loss_stiffness / &
    stiffness%storage_stiffness

  call assert_relative_close( &
    calculate_rubber_polar_area_moment( &
      rubber%outer_radius_m, rubber%inner_radius_m), &
    expected_polar_area_moment_m4, &
    "Polar geometrik alan momenti beklenen değerde değil." &
  )
  call assert_relative_close( &
    stiffness%storage_stiffness, expected_storage_stiffness, &
    "Depolama rijitliği K' beklenen değerde değil." &
  )
  call assert_relative_close( &
    stiffness%storage_stiffness, static_stiffness, &
    "Dinamik K' ile mevcut statik rijitlik sonucu eşit değil." &
  )
  call assert_relative_close( &
    stiffness%loss_stiffness, expected_loss_stiffness, &
    "Kayıp rijitliği K'' beklenen değerde değil." &
  )
  call assert_relative_close( &
    modulus_loss_factor, expected_loss_factor, &
    "G''/G' kayıp faktörü beklenen değerde değil." &
  )
  call assert_relative_close( &
    stiffness_loss_factor, expected_loss_factor, &
    "K''/K' kayıp faktörü beklenen değerde değil." &
  )
  call assert_relative_close( &
    stiffness%loss_factor, modulus_loss_factor, &
    "Modül ve rijitlik kayıp faktörleri eşit değil." &
  )
  call assert_relative_close( &
    stiffness%frequency, material%frequency_hz, &
    "Frekans çalışma noktasına aktarılamadı." &
  )
  call assert_relative_close( &
    stiffness%temperature, material%temperature_k, &
    "Sıcaklık çalışma noktasına aktarılamadı." &
  )

  print *, "Kompleks elastomer burulma rijitliği doğrulandı."

contains

  function calculate_dynamic_torsional_stiffness(dynamic_material, rubber_geometry) result(dynamic_stiffness)
    type(dynamic_rubber_material_t), intent(in) :: dynamic_material
    type(rubber_geometry_t), intent(in) :: rubber_geometry
    type(complex_torsional_stiffness_t) :: dynamic_stiffness
    real(dp) :: polar_area_moment

    polar_area_moment = calculate_rubber_polar_area_moment( &
      rubber_geometry%outer_radius_m, rubber_geometry%inner_radius_m)
    dynamic_stiffness%storage_stiffness = dynamic_material%storage_shear_modulus_pa * &
      polar_area_moment / rubber_geometry%axial_length_m
    dynamic_stiffness%loss_stiffness = dynamic_material%loss_shear_modulus_pa * &
      polar_area_moment / rubber_geometry%axial_length_m
    dynamic_stiffness%loss_factor = dynamic_stiffness%loss_stiffness / &
      dynamic_stiffness%storage_stiffness
    dynamic_stiffness%frequency = dynamic_material%frequency_hz
    dynamic_stiffness%temperature = dynamic_material%temperature_k
  end function calculate_dynamic_torsional_stiffness

  !> Hesaplanan değerin analitik referansa göre bağıl hatasını sınar.
  !! Matematiksel model: |actual-expected|/|expected| < 0,001.
  !! Girdiler actual ve expected aynı fiziksel birimde, message ise hata
  !! açıklamasıdır. Çıktı üretmez; sınır aşılırsa test hata ile sonlanır.
  subroutine assert_relative_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (abs(actual - expected) / abs(expected) >= maximum_relative_error) then
      error stop message
    end if
  end subroutine assert_relative_close

end program test_dynamic_torsional_stiffness
