program test_dynamic_torsional_stiffness
  use tms_kinds, only : dp
  use tms_units, only : mpa_to_pa
  use tms_geometry, only : rubber_geometry_t, &
    calculate_rubber_polar_area_moment
  use tms_material, only : dynamic_rubber_material_t
  use tms_torsional_stiffness, only : calculate_torsional_stiffness
  use tms_dynamic_torsional_stiffness, only : &
    complex_torsional_stiffness_t, &
    calculate_dynamic_torsional_stiffness
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-10_dp
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
  character(len=32) :: validation_case

  ! CTest, geçersiz girdilerin üretim yordamları tarafından reddedildiğini
  ! aynı test yürütülebilir dosyasını farklı bir argümanla çağırarak sınar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

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
    stiffness%loss_factor, stiffness_loss_factor, &
    "Sonuç kayıp faktörü K''/K' oranına eşit değil." &
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

  !> Geçersiz fiziksel girdinin ilgili saf üretim yordamında reddedildiğini
  !! sınamak için seçilen hata durumunu çalıştırır.
  !! Model ve birimler: Yarıçaplar ile etkin uzunluk m, G' ve G'' Pa
  !! cinsindedir. Geçersiz durum error stop üretmelidir; yordam normal
  !! dönerse CTest'in WILL_FAIL özelliği regresyon testini başarısız sayar.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(rubber_geometry_t) :: invalid_rubber
    type(dynamic_rubber_material_t) :: invalid_material
    type(complex_torsional_stiffness_t) :: rejected_stiffness
    real(dp) :: rejected_polar_area_moment

    invalid_rubber = rubber_geometry_t( &
      inner_radius_m=0.02_dp, &
      outer_radius_m=0.05_dp, &
      axial_length_m=0.01_dp &
    )
    invalid_material%storage_shear_modulus_pa = mpa_to_pa(1.0_dp)
    invalid_material%loss_shear_modulus_pa = mpa_to_pa(0.1_dp)
    invalid_material%frequency_hz = 100.0_dp
    invalid_material%temperature_k = 293.15_dp
    rejected_polar_area_moment = 0.0_dp

    select case (case_name)
      case ("negative_radius")
        rejected_polar_area_moment = &
          calculate_rubber_polar_area_moment(0.05_dp, -0.01_dp)
      case ("unordered_radii")
        rejected_polar_area_moment = &
          calculate_rubber_polar_area_moment(0.02_dp, 0.05_dp)
      case ("equal_radii")
        rejected_polar_area_moment = &
          calculate_rubber_polar_area_moment(0.05_dp, 0.05_dp)
      case ("zero_length")
        invalid_rubber%axial_length_m = 0.0_dp
        rejected_stiffness = calculate_dynamic_torsional_stiffness( &
          invalid_material, invalid_rubber)
      case ("nonpositive_storage_modulus")
        invalid_material%storage_shear_modulus_pa = 0.0_dp
        rejected_stiffness = calculate_dynamic_torsional_stiffness( &
          invalid_material, invalid_rubber)
      case ("negative_loss_modulus")
        invalid_material%loss_shear_modulus_pa = -mpa_to_pa(0.1_dp)
        rejected_stiffness = calculate_dynamic_torsional_stiffness( &
          invalid_material, invalid_rubber)
      case default
        error stop "Bilinmeyen geçersiz girdi testi istendi."
    end select

    print *, "Geçersiz girdi beklenmedik biçimde kabul edildi: ", &
      case_name, rejected_polar_area_moment, &
      rejected_stiffness%storage_stiffness
  end subroutine exercise_invalid_case

  !> Hesaplanan değerin analitik referansa göre bağıl hatasını sınar.
  !! Matematiksel model: |actual-expected|/|expected| < 1e-10.
  !! Girdiler actual ve expected aynı fiziksel birimde, message ise hata
  !! açıklamasıdır; expected sıfırdan farklı olmalıdır. Çıktı
  !! üretmez; sınır aşılırsa test hata ile sonlanır.
  subroutine assert_relative_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (abs(actual - expected) / abs(expected) >= maximum_relative_error) then
      error stop message
    end if
  end subroutine assert_relative_close

end program test_dynamic_torsional_stiffness
