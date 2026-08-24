program test_dynamic_torsional_stiffness
  use tms_kinds, only : dp
  use tms_units, only : mpa_to_pa
  use tms_geometry, only : rubber_geometry_t, &
    calculate_annular_bush_torsion_geometry_factor
  use tms_material, only : dynamic_rubber_material_t
  use tms_torsional_stiffness, only : calculate_torsional_stiffness
  use tms_dynamic_torsional_stiffness, only : &
    complex_torsional_stiffness_t, &
    calculate_dynamic_torsional_stiffness
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-10_dp
  real(dp), parameter :: maximum_absolute_error = 1.0e-12_dp
  real(dp), parameter :: expected_geometry_factor_m3 = &
    5.98398600683770e-5_dp
  real(dp), parameter :: expected_storage_stiffness = &
    59.8398600683770_dp
  real(dp), parameter :: expected_loss_stiffness = &
    5.98398600683770_dp
  real(dp), parameter :: expected_loss_factor = 0.1_dp
  real(dp), parameter :: near_gap_growth_factor = 100.0_dp

  type(rubber_geometry_t) :: rubber
  type(rubber_geometry_t) :: near_gap_rubber
  type(dynamic_rubber_material_t) :: material
  type(dynamic_rubber_material_t) :: zero_loss_material
  type(complex_torsional_stiffness_t) :: stiffness
  type(complex_torsional_stiffness_t) :: near_gap_stiffness
  type(complex_torsional_stiffness_t) :: zero_loss_stiffness
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

  block
    type(dynamic_rubber_material_t) :: doubled_storage_material
    type(complex_torsional_stiffness_t) :: doubled_storage_stiffness

    ! G' ölçekleme regresyonu üretim formülünü testte yeniden kurmaz.
    ! Geometri ve diğer malzeme alanları sabitken K' ∝ G' davranışını,
    ! iki üretim fonksiyonu sonucunu karşılaştırarak doğrular.
    doubled_storage_material = material
    doubled_storage_material%storage_shear_modulus_pa = mpa_to_pa(2.0_dp)
    doubled_storage_stiffness = calculate_dynamic_torsional_stiffness( &
      doubled_storage_material, rubber)

    call assert_relative_close( &
      doubled_storage_stiffness%storage_stiffness, &
      2.0_dp * stiffness%storage_stiffness, &
      "G' iki katına çıktığında K' iki katına çıkmadı." &
    )
  end block

  block
    type(rubber_geometry_t) :: doubled_width_rubber
    type(complex_torsional_stiffness_t) :: doubled_width_stiffness

    ! Eksenel genişlik regresyonu C_theta denklemini testte yeniden hesaplamaz.
    ! Malzeme ve yarıçaplar sabitken K' ∝ L ve K'' ∝ L davranışlarını,
    ! iki üretim fonksiyonu sonucunu karşılaştırarak doğrular.
    doubled_width_rubber = rubber
    doubled_width_rubber%axial_length_m = 0.02_dp
    doubled_width_stiffness = calculate_dynamic_torsional_stiffness( &
      material, doubled_width_rubber)

    call assert_relative_close( &
      doubled_width_stiffness%storage_stiffness, &
      2.0_dp * stiffness%storage_stiffness, &
      "Eksenel genişlik iki katına çıktığında K' iki katına çıkmadı." &
    )
    call assert_relative_close( &
      doubled_width_stiffness%loss_stiffness, &
      2.0_dp * stiffness%loss_stiffness, &
      "Eksenel genişlik iki katına çıktığında K'' iki katına çıkmadı." &
    )
  end block

  zero_loss_material = material
  zero_loss_material%loss_shear_modulus_pa = 0.0_dp
  zero_loss_stiffness = calculate_dynamic_torsional_stiffness( &
    zero_loss_material, rubber)

  ! Dış yarıçap sabitken ince radyal boşluğa yaklaşmak, ideal tam bağlı
  ! burç modelinde rijitliği kuvvetli biçimde artırmalıdır.
  near_gap_rubber = rubber
  near_gap_rubber%inner_radius_m = 0.049_dp
  near_gap_stiffness = calculate_dynamic_torsional_stiffness( &
    material, near_gap_rubber)

  call assert_relative_close( &
    calculate_annular_bush_torsion_geometry_factor( &
      inner_radius=rubber%inner_radius_m, &
      outer_radius=rubber%outer_radius_m, &
      axial_length=rubber%axial_length_m), &
    expected_geometry_factor_m3, &
    "Annüler kauçuk burç geometri faktörü beklenen değerde değil." &
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
  call assert_absolute_close( &
    zero_loss_stiffness%loss_stiffness, 0.0_dp, &
    "G'' sıfırken K'' sıfır olmadı." &
  )
  call assert_absolute_close( &
    zero_loss_stiffness%loss_factor, 0.0_dp, &
    "G'' sıfırken kayıp faktörü sıfır olmadı." &
  )

  if (near_gap_stiffness%storage_stiffness <= &
      near_gap_growth_factor * stiffness%storage_stiffness) then
    error stop "ri, ro değerine yaklaştığında rijitlik yeterince artmadı."
  end if

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
  !! Model ve birimler: Yarıçaplar ile bağlı eksenel genişlik m, G' ve G'' Pa
  !! cinsindedir. Geçersiz durum error stop üretmelidir; yordam normal
  !! dönerse CTest'in WILL_FAIL özelliği regresyon testini başarısız sayar.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(rubber_geometry_t) :: invalid_rubber
    type(dynamic_rubber_material_t) :: invalid_material
    type(complex_torsional_stiffness_t) :: rejected_stiffness
    real(dp) :: rejected_geometry_factor

    invalid_rubber = rubber_geometry_t( &
      inner_radius_m=0.02_dp, &
      outer_radius_m=0.05_dp, &
      axial_length_m=0.01_dp &
    )
    invalid_material%storage_shear_modulus_pa = mpa_to_pa(1.0_dp)
    invalid_material%loss_shear_modulus_pa = mpa_to_pa(0.1_dp)
    invalid_material%frequency_hz = 100.0_dp
    invalid_material%temperature_k = 293.15_dp
    rejected_geometry_factor = 0.0_dp

    select case (case_name)
      case ("negative_radius")
        rejected_geometry_factor = &
          calculate_annular_bush_torsion_geometry_factor( &
            inner_radius=-0.01_dp, outer_radius=0.05_dp, &
            axial_length=0.01_dp)
      case ("zero_inner_radius")
        rejected_geometry_factor = &
          calculate_annular_bush_torsion_geometry_factor( &
            inner_radius=0.0_dp, outer_radius=0.05_dp, &
            axial_length=0.01_dp)
      case ("unordered_radii")
        rejected_geometry_factor = &
          calculate_annular_bush_torsion_geometry_factor( &
            inner_radius=0.05_dp, outer_radius=0.02_dp, &
            axial_length=0.01_dp)
      case ("equal_radii")
        rejected_geometry_factor = &
          calculate_annular_bush_torsion_geometry_factor( &
            inner_radius=0.05_dp, outer_radius=0.05_dp, &
            axial_length=0.01_dp)
      case ("zero_length")
        invalid_rubber%axial_length_m = 0.0_dp
        rejected_stiffness = calculate_dynamic_torsional_stiffness( &
          invalid_material, invalid_rubber)
      case ("negative_length")
        invalid_rubber%axial_length_m = -0.01_dp
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
      case_name, rejected_geometry_factor, &
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

  !> Sıfır referanslı bir sonucun mutlak hata sınırını denetler.
  !! Girdiler aynı fiziksel birimdedir; |actual-expected| <= 1e-12 koşulu
  !! uygulanır ve sınır aşılırsa test hata ile sonlanır.
  subroutine assert_absolute_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > maximum_absolute_error) then
      error stop message
    end if
  end subroutine assert_absolute_close

end program test_dynamic_torsional_stiffness
