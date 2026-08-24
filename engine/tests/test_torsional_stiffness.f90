program test_torsional_stiffness
  use tms_kinds, only : dp
  use tms_units, only : mm_to_m, mpa_to_pa
  use tms_geometry, only : rubber_geometry_t
  use tms_material, only : dynamic_rubber_material_t
  use tms_torsional_stiffness, only : calculate_torsional_stiffness
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-10_dp
  real(dp), parameter :: expected_stiffness_nm_per_rad = &
    19244.2184986460_dp

  type(rubber_geometry_t) :: rubber
  type(rubber_geometry_t) :: doubled_length_rubber
  type(dynamic_rubber_material_t) :: material
  type(dynamic_rubber_material_t) :: doubled_modulus_material
  real(dp) :: stiffness_nm_per_rad
  real(dp) :: doubled_length_stiffness
  real(dp) :: doubled_modulus_stiffness
  character(len=32) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  ! Analitik örnek, rijit göbek ve dış halkaya tam bağlı 90-110 mm
  ! yarıçaplı lineer elastomer burcu temsil eder.
  rubber = rubber_geometry_t( &
    inner_radius_m=mm_to_m(90.0_dp), &
    outer_radius_m=mm_to_m(110.0_dp), &
    axial_length_m=mm_to_m(50.0_dp) &
  )
  material%name = "Benchmark lineer elastomer"
  material%storage_shear_modulus_pa = mpa_to_pa(1.25_dp)
  material%temperature_k = 293.15_dp
  material%frequency_hz = 25.0_dp

  stiffness_nm_per_rad = calculate_torsional_stiffness(rubber, material)
  doubled_modulus_material = material
  doubled_modulus_material%storage_shear_modulus_pa = &
    2.0_dp * material%storage_shear_modulus_pa
  doubled_modulus_stiffness = calculate_torsional_stiffness( &
    rubber, doubled_modulus_material)

  doubled_length_rubber = rubber
  doubled_length_rubber%axial_length_m = 2.0_dp * rubber%axial_length_m
  doubled_length_stiffness = calculate_torsional_stiffness( &
    doubled_length_rubber, material)

  call assert_relative_close( &
    stiffness_nm_per_rad, expected_stiffness_nm_per_rad, &
    "Bağlı annüler kauçuk burç rijitliği analitik sonuçla uyuşmuyor.")
  call assert_relative_close( &
    doubled_modulus_stiffness, 2.0_dp * stiffness_nm_per_rad, &
    "Kayma modülü iki katına çıktığında rijitlik iki katına çıkmadı.")
  call assert_relative_close( &
    doubled_length_stiffness, 2.0_dp * stiffness_nm_per_rad, &
    "Eksenel genişlik iki katına çıktığında rijitlik iki katına çıkmadı.")

  print *, "Bağlı annüler TVD elastomer rijitliği doğrulandı."

contains

  !> Statik solver'ın fiziksel olmayan kayma modülünü reddettiğini sınar.
  !! Girdi modülü Pa cinsindedir ve G > 0 olmalıdır. Geçersiz değer
  !! error stop üretmelidir; normal dönüş CTest tarafından hata sayılır.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(rubber_geometry_t) :: invalid_rubber
    type(dynamic_rubber_material_t) :: invalid_material
    real(dp) :: rejected_stiffness

    invalid_rubber = rubber_geometry_t( &
      inner_radius_m=0.02_dp, &
      outer_radius_m=0.05_dp, &
      axial_length_m=0.01_dp)
    invalid_material%storage_shear_modulus_pa = mpa_to_pa(1.0_dp)
    rejected_stiffness = 0.0_dp

    select case (case_name)
      case ("nonpositive_modulus")
        invalid_material%storage_shear_modulus_pa = 0.0_dp
        rejected_stiffness = calculate_torsional_stiffness( &
          invalid_rubber, invalid_material)
      case default
        error stop "Bilinmeyen statik geçersiz girdi testi istendi."
    end select

    print *, "Geçersiz statik girdi beklenmedik biçimde kabul edildi: ", &
      case_name, rejected_stiffness
  end subroutine exercise_invalid_case

  !> Hesaplanan rijitliğin bağımsız analitik referansa bağıl hatasını sınar.
  !! Matematiksel model: |actual-expected|/|expected| < 1e-10. Girdiler aynı
  !! fiziksel birimde ve expected sıfırdan farklı olmalıdır.
  subroutine assert_relative_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (abs(actual - expected) / abs(expected) >= maximum_relative_error) then
      error stop message
    end if
  end subroutine assert_relative_close

end program test_torsional_stiffness
