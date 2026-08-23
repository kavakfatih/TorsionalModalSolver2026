program test_dynamic_modulus
  use tms_kinds, only : dp
  use tms_units, only : mpa_to_pa
  use tms_dynamic_modulus, only : dynamic_shear_modulus, calculate_loss_factor
  use tms_material_frequency, only : material_frequency_point
  use tms_material, only : dynamic_rubber_material_t
  implicit none

  type(dynamic_shear_modulus) :: modulus
  type(material_frequency_point) :: point
  type(dynamic_rubber_material_t) :: material
  real(dp), parameter :: relative_tolerance = 64.0_dp * epsilon(1.0_dp)

  modulus%storage_modulus = mpa_to_pa(1.0_dp)
  modulus%loss_modulus = mpa_to_pa(0.1_dp)
  modulus%frequency = 100.0_dp
  modulus%temperature = 293.15_dp

  call assert_close(modulus%storage_modulus, 1.0e6_dp, relative_tolerance, &
    "G' değeri Pa biriminde saklanamadı.")
  call assert_close(modulus%loss_modulus, 1.0e5_dp, relative_tolerance, &
    "G'' değeri Pa biriminde saklanamadı.")
  call assert_close(modulus%frequency, 100.0_dp, relative_tolerance, &
    "Frekans değeri Hz biriminde saklanamadı.")
  call assert_close(modulus%temperature, 293.15_dp, relative_tolerance, &
    "Sıcaklık değeri K biriminde saklanamadı.")
  call assert_close(calculate_loss_factor(modulus), 0.1_dp, &
    relative_tolerance, "Kayıp faktörü beklenen değerde değil.")

  point%storage_modulus = modulus%storage_modulus
  point%loss_modulus = modulus%loss_modulus
  point%frequency = modulus%frequency
  point%temperature = modulus%temperature

  allocate(material%frequency_points(1))
  material%frequency_points(1) = point

  if (.not. allocated(material%frequency_points)) then
    error stop "Dinamik malzeme veri noktaları ayrılamadı."
  end if
  if (size(material%frequency_points) /= 1) then
    error stop "Dinamik malzeme veri noktası sayısı hatalı."
  end if
  call assert_close(material%frequency_points(1)%storage_modulus, 1.0e6_dp, &
    relative_tolerance, "Malzeme içinde G' değeri saklanamadı.")
  call assert_close(material%frequency_points(1)%loss_modulus, 1.0e5_dp, &
    relative_tolerance, "Malzeme içinde G'' değeri saklanamadı.")

contains

  !> İki gerçek sayının ölçeklenmiş bağıl hata sınırında eşitliğini sınar.
  !! Matematiksel model: |actual-expected| <= tolerance*max(1,|expected|).
  !! Girdiler actual ve expected aynı fiziksel birimde, tolerance boyutsuzdur.
  !! Çıktı üretmez; sınır aşılırsa test hata koduyla sonlanır.
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_close

end program test_dynamic_modulus
