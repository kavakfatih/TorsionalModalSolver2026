program test_units
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_units, only : mm_to_m, mpa_to_pa, degree_to_radian
  implicit none

  real(dp), parameter :: tolerance = 64.0_dp * epsilon(1.0_dp)

  ! 1000 mm uzunluğun tam olarak 1 m ölçeğine dönüştüğünü doğrular.
  if (abs(mm_to_m(1000.0_dp) - 1.0_dp) > tolerance) then
    error stop "Milimetre-metre dönüşümü başarısız."
  end if

  ! 2,5 MPa malzeme modülünün 2,5 milyon Pa olduğunu doğrular.
  if (abs(mpa_to_pa(2.5_dp) - 2.5e6_dp) > tolerance * 2.5e6_dp) then
    error stop "MPa-Pa dönüşümü başarısız."
  end if

  ! 180 derecelik açının pi radyana dönüştüğünü doğrular.
  if (abs(degree_to_radian(180.0_dp) - pi) > tolerance) then
    error stop "Derece-radyan dönüşümü başarısız."
  end if

  print *, "Mühendislik birim dönüşümleri doğrulandı."
end program test_units
