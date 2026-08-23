program test_constants
  use tms_kinds, only : dp
  use tms_constants, only : pi, mm_to_m_factor, mpa_to_pa_factor, &
    degree_to_radian_factor
  implicit none

  real(dp), parameter :: tolerance = 64.0_dp * epsilon(1.0_dp)

  ! Sabitlerin değerlerini çift hassasiyet sınırları içinde doğrular.
  if (abs(pi - acos(-1.0_dp)) > tolerance) then
    error stop "pi sabiti beklenen değerde değil."
  end if

  if (abs(mm_to_m_factor - 1.0e-3_dp) > tolerance) then
    error stop "Milimetre-metre dönüşüm çarpanı hatalı."
  end if

  if (abs(mpa_to_pa_factor - 1.0e6_dp) > tolerance * 1.0e6_dp) then
    error stop "MPa-Pa dönüşüm çarpanı hatalı."
  end if

  if (abs(degree_to_radian_factor - pi / 180.0_dp) > tolerance) then
    error stop "Derece-radyan dönüşüm çarpanı hatalı."
  end if

  print *, "Temel matematik ve mühendislik sabitleri doğrulandı."
end program test_constants
