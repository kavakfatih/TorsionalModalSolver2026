program test_frequency_solver
  use tms_kinds, only : dp
  use tms_frequency_solver, only : calculate_natural_frequency
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-3_dp
  real(dp), parameter :: stiffness_nm_per_rad = 19244.2184986460_dp
  real(dp), parameter :: polar_inertia_kg_m2 = 0.220441786591211_dp
  real(dp), parameter :: expected_frequency_hz = 47.0244051144727_dp

  real(dp) :: frequency_hz

  ! Analitik örnek, lineer ve sönümsüz tek serbestlik dereceli sistemi sınar.
  frequency_hz = calculate_natural_frequency( &
    stiffness_nm_per_rad, polar_inertia_kg_m2)

  if (abs(frequency_hz - expected_frequency_hz) / expected_frequency_hz >= &
      maximum_relative_error) then
    error stop "Doğal frekans bağıl hatası yüzde 0,1 sınırını aşıyor."
  end if

  print *, "Tek serbestlik dereceli doğal frekans doğrulandı."
end program test_frequency_solver
