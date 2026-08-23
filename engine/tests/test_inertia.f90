program test_inertia
  use tms_kinds, only : dp
  use tms_units, only : mm_to_m
  use tms_geometry, only : inertia_ring_geometry_t
  use tms_inertia, only : calculate_annular_ring_properties
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-3_dp
  real(dp), parameter :: expected_mass_kg = 12.7422998029602_dp
  real(dp), parameter :: expected_inertia_kg_m2 = 0.220441786591211_dp

  type(inertia_ring_geometry_t) :: ring
  real(dp) :: mass_kg
  real(dp) :: polar_inertia_kg_m2

  ! Analitik örnek, 110-150 mm yarıçaplı homojen çelik halkayı temsil eder.
  ring = inertia_ring_geometry_t( &
    inner_radius_m=mm_to_m(110.0_dp), &
    outer_radius_m=mm_to_m(150.0_dp), &
    axial_length_m=mm_to_m(50.0_dp) &
  )

  call calculate_annular_ring_properties( &
    ring, 7800.0_dp, mass_kg, polar_inertia_kg_m2)

  if (abs(mass_kg - expected_mass_kg) / expected_mass_kg >= &
      maximum_relative_error) then
    error stop "Halka kütlesinin bağıl hatası yüzde 0,1 sınırını aşıyor."
  end if

  if (abs(polar_inertia_kg_m2 - expected_inertia_kg_m2) / &
      expected_inertia_kg_m2 >= maximum_relative_error) then
    error stop "Polar ataletin bağıl hatası yüzde 0,1 sınırını aşıyor."
  end if

  print *, "Katı halka kütlesi ve polar ataleti doğrulandı."
end program test_inertia
