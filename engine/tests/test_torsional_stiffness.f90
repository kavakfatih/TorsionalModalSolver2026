program test_torsional_stiffness
  use tms_kinds, only : dp
  use tms_units, only : mm_to_m, mpa_to_pa
  use tms_geometry, only : rubber_geometry_t
  use tms_material, only : dynamic_rubber_material_t
  use tms_torsional_stiffness, only : calculate_torsional_stiffness
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-3_dp
  real(dp), parameter :: expected_stiffness_nm_per_rad = &
    3173.00858012569_dp

  type(rubber_geometry_t) :: rubber
  type(dynamic_rubber_material_t) :: material
  real(dp) :: stiffness_nm_per_rad

  ! Analitik örnek, 90-110 mm yarıçaplı lineer elastomer halkayı temsil eder.
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

  if (abs(stiffness_nm_per_rad - expected_stiffness_nm_per_rad) / &
      expected_stiffness_nm_per_rad >= maximum_relative_error) then
    error stop "Burulma rijitliği bağıl hatası yüzde 0,1 sınırını aşıyor."
  end if

  print *, "Lineer elastomer burulma rijitliği doğrulandı."
end program test_torsional_stiffness
