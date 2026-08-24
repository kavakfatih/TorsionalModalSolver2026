program test_geometry
  use tms_kinds, only : dp
  use tms_geometry, only : rubber_geometry_t, inertia_ring_geometry_t, &
    hub_geometry_t, tvd_geometry_t, calculate_rubber_polar_area_moment, &
    calculate_annular_bush_torsion_geometry_factor
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-12_dp
  real(dp), parameter :: expected_polar_area_moment_m4 = &
    1.2692034320502767e-4_dp
  real(dp), parameter :: expected_bush_geometry_factor_m3 = &
    1.53953747989168e-2_dp

  type(tvd_geometry_t) :: geometry
  real(dp) :: polar_area_moment_m4
  real(dp) :: bush_geometry_factor_m3

  ! Temel TVD bileşenlerinin SI tabanlı geometrik verilerini oluşturur.
  geometry%rubber = rubber_geometry_t( &
    inner_radius_m=0.090_dp, &
    outer_radius_m=0.110_dp, &
    axial_length_m=0.050_dp &
  )
  geometry%inertia_ring = inertia_ring_geometry_t( &
    inner_radius_m=0.110_dp, &
    outer_radius_m=0.150_dp, &
    axial_length_m=0.050_dp &
  )
  geometry%hub = hub_geometry_t( &
    bore_radius_m=0.025_dp, &
    outer_radius_m=0.090_dp, &
    axial_length_m=0.050_dp &
  )

  if (geometry%rubber%inner_radius_m /= 0.090_dp) then
    error stop "Elastomer geometrisi veri ataması başarısız."
  end if

  if (geometry%inertia_ring%outer_radius_m /= 0.150_dp) then
    error stop "Atalet halkası geometrisi veri ataması başarısız."
  end if

  if (geometry%hub%bore_radius_m /= 0.025_dp) then
    error stop "Göbek geometrisi veri ataması başarısız."
  end if

  ! Annüler kesit için Jp = pi/2*(ro^4-ri^4) analitik sonucunu doğrular.
  polar_area_moment_m4 = calculate_rubber_polar_area_moment( &
    geometry%rubber%outer_radius_m, geometry%rubber%inner_radius_m)

  if (abs(polar_area_moment_m4 - expected_polar_area_moment_m4) / &
      expected_polar_area_moment_m4 >= maximum_relative_error) then
    error stop "Elastomer polar alan momenti analitik sonuçla uyuşmuyor."
  end if

  ! Tam bağlı eş merkezli kauçuk burç için C_theta analitik sonucunu
  ! polar alan momentinden bağımsız olarak doğrular.
  bush_geometry_factor_m3 = &
    calculate_annular_bush_torsion_geometry_factor( &
      inner_radius=geometry%rubber%inner_radius_m, &
      outer_radius=geometry%rubber%outer_radius_m, &
      axial_length=geometry%rubber%axial_length_m)

  if (abs(bush_geometry_factor_m3 - expected_bush_geometry_factor_m3) / &
      expected_bush_geometry_factor_m3 >= maximum_relative_error) then
    error stop "Annüler kauçuk burç geometri faktörü analitik sonuçla uyuşmuyor."
  end if

  print *, "TVD geometri yordamları analitik sonuçlarla doğrulandı."
end program test_geometry
