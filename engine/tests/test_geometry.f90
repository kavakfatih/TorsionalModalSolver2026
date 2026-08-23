program test_geometry
  use tms_kinds, only : dp
  use tms_geometry, only : rubber_geometry_t, inertia_ring_geometry_t, &
    hub_geometry_t, tvd_geometry_t
  implicit none

  type(tvd_geometry_t) :: geometry

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

  print *, "TVD geometri veri türleri doğrulandı."
end program test_geometry
