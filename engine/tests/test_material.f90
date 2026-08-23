program test_material
  use tms_kinds, only : dp
  use tms_material, only : dynamic_rubber_material_t
  implicit none

  type(dynamic_rubber_material_t) :: material

  ! Dinamik elastomerin örnek çalışma noktası verilerini SI birimleriyle atar.
  material%name = "Örnek dinamik elastomer"
  material%density_kg_m3 = 1100.0_dp
  material%storage_shear_modulus_pa = 1.25e6_dp
  material%loss_shear_modulus_pa = 0.30e6_dp
  material%temperature_k = 293.15_dp
  material%frequency_hz = 25.0_dp

  if (trim(material%name) /= "Örnek dinamik elastomer") then
    error stop "Malzeme adı veri ataması başarısız."
  end if

  if (material%density_kg_m3 /= 1100.0_dp) then
    error stop "Malzeme yoğunluğu veri ataması başarısız."
  end if

  if (material%storage_shear_modulus_pa /= 1.25e6_dp) then
    error stop "Storage shear modulus veri ataması başarısız."
  end if

  if (material%loss_shear_modulus_pa /= 0.30e6_dp) then
    error stop "Loss shear modulus veri ataması başarısız."
  end if

  if (material%temperature_k /= 293.15_dp .or. &
      material%frequency_hz /= 25.0_dp) then
    error stop "Malzeme çalışma noktası veri ataması başarısız."
  end if

  print *, "Dinamik elastomer malzeme veri türü doğrulandı."
end program test_material
