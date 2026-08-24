module tms_torsional_stiffness
  use tms_kinds, only : dp
  use tms_geometry, only : rubber_geometry_t, &
    calculate_annular_bush_torsion_geometry_factor
  use tms_material, only : dynamic_rubber_material_t
  implicit none
  private

  public :: calculate_torsional_stiffness

contains

  !> Tam bağlı annüler TVD elastomerinin burulma rijitliğini hesaplar.
  !!
  !! Fiziksel açıklama: Model, bağıl dönme altındaki tam bağlı eş
  !! merkezli silindirik kauçuk tabakayı (bonded concentric cylindrical
  !! rubber layer under relative rotation) temsil eder. Rijit iç göbek ve
  !! rijit dış atalet halkası arasındaki elastomer geri çağırıcı moment
  !! üretir. Malzemenin storage shear modulus G' değeri bu statik lineer
  !! modelde kayma modülü G olarak kullanılır.
  !! Matematiksel açıklama: C_theta =
  !! 4*pi*L*ri^2*ro^2/(ro^2-ri^2) ve k_theta = G'*C_theta uygulanır.
  !! Girdiler: ri ve ro metre (m), bağlı eksenel genişlik L metre (m),
  !! G' paskal (Pa) cinsindendir.
  !! Çıktı: Burulma rijitliği k_theta newton-metre/radyan (N*m/rad) cinsindendir.
  !! Varsayımlar ve geçerlilik: Elastomer homojen, izotrop, ara yüzlerde
  !! kaymasız, küçük deformasyonlu ve lineer elastiktir; ri > 0, ro > ri,
  !! L > 0 ve G' > 0 olmalıdır. Yordam bu önkoşulları doğrular. Kayıp
  !! modülü G'' bu statik modelde kullanılmaz. Ayrıntılar:
  !! docs/mathematics/torsional-physics-core.md.
  pure function calculate_torsional_stiffness(rubber, material) &
      result(stiffness_nm_per_rad)
    type(rubber_geometry_t), intent(in) :: rubber
    type(dynamic_rubber_material_t), intent(in) :: material
    real(dp) :: stiffness_nm_per_rad

    real(dp) :: geometry_factor_m3

    if (material%storage_shear_modulus_pa <= 0.0_dp) then
      error stop "Kayma modülü G sıfırdan büyük olmalıdır."
    end if

    geometry_factor_m3 = calculate_annular_bush_torsion_geometry_factor( &
      inner_radius=rubber%inner_radius_m, &
      outer_radius=rubber%outer_radius_m, &
      axial_length=rubber%axial_length_m)
    stiffness_nm_per_rad = &
      material%storage_shear_modulus_pa * geometry_factor_m3
  end function calculate_torsional_stiffness

end module tms_torsional_stiffness
