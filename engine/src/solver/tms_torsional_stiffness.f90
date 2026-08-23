module tms_torsional_stiffness
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_geometry, only : rubber_geometry_t
  use tms_material, only : dynamic_rubber_material_t
  implicit none
  private

  public :: calculate_torsional_stiffness

contains

  !> Annüler elastomer bölgenin lineer burulma rijitliğini hesaplar.
  !!
  !! Fiziksel açıklama: Burulma rijitliği, elastomer bölgenin göbek ile atalet
  !! halkası arasındaki bağıl dönmeye karşı ürettiği geri çağırıcı momenti
  !! temsil eder. Dinamik malzeme çalışma noktasındaki storage shear modulus
  !! G' değeri, lineer elastik kayma modülü olarak kullanılır.
  !! Matematiksel açıklama: Annüler kesit için Jp = pi/2 * (ro^4 - ri^4) ve
  !! k_theta = G' * Jp / L bağıntıları uygulanır.
  !! Girdiler: ri, ro ve efektif uzunluk L metre (m); G' paskal (Pa) cinsindedir.
  !! Çıktı: Burulma rijitliği k_theta newton-metre/radyan (N*m/rad) cinsindendir.
  !! Varsayımlar ve geçerlilik: Elastomer homojen, izotrop, küçük şekil
  !! değiştirmeli ve lineer elastiktir; 0 <= ri < ro, L > 0 ve G' > 0
  !! olmalıdır. Yordam bu önkoşulları doğrulamaz. Kayıp modülü G'' bu ilk
  !! modelde kullanılmaz. Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure function calculate_torsional_stiffness(rubber, material) &
      result(stiffness_nm_per_rad)
    type(rubber_geometry_t), intent(in) :: rubber
    type(dynamic_rubber_material_t), intent(in) :: material
    real(dp) :: stiffness_nm_per_rad

    real(dp) :: polar_area_moment_m4

    polar_area_moment_m4 = 0.5_dp * pi * &
      (rubber%outer_radius_m**4 - rubber%inner_radius_m**4)
    stiffness_nm_per_rad = material%storage_shear_modulus_pa * &
      polar_area_moment_m4 / rubber%axial_length_m
  end function calculate_torsional_stiffness

end module tms_torsional_stiffness
