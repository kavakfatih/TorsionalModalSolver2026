module tms_inertia
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_geometry, only : inertia_ring_geometry_t
  implicit none
  private

  public :: calculate_annular_ring_properties

contains

  !> Homojen ve eksenel simetrik katı halkanın kütlesini ve polar kütle
  !! atalet momentini hesaplar.
  !!
  !! Fiziksel açıklama: Polar kütle atalet momenti, atalet halkasının dönme
  !! hızındaki değişime karşı gösterdiği direnci temsil eder.
  !! Matematiksel açıklama: Önce V = pi * (ro^2 - ri^2) * b ve m = rho * V,
  !! ardından J = 1/2 * m * (ro^2 + ri^2) bağıntıları uygulanır.
  !! Girdiler: İç yarıçap ri, dış yarıçap ro ve eksenel genişlik b metre (m),
  !! yoğunluk rho kilogram/metreküp (kg/m^3) cinsindendir.
  !! Çıktılar: Kütle m kilogram (kg), polar kütle atalet momenti J
  !! kilogram-metrekare (kg*m^2) cinsindendir.
  !! Varsayımlar ve geçerlilik: Halka homojen ve eksenel simetriktir;
  !! 0 <= ri < ro, b > 0 ve rho > 0 olmalıdır. Yordam bu önkoşulları
  !! doğrulamaz. Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure subroutine calculate_annular_ring_properties( &
      ring, density_kg_m3, mass_kg, polar_inertia_kg_m2)
    type(inertia_ring_geometry_t), intent(in) :: ring
    real(dp), intent(in) :: density_kg_m3
    real(dp), intent(out) :: mass_kg
    real(dp), intent(out) :: polar_inertia_kg_m2

    real(dp) :: volume_m3

    volume_m3 = pi * (ring%outer_radius_m**2 - ring%inner_radius_m**2) * &
      ring%axial_length_m
    mass_kg = density_kg_m3 * volume_m3
    polar_inertia_kg_m2 = 0.5_dp * mass_kg * &
      (ring%outer_radius_m**2 + ring%inner_radius_m**2)
  end subroutine calculate_annular_ring_properties

end module tms_inertia
