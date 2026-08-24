module tms_inertia
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_geometry, only : hub_geometry_t, inertia_ring_geometry_t
  implicit none
  private

  public :: calculate_annular_ring_properties
  public :: calculate_annular_hub_properties

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
  !! doğrular. Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure subroutine calculate_annular_ring_properties( &
      ring, density_kg_m3, mass_kg, polar_inertia_kg_m2)
    type(inertia_ring_geometry_t), intent(in) :: ring
    real(dp), intent(in) :: density_kg_m3
    real(dp), intent(out) :: mass_kg
    real(dp), intent(out) :: polar_inertia_kg_m2

    real(dp) :: unused_volume_m3

    call calculate_homogeneous_annular_properties( &
      inner_radius_m=ring%inner_radius_m, &
      outer_radius_m=ring%outer_radius_m, &
      axial_length_m=ring%axial_length_m, &
      density_kg_m3=density_kg_m3, &
      volume_m3=unused_volume_m3, &
      mass_kg=mass_kg, &
      polar_inertia_kg_m2=polar_inertia_kg_m2)
  end subroutine calculate_annular_ring_properties

  !> Homojen annüler göbeğin hacmini, kütlesini ve polar kütle ataletini
  !! hesaplar.
  !!
  !! Fiziksel açıklama: Polar kütle ataleti, göbeğin dönme hızındaki
  !! değişime karşı gösterdiği direnci; hacim ve kütle ise yoğunlukla birlikte
  !! rijit gövdenin dağıtılmış kütlesini temsil eder.
  !! Matematiksel açıklama: V = pi*(ro^2-ri^2)*L, m = rho*V ve
  !! J_h = 1/2*m*(ro^2+ri^2) bağıntıları uygulanır.
  !! Girdiler: Delik yarıçapı ri, dış yarıçap ro ve eksenel genişlik L metre
  !! (m), yoğunluk rho kilogram/metreküp (kg/m^3) cinsindendir.
  !! Çıktılar: Hacim V metreküp (m^3), kütle m kilogram (kg), polar kütle
  !! ataleti J_h kilogram-metrekare (kg*m^2) cinsindendir.
  !! Varsayımlar ve geçerlilik: Göbek homojen ve eksenel simetriktir;
  !! 0 <= ri < ro, L > 0 ve rho > 0 olmalıdır. Yordam bu önkoşulları
  !! doğrular. Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure subroutine calculate_annular_hub_properties( &
      hub, density_kg_m3, volume_m3, mass_kg, polar_inertia_kg_m2)
    type(hub_geometry_t), intent(in) :: hub
    real(dp), intent(in) :: density_kg_m3
    real(dp), intent(out) :: volume_m3
    real(dp), intent(out) :: mass_kg
    real(dp), intent(out) :: polar_inertia_kg_m2

    call calculate_homogeneous_annular_properties( &
      inner_radius_m=hub%bore_radius_m, &
      outer_radius_m=hub%outer_radius_m, &
      axial_length_m=hub%axial_length_m, &
      density_kg_m3=density_kg_m3, &
      volume_m3=volume_m3, &
      mass_kg=mass_kg, &
      polar_inertia_kg_m2=polar_inertia_kg_m2)
  end subroutine calculate_annular_hub_properties

  !> Homojen annüler rijit bir gövdenin ortak kütle özelliklerini hesaplar.
  !!
  !! Fiziksel açıklama: Aynı eksenel simetrik kütle dağılımı hem göbek hem de
  !! atalet halkası için geçerlidir; bu yardımcı yordam ortak fiziği tek yerde
  !! uygular.
  !! Matematiksel açıklama: V = pi*(ro^2-ri^2)*L, m = rho*V ve
  !! J = 1/2*m*(ro^2+ri^2).
  !! Girdiler: ri, ro ve L metre (m), rho kg/m^3 cinsindedir.
  !! Çıktılar: V m^3, m kg ve J kg*m^2 cinsindendir.
  !! Varsayımlar ve geçerlilik: Gövde homojen, eş merkezli ve rijittir;
  !! 0 <= ri < ro, L > 0 ve rho > 0 olmalıdır. Geçersiz girdiler reddedilir.
  pure subroutine calculate_homogeneous_annular_properties( &
      inner_radius_m, outer_radius_m, axial_length_m, density_kg_m3, &
      volume_m3, mass_kg, polar_inertia_kg_m2)
    real(dp), intent(in) :: inner_radius_m
    real(dp), intent(in) :: outer_radius_m
    real(dp), intent(in) :: axial_length_m
    real(dp), intent(in) :: density_kg_m3
    real(dp), intent(out) :: volume_m3
    real(dp), intent(out) :: mass_kg
    real(dp), intent(out) :: polar_inertia_kg_m2

    if (inner_radius_m < 0.0_dp) then
      error stop "Annüler rijit gövdenin iç yarıçapı negatif olamaz."
    end if

    if (outer_radius_m <= inner_radius_m) then
      error stop "Annüler rijit gövdenin dış yarıçapı iç yarıçapından büyük olmalıdır."
    end if

    if (axial_length_m <= 0.0_dp) then
      error stop "Annüler rijit gövdenin eksenel genişliği pozitif olmalıdır."
    end if

    if (density_kg_m3 <= 0.0_dp) then
      error stop "Annüler rijit gövdenin yoğunluğu pozitif olmalıdır."
    end if

    volume_m3 = pi * (outer_radius_m**2 - inner_radius_m**2) * &
      axial_length_m
    mass_kg = density_kg_m3 * volume_m3
    polar_inertia_kg_m2 = 0.5_dp * mass_kg * &
      (outer_radius_m**2 + inner_radius_m**2)
  end subroutine calculate_homogeneous_annular_properties

end module tms_inertia
