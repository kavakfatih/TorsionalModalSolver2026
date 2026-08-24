module tms_geometry
  use tms_kinds, only : dp
  use tms_constants, only : pi
  implicit none
  private

  !> Eksenel simetrik elastomer halkanın temel geometrik verilerini taşır.
  !! Tüm uzunluklar iç SI birim sözleşmesine göre metre cinsindendir.
  type, public :: rubber_geometry_t
    real(dp) :: inner_radius_m = 0.0_dp
    real(dp) :: outer_radius_m = 0.0_dp
    real(dp) :: axial_length_m = 0.0_dp
  end type rubber_geometry_t

  !> Eksenel simetrik atalet halkasının temel geometrik verilerini taşır.
  !! Tüm uzunluklar iç SI birim sözleşmesine göre metre cinsindendir.
  type, public :: inertia_ring_geometry_t
    real(dp) :: inner_radius_m = 0.0_dp
    real(dp) :: outer_radius_m = 0.0_dp
    real(dp) :: axial_length_m = 0.0_dp
  end type inertia_ring_geometry_t

  !> Eksenel simetrik göbeğin temel geometrik verilerini taşır.
  !! Delik ve dış yarıçap ile eksenel uzunluk metre cinsindendir.
  type, public :: hub_geometry_t
    real(dp) :: bore_radius_m = 0.0_dp
    real(dp) :: outer_radius_m = 0.0_dp
    real(dp) :: axial_length_m = 0.0_dp
  end type hub_geometry_t

  !> Bir TVD bileşiminin göbek, elastomer ve atalet halkası geometrilerini tutar.
  !! Bu tür yalnızca veri taşır; geometrik uygunluk veya hacim hesabı yapmaz.
  type, public :: tvd_geometry_t
    type(hub_geometry_t) :: hub
    type(rubber_geometry_t) :: rubber
    type(inertia_ring_geometry_t) :: inertia_ring
  end type tvd_geometry_t

  public :: calculate_rubber_polar_area_moment
  public :: calculate_annular_bush_torsion_geometry_factor

contains

  !> Annüler elastomer kesitin polar geometrik alan momentini hesaplar.
  !!
  !! Fiziksel açıklama: Polar alan momenti, kesit geometrisinin burulmaya
  !! karşı rijitlik katkısını temsil eder; polar kütle ataleti değildir.
  !! Matematiksel açıklama: Jp = pi/2 * (ro^4 - ri^4).
  !! Girdiler: outer_radius dış yarıçapı, inner_radius iç yarıçapı temsil eder;
  !! iki değer de metre (m) cinsindendir.
  !! Çıktı: jp polar alan momentidir ve metrenin dördüncü kuvveti (m^4)
  !! cinsindendir.
  !! Varsayımlar ve geçerlilik: Kesit eksenel simetrik ve annülerdir;
  !! 0 <= ri < ro olmalıdır. Negatif yarıçap veya ro <= ri durumunda
  !! fiziksel olarak geçersiz geometriyi önlemek için yordam error stop ile
  !! sonlanır.
  !! Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure function calculate_rubber_polar_area_moment(outer_radius, &
      inner_radius) result(jp)
    real(dp), intent(in) :: outer_radius
    real(dp), intent(in) :: inner_radius
    real(dp) :: jp

    if (outer_radius < 0.0_dp .or. inner_radius < 0.0_dp) then
      error stop "Elastomer yarıçapları negatif olamaz."
    end if

    if (outer_radius <= inner_radius) then
      error stop "Elastomer dış yarıçapı iç yarıçapından büyük olmalıdır."
    end if

    jp = 0.5_dp * pi * (outer_radius**4 - inner_radius**4)
  end function calculate_rubber_polar_area_moment

  !> Tam bağlı annüler kauçuk burcun burulma geometri faktörünü hesaplar.
  !!
  !! Fiziksel açıklama: Faktör, rijit iç göbek ile rijit dış halka
  !! arasındaki eş merkezli elastomer tabakanın bağıl dönmeye karşı geometri
  !! katkısıdır. Elastomer iki silindirik yüzeye tam bağlıdır.
  !! Matematiksel açıklama:
  !! C_theta = 4*pi*L*ri^2*ro^2/(ro^2-ri^2).
  !! Girdiler: inner_radius iç yarıçap ri, outer_radius dış yarıçap ro
  !! ve axial_length bağlı eksenel genişlik L'dir; tümü metre (m)
  !! cinsindendir.
  !! Çıktı: geometry_factor_m3, metreküp (m^3) cinsindedir ve kayma
  !! modülüyle çarpıldığında N*m/rad birimli burulma rijitliği verir.
  !! Varsayımlar ve geçerlilik: Silindirler eş merkezli ve rijit, elastomer
  !! homojen, lineer, küçük deformasyon bölgesinde ve ara yüzlerde kaymasızdır.
  !! ri > 0, ro > ri ve L > 0 olmalıdır. ri = 0, bağlı iç silindirik
  !! yüzey bulunmadığından bu burç modelinin kapsamı dışındadır.
  !! Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure function calculate_annular_bush_torsion_geometry_factor( &
      inner_radius, outer_radius, axial_length) result(geometry_factor_m3)
    real(dp), intent(in) :: inner_radius
    real(dp), intent(in) :: outer_radius
    real(dp), intent(in) :: axial_length
    real(dp) :: geometry_factor_m3

    real(dp) :: radius_square_difference_m2

    if (inner_radius < 0.0_dp .or. outer_radius < 0.0_dp) then
      error stop "Elastomer yarıçapları negatif olamaz."
    end if

    if (inner_radius == 0.0_dp) then
      error stop "Annüler kauçuk burç iç yarıçapı pozitif olmalıdır."
    end if

    if (outer_radius <= inner_radius) then
      error stop "Elastomer dış yarıçapı iç yarıçapından büyük olmalıdır."
    end if

    if (axial_length <= 0.0_dp) then
      error stop "Elastomer eksenel genişliği sıfırdan büyük olmalıdır."
    end if

    radius_square_difference_m2 = &
      (outer_radius - inner_radius) * (outer_radius + inner_radius)
    geometry_factor_m3 = 4.0_dp * pi * axial_length * inner_radius**2 * &
      outer_radius**2 / radius_square_difference_m2
  end function calculate_annular_bush_torsion_geometry_factor

end module tms_geometry
