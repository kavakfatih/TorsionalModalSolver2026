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

contains

  !> Annüler elastomer kesitin polar geometrik alan momentini hesaplar.
  !!
  !! Fiziksel açıklama: Polar alan momenti, kesit geometrisinin burulmaya
  !! karşı rijitlik katkısını temsil eder; polar kütle ataleti değildir.
  !! Matematiksel açıklama: Jp = pi/2 * (ro^4 - ri^4).
  !! Girdi: rubber iç ve dış yarıçapları metre (m) cinsindedir. Eksenel
  !! uzunluk bu geometrik kesit hesabında kullanılmaz.
  !! Çıktı: Polar alan momenti metrenin dördüncü kuvveti (m^4) cinsindendir.
  !! Varsayımlar ve geçerlilik: Kesit eksenel simetrik ve annülerdir;
  !! 0 <= ri < ro olmalıdır. Yordam bu önkoşulları doğrulamaz.
  !! Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure elemental function calculate_rubber_polar_area_moment(rubber) &
      result(polar_area_moment_m4)
    type(rubber_geometry_t), intent(in) :: rubber
    real(dp) :: polar_area_moment_m4

    polar_area_moment_m4 = 0.5_dp * pi * &
      (rubber%outer_radius_m**4 - rubber%inner_radius_m**4)
  end function calculate_rubber_polar_area_moment

end module tms_geometry
