module tms_geometry
  use tms_kinds, only : dp
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

end module tms_geometry
