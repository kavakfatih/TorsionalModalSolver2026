module tms_units
  use tms_kinds, only : dp
  use tms_constants, only : mm_to_m_factor, mpa_to_pa_factor, &
    degree_to_radian_factor
  implicit none
  private

  public :: mm_to_m
  public :: mpa_to_pa
  public :: degree_to_radian

contains

  !> Milimetre ile verilen fiziksel uzunluğu SI uzunluk birimi metreye çevirir.
  !!
  !! Matematiksel açıklama: value_m = value_mm * 10^-3 bağıntısı uygulanır.
  !! Girdi birimi: milimetre (mm).
  !! Çıktı birimi: metre (m).
  !! Varsayım ve geçerlilik: Doğrusal ölçek dönüşümüdür; tüm sonlu gerçel
  !! uzunluk değerleri için geçerlidir ve geometrik bir düzeltme uygulamaz.
  pure elemental function mm_to_m(value_mm) result(value_m)
    real(dp), intent(in) :: value_mm
    real(dp) :: value_m

    value_m = value_mm * mm_to_m_factor
  end function mm_to_m

  !> MPa ile verilen basınç, gerilme veya modülü SI birimi paskala çevirir.
  !!
  !! Matematiksel açıklama: value_pa = value_mpa * 10^6 bağıntısı uygulanır.
  !! Girdi birimi: megapaskal (MPa).
  !! Çıktı birimi: paskal (Pa = N/m^2).
  !! Varsayım ve geçerlilik: Yalnızca birim ölçeği değiştirilir; malzeme modeli,
  !! sıcaklık veya frekans bağımlılığı hesaba katılmaz.
  pure elemental function mpa_to_pa(value_mpa) result(value_pa)
    real(dp), intent(in) :: value_mpa
    real(dp) :: value_pa

    value_pa = value_mpa * mpa_to_pa_factor
  end function mpa_to_pa

  !> Derece ile verilen düzlem açısını boyutsuz SI açısı radyana çevirir.
  !!
  !! Matematiksel açıklama: value_rad = value_degree * pi / 180 uygulanır.
  !! Girdi birimi: derece (deg).
  !! Çıktı birimi: radyan (rad, boyutsuz).
  !! Varsayım ve geçerlilik: İşaret ve devir sayısı korunur; açı herhangi bir
  !! temel aralığa indirgenmez.
  pure elemental function degree_to_radian(value_degree) result(value_rad)
    real(dp), intent(in) :: value_degree
    real(dp) :: value_rad

    value_rad = value_degree * degree_to_radian_factor
  end function degree_to_radian

end module tms_units
