module tms_frequency_solver
  use tms_kinds, only : dp
  use tms_constants, only : pi
  implicit none
  private

  public :: calculate_natural_frequency

contains

  !> Tek serbestlik dereceli burulma sisteminin doğal frekansını hesaplar.
  !!
  !! Fiziksel açıklama: Doğal frekans, atalet halkasının elastomer rijitliği
  !! altında dış zorlama olmadan salınım yapma hızını temsil eder.
  !! Matematiksel açıklama: Açısal doğal frekans omega_n = sqrt(k_theta / J),
  !! çevrimsel doğal frekans ise f_n = omega_n / (2*pi) bağıntısıyla bulunur.
  !! Girdiler: Burulma rijitliği k_theta N*m/rad, polar kütle atalet momenti J
  !! kg*m^2 cinsindendir.
  !! Çıktı: Doğal frekans f_n hertz (Hz) cinsindendir.
  !! Varsayımlar ve geçerlilik: Sistem sönümsüz, lineer ve tek serbestlik
  !! derecelidir; k_theta > 0 ve J > 0 olmalıdır. Yordam bu önkoşulları
  !! doğrulamaz. Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure elemental function calculate_natural_frequency( &
      stiffness_nm_per_rad, polar_inertia_kg_m2) result(frequency_hz)
    real(dp), intent(in) :: stiffness_nm_per_rad
    real(dp), intent(in) :: polar_inertia_kg_m2
    real(dp) :: frequency_hz

    frequency_hz = sqrt(stiffness_nm_per_rad / polar_inertia_kg_m2) / &
      (2.0_dp * pi)
  end function calculate_natural_frequency

end module tms_frequency_solver
