module tms_dynamic_modulus
  use tms_kinds, only : dp
  implicit none
  private

  !> Kompleks dinamik kayma modülünün reel ve sanal bileşenlerini,
  !! ölçüm frekansını ve sıcaklığını birlikte saklar.
  !! Matematiksel gösterim G* = G' + iG'' biçimindedir.
  type, public :: dynamic_shear_modulus
    !> Depolama modülü G' [Pa]. Bir çevrimde elastik olarak depolanan
    !! enerjiyle ilişkili, kompleks kayma modülünün reel bileşenidir.
    real(dp) :: storage_modulus = 0.0_dp

    !> Kayıp modülü G'' [Pa]. Bir çevrimde ısıya dönüşerek kaybedilen
    !! enerjiyle ilişkili, kompleks kayma modülünün sanal bileşenidir.
    real(dp) :: loss_modulus = 0.0_dp

    !> Modül değerlerinin geçerli olduğu uyarım frekansı [Hz].
    !! Matematiksel olarak saniye başına çevrim sayısını ifade eder.
    real(dp) :: frequency = 0.0_dp

    !> Modül değerlerinin geçerli olduğu mutlak sıcaklık [K].
    !! Sıcaklığa bağlı elastomer davranışının durum değişkenidir.
    real(dp) :: temperature = 0.0_dp
  end type dynamic_shear_modulus

  public :: calculate_loss_factor

contains

  !> Kompleks kayma modülünün boyutsuz kayıp faktörünü hesaplar.
  !!
  !! Fiziksel model: Kayıp faktörü, elastomerin bir çevrimde kaybettiği
  !! enerji ile depoladığı elastik enerji arasındaki sönüm göstergesidir.
  !! Matematiksel model: tan(delta) = G'' / G'.
  !! Girdi: modulus; G' ve G'' bileşenleri Pa, frekans Hz ve sıcaklık K.
  !! Çıktı: loss_factor; boyutsuzdur.
  !! Varsayım ve geçerlilik sınırı: G' sıfırdan büyük olmalıdır. Fonksiyon
  !! interpolasyon, eğri uydurma veya sıcaklık-frekans dönüşümü yapmaz.
  !! Ayrıntı: docs/mathematics/dynamic_elastomer_model.md
  pure elemental function calculate_loss_factor(modulus) result(loss_factor)
    type(dynamic_shear_modulus), intent(in) :: modulus
    real(dp) :: loss_factor

    loss_factor = modulus%loss_modulus / modulus%storage_modulus
  end function calculate_loss_factor

end module tms_dynamic_modulus
