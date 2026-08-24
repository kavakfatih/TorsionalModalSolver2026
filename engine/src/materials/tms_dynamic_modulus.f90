module tms_dynamic_modulus
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
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
  !! Varsayım ve geçerlilik sınırı: Tüm alanlar sonlu; G' > 0, G'' >= 0,
  !! frekans >= 0 Hz ve mutlak sıcaklık > 0 K olmalıdır. Fonksiyon
  !! interpolasyon, eğri uydurma veya sıcaklık-frekans dönüşümü yapmaz.
  !! Ayrıntı: docs/mathematics/dynamic_elastomer_model.md
  pure elemental function calculate_loss_factor(modulus) result(loss_factor)
    type(dynamic_shear_modulus), intent(in) :: modulus
    real(dp) :: loss_factor

    if (.not. ieee_is_finite(modulus%storage_modulus) .or. &
        modulus%storage_modulus <= 0.0_dp) then
      error stop "Depolama modülü G' sonlu ve pozitif olmalıdır."
    end if

    if (.not. ieee_is_finite(modulus%loss_modulus) .or. &
        modulus%loss_modulus < 0.0_dp) then
      error stop "Kayıp modülü G'' sonlu ve negatif olmayan bir değer olmalıdır."
    end if

    if (.not. ieee_is_finite(modulus%frequency) .or. &
        modulus%frequency < 0.0_dp) then
      error stop "Dinamik modül frekansı sonlu ve negatif olmayan bir değer olmalıdır."
    end if

    if (.not. ieee_is_finite(modulus%temperature) .or. &
        modulus%temperature <= 0.0_dp) then
      error stop "Dinamik modül sıcaklığı sonlu ve pozitif olmalıdır."
    end if

    loss_factor = modulus%loss_modulus / modulus%storage_modulus

    if (.not. ieee_is_finite(loss_factor) .or. loss_factor < 0.0_dp) then
      error stop "Kayıp faktörü sonlu ve negatif olmayan bir değer olmalıdır."
    end if
  end function calculate_loss_factor

end module tms_dynamic_modulus
