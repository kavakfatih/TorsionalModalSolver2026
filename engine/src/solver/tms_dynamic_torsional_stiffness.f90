module tms_dynamic_torsional_stiffness
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_geometry, only : rubber_geometry_t, &
    calculate_annular_bush_torsion_geometry_factor
  use tms_material, only : dynamic_rubber_material_t
  implicit none
  private

  !> Kompleks burulma rijitliğinin reel ve sanal bileşenlerini ve bunların
  !! geçerli olduğu dinamik malzeme çalışma noktasını birlikte taşır.
  type, public :: complex_torsional_stiffness_t
    !> Depolama rijitliği K' [N*m/rad]. Elastomerin elastik enerji depolayan
    !! davranışını temsil eden kompleks burulma rijitliğinin reel bileşenidir.
    real(dp) :: storage_stiffness = 0.0_dp

    !> Kayıp rijitliği K'' [N*m/rad]. Elastomerin çevrimsel enerji kaybını
    !! temsil eden kompleks burulma rijitliğinin sanal bileşenidir.
    real(dp) :: loss_stiffness = 0.0_dp

    !> Kayıp faktörü tan(delta) = K''/K' = G''/G' [-]. Kompleks rijitliğin
    !! sönüm bileşeninin depolama bileşenine boyutsuz oranıdır.
    real(dp) :: loss_factor = 0.0_dp

    !> Rijitlik bileşenlerinin geçerli olduğu uyarım frekansı f [Hz].
    !! Matematiksel olarak saniye başına çevrim sayısını ifade eder.
    real(dp) :: frequency = 0.0_dp

    !> Rijitlik bileşenlerinin geçerli olduğu mutlak sıcaklık T [K].
    !! Elastomerin sıcaklığa bağlı dinamik davranışının durum değişkenidir.
    real(dp) :: temperature = 0.0_dp
  end type complex_torsional_stiffness_t

  public :: calculate_dynamic_torsional_stiffness

contains

  !> Tam bağlı annüler TVD elastomerinin kompleks rijitliğini hesaplar.
  !!
  !! Fiziksel açıklama: Rijit iç göbek ve rijit dış atalet halkasına
  !! tam bağlı eş merkezli silindirik elastomer tabakasında K' enerji
  !! depolayan, K'' ise enerji kaybettiren moment tepkisidir. Frekans ve
  !! sıcaklık malzemenin çalışma noktasından aktarılır.
  !! Matematiksel açıklama: C_theta =
  !! 4*pi*L*ri^2*ro^2/(ro^2-ri^2), K' = G'*C_theta,
  !! K'' = G''*C_theta ve tan(delta) = K''/K' bağıntıları uygulanır.
  !! Girdiler: material içindeki G' ve G'' Pa, frekans Hz, sıcaklık K;
  !! rubber içindeki ri, ro ve bağlı eksenel genişlik L metre (m) cinsindendir.
  !! Çıktı: K' ve K'' N*m/rad, kayıp faktörü boyutsuz, frekans Hz ve
  !! sıcaklık K alanlarından oluşan complex_torsional_stiffness_t değeridir.
  !! Varsayımlar ve geçerlilik: Elastomer homojen, lineer viskoelastik ve
  !! küçük deformasyon bölgesindedir; 0 < ri < ro, L > 0, G' > 0 ve
  !! G'' >= 0, frekans >= 0 Hz ve mutlak sıcaklık > 0 K olmalıdır. Tüm
  !! girdiler ve hesaplanan çıktılar sonlu olmalıdır. Yoğunluk bu denklemde
  !! kullanılmaz; malzeme bütünlüğü için sonlu olduğu doğrulanır. Yordam bu
  !! fiziksel önkoşulları doğrular ve geçersiz girdide error stop ile
  !! sonlanır. frequency_points dizisinde seçim veya interpolasyon yapmaz;
  !! mevcut tek çalışma noktası alanlarını kullanır. Ayrıntılar:
  !! docs/mathematics/dynamic_elastomer_model.md ve
  !! docs/physics/complex_torsional_stiffness.md.
  pure function calculate_dynamic_torsional_stiffness(material, rubber) &
      result(stiffness)
    type(dynamic_rubber_material_t), intent(in) :: material
    type(rubber_geometry_t), intent(in) :: rubber
    type(complex_torsional_stiffness_t) :: stiffness

    real(dp) :: geometry_factor_m3

    if (.not. ieee_is_finite(material%storage_shear_modulus_pa) .or. &
        material%storage_shear_modulus_pa <= 0.0_dp) then
      error stop "Depolama kayma modülü G' sonlu ve pozitif olmalıdır."
    end if

    if (.not. ieee_is_finite(material%loss_shear_modulus_pa) .or. &
        material%loss_shear_modulus_pa < 0.0_dp) then
      error stop "Kayıp kayma modülü G'' sonlu ve negatif olmayan bir değer olmalıdır."
    end if

    if (.not. ieee_is_finite(material%frequency_hz) .or. &
        material%frequency_hz < 0.0_dp) then
      error stop "Malzeme frekansı sonlu ve negatif olmayan bir değer olmalıdır."
    end if

    if (.not. ieee_is_finite(material%temperature_k) .or. &
        material%temperature_k <= 0.0_dp) then
      error stop "Malzeme sıcaklığı sonlu ve pozitif olmalıdır."
    end if

    if (.not. ieee_is_finite(material%density_kg_m3)) then
      error stop "Elastomer yoğunluğu sonlu olmalıdır."
    end if

    geometry_factor_m3 = calculate_annular_bush_torsion_geometry_factor( &
      inner_radius=rubber%inner_radius_m, &
      outer_radius=rubber%outer_radius_m, &
      axial_length=rubber%axial_length_m)
    stiffness%storage_stiffness = material%storage_shear_modulus_pa * &
      geometry_factor_m3

    if (.not. ieee_is_finite(stiffness%storage_stiffness) .or. &
        stiffness%storage_stiffness <= 0.0_dp) then
      error stop "Depolama rijitliği K' sonlu ve pozitif olmalıdır."
    end if

    stiffness%loss_stiffness = material%loss_shear_modulus_pa * &
      geometry_factor_m3

    if (.not. ieee_is_finite(stiffness%loss_stiffness) .or. &
        stiffness%loss_stiffness < 0.0_dp) then
      error stop "Kayıp rijitliği K'' sonlu ve negatif olmayan bir değer olmalıdır."
    end if

    stiffness%loss_factor = stiffness%loss_stiffness / &
      stiffness%storage_stiffness

    if (.not. ieee_is_finite(stiffness%loss_factor) .or. &
        stiffness%loss_factor < 0.0_dp) then
      error stop "Kayıp faktörü sonlu ve negatif olmayan bir değer olmalıdır."
    end if

    stiffness%frequency = material%frequency_hz
    stiffness%temperature = material%temperature_k
  end function calculate_dynamic_torsional_stiffness

end module tms_dynamic_torsional_stiffness
