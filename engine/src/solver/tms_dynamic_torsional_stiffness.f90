module tms_dynamic_torsional_stiffness
  use tms_kinds, only : dp
  use tms_geometry, only : rubber_geometry_t, &
    calculate_rubber_polar_area_moment
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

  !> Annüler elastomer bölgenin kompleks burulma rijitliğini hesaplar.
  !!
  !! Fiziksel açıklama: K' elastomerin bağıl dönmeye karşı elastik ve enerji
  !! depolayan tepkisini, K'' ise faz dışı ve enerji kaybettiren tepkisini
  !! temsil eder. Frekans ve sıcaklık malzemenin çalışma noktasından aktarılır.
  !! Matematiksel açıklama: Jp = pi/2*(ro^4-ri^4), K' = G'*Jp/L,
  !! K'' = G''*Jp/L ve tan(delta) = K''/K' bağıntıları uygulanır.
  !! Girdiler: material içindeki G' ve G'' Pa, frekans Hz, sıcaklık K;
  !! rubber içindeki ri, ro ve etkin uzunluk L metre (m) cinsindendir.
  !! Çıktı: K' ve K'' N*m/rad, kayıp faktörü boyutsuz, frekans Hz ve
  !! sıcaklık K alanlarından oluşan complex_torsional_stiffness_t değeridir.
  !! Varsayımlar ve geçerlilik: Elastomer homojen, lineer viskoelastik ve
  !! küçük deformasyon bölgesindedir; 0 <= ri < ro, L > 0, G' > 0 ve
  !! G'' >= 0 olmalıdır. Yordam bu fiziksel önkoşulları doğrular ve
  !! geçersiz girdide error stop ile sonlanır. frequency_points dizisinde
  !! seçim veya interpolasyon yapmaz; mevcut tek çalışma noktası alanlarını
  !! kullanır. Ayrıntılar:
  !! docs/mathematics/dynamic_elastomer_model.md ve
  !! docs/physics/complex_torsional_stiffness.md.
  pure function calculate_dynamic_torsional_stiffness(material, rubber) &
      result(stiffness)
    type(dynamic_rubber_material_t), intent(in) :: material
    type(rubber_geometry_t), intent(in) :: rubber
    type(complex_torsional_stiffness_t) :: stiffness

    real(dp) :: polar_area_moment_m4

    if (rubber%axial_length_m <= 0.0_dp) then
      error stop "Elastomer etkin uzunluğu sıfırdan büyük olmalıdır."
    end if

    if (material%storage_shear_modulus_pa <= 0.0_dp) then
      error stop "Depolama kayma modülü G' sıfırdan büyük olmalıdır."
    end if

    if (material%loss_shear_modulus_pa < 0.0_dp) then
      error stop "Kayıp kayma modülü G'' negatif olamaz."
    end if

    polar_area_moment_m4 = calculate_rubber_polar_area_moment( &
      rubber%outer_radius_m, rubber%inner_radius_m)
    stiffness%storage_stiffness = material%storage_shear_modulus_pa * &
      polar_area_moment_m4 / rubber%axial_length_m
    stiffness%loss_stiffness = material%loss_shear_modulus_pa * &
      polar_area_moment_m4 / rubber%axial_length_m
    stiffness%loss_factor = stiffness%loss_stiffness / &
      stiffness%storage_stiffness
    stiffness%frequency = material%frequency_hz
    stiffness%temperature = material%temperature_k
  end function calculate_dynamic_torsional_stiffness

end module tms_dynamic_torsional_stiffness
