module tms_torsional_system
  use tms_kinds, only : dp
  use tms_geometry, only : tvd_geometry_t
  use tms_material, only : dynamic_rubber_material_t
  use tms_inertia, only : calculate_annular_hub_properties, &
    calculate_annular_ring_properties
  use tms_dynamic_torsional_stiffness, only : &
    complex_torsional_stiffness_t, calculate_dynamic_torsional_stiffness
  use tms_frequency_solver, only : calculate_natural_frequency
  implicit none
  private

  !> Göbek, bağlı annüler elastomer ve atalet halkasından oluşan ayrık
  !! parametreli TVD sisteminin modal hesap için gerekli fiziksel verileridir.
  !! Elastomer kütlesi ve polar ataleti bu iki rijit gövdeli modelde ihmal
  !! edilir.
  type, public :: two_inertia_tvd_system_t
    !> Göbeğin dönme ekseni çevresindeki polar kütle ataleti J_h [kg*m^2].
    real(dp) :: hub_polar_inertia_kg_m2 = 0.0_dp

    !> Atalet halkasının dönme ekseni çevresindeki polar kütle ataleti
    !! J_r [kg*m^2].
    real(dp) :: ring_polar_inertia_kg_m2 = 0.0_dp

    !> Elastomerin enerji depolayan burulma rijitliği K' [N*m/rad].
    real(dp) :: storage_stiffness_nm_per_rad = 0.0_dp

    !> Elastomerin çevrimsel enerji kaybını temsil eden burulma rijitliği
    !! K'' [N*m/rad]. V0.2.2 sönümsüz modal hesabında kullanılmaz.
    real(dp) :: loss_stiffness_nm_per_rad = 0.0_dp

    !> Kompleks rijitliğin boyutsuz kayıp faktörü tan(delta) = K''/K' [-].
    !! V0.2.2 bu değerden sönüm oranı hesaplamaz.
    real(dp) :: loss_factor = 0.0_dp

    !> G' ve G'' değerlerinin alındığı malzeme referans frekansı [Hz].
    !! Modal sonuç farklı olsa bile bu sürümde frekans iterasyonu yapılmaz.
    real(dp) :: material_reference_frequency_hz = 0.0_dp

    !> Dinamik malzeme özelliklerinin geçerli olduğu mutlak sıcaklık [K].
    real(dp) :: material_temperature_k = 0.0_dp
  end type two_inertia_tvd_system_t

  !> Serbest-serbest iki ataletli sistemin analitik modal sonucunu taşır.
  !! Mod şekillerinin DOF sırası [göbek, atalet halkası] olarak sabittir.
  type, public :: two_inertia_modal_result_t
    !> Birlikte dönme rijit-cisim modunun doğal frekansı f_0 [Hz].
    real(dp) :: rigid_body_frequency_hz = 0.0_dp

    !> Göbek ile halkanın zıt yönde döndüğü elastik mod frekansı f_e [Hz].
    real(dp) :: elastic_frequency_hz = 0.0_dp

    !> Rijit-cisim modunda birim değere normalize edilmiş göbek genliği [-].
    real(dp) :: rigid_body_mode_hub = 0.0_dp

    !> Rijit-cisim modunda birim değere normalize edilmiş halka genliği [-].
    real(dp) :: rigid_body_mode_ring = 0.0_dp

    !> Elastik modda birim değere normalize edilmiş göbek genliği [-].
    real(dp) :: elastic_mode_hub = 0.0_dp

    !> Elastik modda kütle-ortogonalliğinden bulunan halka genliği [-].
    real(dp) :: elastic_mode_ring = 0.0_dp
  end type two_inertia_modal_result_t

  public :: build_two_inertia_tvd_system
  public :: calculate_fixed_hub_natural_frequency
  public :: solve_free_free_two_inertia_modes

contains

  !> Bileşen geometrisi ve dinamik malzemeden iki ataletli TVD sistemi kurar.
  !!
  !! Fiziksel açıklama: Homojen rijit göbek ve halkanın polar ataletleri ile
  !! tam bağlı annüler elastomerin K', K'' ve kayıp faktörü tek bir sistem
  !! durumunda birleştirilir.
  !! Matematiksel açıklama: Bu yordam atalet veya rijitlik denklemlerini
  !! yinelemez; mevcut bileşen üretim yordamlarının sonuçlarını aktarır.
  !! Girdiler: Geometri uzunlukları m, göbek ve halka yoğunlukları kg/m^3,
  !! malzeme modülleri Pa, referans frekansı Hz ve sıcaklığı K cinsindendir.
  !! Çıktı: Ataletleri kg*m^2, rijitlikleri N*m/rad ve çalışma noktası
  !! bilgilerini SI birimleriyle taşıyan two_inertia_tvd_system_t değeridir.
  !! Varsayımlar ve geçerlilik: Göbek ve halka homojen ve rijit, elastomer
  !! lineer, küçük deformasyonlu ve yüzeylere tam bağlıdır. Elastomer kütlesi
  !! ve polar ataleti M matrisine eklenmez. Yoğunluklar, geometriler ve dinamik
  !! modüller ilgili üretim yordamlarınca doğrulanır.
  !! Ayrıntılar: docs/physics/two_inertia_torsional_system.md.
  pure function build_two_inertia_tvd_system( &
      geometry, hub_density_kg_m3, ring_density_kg_m3, material) result(system)
    type(tvd_geometry_t), intent(in) :: geometry
    real(dp), intent(in) :: hub_density_kg_m3
    real(dp), intent(in) :: ring_density_kg_m3
    type(dynamic_rubber_material_t), intent(in) :: material
    type(two_inertia_tvd_system_t) :: system

    real(dp) :: unused_hub_volume_m3
    real(dp) :: unused_hub_mass_kg
    real(dp) :: unused_ring_mass_kg
    type(complex_torsional_stiffness_t) :: stiffness

    call calculate_annular_hub_properties( &
      geometry%hub, hub_density_kg_m3, unused_hub_volume_m3, &
      unused_hub_mass_kg, system%hub_polar_inertia_kg_m2)
    call calculate_annular_ring_properties( &
      geometry%inertia_ring, ring_density_kg_m3, unused_ring_mass_kg, &
      system%ring_polar_inertia_kg_m2)
    stiffness = calculate_dynamic_torsional_stiffness( &
      material, geometry%rubber)

    system%storage_stiffness_nm_per_rad = stiffness%storage_stiffness
    system%loss_stiffness_nm_per_rad = stiffness%loss_stiffness
    system%loss_factor = stiffness%loss_factor
    system%material_reference_frequency_hz = stiffness%frequency
    system%material_temperature_k = stiffness%temperature
  end function build_two_inertia_tvd_system

  !> Göbeği sabitlenmiş TVD sisteminin sönümsüz doğal frekansını hesaplar.
  !!
  !! Fiziksel açıklama: Göbek hareket etmez; atalet halkası K' depolama
  !! rijitliği altında tek dönel serbestlik derecesiyle salınır. K'' modal
  !! denkleme dahil edilmez.
  !! Matematiksel açıklama: J_r*theta_r'' + K'*theta_r = 0 ve
  !! f_n = sqrt(K'/J_r)/(2*pi). Uygulama mevcut calculate_natural_frequency
  !! yordamını yeniden kullanır.
  !! Girdi: Sistem ataletleri kg*m^2, K' N*m/rad cinsindedir.
  !! Çıktı: Sönümsüz doğal frekans f_n hertz (Hz) cinsindendir.
  !! Varsayımlar ve geçerlilik: J_h > 0, J_r > 0 ve K' > 0 olmalıdır.
  !! G' referans frekansı ve sıcaklığında sabitlenir; modal frekansla
  !! öz-tutarlı interpolasyon veya iterasyon yapılmaz.
  pure function calculate_fixed_hub_natural_frequency(system) &
      result(frequency_hz)
    type(two_inertia_tvd_system_t), intent(in) :: system
    real(dp) :: frequency_hz

    call validate_modal_system(system)
    frequency_hz = calculate_natural_frequency( &
      system%storage_stiffness_nm_per_rad, &
      system%ring_polar_inertia_kg_m2)
  end function calculate_fixed_hub_natural_frequency

  !> Serbest-serbest iki ataletli TVD sisteminin iki analitik modunu çözer.
  !!
  !! Fiziksel açıklama: Rijit-cisim modunda göbek ve halka birlikte döner ve
  !! elastomer şekil değiştirmez. Elastik modda iki rijit gövde zıt yönde
  !! dönerek elastomeri burar.
  !! Matematiksel açıklama: omega_0 = 0 ve
  !! omega_e^2 = K'*(1/J_h + 1/J_r). Eşdeğer atalet
  !! J_eq = 1/(1/J_h + 1/J_r) ile mevcut doğal frekans yordamı kullanılır.
  !! DOF sırası [göbek, halka] olup modlar [1,1] ve [1,-J_h/J_r] biçiminde
  !! göbek genliğine göre normalize edilir.
  !! Girdi: J_h ve J_r kg*m^2, K' N*m/rad cinsindedir.
  !! Çıktı: Frekansları Hz ve boyutsuz normalize mod genliklerini taşıyan
  !! two_inertia_modal_result_t değeridir.
  !! Varsayımlar ve geçerlilik: Sistem lineer ve sönümsüzdür; J_h > 0,
  !! J_r > 0 ve K' > 0 olmalıdır. Elastomer kütlesi ve polar ataleti ihmal
  !! edilir. K'' saklanır ancak bu reel analitik özdeğer tahmininde kullanılmaz.
  !! Genel özdeğer çözümü yapılmaz.
  pure function solve_free_free_two_inertia_modes(system) result(modal_result)
    type(two_inertia_tvd_system_t), intent(in) :: system
    type(two_inertia_modal_result_t) :: modal_result

    real(dp) :: equivalent_inertia_kg_m2

    call validate_modal_system(system)

    equivalent_inertia_kg_m2 = 1.0_dp / &
      (1.0_dp / system%hub_polar_inertia_kg_m2 + &
       1.0_dp / system%ring_polar_inertia_kg_m2)

    modal_result%rigid_body_frequency_hz = 0.0_dp
    modal_result%elastic_frequency_hz = calculate_natural_frequency( &
      system%storage_stiffness_nm_per_rad, equivalent_inertia_kg_m2)
    modal_result%rigid_body_mode_hub = 1.0_dp
    modal_result%rigid_body_mode_ring = 1.0_dp
    modal_result%elastic_mode_hub = 1.0_dp
    modal_result%elastic_mode_ring = &
      -system%hub_polar_inertia_kg_m2 / system%ring_polar_inertia_kg_m2
  end function solve_free_free_two_inertia_modes

  !> Modal sistemin pozitif atalet ve depolama rijitliği koşullarını doğrular.
  !!
  !! Fiziksel açıklama: Sıfır veya negatif polar atalet rijit gövde modelini,
  !! pozitif olmayan K' ise kararlı elastik geri çağırıcı moment varsayımını
  !! ihlal eder.
  !! Matematiksel açıklama: J_h > 0, J_r > 0 ve K' > 0 koşulları, karekök
  !! içinin pozitif ve modal frekansın reel olması için zorunludur.
  !! Girdi: J_h ve J_r kg*m^2, K' N*m/rad cinsindedir. Çıktı üretmez;
  !! geçersiz sistem error stop ile reddedilir.
  pure subroutine validate_modal_system(system)
    type(two_inertia_tvd_system_t), intent(in) :: system

    if (system%hub_polar_inertia_kg_m2 <= 0.0_dp) then
      error stop "Göbek polar kütle ataleti sıfırdan büyük olmalıdır."
    end if

    if (system%ring_polar_inertia_kg_m2 <= 0.0_dp) then
      error stop "Atalet halkası polar kütle ataleti sıfırdan büyük olmalıdır."
    end if

    if (system%storage_stiffness_nm_per_rad <= 0.0_dp) then
      error stop "Depolama burulma rijitliği K' sıfırdan büyük olmalıdır."
    end if
  end subroutine validate_modal_system

end module tms_torsional_system
