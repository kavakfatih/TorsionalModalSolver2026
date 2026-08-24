program test_torsional_system
  use tms_kinds, only : dp
  use tms_geometry, only : hub_geometry_t, inertia_ring_geometry_t, &
    rubber_geometry_t, tvd_geometry_t
  use tms_material, only : dynamic_rubber_material_t
  use tms_inertia, only : calculate_annular_hub_properties, &
    calculate_annular_ring_properties
  use tms_dynamic_torsional_stiffness, only : &
    complex_torsional_stiffness_t, calculate_dynamic_torsional_stiffness
  use tms_frequency_solver, only : calculate_natural_frequency
  use tms_torsional_system, only : two_inertia_tvd_system_t, &
    two_inertia_modal_result_t, build_two_inertia_tvd_system, &
    calculate_fixed_hub_natural_frequency, &
    solve_free_free_two_inertia_modes
  implicit none

  real(dp), parameter :: maximum_relative_error = 1.0e-10_dp
  real(dp), parameter :: maximum_absolute_error = 1.0e-12_dp
  real(dp), parameter :: expected_fixed_frequency_hz = &
    11.253953951963826_dp
  real(dp), parameter :: expected_elastic_frequency_hz = &
    19.492420030841906_dp

  type(tvd_geometry_t) :: geometry
  type(dynamic_rubber_material_t) :: material
  type(two_inertia_tvd_system_t) :: built_system
  type(two_inertia_tvd_system_t) :: reference_system
  type(two_inertia_tvd_system_t) :: scaled_system
  type(two_inertia_tvd_system_t) :: limit_system
  type(two_inertia_modal_result_t) :: modal_result
  type(two_inertia_modal_result_t) :: scaled_modal_result
  type(two_inertia_modal_result_t) :: limit_modal_result
  type(complex_torsional_stiffness_t) :: direct_stiffness
  real(dp) :: direct_hub_volume_m3
  real(dp) :: direct_hub_mass_kg
  real(dp) :: direct_hub_inertia_kg_m2
  real(dp) :: direct_ring_mass_kg
  real(dp) :: direct_ring_inertia_kg_m2
  real(dp) :: fixed_frequency_hz
  real(dp) :: direct_frequency_hz
  real(dp) :: limit_fixed_frequency_hz
  real(dp) :: loss_metadata_fixed_frequency_hz
  character(len=40) :: validation_case

  ! CTest geçersiz sistem ve builder girdilerini ayrı WILL_FAIL vakaları
  ! olarak aynı üretim arayüzleri üzerinden sınar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call initialize_builder_inputs(geometry, material)
  built_system = build_two_inertia_tvd_system( &
    geometry, 7800.0_dp, 7800.0_dp, material)

  ! Builder testinde fizik denklemleri kopyalanmaz; aynı girdilerle çağrılan
  ! mevcut bileşen üretim yordamlarının sonuçları sistem alanlarıyla kıyaslanır.
  call calculate_annular_hub_properties( &
    geometry%hub, 7800.0_dp, direct_hub_volume_m3, &
    direct_hub_mass_kg, direct_hub_inertia_kg_m2)
  call calculate_annular_ring_properties( &
    geometry%inertia_ring, 7800.0_dp, direct_ring_mass_kg, &
    direct_ring_inertia_kg_m2)
  direct_stiffness = calculate_dynamic_torsional_stiffness( &
    material, geometry%rubber)

  call assert_relative_close( &
    built_system%hub_polar_inertia_kg_m2, direct_hub_inertia_kg_m2, &
    "Builder göbek ataletini üretim yordamından aktaramadı.")
  call assert_relative_close( &
    built_system%ring_polar_inertia_kg_m2, direct_ring_inertia_kg_m2, &
    "Builder halka ataletini üretim yordamından aktaramadı.")
  call assert_relative_close( &
    built_system%storage_stiffness_nm_per_rad, &
    direct_stiffness%storage_stiffness, &
    "Builder K' depolama rijitliğini aktaramadı.")
  call assert_relative_close( &
    built_system%loss_stiffness_nm_per_rad, &
    direct_stiffness%loss_stiffness, &
    "Builder K'' kayıp rijitliğini aktaramadı.")
  call assert_relative_close( &
    built_system%loss_factor, direct_stiffness%loss_factor, &
    "Builder kayıp faktörünü aktaramadı.")
  call assert_relative_close( &
    built_system%material_reference_frequency_hz, &
    material%frequency_hz, &
    "Builder malzeme referans frekansını aktaramadı.")
  call assert_relative_close( &
    built_system%material_temperature_k, material%temperature_k, &
    "Builder malzeme sıcaklığını aktaramadı.")

  reference_system = two_inertia_tvd_system_t( &
    hub_polar_inertia_kg_m2=0.10_dp, &
    ring_polar_inertia_kg_m2=0.20_dp, &
    storage_stiffness_nm_per_rad=1000.0_dp, &
    loss_stiffness_nm_per_rad=100.0_dp, &
    loss_factor=0.1_dp, &
    material_reference_frequency_hz=100.0_dp, &
    material_temperature_k=293.15_dp)

  fixed_frequency_hz = &
    calculate_fixed_hub_natural_frequency(reference_system)
  direct_frequency_hz = calculate_natural_frequency( &
    reference_system%storage_stiffness_nm_per_rad, &
    reference_system%ring_polar_inertia_kg_m2)
  modal_result = solve_free_free_two_inertia_modes(reference_system)

  ! Fixed-hub sistem arayüzü aynı fiziğin ikinci bir uygulamasını kurmak
  ! yerine mevcut tek-DOF frekans yordamıyla bire bir aynı sonucu vermelidir.
  call assert_relative_close( &
    fixed_frequency_hz, direct_frequency_hz, &
    "Fixed-hub sonucu mevcut doğal frekans yordamıyla aynı değil.")
  call assert_relative_close( &
    fixed_frequency_hz, expected_fixed_frequency_hz, &
    "Fixed-hub benchmark frekansı bağımsız referansla uyuşmuyor.")

  call assert_absolute_close( &
    modal_result%rigid_body_frequency_hz, 0.0_dp, &
    "Serbest-serbest rijit-cisim modu sıfır frekanslı değil.")
  call assert_relative_close( &
    modal_result%elastic_frequency_hz, expected_elastic_frequency_hz, &
    "Serbest-serbest elastik mod frekansı referansla uyuşmuyor.")
  call assert_absolute_close( &
    modal_result%rigid_body_mode_hub, 1.0_dp, &
    "Rijit-cisim modunun göbek genliği bir değil.")
  call assert_absolute_close( &
    modal_result%rigid_body_mode_ring, 1.0_dp, &
    "Rijit-cisim modunun halka genliği bir değil.")
  call assert_absolute_close( &
    modal_result%elastic_mode_hub, 1.0_dp, &
    "Elastik modun normalize göbek genliği bir değil.")
  call assert_absolute_close( &
    modal_result%elastic_mode_ring, -0.5_dp, &
    "Elastik modun normalize halka genliği -0,5 değil.")

  ! Sönümsüz V0.2.2 modeli yalnız K' kullanır. K'' ile kayıp faktörünü
  ! değiştirmek, üretim formülünü testte kopyalamadan fixed-hub veya free-free
  ! frekanslarını ve mod şekillerini değiştirmemelidir.
  scaled_system = reference_system
  scaled_system%loss_stiffness_nm_per_rad = 900.0_dp
  scaled_system%loss_factor = 0.9_dp
  scaled_modal_result = solve_free_free_two_inertia_modes(scaled_system)
  loss_metadata_fixed_frequency_hz = &
    calculate_fixed_hub_natural_frequency(scaled_system)
  call assert_relative_close( &
    loss_metadata_fixed_frequency_hz, fixed_frequency_hz, &
    "K'' değişimi fixed-hub sönümsüz frekansını etkiledi.")
  call assert_relative_close( &
    scaled_modal_result%elastic_frequency_hz, &
    modal_result%elastic_frequency_hz, &
    "K'' değişimi serbest-serbest elastik frekansı etkiledi.")
  call assert_absolute_close( &
    scaled_modal_result%rigid_body_frequency_hz, &
    modal_result%rigid_body_frequency_hz, &
    "K'' değişimi rijit-cisim frekansını etkiledi.")
  call assert_absolute_close( &
    scaled_modal_result%elastic_mode_ring, modal_result%elastic_mode_ring, &
    "K'' değişimi elastik mod şeklini etkiledi.")

  ! Ölçekleme regresyonu üretim formülünü testte yeniden hesaplamaz.
  ! K' dört katına çıktığında f ∝ sqrt(K') gereği elastik frekans iki katına
  ! çıkmalıdır; iki üretim solver sonucu doğrudan karşılaştırılır.
  scaled_system = reference_system
  scaled_system%storage_stiffness_nm_per_rad = &
    4.0_dp * reference_system%storage_stiffness_nm_per_rad
  scaled_modal_result = solve_free_free_two_inertia_modes(scaled_system)
  call assert_relative_close( &
    scaled_modal_result%elastic_frequency_hz, &
    2.0_dp * modal_result%elastic_frequency_hz, &
    "K' dört katına çıktığında elastik frekans iki katına çıkmadı.")

  ! Atalet ölçekleme regresyonu formülü kopyalamaz. J_h ve J_r birlikte dört
  ! katına çıktığında f ∝ 1/sqrt(J) gereği elastik frekans yarıya inmelidir.
  scaled_system = reference_system
  scaled_system%hub_polar_inertia_kg_m2 = &
    4.0_dp * reference_system%hub_polar_inertia_kg_m2
  scaled_system%ring_polar_inertia_kg_m2 = &
    4.0_dp * reference_system%ring_polar_inertia_kg_m2
  scaled_modal_result = solve_free_free_two_inertia_modes(scaled_system)
  call assert_relative_close( &
    scaled_modal_result%elastic_frequency_hz, &
    0.5_dp * modal_result%elastic_frequency_hz, &
    "Ataletler dört katına çıktığında elastik frekans yarıya inmedi.")

  ! J_h çok büyüdüğünde göbek hareketi kaybolur; serbest-serbest elastik mod
  ! aynı sistemin fixed-hub frekansına yaklaşmalıdır. Bu bir limit regresyonu
  ! olup denklem test içinde yeniden kurulmaz.
  limit_system = reference_system
  limit_system%hub_polar_inertia_kg_m2 = 1.0e12_dp
  limit_modal_result = solve_free_free_two_inertia_modes(limit_system)
  limit_fixed_frequency_hz = &
    calculate_fixed_hub_natural_frequency(limit_system)
  call assert_relative_close( &
    limit_modal_result%elastic_frequency_hz, limit_fixed_frequency_hz, &
    "Büyük göbek ataleti limitinde free-free sonuç fixed-hub'a yaklaşmadı.")

  print *, "Fixed-hub ve serbest-serbest iki ataletli TVD sistemi doğrulandı."

contains

  !> Builder entegrasyon testinde kullanılan geçerli TVD girdilerini kurar.
  !! Fiziksel model: Eş merkezli homojen göbek ve halka, aralarında tam bağlı
  !! annüler elastomer. Uzunluklar m, yoğunluklar çağrıda kg/m^3, modüller Pa,
  !! frekans Hz ve sıcaklık K cinsindedir. Çıktılar geometry ve material'dır.
  subroutine initialize_builder_inputs(test_geometry, test_material)
    type(tvd_geometry_t), intent(out) :: test_geometry
    type(dynamic_rubber_material_t), intent(out) :: test_material

    test_geometry%hub = hub_geometry_t( &
      bore_radius_m=0.01_dp, &
      outer_radius_m=0.03_dp, &
      axial_length_m=0.02_dp)
    test_geometry%rubber = rubber_geometry_t( &
      inner_radius_m=0.03_dp, &
      outer_radius_m=0.04_dp, &
      axial_length_m=0.02_dp)
    test_geometry%inertia_ring = inertia_ring_geometry_t( &
      inner_radius_m=0.04_dp, &
      outer_radius_m=0.06_dp, &
      axial_length_m=0.02_dp)

    test_material%name = "EPDM sistem örneği"
    test_material%density_kg_m3 = 1100.0_dp
    test_material%storage_shear_modulus_pa = 1.0e6_dp
    test_material%loss_shear_modulus_pa = 0.1e6_dp
    test_material%frequency_hz = 100.0_dp
    test_material%temperature_k = 293.15_dp
  end subroutine initialize_builder_inputs

  !> Geçersiz fiziksel sistem veya builder girdisinin üretim yordamınca
  !! reddedildiğini sınayan seçili CTest vakasını çalıştırır.
  !! J değerleri kg*m^2, K' N*m/rad, yoğunluk kg/m^3 ve uzunluklar m
  !! cinsindedir. Çıktı üretmez; geçersiz girdi error stop oluşturmalıdır.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(tvd_geometry_t) :: invalid_geometry
    type(dynamic_rubber_material_t) :: valid_material
    type(two_inertia_tvd_system_t) :: invalid_system
    type(two_inertia_tvd_system_t) :: rejected_system
    type(two_inertia_modal_result_t) :: rejected_modal_result

    call initialize_builder_inputs(invalid_geometry, valid_material)
    invalid_system = two_inertia_tvd_system_t( &
      hub_polar_inertia_kg_m2=0.10_dp, &
      ring_polar_inertia_kg_m2=0.20_dp, &
      storage_stiffness_nm_per_rad=1000.0_dp)

    select case (case_name)
      case ("nonpositive_hub_inertia")
        invalid_system%hub_polar_inertia_kg_m2 = 0.0_dp
        rejected_modal_result = &
          solve_free_free_two_inertia_modes(invalid_system)
      case ("nonpositive_ring_inertia")
        invalid_system%ring_polar_inertia_kg_m2 = 0.0_dp
        rejected_modal_result = &
          solve_free_free_two_inertia_modes(invalid_system)
      case ("nonpositive_storage_stiffness")
        invalid_system%storage_stiffness_nm_per_rad = 0.0_dp
        rejected_modal_result = &
          solve_free_free_two_inertia_modes(invalid_system)
      case ("nonpositive_hub_density")
        rejected_system = build_two_inertia_tvd_system( &
          invalid_geometry, 0.0_dp, 7800.0_dp, valid_material)
      case ("nonpositive_ring_density")
        rejected_system = build_two_inertia_tvd_system( &
          invalid_geometry, 7800.0_dp, 0.0_dp, valid_material)
      case ("invalid_hub_geometry")
        invalid_geometry%hub%outer_radius_m = &
          invalid_geometry%hub%bore_radius_m
        rejected_system = build_two_inertia_tvd_system( &
          invalid_geometry, 7800.0_dp, 7800.0_dp, valid_material)
      case ("invalid_ring_geometry")
        invalid_geometry%inertia_ring%outer_radius_m = &
          invalid_geometry%inertia_ring%inner_radius_m
        rejected_system = build_two_inertia_tvd_system( &
          invalid_geometry, 7800.0_dp, 7800.0_dp, valid_material)
      case default
        error stop "Bilinmeyen geçersiz torsional sistem testi istendi."
    end select

    print *, "Geçersiz sistem beklenmedik biçimde kabul edildi: ", &
      case_name, rejected_system%hub_polar_inertia_kg_m2, &
      rejected_modal_result%elastic_frequency_hz
  end subroutine exercise_invalid_case

  !> Pozitif bir fiziksel sonucun analitik referansa bağıl yakınlığını sınar.
  !! Matematiksel açıklama: |actual-expected|/|expected| < 1e-10. Girdiler
  !! aynı SI biriminde olmalıdır; koşul sağlanmazsa test error stop üretir.
  subroutine assert_relative_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (abs(actual - expected) / abs(expected) >= &
        maximum_relative_error) then
      error stop message
    end if
  end subroutine assert_relative_close

  !> Sıfır veya boyutsuz mode-shape referansını mutlak toleransla sınar.
  !! Matematiksel açıklama: |actual-expected| <= 1e-12. Girdiler aynı birimde
  !! veya boyutsuz olmalıdır; koşul sağlanmazsa test error stop üretir.
  subroutine assert_absolute_close(actual, expected, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > maximum_absolute_error) then
      error stop message
    end if
  end subroutine assert_absolute_close

end program test_torsional_system
