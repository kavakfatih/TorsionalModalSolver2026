program test_numerical_hardening
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_geometry, only : rubber_geometry_t, inertia_ring_geometry_t, &
    hub_geometry_t, tvd_geometry_t, validate_tvd_geometry, &
    calculate_rubber_polar_area_moment, &
    calculate_annular_bush_torsion_geometry_factor
  use tms_material, only : dynamic_rubber_material_t
  use tms_dynamic_modulus, only : dynamic_shear_modulus, &
    calculate_loss_factor
  use tms_inertia, only : calculate_annular_ring_properties
  use tms_torsional_stiffness, only : calculate_torsional_stiffness
  use tms_dynamic_torsional_stiffness, only : &
    complex_torsional_stiffness_t, calculate_dynamic_torsional_stiffness
  use tms_frequency_solver, only : calculate_natural_frequency
  use tms_torsional_system, only : two_inertia_tvd_system_t, &
    two_inertia_modal_result_t, build_two_inertia_tvd_system, &
    solve_free_free_two_inertia_modes
  implicit none

  real(dp), parameter :: relative_tolerance = 1.0e-12_dp
  character(len=64) :: validation_case

  ! CTest her geçersiz girdiyi ayrı WILL_FAIL sürecinde aynı üretim API'lerine
  ! gönderir. API error stop üretmezse dispatcher normal çıkar ve test kırmızı
  ! kalarak sessiz kabulü görünür yapar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_interface_tolerance_acceptance()
  call test_extreme_scale_natural_frequency()
  call test_tiny_scale_two_inertia_frequency()

  print *, "V0.3.2 sayısal sağlamlaştırma regresyonları geçti."

contains

  !> Sonlu ve pozitif bir sonucu bağıl toleransla bağımsız referansa bağlar.
  !! Girdiler aynı fiziksel birimde, tolerance boyutsuzdur. Matematiksel
  !! koşul |actual-expected|/|expected| <= tolerance biçimindedir.
  subroutine assert_relative_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(tolerance) .or. &
        actual <= 0.0_dp .or. expected <= 0.0_dp .or. &
        tolerance < 0.0_dp) then
      error stop "Bağıl sayısal doğrulama girdileri sonlu ve pozitif olmalıdır."
    end if

    if (abs(actual - expected) / abs(expected) > tolerance) then
      error stop message
    end if
  end subroutine assert_relative_close

  !> TVD temas yarıçaplarının karma mutlak/bağıl tolerans içindeki küçük
  !! farklarla kabul edildiğini üretim validator ve builder üzerinden sınar.
  !! Tüm uzunluklar m, yoğunluklar kg/m^3, modüller Pa cinsindedir. Bu test
  !! toleransı yeniden hesaplamaz; nominal yüzeylere 1e-12 m fark uygular.
  subroutine test_interface_tolerance_acceptance()
    type(tvd_geometry_t) :: geometry
    type(dynamic_rubber_material_t) :: material
    type(two_inertia_tvd_system_t) :: system

    call initialize_valid_tvd_geometry(geometry)
    call initialize_valid_material(material)
    geometry%rubber%inner_radius_m = &
      geometry%hub%outer_radius_m + 1.0e-12_dp
    geometry%inertia_ring%inner_radius_m = &
      geometry%rubber%outer_radius_m - 1.0e-12_dp

    call validate_tvd_geometry(geometry)
    system = build_two_inertia_tvd_system( &
      geometry, 7800.0_dp, 7800.0_dp, material)

    if (.not. ieee_is_finite(system%hub_polar_inertia_kg_m2) .or. &
        .not. ieee_is_finite(system%ring_polar_inertia_kg_m2) .or. &
        .not. ieee_is_finite(system%storage_stiffness_nm_per_rad) .or. &
        system%hub_polar_inertia_kg_m2 <= 0.0_dp .or. &
        system%ring_polar_inertia_kg_m2 <= 0.0_dp .or. &
        system%storage_stiffness_nm_per_rad <= 0.0_dp) then
      error stop "Tolerans içindeki TVD ara yüzleri sonlu sistem üretmedi."
    end if
  end subroutine test_interface_tolerance_acceptance

  !> sqrt(K/J) ara oranı temsil aralığını aşsa bile matematiksel frekansın
  !! temsil edilebilir olduğu iki uç ölçeği sınar. K [N*m/rad], J [kg*m^2]
  !! ve çıktı Hz cinsindedir. Referanslar karekök ölçekleri ayrılarak önceden
  !! sadeleştirilmiştir; testte taşmaya açık K/J oranı kurulmaz.
  subroutine test_extreme_scale_natural_frequency()
    real(dp) :: high_ratio_frequency_hz
    real(dp) :: low_ratio_frequency_hz
    real(dp) :: expected_high_ratio_frequency_hz
    real(dp) :: expected_low_ratio_frequency_hz

    high_ratio_frequency_hz = calculate_natural_frequency( &
      1.0e300_dp, 1.0e-100_dp)
    low_ratio_frequency_hz = calculate_natural_frequency( &
      1.0e-300_dp, 1.0e200_dp)
    expected_high_ratio_frequency_hz = 1.0e200_dp / (2.0_dp * pi)
    expected_low_ratio_frequency_hz = 1.0e-250_dp / (2.0_dp * pi)

    call assert_relative_close( &
      high_ratio_frequency_hz, expected_high_ratio_frequency_hz, &
      relative_tolerance, &
      "Yüksek K/J ölçeğinde doğal frekans sonlu referansla uyuşmuyor.")
    call assert_relative_close( &
      low_ratio_frequency_hz, expected_low_ratio_frequency_hz, &
      relative_tolerance, &
      "Düşük K/J ölçeğinde doğal frekans pozitif referansla uyuşmuyor.")
  end subroutine test_extreme_scale_natural_frequency

  !> İki eşit, çok küçük ataletin harmonik birleşiminde ters alma/toplama
  !! taşması oluşmadığını sınar. J_h=J_r=K=s ve s=0.25*tiny(dp) için
  !! J_eq=s/2, K/J_eq=2 ve f=sqrt(2)/(2*pi) olur. J değerleri kg*m^2,
  !! K N*m/rad ve sonuç Hz cinsindedir; sistem lineer ve sönümsüzdür.
  subroutine test_tiny_scale_two_inertia_frequency()
    real(dp) :: small_scale
    real(dp) :: expected_frequency_hz
    type(two_inertia_tvd_system_t) :: system
    type(two_inertia_modal_result_t) :: modal_result

    small_scale = 0.25_dp * tiny(1.0_dp)
    if (small_scale <= 0.0_dp .or. .not. ieee_is_finite(small_scale)) then
      error stop "IEEE küçük ölçek test değeri pozitif ve sonlu üretilemedi."
    end if

    system = two_inertia_tvd_system_t( &
      hub_polar_inertia_kg_m2=small_scale, &
      ring_polar_inertia_kg_m2=small_scale, &
      storage_stiffness_nm_per_rad=small_scale)
    modal_result = solve_free_free_two_inertia_modes(system)
    expected_frequency_hz = sqrt(2.0_dp) / (2.0_dp * pi)

    call assert_relative_close( &
      modal_result%elastic_frequency_hz, expected_frequency_hz, &
      relative_tolerance, &
      "Küçük ölçek iki-atalet frekansı harmonik referansla uyuşmuyor.")
  end subroutine test_tiny_scale_two_inertia_frequency

  !> NaN ve artı/eksi sonsuz değerleri fiziksel büyüklük kategorilerine
  !! yönlendirir; ayrıca geometri arayüzü, sıralama ve temsil aralığı aşımı
  !! vakalarını doğrudan üretim API'lerinde çalıştırır. Girdi selector
  !! boyutsuzdur. Geçersiz değer kabul edilirse normal çıkış, CTest WILL_FAIL
  !! vakasını başarısız yapar.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    character(len=64) :: quantity_name
    logical :: is_nonfinite_case
    real(dp) :: invalid_value

    call parse_nonfinite_case( &
      case_name, is_nonfinite_case, quantity_name, invalid_value)
    if (is_nonfinite_case) then
      call exercise_nonfinite_quantity(trim(quantity_name), invalid_value)
      return
    end if

    select case (case_name)
      case ("hub_rubber_interface_mismatch")
        call exercise_hub_rubber_interface_mismatch()
      case ("rubber_ring_interface_mismatch")
        call exercise_rubber_ring_interface_mismatch()
      case ("zero_length")
        call exercise_zero_length()
      case ("reversed_annular_geometry")
        call exercise_reversed_annular_geometry()
      case ("output_overflow_geometry")
        call exercise_geometry_output_overflow()
      case ("output_overflow_inertia")
        call exercise_inertia_output_overflow()
      case ("output_overflow_stiffness")
        call exercise_stiffness_output_overflow()
      case ("output_overflow_loss_factor")
        call exercise_loss_factor_output_overflow()
      case default
        print *, "Bilinmeyen sayısal sağlamlaştırma selector'ı: ", case_name
        return
    end select
  end subroutine exercise_invalid_case

  !> Selector önekinden IEEE özel değerini ve fiziksel büyüklük adını ayırır.
  !! NaN, +Inf ve -Inf değerleri ieee_value ile çalışma zamanında üretilir.
  subroutine parse_nonfinite_case( &
      case_name, is_nonfinite_case, quantity_name, invalid_value)
    character(len=*), intent(in) :: case_name
    logical, intent(out) :: is_nonfinite_case
    character(len=*), intent(out) :: quantity_name
    real(dp), intent(out) :: invalid_value

    quantity_name = ""
    invalid_value = 0.0_dp
    is_nonfinite_case = .true.

    if (index(case_name, "nan_") == 1) then
      invalid_value = ieee_value(0.0_dp, ieee_quiet_nan)
      quantity_name = case_name(5:)
    else if (index(case_name, "positive_infinity_") == 1) then
      invalid_value = ieee_value(0.0_dp, ieee_positive_inf)
      quantity_name = case_name(19:)
    else if (index(case_name, "negative_infinity_") == 1) then
      invalid_value = ieee_value(0.0_dp, ieee_negative_inf)
      quantity_name = case_name(19:)
    else
      is_nonfinite_case = .false.
    end if
  end subroutine parse_nonfinite_case

  !> Bir IEEE özel değerini seçilen SI fiziksel alanına yerleştirip ilgili
  !! üretim API'sini çağırır. Her kategori için diğer alanlar geçerli bırakılır.
  subroutine exercise_nonfinite_quantity(quantity_name, invalid_value)
    character(len=*), intent(in) :: quantity_name
    real(dp), intent(in) :: invalid_value

    type(inertia_ring_geometry_t) :: ring
    type(rubber_geometry_t) :: rubber
    type(dynamic_rubber_material_t) :: material
    type(complex_torsional_stiffness_t) :: dynamic_stiffness
    type(two_inertia_tvd_system_t) :: modal_system
    type(two_inertia_modal_result_t) :: modal_result
    real(dp) :: scalar_result
    real(dp) :: mass_kg
    real(dp) :: polar_inertia_kg_m2

    ring = inertia_ring_geometry_t( &
      inner_radius_m=0.04_dp, outer_radius_m=0.06_dp, &
      axial_length_m=0.02_dp)
    rubber = rubber_geometry_t( &
      inner_radius_m=0.03_dp, outer_radius_m=0.04_dp, &
      axial_length_m=0.02_dp)
    call initialize_valid_material(material)

    select case (quantity_name)
      case ("polar_inertia")
        modal_system = two_inertia_tvd_system_t( &
          hub_polar_inertia_kg_m2=invalid_value, &
          ring_polar_inertia_kg_m2=0.20_dp, &
          storage_stiffness_nm_per_rad=1000.0_dp)
        modal_result = solve_free_free_two_inertia_modes(modal_system)
        print *, "Geçersiz polar atalet kabul edildi: ", &
          modal_result%elastic_frequency_hz
      case ("density")
        call calculate_annular_ring_properties( &
          ring, invalid_value, mass_kg, polar_inertia_kg_m2)
        print *, "Geçersiz yoğunluk kabul edildi: ", mass_kg
      case ("radius")
        scalar_result = calculate_rubber_polar_area_moment( &
          invalid_value, 0.01_dp)
        print *, "Geçersiz yarıçap kabul edildi: ", scalar_result
      case ("length")
        scalar_result = calculate_annular_bush_torsion_geometry_factor( &
          0.03_dp, 0.04_dp, invalid_value)
        print *, "Geçersiz uzunluk kabul edildi: ", scalar_result
      case ("stiffness")
        scalar_result = calculate_natural_frequency(invalid_value, 0.20_dp)
        print *, "Geçersiz rijitlik kabul edildi: ", scalar_result
      case ("modulus")
        material%storage_shear_modulus_pa = invalid_value
        dynamic_stiffness = calculate_dynamic_torsional_stiffness( &
          material, rubber)
        print *, "Geçersiz modül kabul edildi: ", &
          dynamic_stiffness%storage_stiffness
      case ("frequency")
        material%frequency_hz = invalid_value
        dynamic_stiffness = calculate_dynamic_torsional_stiffness( &
          material, rubber)
        print *, "Geçersiz frekans kabul edildi: ", dynamic_stiffness%frequency
      case ("temperature")
        material%temperature_k = invalid_value
        dynamic_stiffness = calculate_dynamic_torsional_stiffness( &
          material, rubber)
        print *, "Geçersiz sıcaklık kabul edildi: ", &
          dynamic_stiffness%temperature
      case default
        print *, "Bilinmeyen sonlu olmayan büyüklük selector'ı: ", &
          quantity_name
    end select
  end subroutine exercise_nonfinite_quantity

  !> Göbek dış yarıçapı ile elastomer iç yarıçapı tolerans dışındayken
  !! birleşik geometri validator'ının modeli reddetmesini bekler.
  subroutine exercise_hub_rubber_interface_mismatch()
    type(tvd_geometry_t) :: geometry

    call initialize_valid_tvd_geometry(geometry)
    geometry%rubber%inner_radius_m = 0.031_dp
    call validate_tvd_geometry(geometry)
    print *, "Göbek-elastomer arayüz uyumsuzluğu kabul edildi."
  end subroutine exercise_hub_rubber_interface_mismatch

  !> Elastomer dış yarıçapı ile halka iç yarıçapı tolerans dışındayken
  !! birleşik geometri validator'ının modeli reddetmesini bekler.
  subroutine exercise_rubber_ring_interface_mismatch()
    type(tvd_geometry_t) :: geometry

    call initialize_valid_tvd_geometry(geometry)
    geometry%inertia_ring%inner_radius_m = 0.041_dp
    call validate_tvd_geometry(geometry)
    print *, "Elastomer-halka arayüz uyumsuzluğu kabul edildi."
  end subroutine exercise_rubber_ring_interface_mismatch

  !> Sıfır eksenel uzunluğun hacimsiz TVD bileşeni üretemeyeceğini sınar.
  subroutine exercise_zero_length()
    type(tvd_geometry_t) :: geometry

    call initialize_valid_tvd_geometry(geometry)
    geometry%rubber%axial_length_m = 0.0_dp
    call validate_tvd_geometry(geometry)
    print *, "Sıfır elastomer eksenel uzunluğu kabul edildi."
  end subroutine exercise_zero_length

  !> İç yarıçapı dış yarıçaptan büyük annüler halkanın reddini sınar.
  subroutine exercise_reversed_annular_geometry()
    type(tvd_geometry_t) :: geometry

    call initialize_valid_tvd_geometry(geometry)
    geometry%inertia_ring%inner_radius_m = 0.07_dp
    geometry%inertia_ring%outer_radius_m = 0.06_dp
    call validate_tvd_geometry(geometry)
    print *, "Ters sıralı annüler geometri kabul edildi."
  end subroutine exercise_reversed_annular_geometry

  !> Sonlu geometri girdilerinin temsil edilemeyen m^3 çıktı üretmesi halinde
  !! üretim geometri yordamının sessiz Inf döndürmemesini sınar.
  subroutine exercise_geometry_output_overflow()
    real(dp) :: geometry_factor_m3

    geometry_factor_m3 = &
      calculate_annular_bush_torsion_geometry_factor( &
        1.0e150_dp, 2.0e150_dp, 1.0e150_dp)
    print *, "Taşan geometri faktörü kabul edildi: ", geometry_factor_m3
  end subroutine exercise_geometry_output_overflow

  !> Sonlu annüler geometri ve yoğunluğun temsil edilemeyen kütle/atalet
  !! üretmesi halinde üretim inertia yordamının sessiz Inf döndürmemesini sınar.
  subroutine exercise_inertia_output_overflow()
    type(inertia_ring_geometry_t) :: ring
    real(dp) :: mass_kg
    real(dp) :: polar_inertia_kg_m2

    ring = inertia_ring_geometry_t( &
      inner_radius_m=1.0e100_dp, outer_radius_m=2.0e100_dp, &
      axial_length_m=1.0e100_dp)
    call calculate_annular_ring_properties( &
      ring, 1.0e100_dp, mass_kg, polar_inertia_kg_m2)
    print *, "Taşan atalet çıktısı kabul edildi: ", &
      mass_kg, polar_inertia_kg_m2
  end subroutine exercise_inertia_output_overflow

  !> Sonlu G' ve geometrinin temsil edilemeyen N*m/rad çıktı üretmesi halinde
  !! statik rijitlik yordamının sessiz Inf döndürmemesini sınar.
  subroutine exercise_stiffness_output_overflow()
    type(rubber_geometry_t) :: rubber
    type(dynamic_rubber_material_t) :: material
    real(dp) :: stiffness_nm_per_rad

    rubber = rubber_geometry_t( &
      inner_radius_m=1.0e3_dp, outer_radius_m=2.0e3_dp, &
      axial_length_m=1.0e3_dp)
    call initialize_valid_material(material)
    material%storage_shear_modulus_pa = 1.0e300_dp
    stiffness_nm_per_rad = calculate_torsional_stiffness(rubber, material)
    print *, "Taşan torsional rijitlik kabul edildi: ", &
      stiffness_nm_per_rad
  end subroutine exercise_stiffness_output_overflow

  !> Sonlu G' ve G'' girdilerinin temsil edilemeyen G''/G' oranı üretmesi
  !! halinde kayıp faktörü yordamının sessiz Inf döndürmemesini sınar.
  subroutine exercise_loss_factor_output_overflow()
    type(dynamic_shear_modulus) :: modulus
    real(dp) :: loss_factor

    modulus = dynamic_shear_modulus( &
      storage_modulus=1.0e-300_dp, loss_modulus=1.0e300_dp, &
      frequency=100.0_dp, temperature=293.15_dp)
    loss_factor = calculate_loss_factor(modulus)
    print *, "Taşan kayıp faktörü kabul edildi: ", loss_factor
  end subroutine exercise_loss_factor_output_overflow

  !> Testlerde kullanılan fiziksel ve ara yüzleri tam eşleşen TVD geometrisini
  !! metre (m) cinsinden kurar.
  subroutine initialize_valid_tvd_geometry(geometry)
    type(tvd_geometry_t), intent(out) :: geometry

    geometry%hub = hub_geometry_t( &
      bore_radius_m=0.01_dp, outer_radius_m=0.03_dp, &
      axial_length_m=0.02_dp)
    geometry%rubber = rubber_geometry_t( &
      inner_radius_m=0.03_dp, outer_radius_m=0.04_dp, &
      axial_length_m=0.02_dp)
    geometry%inertia_ring = inertia_ring_geometry_t( &
      inner_radius_m=0.04_dp, outer_radius_m=0.06_dp, &
      axial_length_m=0.02_dp)
  end subroutine initialize_valid_tvd_geometry

  !> Testlerin geçerli EPDM benzeri malzeme çalışma noktasını SI birimleriyle
  !! kurar: yoğunluk kg/m^3, modüller Pa, frekans Hz ve sıcaklık K.
  subroutine initialize_valid_material(material)
    type(dynamic_rubber_material_t), intent(out) :: material

    material%name = "EPDM sayısal sağlamlaştırma örneği"
    material%density_kg_m3 = 1100.0_dp
    material%storage_shear_modulus_pa = 1.0e6_dp
    material%loss_shear_modulus_pa = 0.1e6_dp
    material%frequency_hz = 100.0_dp
    material%temperature_k = 293.15_dp
  end subroutine initialize_valid_material

end program test_numerical_hardening
