module tms_geometry
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : pi
  implicit none
  private

  !> TVD bileşenlerinin eşleşen yarıçapları için mutlak tolerans [m].
  real(dp), parameter :: interface_absolute_tolerance_m = 1.0e-12_dp

  !> TVD bileşenlerinin eşleşen yarıçapları için boyutsuz bağıl tolerans.
  real(dp), parameter :: interface_relative_tolerance = 1.0e-9_dp

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
  !! Bu tür yalnızca veri taşır; geometrik uygunluk validate_tvd_geometry
  !! yordamıyla ayrıca doğrulanır.
  type, public :: tvd_geometry_t
    type(hub_geometry_t) :: hub
    type(rubber_geometry_t) :: rubber
    type(inertia_ring_geometry_t) :: inertia_ring
  end type tvd_geometry_t

  public :: calculate_rubber_polar_area_moment
  public :: calculate_annular_bush_torsion_geometry_factor
  public :: validate_tvd_geometry

contains

  !> Annüler elastomer kesitin polar geometrik alan momentini hesaplar.
  !!
  !! Fiziksel açıklama: Polar alan momenti, kesit geometrisinin burulmaya
  !! karşı rijitlik katkısını temsil eder; polar kütle ataleti değildir.
  !! Matematiksel açıklama: Jp = pi/2 * (ro^4 - ri^4).
  !! Girdiler: outer_radius dış yarıçapı, inner_radius iç yarıçapı temsil eder;
  !! iki değer de metre (m) cinsindendir.
  !! Çıktı: jp polar alan momentidir ve metrenin dördüncü kuvveti (m^4)
  !! cinsindendir.
  !! Varsayımlar ve geçerlilik: Kesit eksenel simetrik ve annülerdir;
  !! 0 <= ri < ro olmalı; girdiler ve çıktı sonlu olmalıdır. Negatif,
  !! sonlu olmayan veya yanlış sıralanmış yarıçaplarda ve sonlu-pozitif bir
  !! çıktı üretilemediğinde yordam error stop ile sonlanır.
  !! Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure function calculate_rubber_polar_area_moment(outer_radius, &
      inner_radius) result(jp)
    real(dp), intent(in) :: outer_radius
    real(dp), intent(in) :: inner_radius
    real(dp) :: jp

    real(dp) :: radius_square_difference_m2
    real(dp) :: radius_square_sum_m2

    if (.not. ieee_is_finite(outer_radius) .or. &
        .not. ieee_is_finite(inner_radius)) then
      error stop "Elastomer yarıçapları sonlu olmalıdır."
    end if

    if (outer_radius < 0.0_dp .or. inner_radius < 0.0_dp) then
      error stop "Elastomer yarıçapları negatif olamaz."
    end if

    if (outer_radius <= inner_radius) then
      error stop "Elastomer dış yarıçapı iç yarıçapından büyük olmalıdır."
    end if

    radius_square_difference_m2 = &
      (outer_radius - inner_radius) * (outer_radius + inner_radius)
    radius_square_sum_m2 = outer_radius**2 + inner_radius**2

    if (.not. ieee_is_finite(radius_square_difference_m2) .or. &
        radius_square_difference_m2 <= 0.0_dp .or. &
        .not. ieee_is_finite(radius_square_sum_m2) .or. &
        radius_square_sum_m2 <= 0.0_dp) then
      error stop "Polar alan momenti ara değerleri sonlu ve pozitif olmalıdır."
    end if

    jp = 0.5_dp * pi * radius_square_difference_m2 * radius_square_sum_m2

    if (.not. ieee_is_finite(jp) .or. jp <= 0.0_dp) then
      error stop "Polar alan momenti sonlu ve pozitif olmalıdır."
    end if
  end function calculate_rubber_polar_area_moment

  !> Tam bağlı annüler kauçuk burcun burulma geometri faktörünü hesaplar.
  !!
  !! Fiziksel açıklama: Faktör, rijit iç göbek ile rijit dış halka
  !! arasındaki eş merkezli elastomer tabakanın bağıl dönmeye karşı geometri
  !! katkısıdır. Elastomer iki silindirik yüzeye tam bağlıdır.
  !! Matematiksel açıklama:
  !! C_theta = 4*pi*L*ri^2*ro^2/(ro^2-ri^2).
  !! Girdiler: inner_radius iç yarıçap ri, outer_radius dış yarıçap ro
  !! ve axial_length bağlı eksenel genişlik L'dir; tümü metre (m)
  !! cinsindendir.
  !! Çıktı: geometry_factor_m3, metreküp (m^3) cinsindedir ve kayma
  !! modülüyle çarpıldığında N*m/rad birimli burulma rijitliği verir.
  !! Varsayımlar ve geçerlilik: Silindirler eş merkezli ve rijit, elastomer
  !! homojen, lineer, küçük deformasyon bölgesinde ve ara yüzlerde kaymasızdır.
  !! ri > 0, ro > ri ve L > 0 olmalı; girdiler ve çıktı sonlu olmalıdır.
  !! ri = 0, bağlı iç silindirik yüzey bulunmadığından bu burç modelinin
  !! kapsamı dışındadır.
  !! Ayrıntılar: docs/mathematics/torsional-physics-core.md.
  pure function calculate_annular_bush_torsion_geometry_factor( &
      inner_radius, outer_radius, axial_length) result(geometry_factor_m3)
    real(dp), intent(in) :: inner_radius
    real(dp), intent(in) :: outer_radius
    real(dp), intent(in) :: axial_length
    real(dp) :: geometry_factor_m3

    real(dp) :: radius_square_difference_m2

    if (.not. ieee_is_finite(inner_radius) .or. &
        .not. ieee_is_finite(outer_radius) .or. &
        .not. ieee_is_finite(axial_length)) then
      error stop "Elastomer geometri girdileri sonlu olmalıdır."
    end if

    if (inner_radius < 0.0_dp .or. outer_radius < 0.0_dp) then
      error stop "Elastomer yarıçapları negatif olamaz."
    end if

    if (inner_radius == 0.0_dp) then
      error stop "Annüler kauçuk burç iç yarıçapı pozitif olmalıdır."
    end if

    if (outer_radius <= inner_radius) then
      error stop "Elastomer dış yarıçapı iç yarıçapından büyük olmalıdır."
    end if

    if (axial_length <= 0.0_dp) then
      error stop "Elastomer eksenel genişliği sıfırdan büyük olmalıdır."
    end if

    radius_square_difference_m2 = &
      (outer_radius - inner_radius) * (outer_radius + inner_radius)

    if (.not. ieee_is_finite(radius_square_difference_m2) .or. &
        radius_square_difference_m2 <= 0.0_dp) then
      error stop "Yarıçap kareleri farkı sonlu ve pozitif olmalıdır."
    end if

    geometry_factor_m3 = 4.0_dp * pi * axial_length * inner_radius**2 * &
      outer_radius**2 / radius_square_difference_m2

    if (.not. ieee_is_finite(geometry_factor_m3) .or. &
        geometry_factor_m3 <= 0.0_dp) then
      error stop "Burulma geometri faktörü sonlu ve pozitif olmalıdır."
    end if
  end function calculate_annular_bush_torsion_geometry_factor

  !> Birleşik TVD geometrisinin bileşen ve ara yüz uygunluğunu doğrular.
  !!
  !! Fiziksel açıklama: Göbek, bağlı elastomer ve atalet halkası eş merkezli
  !! ve birbiriyle temas eden silindirik bölgeler olarak modellenir.
  !! Matematiksel açıklama: Her bileşenin yarıçap sırası ve pozitif eksenel
  !! uzunluğu doğrulanır. Temas yarıçapları |a-b| <= atol + rtol*max(|a|,|b|)
  !! koşuluyla karşılaştırılır; atol = 1e-12 m ve rtol = 1e-9 değerleridir.
  !! Girdi: geometry içindeki tüm yarıçap ve uzunluklar metre (m) cinsindedir.
  !! Çıktı: Değer döndürmez; geçersiz veya sonlu olmayan geometride error stop
  !! ile sonlanır.
  !! Varsayımlar ve geçerlilik: Göbek için 0 <= rb < rh ve Lh > 0;
  !! elastomer için 0 < ri < ro ve Le > 0; atalet halkası için
  !! 0 <= rr_i < rr_o ve Lr > 0 olmalıdır. Göbek dış yarıçapı elastomer iç
  !! yarıçapıyla, elastomer dış yarıçapı da halka iç yarıçapıyla eşleşmelidir.
  pure subroutine validate_tvd_geometry(geometry)
    type(tvd_geometry_t), intent(in) :: geometry

    call validate_hub_geometry(geometry%hub)
    call validate_rubber_geometry(geometry%rubber)
    call validate_inertia_ring_geometry(geometry%inertia_ring)

    if (.not. are_interface_radii_compatible( &
        geometry%hub%outer_radius_m, geometry%rubber%inner_radius_m)) then
      error stop "Göbek-elastomer ara yüz yarıçapları eşleşmelidir."
    end if

    if (.not. are_interface_radii_compatible( &
        geometry%rubber%outer_radius_m, &
        geometry%inertia_ring%inner_radius_m)) then
      error stop "Elastomer-halka ara yüz yarıçapları eşleşmelidir."
    end if
  end subroutine validate_tvd_geometry

  !> Göbek geometrisinin sonlu ve fiziksel olarak sıralı olduğunu doğrular.
  !! Girdiler metre (m) cinsindedir; delik yarıçapı sıfır olabilir, dış
  !! yarıçap ve eksenel uzunluk fiziksel hacim için pozitif olmalıdır.
  pure subroutine validate_hub_geometry(hub)
    type(hub_geometry_t), intent(in) :: hub

    if (.not. ieee_is_finite(hub%bore_radius_m) .or. &
        .not. ieee_is_finite(hub%outer_radius_m) .or. &
        .not. ieee_is_finite(hub%axial_length_m)) then
      error stop "Göbek geometri değerleri sonlu olmalıdır."
    end if

    if (hub%bore_radius_m < 0.0_dp) then
      error stop "Göbek delik yarıçapı negatif olamaz."
    end if

    if (hub%outer_radius_m <= hub%bore_radius_m) then
      error stop "Göbek dış yarıçapı delik yarıçapından büyük olmalıdır."
    end if

    if (hub%axial_length_m <= 0.0_dp) then
      error stop "Göbek eksenel uzunluğu pozitif olmalıdır."
    end if
  end subroutine validate_hub_geometry

  !> Elastomer geometrisinin sonlu ve bonded burç modeline uygun olduğunu
  !! doğrular. Tüm girdiler metre (m) cinsindedir; 0 < ri < ro ve L > 0
  !! koşulları aranır.
  pure subroutine validate_rubber_geometry(rubber)
    type(rubber_geometry_t), intent(in) :: rubber

    if (.not. ieee_is_finite(rubber%inner_radius_m) .or. &
        .not. ieee_is_finite(rubber%outer_radius_m) .or. &
        .not. ieee_is_finite(rubber%axial_length_m)) then
      error stop "Elastomer geometri değerleri sonlu olmalıdır."
    end if

    if (rubber%inner_radius_m <= 0.0_dp) then
      error stop "Elastomer iç yarıçapı pozitif olmalıdır."
    end if

    if (rubber%outer_radius_m <= rubber%inner_radius_m) then
      error stop "Elastomer dış yarıçapı iç yarıçapından büyük olmalıdır."
    end if

    if (rubber%axial_length_m <= 0.0_dp) then
      error stop "Elastomer eksenel uzunluğu pozitif olmalıdır."
    end if
  end subroutine validate_rubber_geometry

  !> Atalet halkası geometrisinin sonlu ve fiziksel olarak sıralı olduğunu
  !! doğrular. Tüm girdiler metre (m) cinsindedir; iç yarıçap sıfır olabilir,
  !! dış yarıçap iç yarıçaptan ve eksenel uzunluk sıfırdan büyük olmalıdır.
  pure subroutine validate_inertia_ring_geometry(ring)
    type(inertia_ring_geometry_t), intent(in) :: ring

    if (.not. ieee_is_finite(ring%inner_radius_m) .or. &
        .not. ieee_is_finite(ring%outer_radius_m) .or. &
        .not. ieee_is_finite(ring%axial_length_m)) then
      error stop "Atalet halkası geometri değerleri sonlu olmalıdır."
    end if

    if (ring%inner_radius_m < 0.0_dp) then
      error stop "Atalet halkası iç yarıçapı negatif olamaz."
    end if

    if (ring%outer_radius_m <= ring%inner_radius_m) then
      error stop "Atalet halkası dış yarıçapı iç yarıçapından büyük olmalıdır."
    end if

    if (ring%axial_length_m <= 0.0_dp) then
      error stop "Atalet halkası eksenel uzunluğu pozitif olmalıdır."
    end if
  end subroutine validate_inertia_ring_geometry

  !> İki sonlu ara yüz yarıçapının karma mutlak/bağıl tolerans içinde
  !! eşleşip eşleşmediğini belirler. Girdiler metre (m), çıktı mantıksaldır.
  pure function are_interface_radii_compatible(first_radius_m, &
      second_radius_m) result(are_compatible)
    real(dp), intent(in) :: first_radius_m
    real(dp), intent(in) :: second_radius_m
    logical :: are_compatible

    real(dp) :: tolerance_m

    if (.not. ieee_is_finite(first_radius_m) .or. &
        .not. ieee_is_finite(second_radius_m)) then
      are_compatible = .false.
      return
    end if

    tolerance_m = interface_absolute_tolerance_m + &
      interface_relative_tolerance * &
      max(abs(first_radius_m), abs(second_radius_m))
    are_compatible = abs(first_radius_m - second_radius_m) <= tolerance_m
  end function are_interface_radii_compatible

end module tms_geometry
