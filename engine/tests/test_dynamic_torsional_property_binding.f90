program test_dynamic_torsional_property_binding
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_geometry, only : rubber_geometry_t
  use tms_material_frequency, only : material_frequency_point
  use tms_dynamic_material_metadata, only : dynamic_material_metadata_t, &
    DYNAMIC_DEFORMATION_MODE_SHEAR, DYNAMIC_DEFORMATION_MODE_TENSILE
  use tms_dynamic_modulus_provider, only : LINEAR_FREQUENCY
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t, &
    create_tabulated_dynamic_modulus_provider
  use tms_dynamic_torsional_property_binding, only : &
    dynamic_torsional_property_binding_t, &
    dynamic_torsional_property_state_t, &
    create_dynamic_torsional_property_binding, &
    evaluate_dynamic_torsional_property, &
    get_dynamic_binding_geometry_factor, get_dynamic_binding_metadata
  implicit none

  real(dp), parameter :: tolerance = 5.0e-13_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_annular_mapping_and_passivity()
  call test_binding_input_immutability()
  print *, "V0.7 dynamic torsional property binding doğrulandı."

contains

  !> Known C_theta ve interpolated G'/G'' ile K*=C_theta*G* mapping'ini
  !! analitik doğrular. G'>0, G''>=0 sonucunda K'>0, K''>=0 ve
  !! K''/K'=G''/G'=tan(delta) passivity eşitliği korunmalıdır.
  subroutine test_annular_mapping_and_passivity()
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: binding
    type(dynamic_torsional_property_state_t) :: state
    type(dynamic_material_metadata_t) :: metadata
    type(rubber_geometry_t) :: rubber
    real(dp) :: expected_geometry_factor

    rubber = rubber_geometry_t( &
      inner_radius_m=0.02_dp, outer_radius_m=0.05_dp, &
      axial_length_m=0.01_dp)
    provider = make_provider(DYNAMIC_DEFORMATION_MODE_SHEAR)
    binding = create_dynamic_torsional_property_binding( &
      10, provider, rubber)
    state = evaluate_dynamic_torsional_property( &
      binding, 15.0_dp, 293.15_dp)

    expected_geometry_factor = &
      4.0_dp*pi*rubber%axial_length_m*rubber%inner_radius_m**2 * &
      rubber%outer_radius_m**2 / &
      (rubber%outer_radius_m**2-rubber%inner_radius_m**2)
    call assert_close(get_dynamic_binding_geometry_factor(binding), &
      expected_geometry_factor, tolerance, "Binding C_theta hatalı.")
    call assert_close(state%storage_modulus_pa, 1.5e6_dp, tolerance, &
      "Binding interpolated G' hatalı.")
    call assert_close(state%loss_modulus_pa, 0.15e6_dp, tolerance, &
      "Binding interpolated G'' hatalı.")
    call assert_close(state%storage_stiffness_nm_per_rad, &
      expected_geometry_factor*1.5e6_dp, tolerance, &
      "Binding mapped K' hatalı.")
    call assert_close(state%loss_stiffness_nm_per_rad, &
      expected_geometry_factor*0.15e6_dp, tolerance, &
      "Binding mapped K'' hatalı.")
    call assert_close(state%loss_factor, 0.1_dp, tolerance, &
      "Binding tan(delta) hatalı.")
    call assert_close( &
      state%loss_stiffness_nm_per_rad / &
        state%storage_stiffness_nm_per_rad, &
      state%loss_modulus_pa/state%storage_modulus_pa, tolerance, &
      "K''/K' ile G''/G' eşit değil.")
    if (state%storage_stiffness_nm_per_rad <= 0.0_dp .or. &
        state%loss_stiffness_nm_per_rad < 0.0_dp) then
      error stop "Binding passive K*/energy koşulunu korumadı."
    end if
    metadata = get_dynamic_binding_metadata(binding)
    if (trim(metadata%dataset_identifier) /= "BINDING-EPDM") then
      error stop "Binding dataset trace metadata'sını korumadı."
    end if
  end subroutine test_annular_mapping_and_passivity

  !> Constructor'ın geometry ve provider'ı independent değer olarak
  !! sakladığını doğrular. Caller rubber yapısını değiştirdikten sonra aynı
  !! binding C_theta ve K*(f) değerlerini üretmeye devam etmelidir.
  subroutine test_binding_input_immutability()
    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: binding
    type(dynamic_torsional_property_state_t) :: before_state
    type(dynamic_torsional_property_state_t) :: after_state
    type(rubber_geometry_t) :: rubber

    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    provider = make_provider(DYNAMIC_DEFORMATION_MODE_SHEAR)
    binding = create_dynamic_torsional_property_binding(1, provider, rubber)
    before_state = evaluate_dynamic_torsional_property( &
      binding, 15.0_dp, 293.15_dp)
    rubber%outer_radius_m = 0.5_dp
    after_state = evaluate_dynamic_torsional_property( &
      binding, 15.0_dp, 293.15_dp)
    call assert_close(after_state%storage_stiffness_nm_per_rad, &
      before_state%storage_stiffness_nm_per_rad, 0.0_dp, &
      "Caller geometry mutation binding'i değiştirdi.")
  end subroutine test_binding_input_immutability

  !> Invalid ID/geometri/provider ve direct-SHEAR contract hata yollarını
  !! ayrı CTest süreçlerinde tetikler. TENSILE veri için otomatik E->G yoktur.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(tabulated_dynamic_modulus_provider_t) :: provider
    type(dynamic_torsional_property_binding_t) :: binding
    type(rubber_geometry_t) :: rubber

    rubber = rubber_geometry_t(0.02_dp, 0.05_dp, 0.01_dp)
    select case (case_name)
    case ("nonpositive_element_id")
      provider = make_provider(DYNAMIC_DEFORMATION_MODE_SHEAR)
      binding = create_dynamic_torsional_property_binding(0, provider, rubber)
    case ("negative_radius")
      provider = make_provider(DYNAMIC_DEFORMATION_MODE_SHEAR)
      rubber%inner_radius_m = -0.02_dp
      binding = create_dynamic_torsional_property_binding(1, provider, rubber)
    case ("unordered_radii")
      provider = make_provider(DYNAMIC_DEFORMATION_MODE_SHEAR)
      rubber%inner_radius_m = rubber%outer_radius_m
      binding = create_dynamic_torsional_property_binding(1, provider, rubber)
    case ("zero_length")
      provider = make_provider(DYNAMIC_DEFORMATION_MODE_SHEAR)
      rubber%axial_length_m = 0.0_dp
      binding = create_dynamic_torsional_property_binding(1, provider, rubber)
    case ("uninitialized_provider")
      binding = create_dynamic_torsional_property_binding(1, provider, rubber)
    case ("tensile_mode")
      provider = make_provider(DYNAMIC_DEFORMATION_MODE_TENSILE)
      binding = create_dynamic_torsional_property_binding(1, provider, rubber)
    case default
      error stop "Bilinmeyen dynamic binding validation selector."
    end select
  end subroutine exercise_invalid_case

  pure function make_provider(mode) result(provider)
    integer, intent(in) :: mode
    type(tabulated_dynamic_modulus_provider_t) :: provider

    type(material_frequency_point) :: points(2)
    type(dynamic_material_metadata_t) :: metadata

    points(1)%frequency = 10.0_dp
    points(1)%temperature = 293.15_dp
    points(1)%storage_modulus = 1.0e6_dp
    points(1)%loss_modulus = 0.1e6_dp
    points(2)%frequency = 20.0_dp
    points(2)%temperature = 293.15_dp
    points(2)%storage_modulus = 2.0e6_dp
    points(2)%loss_modulus = 0.2e6_dp
    metadata%dataset_identifier = "BINDING-EPDM"
    metadata%material_identifier = "EPDM"
    metadata%dataset_temperature_k = 293.15_dp
    metadata%deformation_mode = mode
    provider = create_tabulated_dynamic_modulus_provider( &
      points, metadata, LINEAR_FREQUENCY)
  end function make_provider

  subroutine assert_close(actual, expected, relative_tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: relative_tolerance
    character(len=*), intent(in) :: message

    if (abs(actual - expected) > &
        relative_tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_close

end program test_dynamic_torsional_property_binding
