module tms_tts_identification
  use tms_tts_types, only : tts_material_family_t, &
    tts_pair_shift_configuration_t, tts_shift_chain_result_t, &
    tts_identification_result_t, tts_validation_result_t, &
    TTS_IDENTIFICATION_SUCCESS, TTS_IDENTIFICATION_REFERENCE_NOT_FOUND, &
    TTS_IDENTIFICATION_CHAIN_BROKEN, &
    TTS_IDENTIFICATION_MASTER_CURVE_FAILED, &
    TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED, PAIR_SHIFT_SUCCESS, &
    validate_tts_material_family
  use tms_tts_shift_chain, only : build_tts_shift_chain
  use tms_tts_master_curve, only : build_tts_master_experimental_cloud, &
    stitch_tts_runtime_master_table
  use tms_tts_diagnostics, only : build_tts_vgp_cloud, &
    build_tts_cole_cole_cloud
  implicit none
  private

  public :: identify_tts_master_curve

contains

  !> Measured dynamic-shear isotherm family'den explicit measured reference
  !! kullanarak V0.8.1 identification sonucunu üretir. Akış validation,
  !! adjacent-pair shift chain, provenance-preserving experimental cloud,
  !! reference-centered runtime stitching ve VGP/Cole-Cole evidence'tır.
  !! Girdiler SI birimlerindedir; input family deep-copy edilir. Failure'lar
  !! programmer error değildir ve status ile döner. SUCCESS yalnız numerical
  !! construction başarısıdır; material için universal TRS PASS değildir.
  pure function identify_tts_master_curve( &
      family, reference_isotherm_identifier, configuration) result(result)
    type(tts_material_family_t), intent(in) :: family
    character(len=*), intent(in) :: reference_isotherm_identifier
    type(tts_pair_shift_configuration_t), intent(in), optional :: configuration
    type(tts_identification_result_t) :: result

    type(tts_shift_chain_result_t) :: chain
    type(tts_validation_result_t) :: validation
    logical :: master_success
    logical :: runtime_success

    result%source_family = family
    validation = validate_tts_material_family(family)
    if (.not. validation%valid) then
      result%status = validation%status
      result%message = validation%message
      return
    end if

    if (present(configuration)) then
      chain = build_tts_shift_chain(family, reference_isotherm_identifier, &
        configuration)
    else
      chain = build_tts_shift_chain(family, reference_isotherm_identifier)
    end if
    result%reference_isotherm_index = chain%reference_isotherm_index
    if (chain%reference_isotherm_index > 0) then
      result%reference_isotherm_identifier = family%isotherms( &
        chain%reference_isotherm_index)%isotherm_identifier
      result%reference_temperature_k = family%isotherms( &
        chain%reference_isotherm_index)%temperature_k
    end if
    if (allocated(chain%pair_shift_results)) then
      result%pair_shift_results = chain%pair_shift_results
    end if
    if (.not. chain%available) then
      result%status = chain%status
      if (chain%status == TTS_IDENTIFICATION_REFERENCE_NOT_FOUND) then
        result%message = "Explicit measured reference isotherm bulunamadı."
      else if (chain%status == TTS_IDENTIFICATION_CHAIN_BROKEN) then
        result%message = "Zorunlu adjacent pair çözülemedi; shift chain kırıldı."
      else
        result%message = "Experimental shift chain oluşturulamadı."
      end if
      return
    end if
    result%shift_chain_available = .true.
    result%empirical_shifts = chain%empirical_shifts

    call build_tts_master_experimental_cloud(family, &
      result%empirical_shifts, result%master_cloud, master_success)
    if (.not. master_success) then
      result%status = TTS_IDENTIFICATION_MASTER_CURVE_FAILED
      result%message = "Experimental master cloud oluşturulamadı."
      return
    end if
    result%master_cloud_available = .true.
    result%diagnostics%vgp_points = build_tts_vgp_cloud(family)
    result%diagnostics%cole_cole_points = build_tts_cole_cole_cloud(family)

    call stitch_tts_runtime_master_table(family, &
      result%reference_isotherm_index, result%empirical_shifts, &
      result%master_cloud, result%runtime_master_table, &
      result%diagnostics%boundaries, runtime_success)
    if (.not. runtime_success) then
      result%status = TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED
      result%message = "Solver-ready runtime master table oluşturulamadı."
      return
    end if

    result%diagnostics%full_complex_pair_support = &
      all_pair_shifts_have_joint_support(result)
    result%runtime_export_ready = .true.
    result%status = TTS_IDENTIFICATION_SUCCESS
    result%message = &
      "Experimental master curve ve runtime export başarıyla oluşturuldu."
  end function identify_tts_master_curve

  pure function all_pair_shifts_have_joint_support(result) result(supported)
    type(tts_identification_result_t), intent(in) :: result
    logical :: supported
    integer :: i

    supported = allocated(result%pair_shift_results)
    if (.not. supported) return
    supported = size(result%pair_shift_results) > 0
    if (.not. supported) return
    do i = 1, size(result%pair_shift_results)
      if (result%pair_shift_results(i)%status /= PAIR_SHIFT_SUCCESS .or. &
          .not. result%pair_shift_results(i)%joint_shift_available) then
        supported = .false.
        return
      end if
    end do
  end function all_pair_shifts_have_joint_support

end module tms_tts_identification
