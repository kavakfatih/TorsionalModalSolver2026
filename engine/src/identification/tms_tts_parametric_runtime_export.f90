module tms_tts_parametric_runtime_export
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_identification_result_t, tts_text_length
  use tms_tts_shift_law_types, only : tts_arrhenius_fit_result_t, &
    tts_wlf_fit_result_t, SHIFT_LAW_FIT_SUCCESS
  use tms_tts_runtime_export, only : tts_runtime_export_t, &
    create_tts_runtime_export
  use tms_tabulated_dynamic_modulus_provider, only : &
    tabulated_dynamic_modulus_provider_t
  use tms_arrhenius_temperature_shift, only : &
    arrhenius_temperature_shift_provider_t, &
    create_arrhenius_temperature_shift_provider
  use tms_wlf_temperature_shift, only : wlf_temperature_shift_provider_t, &
    create_wlf_temperature_shift_provider
  use tms_thermorheological_dynamic_modulus_provider, only : &
    thermorheological_dynamic_modulus_provider_t, &
    create_thermorheological_dynamic_modulus_provider
  implicit none
  private

  !> V0.8.1 authoritative master table ile identified Arrhenius shift
  !! provider'ını mevcut V0.8.0 thermorheological runtime API'sinde taşır.
  !! predictive_validation_available metadata'dır; automatic model choice
  !! veya experimental empirical table replacement anlamına gelmez.
  type, public :: tts_arrhenius_runtime_export_t
    logical :: available = .false.
    logical :: predictive_validation_available = .false.
    character(len=tts_text_length) :: message = ""
    type(tabulated_dynamic_modulus_provider_t) :: master_curve_provider
    type(arrhenius_temperature_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: &
      thermorheological_provider
  end type tts_arrhenius_runtime_export_t

  !> V0.8.1 authoritative master table ile identifiable WLF shift provider'ını
  !! mevcut V0.8.0 runtime API'sinde taşır. parameter_identifiable=false olan
  !! large-C2 limitleri güvenli runtime modeli olarak export edilmez.
  type, public :: tts_wlf_runtime_export_t
    logical :: available = .false.
    logical :: predictive_validation_available = .false.
    character(len=tts_text_length) :: message = ""
    type(tabulated_dynamic_modulus_provider_t) :: master_curve_provider
    type(wlf_temperature_shift_provider_t) :: shift_provider
    type(thermorheological_dynamic_modulus_provider_t) :: &
      thermorheological_provider
  end type tts_wlf_runtime_export_t

  public :: create_tts_arrhenius_runtime_export
  public :: create_tts_wlf_runtime_export

contains

  !> Başarılı adjacent-pair Ea_app fit'ini mevcut Arrhenius runtime provider'a
  !! bağlar. Master G'/G'' [Pa]-f_r [Hz] tablosu V0.8.1 empirical export'tan
  !! aynen gelir; yalnız T [K] -> s=log10(a_T) provider'ı parametrik olur.
  !! Domain measured [Tmin,Tmax] ile sınırlıdır ve extrapolation yapılmaz.
  function create_tts_arrhenius_runtime_export( &
      identification, fit) result(runtime_export)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_arrhenius_fit_result_t), intent(in) :: fit
    type(tts_arrhenius_runtime_export_t) :: runtime_export

    type(tts_runtime_export_t) :: empirical_export

    if (.not. valid_arrhenius_export_request(identification, fit)) then
      runtime_export%message = &
        "Arrhenius fit/runtime reference veya calibrated domain geçersiz."
      return
    end if
    empirical_export = create_tts_runtime_export(identification)
    if (.not. empirical_export%available) then
      runtime_export%message = &
        "Authoritative V0.8.1 master runtime export hazır değil."
      return
    end if

    runtime_export%master_curve_provider = &
      empirical_export%master_curve_provider
    runtime_export%shift_provider = &
      create_arrhenius_temperature_shift_provider( &
        fit%apparent_activation_energy_j_per_mol, &
        fit%reference_temperature_k, fit%minimum_temperature_k, &
        fit%maximum_temperature_k)
    runtime_export%thermorheological_provider = &
      create_thermorheological_dynamic_modulus_provider( &
        runtime_export%master_curve_provider, runtime_export%shift_provider)
    runtime_export%predictive_validation_available = &
      fit%predictive_validation_available
    runtime_export%available = .true.
    runtime_export%message = &
      "Measured-domain Arrhenius parametric runtime export hazır."
  end function create_tts_arrhenius_runtime_export

  !> Identifiable profiled WLF fit'ini mevcut WLF runtime provider'a bağlar.
  !! Master table V0.8.1 empirical representation'dır; C1/C2 yalnız horizontal
  !! temperature shifting yapar. Calibrated [Tmin,Tmax] [K] dışında
  !! extrapolation ve poorly-identified WLF export kesinlikle yapılmaz.
  function create_tts_wlf_runtime_export( &
      identification, fit) result(runtime_export)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_wlf_fit_result_t), intent(in) :: fit
    type(tts_wlf_runtime_export_t) :: runtime_export

    type(tts_runtime_export_t) :: empirical_export

    if (.not. valid_wlf_export_request(identification, fit)) then
      runtime_export%message = &
        "WLF fit identifiable değil veya calibrated pole-safe domain geçersiz."
      return
    end if
    empirical_export = create_tts_runtime_export(identification)
    if (.not. empirical_export%available) then
      runtime_export%message = &
        "Authoritative V0.8.1 master runtime export hazır değil."
      return
    end if

    runtime_export%master_curve_provider = &
      empirical_export%master_curve_provider
    runtime_export%shift_provider = create_wlf_temperature_shift_provider( &
      fit%c1, fit%c2_k, fit%reference_temperature_k, &
      fit%minimum_temperature_k, fit%maximum_temperature_k)
    runtime_export%thermorheological_provider = &
      create_thermorheological_dynamic_modulus_provider( &
        runtime_export%master_curve_provider, runtime_export%shift_provider)
    runtime_export%predictive_validation_available = &
      fit%predictive_validation_available
    runtime_export%available = .true.
    runtime_export%message = &
      "Measured-domain identifiable WLF parametric runtime export hazır."
  end function create_tts_wlf_runtime_export

  pure function valid_arrhenius_export_request( &
      identification, fit) result(valid)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_arrhenius_fit_result_t), intent(in) :: fit
    logical :: valid

    valid = identification%runtime_export_ready .and. &
      fit%status == SHIFT_LAW_FIT_SUCCESS .and. fit%fit_available .and. &
      ieee_is_finite(fit%apparent_activation_energy_j_per_mol) .and. &
      fit%apparent_activation_energy_j_per_mol > 0.0_dp .and. &
      valid_common_export_domain(identification, &
        fit%reference_temperature_k, fit%minimum_temperature_k, &
        fit%maximum_temperature_k)
  end function valid_arrhenius_export_request

  pure function valid_wlf_export_request(identification, fit) result(valid)
    type(tts_identification_result_t), intent(in) :: identification
    type(tts_wlf_fit_result_t), intent(in) :: fit
    logical :: valid
    real(dp) :: minimum_denominator_k
    real(dp) :: pole_scale_k

    valid = identification%runtime_export_ready .and. &
      fit%status == SHIFT_LAW_FIT_SUCCESS .and. fit%fit_available .and. &
      fit%parameter_identifiable .and. ieee_is_finite(fit%c1) .and. &
      fit%c1 > 0.0_dp .and. ieee_is_finite(fit%c2_k) .and. &
      fit%c2_k > 0.0_dp .and. ieee_is_finite(fit%p_c1_over_c2_per_k) .and. &
      fit%p_c1_over_c2_per_k > 0.0_dp .and. &
      ieee_is_finite(fit%q_inverse_c2_per_k) .and. &
      fit%q_inverse_c2_per_k > 0.0_dp .and. &
      valid_common_export_domain(identification, &
        fit%reference_temperature_k, fit%minimum_temperature_k, &
        fit%maximum_temperature_k)
    if (.not. valid) return

    minimum_denominator_k = fit%c2_k + fit%minimum_temperature_k - &
      fit%reference_temperature_k
    pole_scale_k = max(1.0_dp, abs(fit%c2_k), &
      abs(fit%minimum_temperature_k-fit%reference_temperature_k))
    valid = ieee_is_finite(minimum_denominator_k) .and. &
      minimum_denominator_k > 64.0_dp*epsilon(1.0_dp)*pole_scale_k
  end function valid_wlf_export_request

  pure function valid_common_export_domain( &
      identification, reference_temperature_k, minimum_temperature_k, &
      maximum_temperature_k) result(valid)
    type(tts_identification_result_t), intent(in) :: identification
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: minimum_temperature_k
    real(dp), intent(in) :: maximum_temperature_k
    logical :: valid

    valid = ieee_is_finite(reference_temperature_k) .and. &
      ieee_is_finite(minimum_temperature_k) .and. &
      ieee_is_finite(maximum_temperature_k) .and. &
      minimum_temperature_k > 0.0_dp .and. &
      maximum_temperature_k > minimum_temperature_k .and. &
      (reference_temperature_k > minimum_temperature_k .or. &
        temperatures_are_machine_equivalent(reference_temperature_k, &
          minimum_temperature_k)) .and. &
      (reference_temperature_k < maximum_temperature_k .or. &
        temperatures_are_machine_equivalent(reference_temperature_k, &
          maximum_temperature_k)) .and. &
      temperatures_are_machine_equivalent(reference_temperature_k, &
        identification%reference_temperature_k)
  end function valid_common_export_domain

  pure elemental function temperatures_are_machine_equivalent(a, b) &
      result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    equivalent = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
      abs(a-b) <= 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(a), abs(b))
  end function temperatures_are_machine_equivalent

end module tms_tts_parametric_runtime_export
