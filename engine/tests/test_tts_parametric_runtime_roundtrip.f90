program test_tts_parametric_runtime_roundtrip
  use tms_kinds, only : dp
  use tms_constants, only : universal_gas_constant_j_per_mol_k
  use tms_tts_types, only : tts_material_family_t, &
    tts_identification_result_t, TTS_IDENTIFICATION_SUCCESS
  use tms_tts_shift_law_types, only : &
    tts_shift_law_identification_result_t, tts_arrhenius_fit_result_t, &
    tts_wlf_fit_result_t, SHIFT_LAW_FIT_SUCCESS
  use tms_tts_identification, only : identify_tts_master_curve
  use tms_tts_shift_law_validation, only : fit_tts_shift_laws
  use tms_tts_parametric_runtime_export, only : &
    tts_arrhenius_runtime_export_t, tts_wlf_runtime_export_t, &
    create_tts_arrhenius_runtime_export, create_tts_wlf_runtime_export
  use tms_temperature_shift_types, only : temperature_shift_domain_t, &
    ARRHENIUS_TEMPERATURE_SHIFT, WLF_TEMPERATURE_SHIFT
  use tms_temperature_shift_provider, only : get_temperature_shift_domain
  use tms_dynamic_modulus_provider, only : dynamic_modulus_evaluation_t, &
    evaluate_dynamic_shear_modulus
  use tms_tts_test_support, only : make_generalized_maxwell_trs_family, &
    assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures_k(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: reference_temperature_k = 293.15_dp
  real(dp), parameter :: apparent_activation_energy_j_per_mol = 12000.0_dp
  real(dp), parameter :: known_wlf_c1 = 2.0_dp
  real(dp), parameter :: known_wlf_c2_k = 120.0_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_outside_domain(trim(validation_case))
    stop 0
  end if

  call test_arrhenius_runtime_roundtrip()
  call test_wlf_runtime_roundtrip()
  call test_invalid_fit_exports()

  print *, "V0.8.2 parametric runtime provider round-trip doğrulandı."

contains

  !> Exact Arrhenius horizontal shifts ile üretilen measured family'yi önce
  !! V0.8.1 empirical master'a, sonra adjacent-pair Ea_app fit'ine taşır.
  !! Off-reference physical f [Hz], fitted shift ile aynı f_r master noktasına
  !! gitmeli; G'/G'' [Pa] ve trace alanları korunmalıdır.
  subroutine test_arrhenius_runtime_roundtrip()
    type(tts_identification_result_t) :: identification
    type(tts_shift_law_identification_result_t) :: laws
    type(tts_arrhenius_runtime_export_t) :: runtime_export
    type(dynamic_modulus_evaluation_t) :: evaluation
    type(temperature_shift_domain_t) :: domain
    real(dp) :: beta_k
    real(dp) :: expected_shift
    real(dp) :: physical_frequency_hz
    real(dp) :: reduced_frequency_hz
    integer :: point_index

    call build_arrhenius_identification(identification, laws)
    runtime_export = create_tts_arrhenius_runtime_export( &
      identification, laws%arrhenius)
    call assert_true(runtime_export%available, &
      "Identified Arrhenius runtime export oluşturulamadı.")
    domain = get_temperature_shift_domain(runtime_export%shift_provider)
    call assert_close(domain%minimum_temperature_k, temperatures_k(1), &
      0.0_dp, "Arrhenius runtime Tmin measured domain değil.")
    call assert_close(domain%maximum_temperature_k, temperatures_k(5), &
      0.0_dp, "Arrhenius runtime Tmax measured domain değil.")

    point_index = size(identification%runtime_master_table)/2
    reduced_frequency_hz = identification%runtime_master_table(point_index) &
      %reduced_frequency_hz
    beta_k = apparent_activation_energy_j_per_mol / &
      (universal_gas_constant_j_per_mol_k*log(10.0_dp))
    expected_shift = beta_k*(1.0_dp/temperatures_k(1) - &
      1.0_dp/reference_temperature_k)
    physical_frequency_hz = reduced_frequency_hz/(10.0_dp**expected_shift)
    evaluation = evaluate_dynamic_shear_modulus( &
      runtime_export%thermorheological_provider, physical_frequency_hz, &
      temperatures_k(1))
    call assert_true(evaluation%shift_model_kind == &
      ARRHENIUS_TEMPERATURE_SHIFT, &
      "Arrhenius parametric runtime model trace'i hatalı.")
    call assert_close(evaluation%log10_a_t, expected_shift, 1.0e-5_dp, &
      "Arrhenius runtime fitted shift analytical truth ile uyuşmuyor.")
    call assert_close(evaluation%lookup_frequency_hz, reduced_frequency_hz, &
      3.0e-5_dp, "Arrhenius runtime reduced-frequency round-trip hatalı.")
    call assert_close(evaluation%modulus%storage_modulus, &
      identification%runtime_master_table(point_index)%storage_modulus_pa, &
      3.0e-5_dp, "Arrhenius runtime G' master-table round-trip hatalı.")
    call assert_close(evaluation%modulus%loss_modulus, &
      identification%runtime_master_table(point_index)%loss_modulus_pa, &
      3.0e-5_dp, "Arrhenius runtime G'' master-table round-trip hatalı.")
  end subroutine test_arrhenius_runtime_roundtrip

  !> Exact WLF horizontal shifts ile V0.8.1 master table ve identifiable
  !! profiled C1/C2 fit'i oluşturur. Existing V0.8.0 WLF provider üzerinden
  !! off-reference f_r=a_T*f ve G'/G'' [Pa] round-trip doğrulanır.
  subroutine test_wlf_runtime_roundtrip()
    type(tts_identification_result_t) :: identification
    type(tts_shift_law_identification_result_t) :: laws
    type(tts_wlf_runtime_export_t) :: runtime_export
    type(dynamic_modulus_evaluation_t) :: evaluation
    type(temperature_shift_domain_t) :: domain
    real(dp) :: expected_shift
    real(dp) :: physical_frequency_hz
    real(dp) :: reduced_frequency_hz
    integer :: point_index

    call build_wlf_identification(identification, laws)
    runtime_export = create_tts_wlf_runtime_export( &
      identification, laws%wlf)
    call assert_true(runtime_export%available, &
      "Identified WLF runtime export oluşturulamadı.")
    domain = get_temperature_shift_domain(runtime_export%shift_provider)
    call assert_close(domain%minimum_temperature_k, temperatures_k(1), &
      0.0_dp, "WLF runtime Tmin measured domain değil.")
    call assert_close(domain%maximum_temperature_k, temperatures_k(5), &
      0.0_dp, "WLF runtime Tmax measured domain değil.")

    point_index = size(identification%runtime_master_table)/2
    reduced_frequency_hz = identification%runtime_master_table(point_index) &
      %reduced_frequency_hz
    expected_shift = -known_wlf_c1 * &
      (temperatures_k(5)-reference_temperature_k) / &
      (known_wlf_c2_k+temperatures_k(5)-reference_temperature_k)
    physical_frequency_hz = reduced_frequency_hz/(10.0_dp**expected_shift)
    evaluation = evaluate_dynamic_shear_modulus( &
      runtime_export%thermorheological_provider, physical_frequency_hz, &
      temperatures_k(5))
    call assert_true(evaluation%shift_model_kind == WLF_TEMPERATURE_SHIFT, &
      "WLF parametric runtime model trace'i hatalı.")
    call assert_close(evaluation%log10_a_t, expected_shift, 5.0e-6_dp, &
      "WLF runtime fitted shift analytical truth ile uyuşmuyor.")
    call assert_close(evaluation%lookup_frequency_hz, reduced_frequency_hz, &
      2.0e-5_dp, "WLF runtime reduced-frequency round-trip hatalı.")
    call assert_close(evaluation%modulus%storage_modulus, &
      identification%runtime_master_table(point_index)%storage_modulus_pa, &
      2.0e-5_dp, "WLF runtime G' master-table round-trip hatalı.")
    call assert_close(evaluation%modulus%loss_modulus, &
      identification%runtime_master_table(point_index)%loss_modulus_pa, &
      2.0e-5_dp, "WLF runtime G'' master-table round-trip hatalı.")
  end subroutine test_wlf_runtime_roundtrip

  subroutine test_invalid_fit_exports()
    type(tts_identification_result_t) :: identification
    type(tts_shift_law_identification_result_t) :: laws
    type(tts_arrhenius_fit_result_t) :: invalid_arrhenius
    type(tts_wlf_fit_result_t) :: invalid_wlf
    type(tts_arrhenius_runtime_export_t) :: arrhenius_export
    type(tts_wlf_runtime_export_t) :: wlf_export

    call build_arrhenius_identification(identification, laws)
    invalid_arrhenius = laws%arrhenius
    invalid_arrhenius%fit_available = .false.
    arrhenius_export = create_tts_arrhenius_runtime_export( &
      identification, invalid_arrhenius)
    call assert_true(.not. arrhenius_export%available, &
      "Invalid Arrhenius fit runtime'a açıldı.")

    call build_wlf_identification(identification, laws)
    invalid_wlf = laws%wlf
    invalid_wlf%parameter_identifiable = .false.
    wlf_export = create_tts_wlf_runtime_export(identification, invalid_wlf)
    call assert_true(.not. wlf_export%available, &
      "Poorly identified WLF runtime'a açıldı.")

    invalid_wlf = laws%wlf
    invalid_wlf%c1 = -1.0_dp
    wlf_export = create_tts_wlf_runtime_export(identification, invalid_wlf)
    call assert_true(.not. wlf_export%available, &
      "Invalid negative WLF C1 runtime'a açıldı.")

    invalid_wlf = laws%wlf
    invalid_wlf%c2_k = 10.0_dp
    invalid_wlf%p_c1_over_c2_per_k = invalid_wlf%c1/invalid_wlf%c2_k
    invalid_wlf%q_inverse_c2_per_k = 1.0_dp/invalid_wlf%c2_k
    wlf_export = create_tts_wlf_runtime_export(identification, invalid_wlf)
    call assert_true(.not. wlf_export%available, &
      "Calibrated domain içinde pole oluşturan WLF C2 runtime'a açıldı.")
  end subroutine test_invalid_fit_exports

  subroutine exercise_outside_domain(case_name)
    character(len=*), intent(in) :: case_name
    type(tts_identification_result_t) :: identification
    type(tts_shift_law_identification_result_t) :: laws
    type(tts_arrhenius_runtime_export_t) :: arrhenius_export
    type(tts_wlf_runtime_export_t) :: wlf_export
    type(dynamic_modulus_evaluation_t) :: evaluation

    select case (case_name)
    case ("arrhenius_outside_domain")
      call build_arrhenius_identification(identification, laws)
      arrhenius_export = create_tts_arrhenius_runtime_export( &
        identification, laws%arrhenius)
      evaluation = evaluate_dynamic_shear_modulus( &
        arrhenius_export%thermorheological_provider, 1.0_dp, 250.0_dp)
    case ("wlf_outside_domain")
      call build_wlf_identification(identification, laws)
      wlf_export = create_tts_wlf_runtime_export(identification, laws%wlf)
      evaluation = evaluate_dynamic_shear_modulus( &
        wlf_export%thermorheological_provider, 1.0_dp, 340.0_dp)
    case default
      error stop "Bilinmeyen parametric runtime validation selector."
    end select
  end subroutine exercise_outside_domain

  subroutine build_arrhenius_identification(identification, laws)
    type(tts_identification_result_t), intent(out) :: identification
    type(tts_shift_law_identification_result_t), intent(out) :: laws
    type(tts_material_family_t) :: family
    real(dp) :: beta_k
    real(dp) :: shifts(size(temperatures_k))
    integer :: i

    beta_k = apparent_activation_energy_j_per_mol / &
      (universal_gas_constant_j_per_mol_k*log(10.0_dp))
    do i = 1, size(temperatures_k)
      shifts(i) = beta_k*(1.0_dp/temperatures_k(i) - &
        1.0_dp/reference_temperature_k)
    end do
    family = make_generalized_maxwell_trs_family(temperatures_k, shifts)
    identification = identify_tts_master_curve(family, "ISO-3")
    call assert_true(identification%status == TTS_IDENTIFICATION_SUCCESS .and. &
      identification%runtime_export_ready, &
      "Arrhenius synthetic V0.8.1 empirical identification başarısız.")
    laws = fit_tts_shift_laws(identification)
    call assert_true(laws%arrhenius%status == SHIFT_LAW_FIT_SUCCESS .and. &
      laws%arrhenius%fit_available, &
      "Arrhenius synthetic adjacent-pair law fit başarısız.")
  end subroutine build_arrhenius_identification

  subroutine build_wlf_identification(identification, laws)
    type(tts_identification_result_t), intent(out) :: identification
    type(tts_shift_law_identification_result_t), intent(out) :: laws
    type(tts_material_family_t) :: family
    real(dp) :: shifts(size(temperatures_k))
    integer :: i

    do i = 1, size(temperatures_k)
      shifts(i) = -known_wlf_c1 * &
        (temperatures_k(i)-reference_temperature_k) / &
        (known_wlf_c2_k+temperatures_k(i)-reference_temperature_k)
    end do
    family = make_generalized_maxwell_trs_family(temperatures_k, shifts)
    identification = identify_tts_master_curve(family, "ISO-3")
    call assert_true(identification%status == TTS_IDENTIFICATION_SUCCESS .and. &
      identification%runtime_export_ready, &
      "WLF synthetic V0.8.1 empirical identification başarısız.")
    laws = fit_tts_shift_laws(identification)
    call assert_true(laws%wlf%status == SHIFT_LAW_FIT_SUCCESS .and. &
      laws%wlf%fit_available .and. laws%wlf%parameter_identifiable, &
      "WLF synthetic adjacent-pair law fit başarısız.")
  end subroutine build_wlf_identification

end program test_tts_parametric_runtime_roundtrip
