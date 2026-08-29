module tms_tts_repeatability_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_sample_statistics, only : calculate_sample_statistics
  use tms_bootstrap, only : cluster_bootstrap_plan_t, &
    bootstrap_mean_result_t, BOOTSTRAP_SUCCESS, &
    BOOTSTRAP_INSUFFICIENT_CAMPAIGNS, is_valid_bootstrap_configuration, &
    create_cluster_bootstrap_plan, calculate_cluster_bootstrap_mean
  use tms_tts_types, only : tts_common_test_state_t, &
    tts_identification_result_t, tts_validation_result_t, &
    TTS_IDENTIFICATION_SUCCESS, PAIR_SHIFT_SUCCESS, &
    PAIR_SHIFT_STORAGE_ONLY, &
    validate_tts_material_family, are_tts_values_machine_equivalent
  use tms_tts_shift_law_types, only : &
    tts_shift_law_identification_result_t, &
    SHIFT_LAW_FIT_INSUFFICIENT_DATA, WLF_FIT_POORLY_IDENTIFIED
  use tms_tts_shift_law_validation, only : fit_tts_shift_laws
  use tms_tts_repeatability_types, only : &
    tts_repeatability_campaign_t, tts_repeatability_study_configuration_t, &
    tts_repeatability_study_result_t, REPLICATE_BASIS_UNSPECIFIED, &
    INDEPENDENT_SPECIMEN_CAMPAIGN, SAME_SPECIMEN_RERUN, &
    TTS_REPEATABILITY_SUCCESS, TTS_REPEATABILITY_INVALID_INPUT, &
    TTS_REPEATABILITY_INCOMPATIBLE_STATE, &
    TTS_REPEATABILITY_TEMPERATURE_SET_MISMATCH, &
    TTS_REPEATABILITY_REFERENCE_NOT_FOUND, &
    TTS_REPEATABILITY_INSUFFICIENT_REPLICATES, &
    TTS_REPEATABILITY_NO_INDEPENDENT_CAMPAIGNS, &
    TTS_REPEATABILITY_BOOTSTRAP_UNAVAILABLE, &
    TTS_REPEATABILITY_NONFINITE_DATA
  implicit none
  private

  type :: canonical_campaign_data_t
    real(dp), allocatable :: temperatures_k(:)
    real(dp), allocatable :: normalized_shifts(:)
    real(dp), allocatable :: adjacent_pair_shifts(:)
    type(tts_shift_law_identification_result_t) :: shift_laws
  end type canonical_campaign_data_t

  public :: analyze_tts_repeatability

contains

  !> Birden fazla complete V0.8.1 campaign'i canonical physical temperature
  !! keys ile eşler, s_common(T)=s_original(T)-s_original(T_common)
  !! normalizasyonunu uygular ve campaign-level descriptive/cluster bootstrap
  !! evidence üretir. Sıcaklık [K], s ve delta_s boyutsuzdur. V0.8.2
  !! fit_tts_shift_laws aynen reuse edilir; input campaign/identification
  !! sonuçları intent(in) kalır. Complete campaign statistical sample unit'tir;
  !! frequency, isotherm ve pair noktaları replicate sayılmaz. Ayrıntılı model
  !! docs/mathematics/tts_repeatability_bootstrap.md belgesindedir.
  pure function analyze_tts_repeatability( &
      campaigns, common_reference_temperature_k, configuration) result(study)
    type(tts_repeatability_campaign_t), intent(in) :: campaigns(:)
    real(dp), intent(in) :: common_reference_temperature_k
    type(tts_repeatability_study_configuration_t), intent(in), optional :: &
      configuration
    type(tts_repeatability_study_result_t) :: study

    type(tts_repeatability_study_configuration_t) :: settings
    type(canonical_campaign_data_t), allocatable :: canonical_data(:)
    type(cluster_bootstrap_plan_t) :: bootstrap_plan
    type(tts_common_test_state_t) :: reference_common_state
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: descriptive_values(:)
    integer, allocatable :: all_campaign_indices(:)
    integer, allocatable :: independent_campaign_indices(:)
    logical, allocatable :: descriptive_mask(:)
    logical, allocatable :: independent_mask(:)
    logical, allocatable :: usable(:)
    integer :: campaign_index
    integer :: extraction_status
    integer :: pair_index
    integer :: temperature_index
    character(len=256) :: extraction_message

    if (present(configuration)) settings = configuration
    study%study_identifier = settings%study_identifier
    study%bootstrap_configuration = settings%bootstrap
    study%common_reference_temperature_k = common_reference_temperature_k
    study%total_campaign_count = size(campaigns)

    if (size(campaigns) == 0) then
      study%status = TTS_REPEATABILITY_INSUFFICIENT_REPLICATES
      study%message = "Repeatability study en az bir complete campaign ister."
      return
    end if
    if (.not. ieee_is_finite(common_reference_temperature_k)) then
      study%status = TTS_REPEATABILITY_NONFINITE_DATA
      study%message = "Common reference temperature sonlu olmalıdır."
      return
    end if
    if (common_reference_temperature_k <= 0.0_dp) then
      study%status = TTS_REPEATABILITY_INVALID_INPUT
      study%message = "Common reference temperature pozitif Kelvin olmalıdır."
      return
    end if
    if (.not. is_valid_bootstrap_configuration(settings%bootstrap)) then
      study%status = TTS_REPEATABILITY_INVALID_INPUT
      study%message = "Bootstrap draw count/confidence ayarları geçersiz."
      return
    end if
    if (.not. campaign_identifiers_are_valid(campaigns)) then
      study%status = TTS_REPEATABILITY_INVALID_INPUT
      study%message = "Campaign kimlikleri boş olamaz ve unique olmalıdır."
      return
    end if
    if (.not. replicate_bases_are_known(campaigns)) then
      study%status = TTS_REPEATABILITY_INVALID_INPUT
      study%message = "Tanımsız replicate basis değeri bulundu."
      return
    end if

    allocate(canonical_data(size(campaigns)))
    do campaign_index = 1, size(campaigns)
      if (campaign_index == 1) then
        reference_common_state = &
          campaigns(campaign_index)%identification%source_family%common_state
      else
        if (.not. common_states_are_compatible(reference_common_state, &
            campaigns(campaign_index)%identification%source_family &
              %common_state)) then
          study%status = TTS_REPEATABILITY_INCOMPATIBLE_STATE
          study%message = &
            "Campaign physics-critical common test states uyumsuz."
          return
        end if
      end if
      call extract_canonical_campaign_data(campaigns(campaign_index), &
        common_reference_temperature_k, canonical_data(campaign_index), &
        extraction_status, extraction_message)
      if (extraction_status /= TTS_REPEATABILITY_SUCCESS) then
        study%status = extraction_status
        study%message = extraction_message
        return
      end if
      if (campaign_index > 1) then
        if (.not. temperature_sets_are_equivalent( &
            canonical_data(1)%temperatures_k, &
            canonical_data(campaign_index)%temperatures_k)) then
          study%status = TTS_REPEATABILITY_TEMPERATURE_SET_MISMATCH
          study%message = &
            "Campaign canonical measured temperature setleri eşleşmiyor."
          return
        end if
      end if
      canonical_data(campaign_index)%shift_laws = &
        fit_tts_shift_laws(campaigns(campaign_index)%identification)
    end do

    allocate(descriptive_mask(size(campaigns)), &
      independent_mask(size(campaigns)), all_campaign_indices(size(campaigns)))
    descriptive_mask = .true.
    independent_mask = .false.
    do campaign_index = 1, size(campaigns)
      all_campaign_indices(campaign_index) = campaign_index
      select case (campaigns(campaign_index)%replicate_basis)
      case (INDEPENDENT_SPECIMEN_CAMPAIGN)
        independent_mask(campaign_index) = .true.
        study%independent_campaign_count = &
          study%independent_campaign_count + 1
      case (SAME_SPECIMEN_RERUN)
        study%same_specimen_rerun_count = &
          study%same_specimen_rerun_count + 1
        if (.not. settings%include_same_specimen_reruns_in_descriptive) &
          descriptive_mask(campaign_index) = .false.
      case (REPLICATE_BASIS_UNSPECIFIED)
        study%unspecified_replicate_basis_count = &
          study%unspecified_replicate_basis_count + 1
      end select
    end do
    study%descriptive_campaign_count = count(descriptive_mask)
    if (study%descriptive_campaign_count == 0) then
      study%status = TTS_REPEATABILITY_INSUFFICIENT_REPLICATES
      study%message = "Study policy descriptive campaign bırakmadı."
      return
    end if
    study%descriptive_statistics_available = .true.
    study%intralaboratory_context_explicit = &
      intralaboratory_context_is_explicit(campaigns)

    allocate(study%campaign_provenance(size(campaigns)))
    do campaign_index = 1, size(campaigns)
      call copy_campaign_provenance(campaigns(campaign_index), &
        descriptive_mask(campaign_index), independent_mask(campaign_index), &
        study%campaign_provenance(campaign_index))
    end do

    allocate(values(size(campaigns)), usable(size(campaigns)))
    usable = .true.
    allocate(study%pair_results( &
      size(canonical_data(1)%adjacent_pair_shifts)))
    do pair_index = 1, size(study%pair_results)
      study%pair_results(pair_index)%lower_temperature_k = &
        canonical_data(1)%temperatures_k(pair_index)
      study%pair_results(pair_index)%upper_temperature_k = &
        canonical_data(1)%temperatures_k(pair_index + 1)
      do campaign_index = 1, size(campaigns)
        values(campaign_index) = &
          canonical_data(campaign_index)%adjacent_pair_shifts(pair_index)
      end do
      descriptive_values = pack(values, descriptive_mask)
      study%pair_results(pair_index)%delta_s_statistics = &
        calculate_sample_statistics(descriptive_values)
    end do

    allocate(study%absolute_shift_results( &
      size(canonical_data(1)%temperatures_k)))
    do temperature_index = 1, size(study%absolute_shift_results)
      study%absolute_shift_results(temperature_index)%temperature_k = &
        canonical_data(1)%temperatures_k(temperature_index)
      study%absolute_shift_results(temperature_index)%is_reference_anchor = &
        are_tts_values_machine_equivalent( &
          canonical_data(1)%temperatures_k(temperature_index), &
          common_reference_temperature_k)
      study%absolute_shift_results(temperature_index) &
        %uncertainty_informative = .not. study%absolute_shift_results( &
          temperature_index)%is_reference_anchor
      do campaign_index = 1, size(campaigns)
        values(campaign_index) = &
          canonical_data(campaign_index)%normalized_shifts(temperature_index)
      end do
      descriptive_values = pack(values, descriptive_mask)
      study%absolute_shift_results(temperature_index) &
        %log10_a_t_statistics = &
          calculate_sample_statistics(descriptive_values)
    end do

    call build_arrhenius_cohort(canonical_data, descriptive_mask, &
      study%arrhenius)
    call build_wlf_cohort(canonical_data, descriptive_mask, study%wlf)

    independent_campaign_indices = pack(all_campaign_indices, independent_mask)
    if (size(independent_campaign_indices) == 0) then
      study%bootstrap_status = TTS_REPEATABILITY_NO_INDEPENDENT_CAMPAIGNS
      study%bootstrap_message = &
        "Independent specimen campaign yok; cluster bootstrap unavailable."
    else if (size(independent_campaign_indices) == 1) then
      study%bootstrap_status = TTS_REPEATABILITY_BOOTSTRAP_UNAVAILABLE
      study%bootstrap_message = &
        "Tek independent campaign repeatability CI üretmez."
    else
      bootstrap_plan = create_cluster_bootstrap_plan( &
        independent_campaign_indices, settings%bootstrap)
      if (bootstrap_plan%status == BOOTSTRAP_SUCCESS) then
        study%independent_cluster_bootstrap_available = .true.
        study%bootstrap_status = TTS_REPEATABILITY_SUCCESS
        study%bootstrap_message = &
          "Independent complete-campaign cluster bootstrap hazır."
        allocate(study%bootstrap_population_campaign_identifiers( &
          size(independent_campaign_indices)))
        do campaign_index = 1, size(independent_campaign_indices)
          study%bootstrap_population_campaign_identifiers(campaign_index) = &
            campaigns(independent_campaign_indices(campaign_index)) &
              %campaign_identifier
        end do
        call attach_bootstrap_intervals( &
          canonical_data, bootstrap_plan, study)
      else if (bootstrap_plan%status == &
          BOOTSTRAP_INSUFFICIENT_CAMPAIGNS) then
        study%bootstrap_status = TTS_REPEATABILITY_BOOTSTRAP_UNAVAILABLE
        study%bootstrap_message = &
          "Independent campaign population bootstrap için yetersiz."
      else
        study%bootstrap_status = TTS_REPEATABILITY_BOOTSTRAP_UNAVAILABLE
        study%bootstrap_message = "Cluster bootstrap planı oluşturulamadı."
      end if
    end if

    study%status = TTS_REPEATABILITY_SUCCESS
    study%message = &
      "Campaign-level repeatability ve uncertainty evidence üretildi."
  end function analyze_tts_repeatability

  !> Complete identification'ı canonical ascending T key'lerine dönüştürür.
  !! Empirical shift ve adjacent pair eşlemesi array index sırasına değil
  !! physical (T_i,T_j) anahtarına dayanır; ters orientation delta_s işaretini
  !! çevirir. Common-reference normalizasyonu input'u mutate etmez.
  pure subroutine extract_canonical_campaign_data( &
      campaign, common_reference_temperature_k, data, status, message)
    type(tts_repeatability_campaign_t), intent(in) :: campaign
    real(dp), intent(in) :: common_reference_temperature_k
    type(canonical_campaign_data_t), intent(out) :: data
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    type(tts_identification_result_t) :: identification
    type(tts_validation_result_t) :: validation
    integer, allocatable :: sorted_source_indices(:)
    real(dp), allocatable :: original_shifts(:)
    integer :: canonical_index
    integer :: match_count
    integer :: reference_position

    status = TTS_REPEATABILITY_INVALID_INPUT
    message = ""
    identification = campaign%identification
    if (identification%status /= TTS_IDENTIFICATION_SUCCESS .or. &
        .not. identification%shift_chain_available .or. &
        .not. identification%master_cloud_available .or. &
        .not. identification%runtime_export_ready) then
      message = "Campaign complete/usable V0.8.1 identification içermiyor."
      return
    end if
    validation = validate_tts_material_family(identification%source_family)
    if (.not. validation%valid) then
      message = "Campaign source family V0.8.1 validation'dan geçmedi."
      return
    end if
    if (.not. allocated(identification%empirical_shifts)) then
      message = "Campaign empirical shift table içermiyor."
      return
    end if
    if (.not. allocated(identification%pair_shift_results)) then
      message = "Campaign adjacent pair shift results içermiyor."
      return
    end if
    if (size(identification%empirical_shifts) /= &
        size(identification%source_family%isotherms)) then
      message = "Campaign empirical shift/temperature boyutları tutarsız."
      return
    end if
    if (size(identification%pair_shift_results) /= &
        size(identification%source_family%isotherms) - 1) then
      message = "Campaign complete adjacent pair zinciri içermiyor."
      return
    end if
    if (.not. ieee_is_finite(identification%reference_temperature_k)) then
      status = TTS_REPEATABILITY_NONFINITE_DATA
      message = "Campaign original reference temperature sonlu değil."
      return
    end if
    if (identification%reference_isotherm_index < 1 .or. &
        identification%reference_isotherm_index > &
          size(identification%source_family%isotherms)) then
      message = "Campaign original reference provenance geçersiz."
      return
    end if
    if (.not. are_tts_values_machine_equivalent( &
        identification%reference_temperature_k, &
        identification%source_family%isotherms( &
          identification%reference_isotherm_index)%temperature_k)) then
      message = "Campaign original reference provenance tutarsız."
      return
    end if

    sorted_source_indices = canonical_temperature_indices( &
      identification%source_family%isotherms%temperature_k)
    allocate(data%temperatures_k(size(sorted_source_indices)), &
      original_shifts(size(sorted_source_indices)), &
      data%normalized_shifts(size(sorted_source_indices)), &
      data%adjacent_pair_shifts(size(sorted_source_indices) - 1))
    reference_position = 0
    do canonical_index = 1, size(sorted_source_indices)
      data%temperatures_k(canonical_index) = &
        identification%source_family%isotherms( &
          sorted_source_indices(canonical_index))%temperature_k
      call find_empirical_shift_by_temperature(identification, &
        data%temperatures_k(canonical_index), original_shifts(canonical_index), &
        match_count, status)
      if (status /= TTS_REPEATABILITY_SUCCESS) then
        message = "Campaign empirical shift value/provenance geçersiz."
        return
      end if
      if (match_count /= 1) then
        message = "Physical temperature key unique empirical shift bulamadı."
        return
      end if
      if (are_tts_values_machine_equivalent( &
          data%temperatures_k(canonical_index), &
          common_reference_temperature_k)) reference_position = &
            canonical_index
    end do
    if (reference_position == 0) then
      status = TTS_REPEATABILITY_REFERENCE_NOT_FOUND
      message = "Common reference measured temperature setinde bulunamadı."
      return
    end if
    data%normalized_shifts = original_shifts - &
      original_shifts(reference_position)
    data%normalized_shifts(reference_position) = 0.0_dp
    if (.not. all(ieee_is_finite(data%normalized_shifts))) then
      status = TTS_REPEATABILITY_NONFINITE_DATA
      message = "Common-reference normalization nonfinite sonuç üretti."
      return
    end if

    do canonical_index = 1, size(data%adjacent_pair_shifts)
      call find_canonical_pair_shift(identification, &
        data%temperatures_k(canonical_index), &
        data%temperatures_k(canonical_index + 1), &
        data%adjacent_pair_shifts(canonical_index), match_count, status)
      if (status /= TTS_REPEATABILITY_SUCCESS) then
        message = "Campaign adjacent pair value/provenance geçersiz."
        return
      end if
      if (match_count /= 1) then
        message = "Canonical physical pair key unique shift bulamadı."
        return
      end if
      if (.not. are_tts_values_machine_equivalent( &
          data%adjacent_pair_shifts(canonical_index), &
          original_shifts(canonical_index + 1) - &
            original_shifts(canonical_index))) then
        message = "Empirical shift chain ile adjacent pair delta_s tutarsız."
        return
      end if
    end do
    status = TTS_REPEATABILITY_SUCCESS
    message = "Campaign canonical temperature map hazır."
  end subroutine extract_canonical_campaign_data

  pure subroutine find_empirical_shift_by_temperature( &
      identification, temperature_k, shift, match_count, status)
    type(tts_identification_result_t), intent(in) :: identification
    real(dp), intent(in) :: temperature_k
    real(dp), intent(out) :: shift
    integer, intent(out) :: match_count
    integer, intent(out) :: status
    integer :: i

    shift = 0.0_dp
    match_count = 0
    status = TTS_REPEATABILITY_INVALID_INPUT
    do i = 1, size(identification%empirical_shifts)
      if (.not. ieee_is_finite( &
          identification%empirical_shifts(i)%temperature_k) .or. &
          .not. ieee_is_finite( &
            identification%empirical_shifts(i)%log10_a_t) .or. &
          .not. ieee_is_finite(identification%empirical_shifts(i)%a_t)) then
        status = TTS_REPEATABILITY_NONFINITE_DATA
        return
      end if
      if (identification%empirical_shifts(i)%a_t <= 0.0_dp) return
      if (are_tts_values_machine_equivalent( &
          identification%empirical_shifts(i)%temperature_k, &
          temperature_k)) then
        match_count = match_count + 1
        shift = identification%empirical_shifts(i)%log10_a_t
      end if
    end do
    status = TTS_REPEATABILITY_SUCCESS
  end subroutine find_empirical_shift_by_temperature

  pure subroutine find_canonical_pair_shift( &
      identification, lower_temperature_k, upper_temperature_k, shift, &
      match_count, status)
    type(tts_identification_result_t), intent(in) :: identification
    real(dp), intent(in) :: lower_temperature_k
    real(dp), intent(in) :: upper_temperature_k
    real(dp), intent(out) :: shift
    integer, intent(out) :: match_count
    integer, intent(out) :: status

    real(dp) :: moving_temperature_k
    real(dp) :: reference_temperature_k
    integer :: moving_index
    integer :: pair_index
    integer :: reference_index

    shift = 0.0_dp
    match_count = 0
    status = TTS_REPEATABILITY_INVALID_INPUT
    do pair_index = 1, size(identification%pair_shift_results)
      reference_index = identification%pair_shift_results(pair_index) &
        %reference_isotherm_index
      moving_index = identification%pair_shift_results(pair_index) &
        %moving_isotherm_index
      if (reference_index < 1 .or. reference_index > &
          size(identification%source_family%isotherms) .or. &
          moving_index < 1 .or. moving_index > &
            size(identification%source_family%isotherms) .or. &
          reference_index == moving_index) return
      if (.not. identification%pair_shift_results(pair_index) &
          %shift_available .or. &
          (identification%pair_shift_results(pair_index)%status /= &
            PAIR_SHIFT_SUCCESS .and. &
           identification%pair_shift_results(pair_index)%status /= &
            PAIR_SHIFT_STORAGE_ONLY)) return
      if (.not. ieee_is_finite( &
          identification%pair_shift_results(pair_index)%delta_s)) then
        status = TTS_REPEATABILITY_NONFINITE_DATA
        return
      end if
      reference_temperature_k = identification%source_family%isotherms( &
        reference_index)%temperature_k
      moving_temperature_k = identification%source_family%isotherms( &
        moving_index)%temperature_k
      if (are_tts_values_machine_equivalent(reference_temperature_k, &
          lower_temperature_k) .and. &
          are_tts_values_machine_equivalent(moving_temperature_k, &
            upper_temperature_k)) then
        match_count = match_count + 1
        shift = identification%pair_shift_results(pair_index)%delta_s
      else if (are_tts_values_machine_equivalent(reference_temperature_k, &
          upper_temperature_k) .and. &
          are_tts_values_machine_equivalent(moving_temperature_k, &
            lower_temperature_k)) then
        match_count = match_count + 1
        shift = -identification%pair_shift_results(pair_index)%delta_s
      end if
    end do
    status = TTS_REPEATABILITY_SUCCESS
  end subroutine find_canonical_pair_shift

  pure subroutine build_arrhenius_cohort(data, descriptive_mask, cohort)
    use tms_tts_repeatability_types, only : &
      tts_arrhenius_repeatability_result_t
    type(canonical_campaign_data_t), intent(in) :: data(:)
    logical, intent(in) :: descriptive_mask(:)
    type(tts_arrhenius_repeatability_result_t), intent(out) :: cohort

    real(dp), allocatable :: beta_values(:)
    real(dp), allocatable :: energy_values(:)
    logical, allocatable :: usable(:)
    integer :: i

    cohort%total_campaign_count = size(data)
    allocate(beta_values(size(data)), energy_values(size(data)), &
      usable(size(data)))
    beta_values = 0.0_dp
    energy_values = 0.0_dp
    usable = .false.
    do i = 1, size(data)
      if (data(i)%shift_laws%arrhenius%fit_available) then
        cohort%fit_available_count = cohort%fit_available_count + 1
        beta_values(i) = data(i)%shift_laws%arrhenius%beta_k
        energy_values(i) = data(i)%shift_laws%arrhenius &
          %apparent_activation_energy_j_per_mol
        usable(i) = descriptive_mask(i)
      else
        cohort%fit_unavailable_count = cohort%fit_unavailable_count + 1
        if (data(i)%shift_laws%arrhenius%status /= &
            SHIFT_LAW_FIT_INSUFFICIENT_DATA) then
          cohort%invalid_or_nonphysical_fit_count = &
            cohort%invalid_or_nonphysical_fit_count + 1
        end if
      end if
    end do
    cohort%beta_k_statistics = &
      calculate_sample_statistics(pack(beta_values, usable))
    cohort%apparent_activation_energy_statistics = &
      calculate_sample_statistics(pack(energy_values, usable))
  end subroutine build_arrhenius_cohort

  pure subroutine build_wlf_cohort(data, descriptive_mask, cohort)
    use tms_tts_repeatability_types, only : tts_wlf_repeatability_result_t
    type(canonical_campaign_data_t), intent(in) :: data(:)
    logical, intent(in) :: descriptive_mask(:)
    type(tts_wlf_repeatability_result_t), intent(out) :: cohort

    real(dp), allocatable :: c1_values(:)
    real(dp), allocatable :: c2_values(:)
    real(dp), allocatable :: p_values(:)
    real(dp), allocatable :: q_values(:)
    logical, allocatable :: usable(:)
    integer :: i

    cohort%total_campaign_count = size(data)
    allocate(c1_values(size(data)), c2_values(size(data)), &
      p_values(size(data)), q_values(size(data)), usable(size(data)))
    c1_values = 0.0_dp
    c2_values = 0.0_dp
    p_values = 0.0_dp
    q_values = 0.0_dp
    usable = .false.
    do i = 1, size(data)
      if (data(i)%shift_laws%wlf%fit_available) then
        cohort%fit_available_count = cohort%fit_available_count + 1
      else
        cohort%fit_unavailable_count = cohort%fit_unavailable_count + 1
        cohort%invalid_fit_count = cohort%invalid_fit_count + 1
      end if
      if (data(i)%shift_laws%wlf%status == WLF_FIT_POORLY_IDENTIFIED .or. &
          (data(i)%shift_laws%wlf%fit_available .and. &
            .not. data(i)%shift_laws%wlf%parameter_identifiable)) then
        cohort%poorly_identified_count = cohort%poorly_identified_count + 1
      end if
      if (data(i)%shift_laws%wlf%fit_available .and. &
          data(i)%shift_laws%wlf%parameter_identifiable) then
        cohort%parameter_identifiable_count = &
          cohort%parameter_identifiable_count + 1
        c1_values(i) = data(i)%shift_laws%wlf%c1
        c2_values(i) = data(i)%shift_laws%wlf%c2_k
        p_values(i) = data(i)%shift_laws%wlf%p_c1_over_c2_per_k
        q_values(i) = data(i)%shift_laws%wlf%q_inverse_c2_per_k
        usable(i) = descriptive_mask(i)
      end if
    end do
    cohort%c1_statistics = &
      calculate_sample_statistics(pack(c1_values, usable))
    cohort%c2_k_statistics = &
      calculate_sample_statistics(pack(c2_values, usable))
    cohort%p_c1_over_c2_per_k_statistics = &
      calculate_sample_statistics(pack(p_values, usable))
    cohort%q_inverse_c2_per_k_statistics = &
      calculate_sample_statistics(pack(q_values, usable))
  end subroutine build_wlf_cohort

  !> Tek whole-campaign bootstrap planını bütün empirical ve parametric
  !! quantities için reuse eder. Quantity-specific fit availability plan
  !! üretildikten sonra uygulanır; WLF population önceden filtrelenmez.
  pure subroutine attach_bootstrap_intervals(data, plan, study)
    type(canonical_campaign_data_t), intent(in) :: data(:)
    type(cluster_bootstrap_plan_t), intent(in) :: plan
    type(tts_repeatability_study_result_t), intent(inout) :: study

    type(bootstrap_mean_result_t) :: mean_result
    real(dp), allocatable :: values(:)
    logical, allocatable :: usable(:)
    integer :: i
    integer :: pair_index
    integer :: temperature_index

    allocate(values(size(data)), usable(size(data)))
    usable = .true.
    do pair_index = 1, size(study%pair_results)
      do i = 1, size(data)
        values(i) = data(i)%adjacent_pair_shifts(pair_index)
      end do
      mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
      study%pair_results(pair_index)%mean_bootstrap_interval = &
        mean_result%interval
    end do
    do temperature_index = 1, size(study%absolute_shift_results)
      do i = 1, size(data)
        values(i) = data(i)%normalized_shifts(temperature_index)
      end do
      mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
      study%absolute_shift_results(temperature_index) &
        %mean_bootstrap_interval = mean_result%interval
    end do

    usable = .false.
    values = 0.0_dp
    do i = 1, size(data)
      if (data(i)%shift_laws%arrhenius%fit_available) then
        usable(i) = .true.
        values(i) = data(i)%shift_laws%arrhenius%beta_k
      end if
    end do
    mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
    study%arrhenius%beta_k_mean_bootstrap_interval = mean_result%interval
    do i = 1, size(data)
      if (usable(i)) values(i) = data(i)%shift_laws%arrhenius &
        %apparent_activation_energy_j_per_mol
    end do
    mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
    study%arrhenius%apparent_activation_energy_mean_bootstrap_interval = &
      mean_result%interval

    usable = .false.
    values = 0.0_dp
    do i = 1, size(data)
      if (data(i)%shift_laws%wlf%fit_available .and. &
          data(i)%shift_laws%wlf%parameter_identifiable) then
        usable(i) = .true.
        values(i) = data(i)%shift_laws%wlf%c1
      end if
    end do
    mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
    study%wlf%c1_mean_bootstrap_interval = mean_result%interval
    do i = 1, size(data)
      if (usable(i)) values(i) = data(i)%shift_laws%wlf%c2_k
    end do
    mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
    study%wlf%c2_k_mean_bootstrap_interval = mean_result%interval
    do i = 1, size(data)
      if (usable(i)) values(i) = &
        data(i)%shift_laws%wlf%p_c1_over_c2_per_k
    end do
    mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
    study%wlf%p_mean_bootstrap_interval = mean_result%interval
    do i = 1, size(data)
      if (usable(i)) values(i) = &
        data(i)%shift_laws%wlf%q_inverse_c2_per_k
    end do
    mean_result = calculate_cluster_bootstrap_mean(values, usable, plan)
    study%wlf%q_mean_bootstrap_interval = mean_result%interval
  end subroutine attach_bootstrap_intervals

  pure function canonical_temperature_indices(temperatures_k) result(indices)
    real(dp), intent(in) :: temperatures_k(:)
    integer, allocatable :: indices(:)
    integer :: i
    integer :: j
    integer :: key

    allocate(indices(size(temperatures_k)))
    do i = 1, size(indices)
      indices(i) = i
    end do
    do i = 2, size(indices)
      key = indices(i)
      j = i - 1
      do while (j >= 1)
        if (temperatures_k(indices(j)) <= temperatures_k(key)) exit
        indices(j + 1) = indices(j)
        j = j - 1
      end do
      indices(j + 1) = key
    end do
  end function canonical_temperature_indices

  pure function temperature_sets_are_equivalent(a, b) result(equivalent)
    real(dp), intent(in) :: a(:)
    real(dp), intent(in) :: b(:)
    logical :: equivalent
    integer :: i

    equivalent = size(a) == size(b)
    if (.not. equivalent) return
    do i = 1, size(a)
      if (.not. are_tts_values_machine_equivalent(a(i), b(i))) then
        equivalent = .false.
        return
      end if
    end do
  end function temperature_sets_are_equivalent

  !> Physics-critical common material/test state eşdeğerliğini sınar.
  !! Material, batch, amplitude, prestrain, deformation, conditioning ve test
  !! method gate'tir. Source metadata/specimen/source path provenance'dır ve
  !! farklı independent specimen'i yanlışlıkla engellemez.
  pure function common_states_are_compatible(a, b) result(compatible)
    type(tts_common_test_state_t), intent(in) :: a
    type(tts_common_test_state_t), intent(in) :: b
    logical :: compatible

    compatible = trim(a%material_identifier) == trim(b%material_identifier)
    if (.not. compatible) return
    compatible = trim(a%batch_state_identifier) == &
      trim(b%batch_state_identifier)
    if (.not. compatible) return
    compatible = are_tts_values_machine_equivalent( &
      a%dynamic_strain_amplitude_ratio, b%dynamic_strain_amplitude_ratio)
    if (.not. compatible) return
    compatible = are_tts_values_machine_equivalent( &
      a%static_prestrain_ratio, b%static_prestrain_ratio)
    if (.not. compatible) return
    compatible = a%deformation_mode == b%deformation_mode
    if (.not. compatible) return
    compatible = trim(a%conditioning_description) == &
      trim(b%conditioning_description)
    if (.not. compatible) return
    compatible = trim(a%test_method) == trim(b%test_method)
  end function common_states_are_compatible

  pure function campaign_identifiers_are_valid(campaigns) result(valid)
    type(tts_repeatability_campaign_t), intent(in) :: campaigns(:)
    logical :: valid
    integer :: i
    integer :: j

    valid = .true.
    do i = 1, size(campaigns)
      if (len_trim(campaigns(i)%campaign_identifier) == 0) then
        valid = .false.
        return
      end if
      do j = 1, i - 1
        if (trim(campaigns(i)%campaign_identifier) == &
            trim(campaigns(j)%campaign_identifier)) then
          valid = .false.
          return
        end if
      end do
    end do
  end function campaign_identifiers_are_valid

  pure function replicate_bases_are_known(campaigns) result(known)
    type(tts_repeatability_campaign_t), intent(in) :: campaigns(:)
    logical :: known
    integer :: i

    known = .true.
    do i = 1, size(campaigns)
      if (campaigns(i)%replicate_basis /= REPLICATE_BASIS_UNSPECIFIED .and. &
          campaigns(i)%replicate_basis /= INDEPENDENT_SPECIMEN_CAMPAIGN .and. &
          campaigns(i)%replicate_basis /= SAME_SPECIMEN_RERUN) then
        known = .false.
        return
      end if
    end do
  end function replicate_bases_are_known

  pure function intralaboratory_context_is_explicit(campaigns) &
      result(explicit_context)
    type(tts_repeatability_campaign_t), intent(in) :: campaigns(:)
    logical :: explicit_context
    integer :: i

    explicit_context = &
      len_trim(campaigns(1)%laboratory_identifier) > 0 .and. &
      len_trim(campaigns(1)%instrument_identifier) > 0 .and. &
      len_trim(campaigns(1)%test_protocol_identifier) > 0
    if (.not. explicit_context) return
    do i = 2, size(campaigns)
      explicit_context = &
        trim(campaigns(i)%laboratory_identifier) == &
          trim(campaigns(1)%laboratory_identifier) .and. &
        trim(campaigns(i)%instrument_identifier) == &
          trim(campaigns(1)%instrument_identifier) .and. &
        trim(campaigns(i)%test_protocol_identifier) == &
          trim(campaigns(1)%test_protocol_identifier)
      if (.not. explicit_context) return
    end do
  end function intralaboratory_context_is_explicit

  pure subroutine copy_campaign_provenance( &
      campaign, descriptive, independent, provenance)
    use tms_tts_repeatability_types, only : &
      tts_repeatability_campaign_provenance_t
    type(tts_repeatability_campaign_t), intent(in) :: campaign
    logical, intent(in) :: descriptive
    logical, intent(in) :: independent
    type(tts_repeatability_campaign_provenance_t), intent(out) :: provenance
    integer :: i

    provenance%campaign_identifier = campaign%campaign_identifier
    provenance%replicate_basis = campaign%replicate_basis
    provenance%included_in_descriptive_statistics = descriptive
    provenance%included_in_independent_bootstrap_population = independent
    provenance%laboratory_identifier = campaign%laboratory_identifier
    provenance%operator_identifier = campaign%operator_identifier
    provenance%instrument_identifier = campaign%instrument_identifier
    provenance%test_protocol_identifier = campaign%test_protocol_identifier
    provenance%calibration_reference = campaign%calibration_reference
    provenance%run_identifier = campaign%run_identifier
    provenance%test_date_metadata = campaign%test_date_metadata
    provenance%source_family_identifier = &
      campaign%identification%source_family%family_identifier
    provenance%source_metadata = &
      campaign%identification%source_family%common_state%source_metadata
    allocate(provenance%specimen_identifiers( &
      size(campaign%identification%source_family%isotherms)), &
      provenance%source_identifiers( &
        size(campaign%identification%source_family%isotherms)))
    do i = 1, size(provenance%specimen_identifiers)
      provenance%specimen_identifiers(i) = campaign%identification &
        %source_family%isotherms(i)%specimen_identifier
      provenance%source_identifiers(i) = campaign%identification &
        %source_family%isotherms(i)%source_identifier
    end do
  end subroutine copy_campaign_provenance

end module tms_tts_repeatability_analysis
