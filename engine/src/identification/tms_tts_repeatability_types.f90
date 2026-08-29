module tms_tts_repeatability_types
  use tms_kinds, only : dp
  use tms_sample_statistics, only : sample_statistics_t
  use tms_bootstrap, only : bootstrap_configuration_t, bootstrap_interval_t
  use tms_tts_types, only : tts_identification_result_t, tts_text_length
  implicit none
  private

  integer, parameter, public :: REPLICATE_BASIS_UNSPECIFIED = 0
  integer, parameter, public :: INDEPENDENT_SPECIMEN_CAMPAIGN = 1
  integer, parameter, public :: SAME_SPECIMEN_RERUN = 2

  integer, parameter, public :: TTS_REPEATABILITY_SUCCESS = 0
  integer, parameter, public :: TTS_REPEATABILITY_INVALID_INPUT = 1
  integer, parameter, public :: TTS_REPEATABILITY_INCOMPATIBLE_STATE = 2
  integer, parameter, public :: TTS_REPEATABILITY_TEMPERATURE_SET_MISMATCH = 3
  integer, parameter, public :: TTS_REPEATABILITY_REFERENCE_NOT_FOUND = 4
  integer, parameter, public :: TTS_REPEATABILITY_INSUFFICIENT_REPLICATES = 5
  integer, parameter, public :: TTS_REPEATABILITY_NO_INDEPENDENT_CAMPAIGNS = 6
  integer, parameter, public :: TTS_REPEATABILITY_BOOTSTRAP_UNAVAILABLE = 7
  integer, parameter, public :: TTS_REPEATABILITY_NONFINITE_DATA = 8

  !> Tek complete V0.8.1 DMA/TTS campaign'i ve experimental independence
  !! provenance'ını taşır. Statistical sample unit bu wrapper'dır; içindeki
  !! frequency points, isotherms veya adjacent pairs replicate değildir.
  !! identification authoritative kalır ve analiz API'sine intent(in) verilir.
  type, public :: tts_repeatability_campaign_t
    character(len=tts_text_length) :: campaign_identifier = ""
    integer :: replicate_basis = REPLICATE_BASIS_UNSPECIFIED
    character(len=tts_text_length) :: laboratory_identifier = ""
    character(len=tts_text_length) :: operator_identifier = ""
    character(len=tts_text_length) :: instrument_identifier = ""
    character(len=tts_text_length) :: test_protocol_identifier = ""
    character(len=tts_text_length) :: calibration_reference = ""
    character(len=tts_text_length) :: run_identifier = ""
    character(len=tts_text_length) :: test_date_metadata = ""
    type(tts_identification_result_t) :: identification
  end type tts_repeatability_campaign_t

  !> Study policy ve bootstrap numerical ayarlarıdır. Same-specimen rerun'lar
  !! default descriptive evidence'a dahil edilir fakat hiçbir durumda
  !! independent cluster bootstrap population'ına girmez.
  type, public :: tts_repeatability_study_configuration_t
    character(len=tts_text_length) :: study_identifier = ""
    logical :: include_same_specimen_reruns_in_descriptive = .true.
    type(bootstrap_configuration_t) :: bootstrap
  end type tts_repeatability_study_configuration_t

  !> Campaign ve kaynak izini, full identification sonucunu kopyalamadan
  !! saklar. Specimen/source kimlikleri physical-state equality gate değildir;
  !! farklı specimen bağımsız replicate olabilir.
  type, public :: tts_repeatability_campaign_provenance_t
    character(len=tts_text_length) :: campaign_identifier = ""
    integer :: replicate_basis = REPLICATE_BASIS_UNSPECIFIED
    logical :: included_in_descriptive_statistics = .false.
    logical :: included_in_independent_bootstrap_population = .false.
    character(len=tts_text_length) :: laboratory_identifier = ""
    character(len=tts_text_length) :: operator_identifier = ""
    character(len=tts_text_length) :: instrument_identifier = ""
    character(len=tts_text_length) :: test_protocol_identifier = ""
    character(len=tts_text_length) :: calibration_reference = ""
    character(len=tts_text_length) :: run_identifier = ""
    character(len=tts_text_length) :: test_date_metadata = ""
    character(len=tts_text_length) :: source_family_identifier = ""
    character(len=tts_text_length) :: source_metadata = ""
    character(len=tts_text_length), allocatable :: specimen_identifiers(:)
    character(len=tts_text_length), allocatable :: source_identifiers(:)
  end type tts_repeatability_campaign_provenance_t

  !> Canonical T_low<T_high adjacent pair için delta_s=s(T_high)-s(T_low)
  !! repeatability evidence'ıdır. Shift ve bütün statistics boyutsuzdur.
  type, public :: tts_pair_repeatability_result_t
    real(dp) :: lower_temperature_k = 0.0_dp
    real(dp) :: upper_temperature_k = 0.0_dp
    type(sample_statistics_t) :: delta_s_statistics
    type(bootstrap_interval_t) :: mean_bootstrap_interval
  end type tts_pair_repeatability_result_t

  !> Common-reference normalized empirical s(T)=log10(a_T) marginal
  !! repeatability evidence'ıdır. Reference noktasındaki sıfır yapısal ankordur;
  !! uncertainty_informative=false olması zero physical uncertainty yorumunu
  !! engeller. Shift ve statistics boyutsuz, sıcaklık [K]'dir.
  type, public :: tts_absolute_shift_repeatability_result_t
    real(dp) :: temperature_k = 0.0_dp
    logical :: is_reference_anchor = .false.
    logical :: uncertainty_informative = .true.
    type(sample_statistics_t) :: log10_a_t_statistics
    type(bootstrap_interval_t) :: mean_bootstrap_interval
  end type tts_absolute_shift_repeatability_result_t

  !> Campaign bazlı V0.8.2 Arrhenius fit cohort evidence'ıdır. beta [K],
  !! Ea_app [J/mol] ve statistics/CI aynı birimlerdedir. Invalid fit hiçbir
  !! zaman zero placeholder ile population'a eklenmez.
  type, public :: tts_arrhenius_repeatability_result_t
    integer :: total_campaign_count = 0
    integer :: fit_available_count = 0
    integer :: fit_unavailable_count = 0
    integer :: invalid_or_nonphysical_fit_count = 0
    type(sample_statistics_t) :: beta_k_statistics
    type(sample_statistics_t) :: apparent_activation_energy_statistics
    type(bootstrap_interval_t) :: beta_k_mean_bootstrap_interval
    type(bootstrap_interval_t) :: &
      apparent_activation_energy_mean_bootstrap_interval
  end type tts_arrhenius_repeatability_result_t

  !> Campaign bazlı V0.8.2 WLF cohort evidence'ıdır. C1 boyutsuz, C2 [K],
  !! p=C1/C2 [1/K], q=1/C2 [1/K]'dir. Parameter statistics yalnız fit
  !! available ve identifiable campaigns'ten gelir; poorly identified sayısı
  !! engineering evidence olarak ayrıca korunur.
  type, public :: tts_wlf_repeatability_result_t
    integer :: total_campaign_count = 0
    integer :: fit_available_count = 0
    integer :: fit_unavailable_count = 0
    integer :: parameter_identifiable_count = 0
    integer :: poorly_identified_count = 0
    integer :: invalid_fit_count = 0
    type(sample_statistics_t) :: c1_statistics
    type(sample_statistics_t) :: c2_k_statistics
    type(sample_statistics_t) :: p_c1_over_c2_per_k_statistics
    type(sample_statistics_t) :: q_inverse_c2_per_k_statistics
    type(bootstrap_interval_t) :: c1_mean_bootstrap_interval
    type(bootstrap_interval_t) :: c2_k_mean_bootstrap_interval
    type(bootstrap_interval_t) :: p_mean_bootstrap_interval
    type(bootstrap_interval_t) :: q_mean_bootstrap_interval
  end type tts_wlf_repeatability_result_t

  !> V0.8.3 additive offline repeatability/uncertainty sonucudur. V0.8.1
  !! empirical ve V0.8.2 parametric sonuçları external authoritative inputs
  !! olarak kalır. Availability ve independence sayıları, descriptive evidence
  !! ile independent cluster bootstrap iddiasını birbirinden ayırır.
  type, public :: tts_repeatability_study_result_t
    integer :: status = TTS_REPEATABILITY_INVALID_INPUT
    character(len=tts_text_length) :: message = ""
    character(len=tts_text_length) :: study_identifier = ""
    logical :: descriptive_statistics_available = .false.
    logical :: independent_cluster_bootstrap_available = .false.
    logical :: intralaboratory_context_explicit = .false.
    integer :: total_campaign_count = 0
    integer :: descriptive_campaign_count = 0
    integer :: independent_campaign_count = 0
    integer :: same_specimen_rerun_count = 0
    integer :: unspecified_replicate_basis_count = 0
    real(dp) :: common_reference_temperature_k = 0.0_dp
    integer :: bootstrap_status = TTS_REPEATABILITY_BOOTSTRAP_UNAVAILABLE
    character(len=tts_text_length) :: bootstrap_message = ""
    type(bootstrap_configuration_t) :: bootstrap_configuration
    character(len=tts_text_length), allocatable :: &
      bootstrap_population_campaign_identifiers(:)
    type(tts_repeatability_campaign_provenance_t), allocatable :: &
      campaign_provenance(:)
    type(tts_pair_repeatability_result_t), allocatable :: pair_results(:)
    type(tts_absolute_shift_repeatability_result_t), allocatable :: &
      absolute_shift_results(:)
    type(tts_arrhenius_repeatability_result_t) :: arrhenius
    type(tts_wlf_repeatability_result_t) :: wlf
  end type tts_repeatability_study_result_t

end module tms_tts_repeatability_types
