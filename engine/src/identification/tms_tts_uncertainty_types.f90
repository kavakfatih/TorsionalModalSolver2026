module tms_tts_uncertainty_types
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_text_length, &
    tts_pair_shift_configuration_t
  implicit none
  private

  integer, parameter, public :: TTS_UNCERTAINTY_SUCCESS = 0
  integer, parameter, public :: TTS_UNCERTAINTY_INVALID_INPUT = 1
  integer, parameter, public :: TTS_UNCERTAINTY_DATA_MISMATCH = 2
  integer, parameter, public :: TTS_UNCERTAINTY_NONFINITE_DATA = 3
  integer, parameter, public :: TTS_UNCERTAINTY_NO_SUPPORT = 4

  integer, parameter, public :: WEIGHTED_SHIFT_SUCCESS = 0
  integer, parameter, public :: WEIGHTED_SHIFT_STORAGE_ONLY = 1
  integer, parameter, public :: WEIGHTED_SHIFT_INVALID_INPUT = 2
  integer, parameter, public :: WEIGHTED_SHIFT_NO_OVERLAP = 3
  integer, parameter, public :: WEIGHTED_SHIFT_INSUFFICIENT_SUPPORT = 4
  integer, parameter, public :: WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM = 5
  integer, parameter, public :: WEIGHTED_SHIFT_NUMERICAL_FAILURE = 6

  integer, parameter, public :: HUBER_SHIFT_SUCCESS = 0
  integer, parameter, public :: HUBER_SHIFT_STORAGE_ONLY = 1
  integer, parameter, public :: HUBER_SHIFT_INVALID_INPUT = 2
  integer, parameter, public :: HUBER_SHIFT_NO_OVERLAP = 3
  integer, parameter, public :: HUBER_SHIFT_INSUFFICIENT_SUPPORT = 4
  integer, parameter, public :: HUBER_SHIFT_NO_INTERIOR_MINIMUM = 5
  integer, parameter, public :: HUBER_SHIFT_NUMERICAL_FAILURE = 6

  integer, parameter, public :: TTS_WEIGHTED_L2_OBJECTIVE = 1
  integer, parameter, public :: TTS_STANDARDIZED_HUBER_OBJECTIVE = 2

  integer, parameter, public :: UNCERTAINTY_SOURCE_UNSPECIFIED = 0
  integer, parameter, public :: UNCERTAINTY_SOURCE_TYPE_A = 1
  integer, parameter, public :: UNCERTAINTY_SOURCE_TYPE_B = 2
  integer, parameter, public :: UNCERTAINTY_SOURCE_COMBINED_STANDARD = 3
  integer, parameter, public :: UNCERTAINTY_SOURCE_REPEAT_MEASUREMENT = 4

  real(dp), parameter, public :: DEFAULT_TTS_HUBER_C = 1.345_dp

  !> Tek physical DMA noktasındaki standard uncertainty overlay kaydıdır.
  !! Anahtar (temperature [K], frequency [Hz]); storage/loss u_G değerleri
  !! [Pa] standard uncertainty'dir. Expanded uncertainty doğrudan saklanmaz.
  type, public :: tts_dynamic_modulus_uncertainty_point_t
    real(dp) :: temperature_k = 0.0_dp
    real(dp) :: frequency_hz = 0.0_dp
    real(dp) :: storage_standard_uncertainty_pa = 0.0_dp
    real(dp) :: loss_standard_uncertainty_pa = 0.0_dp
    logical :: storage_uncertainty_available = .false.
    logical :: loss_uncertainty_available = .false.
    integer :: uncertainty_source = UNCERTAINTY_SOURCE_UNSPECIFIED
    character(len=tts_text_length) :: source_metadata = ""
  end type tts_dynamic_modulus_uncertainty_point_t

  !> Aynı physical sıcaklıktaki pointwise uncertainty kayıtlarını taşır.
  !! Dizi sırası eşleme amacıyla kullanılmaz; (T,f) unique physical key'dir.
  type, public :: tts_dynamic_modulus_uncertainty_isotherm_t
    character(len=tts_text_length) :: isotherm_identifier = ""
    real(dp) :: temperature_k = 0.0_dp
    type(tts_dynamic_modulus_uncertainty_point_t), allocatable :: points(:)
  end type tts_dynamic_modulus_uncertainty_isotherm_t

  !> V0.8.1 authoritative family'yi değiştirmeyen additive uncertainty
  !! overlay ailesidir. Provenance, standard uncertainty kaynağını açıklar.
  type, public :: tts_dynamic_modulus_uncertainty_family_t
    character(len=tts_text_length) :: family_identifier = ""
    character(len=tts_text_length) :: provenance = ""
    type(tts_dynamic_modulus_uncertainty_isotherm_t), allocatable :: &
      isotherms(:)
  end type tts_dynamic_modulus_uncertainty_family_t

  !> y=log10(G) için first-order propagated standard uncertainty sonucudur.
  !! u_y=u_G/(G ln(10)); u_y boyutsuzdur.
  type, public :: tts_log_uncertainty_result_t
    integer :: status = TTS_UNCERTAINTY_INVALID_INPUT
    logical :: valid = .false.
    real(dp) :: standard_uncertainty = 0.0_dp
    real(dp) :: variance = 0.0_dp
  end type tts_log_uncertainty_result_t

  !> Uncertainty overlay ile V0.8.1 family arasındaki validation sonucudur.
  type, public :: tts_uncertainty_validation_result_t
    integer :: status = TTS_UNCERTAINTY_INVALID_INPUT
    logical :: valid = .false.
    character(len=tts_text_length) :: message = ""
    integer :: matched_point_count = 0
    integer :: storage_available_point_count = 0
    integer :: loss_available_point_count = 0
  end type tts_uncertainty_validation_result_t

  !> Contiguous weighted support segmentidir. x=log10(f/Hz), y=log10(G/Pa),
  !! variance=u_log10_G^2 boyutsuzdur. Source indeksleri provenance sağlar.
  type, public :: tts_uncertainty_log_segment_t
    integer :: channel = 0
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: variance(:)
    integer, allocatable :: source_point_indices(:)
    integer, allocatable :: uncertainty_point_indices(:)
  end type tts_uncertainty_log_segment_t

  !> Bir shift'teki weighted veya Huber objective ve support tanılarıdır.
  !! Objective boyutsuz, genişlikler log10-frequency decade cinsindedir.
  type, public :: tts_uncertainty_objective_evaluation_t
    logical :: valid = .false.
    logical :: numerical_failure = .false.
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: overlap_width_decades = 0.0_dp
    real(dp) :: overlap_fraction = 0.0_dp
    real(dp) :: quadratic_width_decades = 0.0_dp
    real(dp) :: tail_width_decades = 0.0_dp
    real(dp) :: tail_fraction = 0.0_dp
    real(dp) :: rms_standardized_residual = 0.0_dp
    integer :: interpolation_interval_count = 0
  end type tts_uncertainty_objective_evaluation_t

  !> Weighted/Huber sensitivity analizi numerical ayarlarıdır. Pair scan ve
  !! Brent toleransları V0.8.1 tipiyle aynıdır; Huber c boyutsuzdur.
  type, public :: tts_uncertainty_sensitivity_configuration_t
    type(tts_pair_shift_configuration_t) :: pair_shift
    real(dp) :: huber_c = DEFAULT_TTS_HUBER_C
  end type tts_uncertainty_sensitivity_configuration_t

  !> Tek objective modunun pair çözümünü ve channel diagnostics'ini taşır.
  type, public :: tts_uncertainty_pair_solution_t
    integer :: status = WEIGHTED_SHIFT_INVALID_INPUT
    integer :: objective_mode = 0
    integer :: production_channel = 0
    logical :: shift_available = .false.
    logical :: joint_shift_available = .false.
    logical :: storage_shift_available = .false.
    logical :: loss_shift_available = .false.
    real(dp) :: shift = 0.0_dp
    real(dp) :: joint_shift = 0.0_dp
    real(dp) :: storage_shift = 0.0_dp
    real(dp) :: loss_shift = 0.0_dp
    real(dp) :: objective_minimum = huge(1.0_dp)
    integer :: iteration_count = 0
    integer :: evaluation_count = 0
    type(tts_uncertainty_objective_evaluation_t) :: production_diagnostics
    type(tts_uncertainty_objective_evaluation_t) :: storage_diagnostics
    type(tts_uncertainty_objective_evaluation_t) :: loss_diagnostics
  end type tts_uncertainty_pair_solution_t

  !> Authoritative baseline ile additive weighted/Huber shift evidence'ını
  !! aynı pair kimliği altında karşılaştırır. Delta değerleri boyutsuz s'dir.
  type, public :: tts_uncertainty_shift_sensitivity_t
    integer :: reference_isotherm_index = 0
    integer :: moving_isotherm_index = 0
    character(len=tts_text_length) :: reference_isotherm_identifier = ""
    character(len=tts_text_length) :: moving_isotherm_identifier = ""
    real(dp) :: reference_temperature_k = 0.0_dp
    real(dp) :: moving_temperature_k = 0.0_dp
    integer :: baseline_status = 0
    logical :: baseline_shift_available = .false.
    real(dp) :: baseline_shift = 0.0_dp
    type(tts_uncertainty_pair_solution_t) :: weighted
    type(tts_uncertainty_pair_solution_t) :: huber
    logical :: weighted_delta_available = .false.
    logical :: huber_baseline_delta_available = .false.
    logical :: huber_weighted_delta_available = .false.
    real(dp) :: delta_weighted_vs_baseline = 0.0_dp
    real(dp) :: delta_huber_vs_baseline = 0.0_dp
    real(dp) :: delta_huber_vs_weighted = 0.0_dp
    real(dp) :: huber_c = DEFAULT_TTS_HUBER_C
    character(len=tts_text_length) :: uncertainty_provenance = ""
    logical :: cross_channel_covariance_modeled = .false.
    logical :: cross_isotherm_covariance_modeled = .false.
  end type tts_uncertainty_shift_sensitivity_t

  !> Bütün adjacent pair'ler için offline sensitivity sonucudur. Input
  !! identification ve authoritative runtime yapıları bu sonuçta kopyalanmaz
  !! veya değiştirilmez.
  type, public :: tts_uncertainty_sensitivity_result_t
    integer :: status = TTS_UNCERTAINTY_INVALID_INPUT
    character(len=tts_text_length) :: message = ""
    character(len=tts_text_length) :: uncertainty_family_identifier = ""
    character(len=tts_text_length) :: uncertainty_provenance = ""
    type(tts_uncertainty_sensitivity_configuration_t) :: configuration
    type(tts_uncertainty_shift_sensitivity_t), allocatable :: pair_results(:)
    integer :: baseline_available_count = 0
    integer :: weighted_available_count = 0
    integer :: huber_available_count = 0
    logical :: cross_channel_covariance_modeled = .false.
    logical :: cross_isotherm_covariance_modeled = .false.
  end type tts_uncertainty_sensitivity_result_t

end module tms_tts_uncertainty_types
