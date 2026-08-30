module tms_tts_covariance_types
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_text_length, tts_pair_shift_configuration_t
  use tms_tts_uncertainty_types, only : tts_uncertainty_pair_solution_t
  implicit none
  private

  integer, parameter, public :: TTS_COVARIANCE_SUCCESS = 0
  integer, parameter, public :: TTS_COVARIANCE_INVALID_INPUT = 1
  integer, parameter, public :: TTS_COVARIANCE_DATA_MISMATCH = 2
  integer, parameter, public :: TTS_COVARIANCE_NO_BIVARIATE_SUPPORT = 3
  integer, parameter, public :: TTS_COVARIANCE_INVALID_MATRIX = 4
  integer, parameter, public :: TTS_COVARIANCE_SINGULAR_MATRIX = 5
  integer, parameter, public :: TTS_COVARIANCE_ILL_CONDITIONED = 6
  integer, parameter, public :: TTS_COVARIANCE_NUMERICAL_FAILURE = 7

  integer, parameter, public :: MAHALANOBIS_SHIFT_SUCCESS = 0
  integer, parameter, public :: MAHALANOBIS_SHIFT_NO_OVERLAP = 1
  integer, parameter, public :: MAHALANOBIS_SHIFT_INVALID_COVARIANCE = 2
  integer, parameter, public :: MAHALANOBIS_SHIFT_NO_INTERIOR_MINIMUM = 3
  integer, parameter, public :: MAHALANOBIS_SHIFT_NUMERICAL_FAILURE = 4

  integer, parameter, public :: TTS_MATCHED_DIAGONAL_OBJECTIVE = 1
  integer, parameter, public :: TTS_MAHALANOBIS_OBJECTIVE = 2

  integer, parameter, public :: COVARIANCE_SOURCE_UNSPECIFIED = 0
  integer, parameter, public :: COVARIANCE_SOURCE_DIRECT = 1
  integer, parameter, public :: COVARIANCE_SOURCE_REPEAT_MEASUREMENT = 2
  integer, parameter, public :: &
    COVARIANCE_SOURCE_MAGNITUDE_PHASE_PROPAGATION = 3
  integer, parameter, public :: COVARIANCE_SOURCE_CALIBRATION_MODEL = 4

  !> Aynı fiziksel DMA noktasındaki 2x2 covariance matrisini taşır.
  !! Diyagonal alanlar Var(G') ve Var(G''), çapraz alan Cov(G',G'') olup
  !! tüm alanların SI birimi Pa^2'dir. Matris yapısı gereği simetriktir.
  type, public :: tts_covariance_matrix_2x2_t
    real(dp) :: storage_variance_pa2 = 0.0_dp
    real(dp) :: loss_variance_pa2 = 0.0_dp
    real(dp) :: storage_loss_covariance_pa2 = 0.0_dp
  end type tts_covariance_matrix_2x2_t

  !> Unique physical (T [K], f [Hz]) anahtarındaki point-local measurement
  !! covariance kaydıdır. `covariance_available=.false.` gerçek support gap'tir.
  type, public :: tts_dynamic_modulus_covariance_point_t
    real(dp) :: temperature_k = 0.0_dp
    real(dp) :: frequency_hz = 0.0_dp
    type(tts_covariance_matrix_2x2_t) :: covariance
    logical :: covariance_available = .false.
    integer :: covariance_source_kind = COVARIANCE_SOURCE_UNSPECIFIED
    character(len=tts_text_length) :: measurement_method = ""
    character(len=tts_text_length) :: instrument_identifier = ""
    character(len=tts_text_length) :: calibration_reference = ""
    character(len=tts_text_length) :: source_identifier = ""
    character(len=tts_text_length) :: source_description = ""
  end type tts_dynamic_modulus_covariance_point_t

  !> Aynı physical sıcaklıktaki point-local covariance kayıt koleksiyonudur.
  !! Dizi sırası eşleme semantiği taşımaz; her kayıt (T,f) ile bulunur.
  type, public :: tts_dynamic_modulus_covariance_isotherm_t
    character(len=tts_text_length) :: isotherm_identifier = ""
    real(dp) :: temperature_k = 0.0_dp
    type(tts_dynamic_modulus_covariance_point_t), allocatable :: points(:)
  end type tts_dynamic_modulus_covariance_isotherm_t

  !> V0.8.5 additive covariance overlay ailesidir. Canonical dış gösterim
  !! physical linear G uzayındaki tam 2x2 matristir; log covariance değildir.
  type, public :: tts_dynamic_modulus_covariance_family_t
    character(len=tts_text_length) :: family_identifier = ""
    character(len=tts_text_length) :: provenance = ""
    type(tts_dynamic_modulus_covariance_isotherm_t), allocatable :: &
      isotherms(:)
  end type tts_dynamic_modulus_covariance_family_t

  !> 2x2 SPD kontrolünün boyutsuz numerical conditioning tanılarıdır.
  !! reciprocal_condition_estimate yalnız makine güvenliği içindir; fiziksel
  !! covariance modelinin deneysel doğruluğunu kanıtlamaz.
  type, public :: tts_covariance_matrix_validation_t
    integer :: status = TTS_COVARIANCE_INVALID_MATRIX
    logical :: covariance_valid = .false.
    logical :: covariance_numerically_well_conditioned = .false.
    real(dp) :: determinant = 0.0_dp
    real(dp) :: correlation = 0.0_dp
    real(dp) :: reciprocal_condition_estimate = 0.0_dp
  end type tts_covariance_matrix_validation_t

  !> M [Pa] ve phase delta [rad] girdilerinden G', G'' [Pa] ve Sigma_G [Pa^2]
  !! üreten first-order Jacobian propagation sonucudur.
  type, public :: tts_polar_covariance_propagation_result_t
    integer :: status = TTS_COVARIANCE_INVALID_INPUT
    logical :: valid = .false.
    logical :: first_order_covariance_propagation = .true.
    real(dp) :: storage_modulus_pa = 0.0_dp
    real(dp) :: loss_modulus_pa = 0.0_dp
    type(tts_covariance_matrix_2x2_t) :: covariance
  end type tts_polar_covariance_propagation_result_t

  !> Sigma_y=D*Sigma_G*D^T first-order dönüşüm sonucudur. Log10-modulus
  !! covariance alanları boyutsuzdur; G' ve G'' kesinlikle pozitif olmalıdır.
  type, public :: tts_log_covariance_propagation_result_t
    integer :: status = TTS_COVARIANCE_INVALID_INPUT
    logical :: valid = .false.
    logical :: first_order_covariance_propagation = .true.
    real(dp) :: storage_variance = 0.0_dp
    real(dp) :: loss_variance = 0.0_dp
    real(dp) :: storage_loss_covariance = 0.0_dp
    real(dp) :: correlation = 0.0_dp
  end type tts_log_covariance_propagation_result_t

  !> Covariance overlay, measurement family ve V0.8.4 uncertainty
  !! diyagonallerinin birlikte doğrulanma sonucudur.
  type, public :: tts_covariance_validation_result_t
    integer :: status = TTS_COVARIANCE_INVALID_INPUT
    logical :: valid = .false.
    character(len=tts_text_length) :: message = ""
    integer :: matched_point_count = 0
    integer :: covariance_available_point_count = 0
    integer :: covariance_gap_count = 0
    integer :: ill_conditioned_covariance_count = 0
  end type tts_covariance_validation_result_t

  !> Tek contiguous bivariate support segmentidir. x=log10(f/Hz), y alanları
  !! log10(G/Pa), covariance alanları boyutsuz log-modulus covariance'dır.
  type, public :: tts_bivariate_covariance_log_segment_t
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: storage_y(:)
    real(dp), allocatable :: loss_y(:)
    real(dp), allocatable :: storage_variance(:)
    real(dp), allocatable :: loss_variance(:)
    real(dp), allocatable :: storage_loss_covariance(:)
    integer, allocatable :: source_point_indices(:)
    integer, allocatable :: covariance_point_indices(:)
  end type tts_bivariate_covariance_log_segment_t

  !> Aynı O_B üzerinde matched-diagonal ve Mahalanobis objective tanılarıdır.
  !! Objective boyutsuz, overlap genişliği log10-frequency decade'dir.
  type, public :: tts_covariance_objective_evaluation_t
    logical :: valid = .false.
    logical :: numerical_failure = .false.
    real(dp) :: matched_diagonal_objective = huge(1.0_dp)
    real(dp) :: mahalanobis_objective = huge(1.0_dp)
    real(dp) :: overlap_width_decades = 0.0_dp
    integer :: interpolation_interval_count = 0
  end type tts_covariance_objective_evaluation_t

  !> Bir objective modu için deterministic coarse-scan/Brent pair sonucudur.
  type, public :: tts_covariance_shift_solution_t
    integer :: status = MAHALANOBIS_SHIFT_INVALID_COVARIANCE
    integer :: objective_mode = 0
    logical :: shift_available = .false.
    real(dp) :: shift = 0.0_dp
    real(dp) :: objective_minimum = huge(1.0_dp)
    integer :: iteration_count = 0
    integer :: evaluation_count = 0
    type(tts_covariance_objective_evaluation_t) :: diagnostics
  end type tts_covariance_shift_solution_t

  !> Matched-support diagonal ile Mahalanobis çözümlerini aynı O_B tanımı
  !! altında birlikte taşır; storage-only fallback alanı kasıtlı olarak yoktur.
  type, public :: tts_covariance_pair_solution_t
    type(tts_covariance_shift_solution_t) :: diagonal_matched
    type(tts_covariance_shift_solution_t) :: mahalanobis
  end type tts_covariance_pair_solution_t

  !> V0.8.5 offline sensitivity optimizer ayarlarıdır. V0.8.1'in boyutsuz
  !! s=log10(a_T) coarse scan ve Brent toleransları aynen reuse edilir.
  type, public :: tts_covariance_sensitivity_configuration_t
    type(tts_pair_shift_configuration_t) :: pair_shift
  end type tts_covariance_sensitivity_configuration_t

  !> Dört ayrı shift evidence'ını ve destek/covariance decomposition'ını taşır.
  !! Hiçbir alan preferred veya authoritative shift seçimi değildir.
  type, public :: tts_covariance_pair_sensitivity_t
    integer :: reference_isotherm_index = 0
    integer :: moving_isotherm_index = 0
    character(len=tts_text_length) :: reference_isotherm_identifier = ""
    character(len=tts_text_length) :: moving_isotherm_identifier = ""
    real(dp) :: reference_temperature_k = 0.0_dp
    real(dp) :: moving_temperature_k = 0.0_dp
    integer :: baseline_status = 0
    logical :: baseline_shift_available = .false.
    real(dp) :: baseline_shift = 0.0_dp
    type(tts_uncertainty_pair_solution_t) :: weighted_original
    type(tts_covariance_shift_solution_t) :: diagonal_matched
    type(tts_covariance_shift_solution_t) :: mahalanobis
    logical :: delta_support_available = .false.
    logical :: delta_covariance_available = .false.
    logical :: delta_total_available = .false.
    logical :: delta_mahalanobis_vs_baseline_available = .false.
    real(dp) :: delta_support = 0.0_dp
    real(dp) :: delta_covariance = 0.0_dp
    real(dp) :: delta_total = 0.0_dp
    real(dp) :: delta_mahalanobis_vs_baseline = 0.0_dp
    real(dp) :: bivariate_overlap_width_decades = 0.0_dp
    integer :: covariance_point_count = 0
    integer :: covariance_gap_count = 0
    integer :: ill_conditioned_covariance_count = 0
    logical :: first_order_covariance_propagation = .true.
    logical :: cross_channel_covariance_modeled = .true.
    logical :: cross_isotherm_covariance_modeled = .false.
    logical :: cross_frequency_covariance_modeled = .false.
    logical :: automatic_regularization_used = .false.
  end type tts_covariance_pair_sensitivity_t

  !> Bütün adjacent pair'ler için additive V0.8.5 sensitivity sonucudur.
  !! Identification ve V0.8.4 uncertainty girdileri kopyalanmaz/değiştirilmez.
  type, public :: tts_covariance_sensitivity_result_t
    integer :: status = TTS_COVARIANCE_INVALID_INPUT
    character(len=tts_text_length) :: message = ""
    character(len=tts_text_length) :: covariance_family_identifier = ""
    character(len=tts_text_length) :: covariance_provenance = ""
    type(tts_covariance_sensitivity_configuration_t) :: configuration
    type(tts_covariance_pair_sensitivity_t), allocatable :: pair_results(:)
    integer :: pair_count = 0
    integer :: mahalanobis_available_count = 0
    integer :: matched_diagonal_available_count = 0
    integer :: unsupported_pair_count = 0
    logical :: first_order_covariance_propagation = .true.
    logical :: cross_channel_covariance_modeled = .true.
    logical :: cross_isotherm_covariance_modeled = .false.
    logical :: cross_frequency_covariance_modeled = .false.
    logical :: automatic_regularization_used = .false.
  end type tts_covariance_sensitivity_result_t

end module tms_tts_covariance_types
