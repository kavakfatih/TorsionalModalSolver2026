module tms_tts_shift_law_types
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_text_length
  implicit none
  private

  integer, parameter, public :: SHIFT_LAW_FIT_SUCCESS = 0
  integer, parameter, public :: SHIFT_LAW_FIT_INVALID_INPUT = 1
  integer, parameter, public :: SHIFT_LAW_FIT_INSUFFICIENT_DATA = 2
  integer, parameter, public :: SHIFT_LAW_FIT_NONFINITE_RESULT = 3
  integer, parameter, public :: ARRHENIUS_FIT_INVALID_SLOPE = 10
  integer, parameter, public :: WLF_FIT_NO_VALID_DOMAIN = 20
  integer, parameter, public :: WLF_FIT_NO_INTERIOR_BRACKET = 21
  integer, parameter, public :: WLF_FIT_OPTIMIZATION_FAILED = 22
  integer, parameter, public :: WLF_FIT_POORLY_IDENTIFIED = 23
  integer, parameter, public :: WLF_FIT_INVALID_PARAMETERS = 24

  integer, parameter, public :: DEFAULT_WLF_MAXIMUM_BRACKET_EXPANSIONS = 64
  real(dp), parameter, public :: DEFAULT_WLF_BRACKET_EXPANSION_FACTOR = &
    2.0_dp

  !> Adjacent iki measured isotherm arasındaki primary fit gözlemidir.
  !! Matematiksel tanım delta_s=s(T_j)-s(T_i), s=log10(a_T)'dir. Sıcaklıklar
  !! [K], delta_s boyutsuzdur. Cumulative absolute shift bu tipe girmez;
  !! source indeksleri yalnız provenance içindir.
  type, public :: tts_pair_shift_observation_t
    integer :: isotherm_i_index = 0
    integer :: isotherm_j_index = 0
    real(dp) :: temperature_i_k = 0.0_dp
    real(dp) :: temperature_j_k = 0.0_dp
    real(dp) :: delta_s_j_minus_i = 0.0_dp
  end type tts_pair_shift_observation_t

  !> WLF profile objective için deterministic numerical arama ayarlarıdır.
  !! Expansion factor ve count yalnız interior C2 bracket bulma algoritmasını
  !! sınırlar; material kabul aralığı veya fiziksel C2 limiti değildir.
  type, public :: tts_wlf_fit_configuration_t
    real(dp) :: bracket_expansion_factor = &
      DEFAULT_WLF_BRACKET_EXPANSION_FACTOR
    integer :: maximum_bracket_expansions = &
      DEFAULT_WLF_MAXIMUM_BRACKET_EXPANSIONS
    real(dp) :: absolute_tolerance_k = 8.0_dp*sqrt(epsilon(1.0_dp))
    real(dp) :: relative_tolerance = 8.0_dp*sqrt(epsilon(1.0_dp))
    integer :: maximum_iterations = 200
  end type tts_wlf_fit_configuration_t

  !> Bir Leave-One-Temperature-Out fold'unun predictive tanısını taşır.
  !! Held-out sıcaklık [K], residual değerleri boyutsuz delta_s ölçeğindedir.
  !! available yalnız training fit ve held-out prediction üretilebildiğini
  !! belirtir; material-law kabul kararı değildir.
  type, public :: tts_loto_fold_diagnostic_t
    real(dp) :: held_out_temperature_k = 0.0_dp
    integer :: training_observation_count = 0
    integer :: validation_observation_count = 0
    logical :: available = .false.
    real(dp) :: pair_rmse = 0.0_dp
    real(dp) :: pair_max_abs_residual = 0.0_dp
    real(dp) :: pair_mean_residual = 0.0_dp
  end type tts_loto_fold_diagnostic_t

  !> Equal-weight adjacent-pair Arrhenius fit sonucudur. beta [K], apparent
  !! activation energy [J/mol], sıcaklıklar [K], shift/residual değerleri
  !! boyutsuzdur. fit_available matematiksel parametre çözümünü;
  !! validation flag'leri ise ayrı evidence seviyelerini ifade eder.
  type, public :: tts_arrhenius_fit_result_t
    integer :: status = SHIFT_LAW_FIT_INVALID_INPUT
    character(len=tts_text_length) :: message = ""
    logical :: fit_available = .false.
    logical :: residual_validation_available = .false.
    logical :: predictive_validation_available = .false.
    real(dp) :: apparent_activation_energy_j_per_mol = 0.0_dp
    real(dp) :: beta_k = 0.0_dp
    real(dp) :: reference_temperature_k = 0.0_dp
    real(dp) :: minimum_temperature_k = 0.0_dp
    real(dp) :: maximum_temperature_k = 0.0_dp
    real(dp) :: pair_rmse = 0.0_dp
    real(dp) :: pair_max_abs_residual = 0.0_dp
    real(dp) :: pair_mean_residual = 0.0_dp
    real(dp) :: loto_rmse = 0.0_dp
    real(dp) :: loto_max_abs_residual = 0.0_dp
    integer :: observation_count = 0
    integer :: validation_degree_of_freedom = 0
    real(dp), allocatable :: predicted_pair_shifts(:)
    real(dp), allocatable :: pair_residuals(:)
    type(tts_loto_fold_diagnostic_t), allocatable :: loto_diagnostics(:)
  end type tts_arrhenius_fit_result_t

  !> Equal-weight profiled WLF fit sonucudur. C1 ve p=C1/C2 boyutsuz/K
  !! ilişkisini, C2 [K] ve q=1/C2 [1/K] tanılarını taşır. Low residual ile
  !! parameter identifiability ayrı tutulur; poorly identified sonuç runtime
  !! WLF provider'a otomatik açılmaz.
  type, public :: tts_wlf_fit_result_t
    integer :: status = SHIFT_LAW_FIT_INVALID_INPUT
    character(len=tts_text_length) :: message = ""
    logical :: fit_available = .false.
    logical :: parameter_identifiable = .false.
    logical :: residual_validation_available = .false.
    logical :: predictive_validation_available = .false.
    real(dp) :: reference_temperature_k = 0.0_dp
    real(dp) :: minimum_temperature_k = 0.0_dp
    real(dp) :: maximum_temperature_k = 0.0_dp
    real(dp) :: c1 = 0.0_dp
    real(dp) :: c2_k = 0.0_dp
    real(dp) :: p_c1_over_c2_per_k = 0.0_dp
    real(dp) :: q_inverse_c2_per_k = 0.0_dp
    real(dp) :: pair_rmse = 0.0_dp
    real(dp) :: pair_max_abs_residual = 0.0_dp
    real(dp) :: pair_mean_residual = 0.0_dp
    real(dp) :: loto_rmse = 0.0_dp
    real(dp) :: loto_max_abs_residual = 0.0_dp
    real(dp) :: profile_objective_minimum = 0.0_dp
    real(dp) :: bracket_lower_c2_k = 0.0_dp
    real(dp) :: bracket_middle_c2_k = 0.0_dp
    real(dp) :: bracket_upper_c2_k = 0.0_dp
    integer :: observation_count = 0
    integer :: validation_degree_of_freedom = 0
    integer :: bracket_expansion_count = 0
    integer :: minimizer_iteration_count = 0
    integer :: minimizer_evaluation_count = 0
    real(dp), allocatable :: predicted_pair_shifts(:)
    real(dp), allocatable :: pair_residuals(:)
    type(tts_loto_fold_diagnostic_t), allocatable :: loto_diagnostics(:)
  end type tts_wlf_fit_result_t

  !> V0.8.1 empirical identification'dan türetilen iki bağımsız parametrik
  !! approximation'ı birlikte taşır. Automatic best-model alanı bilinçli
  !! olarak yoktur; empirical shift table bu sonucun dışında authoritative
  !! kalır.
  type, public :: tts_shift_law_identification_result_t
    integer :: status = SHIFT_LAW_FIT_INVALID_INPUT
    character(len=tts_text_length) :: message = ""
    logical :: pair_observations_available = .false.
    integer :: pair_observation_count = 0
    type(tts_arrhenius_fit_result_t) :: arrhenius
    type(tts_wlf_fit_result_t) :: wlf
  end type tts_shift_law_identification_result_t

end module tms_tts_shift_law_types
