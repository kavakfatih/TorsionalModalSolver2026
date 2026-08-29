module tms_tts_shift_law_fit
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : universal_gas_constant_j_per_mol_k
  use tms_scalar_minimizer, only : scalar_minimizer_result_t, &
    SCALAR_MINIMIZER_SUCCESS, minimize_scalar_brent, &
    is_valid_scalar_minimum_bracket
  use tms_tts_shift_law_types, only : tts_pair_shift_observation_t, &
    tts_wlf_fit_configuration_t, tts_arrhenius_fit_result_t, &
    tts_wlf_fit_result_t, SHIFT_LAW_FIT_SUCCESS, &
    SHIFT_LAW_FIT_INVALID_INPUT, SHIFT_LAW_FIT_INSUFFICIENT_DATA, &
    SHIFT_LAW_FIT_NONFINITE_RESULT, ARRHENIUS_FIT_INVALID_SLOPE, &
    WLF_FIT_NO_VALID_DOMAIN, WLF_FIT_NO_INTERIOR_BRACKET, &
    WLF_FIT_OPTIMIZATION_FAILED, WLF_FIT_POORLY_IDENTIFIED, &
    WLF_FIT_INVALID_PARAMETERS
  implicit none
  private

  type :: wlf_profile_evaluation_t
    logical :: valid = .false.
    real(dp) :: c1 = 0.0_dp
    real(dp) :: objective = 0.0_dp
  end type wlf_profile_evaluation_t

  public :: fit_tts_arrhenius_shift_law
  public :: fit_tts_wlf_shift_law
  public :: predict_tts_arrhenius_pair_shift
  public :: predict_tts_wlf_pair_shift
  public :: validate_tts_pair_shift_observations

contains

  !> Adjacent pair gözlemlerinden Arrhenius beta katsayısını equal-weight
  !! analytical least-squares ile hesaplar:
  !! beta=sum(x*delta_s)/sum(x^2), x=1/T_j-1/T_i ve
  !! Ea_app=R*ln(10)*beta. T [K], beta [K], Ea_app [J/mol] ve delta_s
  !! boyutsuzdur. Cumulative absolute shift değerleri kullanılmaz.
  pure function fit_tts_arrhenius_shift_law( &
      observations, reference_temperature_k) result(fit)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: reference_temperature_k
    type(tts_arrhenius_fit_result_t) :: fit

    real(dp), allocatable :: scaled_x(:)
    real(dp), allocatable :: x(:)
    real(dp) :: denominator
    real(dp) :: maximum_abs_x
    real(dp) :: numerator
    integer :: i
    integer :: validation_status
    character(len=256) :: validation_message

    call validate_tts_pair_shift_observations(observations, &
      reference_temperature_k, 1, fit%minimum_temperature_k, &
      fit%maximum_temperature_k, validation_status, validation_message)
    fit%reference_temperature_k = reference_temperature_k
    fit%observation_count = size(observations)
    fit%validation_degree_of_freedom = max(0, size(observations) - 1)
    if (validation_status /= SHIFT_LAW_FIT_SUCCESS) then
      fit%status = validation_status
      fit%message = validation_message
      return
    end if

    allocate(x(size(observations)))
    do i = 1, size(observations)
      x(i) = 1.0_dp/observations(i)%temperature_j_k - &
        1.0_dp/observations(i)%temperature_i_k
    end do
    maximum_abs_x = maxval(abs(x))
    if (.not. ieee_is_finite(maximum_abs_x) .or. &
        maximum_abs_x <= tiny(1.0_dp)) then
      fit%status = SHIFT_LAW_FIT_INSUFFICIENT_DATA
      fit%message = "Arrhenius sum(x^2) için nonzero temperature farkı yok."
      return
    end if

    ! x/max|x| normalizasyonu sum(x^2) underflow riskini azaltır; çözüm
    ! cebirsel olarak aynı equal-weight slope değeridir.
    scaled_x = x/maximum_abs_x
    denominator = sum(scaled_x*scaled_x)
    numerator = sum(scaled_x*observations%delta_s_j_minus_i)
    if (.not. ieee_is_finite(denominator) .or. denominator <= 0.0_dp .or. &
        .not. ieee_is_finite(numerator)) then
      fit%status = SHIFT_LAW_FIT_NONFINITE_RESULT
      fit%message = "Arrhenius normalized least-squares toplamı geçersiz."
      return
    end if

    fit%beta_k = numerator/(denominator*maximum_abs_x)
    fit%apparent_activation_energy_j_per_mol = &
      universal_gas_constant_j_per_mol_k*log(10.0_dp)*fit%beta_k
    if (.not. ieee_is_finite(fit%beta_k) .or. &
        .not. ieee_is_finite(fit%apparent_activation_energy_j_per_mol)) then
      fit%status = SHIFT_LAW_FIT_NONFINITE_RESULT
      fit%message = "Arrhenius beta veya Ea_app sonucu sonlu değil."
      return
    end if
    if (fit%beta_k <= 0.0_dp .or. &
        fit%apparent_activation_energy_j_per_mol <= 0.0_dp) then
      fit%status = ARRHENIUS_FIT_INVALID_SLOPE
      fit%message = "Conventional Arrhenius fit pozitif Ea_app üretmedi."
      return
    end if

    allocate(fit%predicted_pair_shifts(size(observations)))
    do i = 1, size(observations)
      fit%predicted_pair_shifts(i) = predict_tts_arrhenius_pair_shift( &
        fit%beta_k, observations(i))
    end do
    fit%pair_residuals = observations%delta_s_j_minus_i - &
      fit%predicted_pair_shifts
    call calculate_residual_metrics(fit%pair_residuals, fit%pair_rmse, &
      fit%pair_max_abs_residual, fit%pair_mean_residual)
    if (.not. residual_metrics_are_finite(fit%pair_rmse, &
        fit%pair_max_abs_residual, fit%pair_mean_residual)) then
      fit%status = SHIFT_LAW_FIT_NONFINITE_RESULT
      fit%message = "Arrhenius pair residual metrics sonlu değil."
      return
    end if

    fit%fit_available = .true.
    fit%residual_validation_available = &
      fit%validation_degree_of_freedom > 0
    fit%status = SHIFT_LAW_FIT_SUCCESS
    fit%message = "Adjacent-pair Arrhenius Ea_app fit hazır."
  end function fit_tts_arrhenius_shift_law

  !> Fixed C2 [K] için C1'i analytical profile ederek yalnız C2 üzerinde
  !! bounded Brent minimizasyonu yapar. Model
  !! g(T,C2)=-(T-T_ref)/(C2+T-T_ref), delta_s=C1*(g_j-g_i)'dir.
  !! Search lower bound measured Tmin'den pole-safe türetilir; upper bound
  !! temperature span ile başlayıp multiplicatively genişler. Equal weights
  !! kullanılır; covariance veya objective-curvature statistical weight değildir.
  pure function fit_tts_wlf_shift_law( &
      observations, reference_temperature_k, configuration) result(fit)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: reference_temperature_k
    type(tts_wlf_fit_configuration_t), intent(in), optional :: configuration
    type(tts_wlf_fit_result_t) :: fit

    type(tts_wlf_fit_configuration_t) :: settings
    type(wlf_profile_evaluation_t) :: lower_profile
    type(wlf_profile_evaluation_t) :: middle_profile
    type(wlf_profile_evaluation_t) :: upper_profile
    type(wlf_profile_evaluation_t) :: optimum_profile
    type(scalar_minimizer_result_t) :: minimizer
    real(dp) :: c2_lower
    real(dp) :: c2_middle
    real(dp) :: c2_upper
    real(dp) :: domain_scale_k
    real(dp) :: improvement_tolerance
    real(dp) :: linear_objective
    real(dp) :: linear_p
    real(dp) :: pole_distance_k
    real(dp) :: pole_margin_k
    real(dp) :: shift_scale
    logical :: bracket_found
    integer :: expansion
    integer :: validation_status
    character(len=256) :: validation_message

    if (present(configuration)) settings = configuration
    fit%reference_temperature_k = reference_temperature_k
    fit%observation_count = size(observations)
    fit%validation_degree_of_freedom = max(0, size(observations) - 2)
    call validate_tts_pair_shift_observations(observations, &
      reference_temperature_k, 2, fit%minimum_temperature_k, &
      fit%maximum_temperature_k, validation_status, validation_message)
    if (validation_status /= SHIFT_LAW_FIT_SUCCESS) then
      fit%status = validation_status
      fit%message = validation_message
      return
    end if
    if (.not. valid_wlf_configuration(settings)) then
      fit%status = SHIFT_LAW_FIT_INVALID_INPUT
      fit%message = "WLF numerical fit configuration geçersiz."
      return
    end if

    call evaluate_large_c2_limit(observations, linear_p, &
      linear_objective)
    if (.not. ieee_is_finite(linear_p) .or. linear_p <= 0.0_dp .or. &
        .not. ieee_is_finite(linear_objective)) then
      fit%status = WLF_FIT_INVALID_PARAMETERS
      fit%message = "WLF pair yönü positive C1/C2 limitiyle uyumlu değil."
      return
    end if

    pole_distance_k = max(0.0_dp, reference_temperature_k - &
      fit%minimum_temperature_k)
    domain_scale_k = fit%maximum_temperature_k - &
      fit%minimum_temperature_k
    if (.not. ieee_is_finite(domain_scale_k) .or. domain_scale_k <= 0.0_dp) then
      fit%status = WLF_FIT_NO_VALID_DOMAIN
      fit%message = "WLF calibrated temperature span pozitif değil."
      return
    end if
    pole_margin_k = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
      abs(reference_temperature_k), abs(fit%minimum_temperature_k), &
      abs(fit%maximum_temperature_k), domain_scale_k, pole_distance_k)
    c2_lower = pole_distance_k + pole_margin_k
    c2_middle = pole_distance_k + domain_scale_k
    if (c2_middle <= c2_lower) c2_middle = c2_lower + domain_scale_k
    c2_upper = pole_distance_k + &
      settings%bracket_expansion_factor*(c2_middle - pole_distance_k)

    lower_profile = evaluate_wlf_profile(observations, &
      reference_temperature_k, c2_lower)
    middle_profile = evaluate_wlf_profile(observations, &
      reference_temperature_k, c2_middle)
    upper_profile = evaluate_wlf_profile(observations, &
      reference_temperature_k, c2_upper)
    if (.not. lower_profile%valid .or. .not. middle_profile%valid .or. &
        .not. upper_profile%valid) then
      fit%status = WLF_FIT_NO_VALID_DOMAIN
      fit%message = "Pole-safe WLF profile başlangıç domain'i kurulamadı."
      return
    end if

    bracket_found = .false.
    do expansion = 0, settings%maximum_bracket_expansions
      if (is_valid_scalar_minimum_bracket(c2_lower, c2_middle, c2_upper, &
          lower_profile%objective, middle_profile%objective, &
          upper_profile%objective)) then
        bracket_found = .true.
        exit
      end if
      if (middle_profile%objective >= lower_profile%objective .and. &
          upper_profile%objective >= middle_profile%objective) exit
      c2_lower = c2_middle
      lower_profile = middle_profile
      c2_middle = c2_upper
      middle_profile = upper_profile
      c2_upper = pole_distance_k + settings%bracket_expansion_factor * &
        (c2_middle - pole_distance_k)
      if (.not. ieee_is_finite(c2_upper) .or. &
          c2_upper <= c2_middle) exit
      upper_profile = evaluate_wlf_profile(observations, &
        reference_temperature_k, c2_upper)
      if (.not. upper_profile%valid) exit
    end do
    fit%bracket_expansion_count = expansion

    if (.not. bracket_found) then
      if (upper_profile%valid .and. &
          upper_profile%objective <= middle_profile%objective) then
        call assign_poorly_identified_wlf_limit( &
          observations, linear_p, linear_objective, fit)
      else
        fit%status = WLF_FIT_NO_INTERIOR_BRACKET
        fit%message = "WLF C2 profile objective interior minimum üretmedi."
      end if
      return
    end if

    fit%bracket_lower_c2_k = c2_lower
    fit%bracket_middle_c2_k = c2_middle
    fit%bracket_upper_c2_k = c2_upper
    minimizer = minimize_scalar_brent(profile_objective_at_c2, &
      c2_lower, c2_middle, c2_upper, settings%absolute_tolerance_k, &
      settings%relative_tolerance, settings%maximum_iterations)
    fit%minimizer_iteration_count = minimizer%iteration_count
    fit%minimizer_evaluation_count = minimizer%function_evaluation_count
    if (minimizer%status /= SCALAR_MINIMIZER_SUCCESS .or. &
        .not. minimizer%converged) then
      fit%status = WLF_FIT_OPTIMIZATION_FAILED
      fit%message = "WLF profiled C2 Brent refinement başarısız."
      return
    end if

    optimum_profile = evaluate_wlf_profile(observations, &
      reference_temperature_k, minimizer%x_minimum)
    if (.not. optimum_profile%valid) then
      fit%status = SHIFT_LAW_FIT_NONFINITE_RESULT
      fit%message = "WLF optimum profile sonucu sonlu değil."
      return
    end if

    shift_scale = max(sum(observations%delta_s_j_minus_i**2) / &
      real(size(observations), dp), tiny(1.0_dp))
    improvement_tolerance = 2048.0_dp*epsilon(1.0_dp)*shift_scale
    if (linear_objective - optimum_profile%objective <= &
        improvement_tolerance) then
      call assign_poorly_identified_wlf_limit( &
        observations, linear_p, linear_objective, fit)
      return
    end if

    fit%c1 = optimum_profile%c1
    fit%c2_k = minimizer%x_minimum
    fit%p_c1_over_c2_per_k = fit%c1/fit%c2_k
    fit%q_inverse_c2_per_k = 1.0_dp/fit%c2_k
    fit%profile_objective_minimum = optimum_profile%objective
    if (.not. valid_identified_wlf_parameters(fit, pole_distance_k, &
        pole_margin_k)) then
      fit%status = WLF_FIT_INVALID_PARAMETERS
      fit%message = "Identified WLF C1/C2 parameters physical domain'e aykırı."
      return
    end if

    call populate_wlf_predictions(observations, fit)
    if (.not. residual_metrics_are_finite(fit%pair_rmse, &
        fit%pair_max_abs_residual, fit%pair_mean_residual)) then
      fit%status = SHIFT_LAW_FIT_NONFINITE_RESULT
      fit%message = "WLF pair residual metrics sonlu değil."
      return
    end if
    fit%fit_available = .true.
    fit%parameter_identifiable = .true.
    fit%residual_validation_available = &
      fit%validation_degree_of_freedom > 0
    fit%status = SHIFT_LAW_FIT_SUCCESS
    fit%message = "Adjacent-pair profiled WLF fit ve C1/C2 hazır."

  contains

    !> Brent tarafından çağrılan pure profile objective'tir. Girdi C2 [K],
    !! çıktı mean squared pair residual [boyutsuz^2] değeridir.
    pure function profile_objective_at_c2(c2_k) result(objective)
      real(dp), intent(in) :: c2_k
      real(dp) :: objective
      type(wlf_profile_evaluation_t) :: profile

      profile = evaluate_wlf_profile(observations, &
        reference_temperature_k, c2_k)
      if (profile%valid) then
        objective = profile%objective
      else
        objective = huge(1.0_dp)
      end if
    end function profile_objective_at_c2

  end function fit_tts_wlf_shift_law

  !> Arrhenius pair prediction'ını hesaplar:
  !! delta_s=beta*(1/T_j-1/T_i). beta [K], sıcaklıklar [K] ve çıktı
  !! boyutsuzdur. Reference temperature bu relative denklemde yoktur;
  !! dolayısıyla Ea_app reference reparameterization'dan bağımsızdır.
  pure function predict_tts_arrhenius_pair_shift(beta_k, observation) &
      result(delta_s)
    real(dp), intent(in) :: beta_k
    type(tts_pair_shift_observation_t), intent(in) :: observation
    real(dp) :: delta_s

    delta_s = beta_k*(1.0_dp/observation%temperature_j_k - &
      1.0_dp/observation%temperature_i_k)
  end function predict_tts_arrhenius_pair_shift

  !> WLF pair prediction'ını hesaplar. Her endpoint için
  !! g=-(T-T_ref)/(C2+T-T_ref), delta_s=C1*(g_j-g_i) kullanılır. C1
  !! boyutsuz, C2/T [K], çıktı boyutsuzdur; caller pole-safe domain sağlar.
  pure function predict_tts_wlf_pair_shift( &
      c1, c2_k, reference_temperature_k, observation) result(delta_s)
    real(dp), intent(in) :: c1
    real(dp), intent(in) :: c2_k
    real(dp), intent(in) :: reference_temperature_k
    type(tts_pair_shift_observation_t), intent(in) :: observation
    real(dp) :: delta_s
    real(dp) :: g_i
    real(dp) :: g_j

    g_i = -(observation%temperature_i_k-reference_temperature_k) / &
      (c2_k+observation%temperature_i_k-reference_temperature_k)
    g_j = -(observation%temperature_j_k-reference_temperature_k) / &
      (c2_k+observation%temperature_j_k-reference_temperature_k)
    delta_s = c1*(g_j-g_i)
  end function predict_tts_wlf_pair_shift

  !> Pair fit input sözleşmesini doğrular. Her endpoint sonlu/pozitif [K],
  !! aynı pair içinde farklı ve delta_s sonlu olmalıdır. Shared isotherm'ler
  !! adjacent observations arasında beklenir ve duplicate data sayılmaz.
  !! Global reference sonlu/pozitif K olmalıdır. Full identification katmanı
  !! bunun measured reference olduğunu ayrıca doğrular; fit yordamındaki daha
  !! genel sözleşme reference isotherm'in held-out olduğu LOTO fold'larına da
  !! izin verir. Failure error-stop yerine status/message döndürür.
  pure subroutine validate_tts_pair_shift_observations( &
      observations, reference_temperature_k, minimum_observation_count, &
      minimum_temperature_k, maximum_temperature_k, status, message)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: reference_temperature_k
    integer, intent(in) :: minimum_observation_count
    real(dp), intent(out) :: minimum_temperature_k
    real(dp), intent(out) :: maximum_temperature_k
    integer, intent(out) :: status
    character(len=*), intent(out) :: message

    integer :: i

    minimum_temperature_k = 0.0_dp
    maximum_temperature_k = 0.0_dp
    status = SHIFT_LAW_FIT_INVALID_INPUT
    message = ""
    if (minimum_observation_count <= 0) then
      message = "Minimum observation count pozitif olmalıdır."
      return
    end if
    if (size(observations) < minimum_observation_count) then
      status = SHIFT_LAW_FIT_INSUFFICIENT_DATA
      message = "Parametric fit için pair observation sayısı yetersiz."
      return
    end if
    if (.not. ieee_is_finite(reference_temperature_k) .or. &
        reference_temperature_k <= 0.0_dp) then
      message = "Reference temperature sonlu ve pozitif K olmalıdır."
      return
    end if

    minimum_temperature_k = huge(1.0_dp)
    maximum_temperature_k = -huge(1.0_dp)
    do i = 1, size(observations)
      if (.not. ieee_is_finite(observations(i)%temperature_i_k) .or. &
          observations(i)%temperature_i_k <= 0.0_dp .or. &
          .not. ieee_is_finite(observations(i)%temperature_j_k) .or. &
          observations(i)%temperature_j_k <= 0.0_dp) then
        message = "Pair endpoint sıcaklıkları sonlu ve pozitif K olmalıdır."
        return
      end if
      if (temperatures_are_machine_equivalent( &
          observations(i)%temperature_i_k, &
          observations(i)%temperature_j_k)) then
        message = "Bir pair aynı/duplicate endpoint sıcaklığı içeremez."
        return
      end if
      if (.not. ieee_is_finite(observations(i)%delta_s_j_minus_i)) then
        message = "Adjacent pair delta_s gözlemi sonlu olmalıdır."
        return
      end if
      minimum_temperature_k = min(minimum_temperature_k, &
        observations(i)%temperature_i_k, observations(i)%temperature_j_k)
      maximum_temperature_k = max(maximum_temperature_k, &
        observations(i)%temperature_i_k, observations(i)%temperature_j_k)
    end do
    if (maximum_temperature_k <= minimum_temperature_k) then
      status = SHIFT_LAW_FIT_INSUFFICIENT_DATA
      message = "Calibrated temperature domain pozitif span içermiyor."
      return
    end if

    status = SHIFT_LAW_FIT_SUCCESS
    message = "Adjacent pair observations geçerli."
  end subroutine validate_tts_pair_shift_observations

  pure function evaluate_wlf_profile( &
      observations, reference_temperature_k, c2_k) result(profile)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: c2_k
    type(wlf_profile_evaluation_t) :: profile

    real(dp), allocatable :: delta_g(:)
    real(dp), allocatable :: predictions(:)
    real(dp) :: denominator
    real(dp) :: numerator
    real(dp) :: denominator_i
    real(dp) :: denominator_j
    real(dp) :: g_i
    real(dp) :: g_j
    integer :: i

    if (.not. ieee_is_finite(c2_k) .or. c2_k <= 0.0_dp) return
    allocate(delta_g(size(observations)))
    do i = 1, size(observations)
      denominator_i = c2_k + observations(i)%temperature_i_k - &
        reference_temperature_k
      denominator_j = c2_k + observations(i)%temperature_j_k - &
        reference_temperature_k
      if (.not. ieee_is_finite(denominator_i) .or. &
          .not. ieee_is_finite(denominator_j) .or. &
          denominator_i <= 0.0_dp .or. denominator_j <= 0.0_dp) return
      g_i = -(observations(i)%temperature_i_k-reference_temperature_k) / &
        denominator_i
      g_j = -(observations(i)%temperature_j_k-reference_temperature_k) / &
        denominator_j
      delta_g(i) = g_j-g_i
      if (.not. ieee_is_finite(delta_g(i))) return
    end do

    denominator = sum(delta_g*delta_g)
    numerator = sum(delta_g*observations%delta_s_j_minus_i)
    if (.not. ieee_is_finite(denominator) .or. denominator <= tiny(1.0_dp) .or. &
        .not. ieee_is_finite(numerator)) return
    profile%c1 = numerator/denominator
    if (.not. ieee_is_finite(profile%c1) .or. profile%c1 <= 0.0_dp) return
    predictions = profile%c1*delta_g
    profile%objective = sum((observations%delta_s_j_minus_i - &
      predictions)**2)/real(size(observations), dp)
    profile%valid = ieee_is_finite(profile%objective) .and. &
      profile%objective >= 0.0_dp
  end function evaluate_wlf_profile

  pure subroutine evaluate_large_c2_limit(observations, p, objective)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(out) :: p
    real(dp), intent(out) :: objective

    real(dp), allocatable :: delta_temperature(:)
    real(dp), allocatable :: scaled_delta_temperature(:)
    real(dp), allocatable :: predictions(:)
    real(dp) :: denominator
    real(dp) :: maximum_delta_temperature
    real(dp) :: numerator
    integer :: i

    allocate(delta_temperature(size(observations)))
    do i = 1, size(observations)
      delta_temperature(i) = observations(i)%temperature_j_k - &
        observations(i)%temperature_i_k
    end do
    maximum_delta_temperature = maxval(abs(delta_temperature))
    if (.not. ieee_is_finite(maximum_delta_temperature) .or. &
        maximum_delta_temperature <= tiny(1.0_dp)) then
      p = 0.0_dp
      objective = huge(1.0_dp)
      return
    end if
    scaled_delta_temperature = delta_temperature/maximum_delta_temperature
    denominator = sum(scaled_delta_temperature**2)
    numerator = -sum(scaled_delta_temperature * &
      observations%delta_s_j_minus_i)
    p = numerator/(denominator*maximum_delta_temperature)
    predictions = -p*delta_temperature
    objective = sum((observations%delta_s_j_minus_i-predictions)**2) / &
      real(size(observations), dp)
  end subroutine evaluate_large_c2_limit

  pure subroutine assign_poorly_identified_wlf_limit( &
      observations, linear_p, linear_objective, fit)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    real(dp), intent(in) :: linear_p
    real(dp), intent(in) :: linear_objective
    type(tts_wlf_fit_result_t), intent(inout) :: fit

    integer :: i

    fit%status = WLF_FIT_POORLY_IDENTIFIED
    fit%message = "WLF profile large-C2 limitinde; C1 ve C2 ayrı belirlenemiyor."
    fit%fit_available = .true.
    fit%parameter_identifiable = .false.
    fit%residual_validation_available = &
      fit%validation_degree_of_freedom > 0
    fit%c1 = 0.0_dp
    fit%c2_k = 0.0_dp
    fit%p_c1_over_c2_per_k = linear_p
    fit%q_inverse_c2_per_k = 0.0_dp
    fit%profile_objective_minimum = linear_objective
    allocate(fit%predicted_pair_shifts(size(observations)))
    do i = 1, size(observations)
      fit%predicted_pair_shifts(i) = -linear_p * &
        (observations(i)%temperature_j_k-observations(i)%temperature_i_k)
    end do
    fit%pair_residuals = observations%delta_s_j_minus_i - &
      fit%predicted_pair_shifts
    call calculate_residual_metrics(fit%pair_residuals, fit%pair_rmse, &
      fit%pair_max_abs_residual, fit%pair_mean_residual)
  end subroutine assign_poorly_identified_wlf_limit

  pure subroutine populate_wlf_predictions(observations, fit)
    type(tts_pair_shift_observation_t), intent(in) :: observations(:)
    type(tts_wlf_fit_result_t), intent(inout) :: fit
    integer :: i

    allocate(fit%predicted_pair_shifts(size(observations)))
    do i = 1, size(observations)
      fit%predicted_pair_shifts(i) = predict_tts_wlf_pair_shift(fit%c1, &
        fit%c2_k, fit%reference_temperature_k, observations(i))
    end do
    fit%pair_residuals = observations%delta_s_j_minus_i - &
      fit%predicted_pair_shifts
    call calculate_residual_metrics(fit%pair_residuals, fit%pair_rmse, &
      fit%pair_max_abs_residual, fit%pair_mean_residual)
  end subroutine populate_wlf_predictions

  pure subroutine calculate_residual_metrics( &
      residuals, rmse, maximum_absolute_residual, mean_residual)
    real(dp), intent(in) :: residuals(:)
    real(dp), intent(out) :: rmse
    real(dp), intent(out) :: maximum_absolute_residual
    real(dp), intent(out) :: mean_residual

    if (size(residuals) == 0) then
      rmse = 0.0_dp
      maximum_absolute_residual = 0.0_dp
      mean_residual = 0.0_dp
      return
    end if
    rmse = sqrt(sum(residuals**2)/real(size(residuals), dp))
    maximum_absolute_residual = maxval(abs(residuals))
    mean_residual = sum(residuals)/real(size(residuals), dp)
  end subroutine calculate_residual_metrics

  pure elemental function temperatures_are_machine_equivalent(a, b) &
      result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent

    equivalent = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
      abs(a-b) <= 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(a), abs(b))
  end function temperatures_are_machine_equivalent

  pure function valid_wlf_configuration(settings) result(valid)
    type(tts_wlf_fit_configuration_t), intent(in) :: settings
    logical :: valid

    valid = ieee_is_finite(settings%bracket_expansion_factor) .and. &
      settings%bracket_expansion_factor > 1.0_dp .and. &
      settings%maximum_bracket_expansions > 0 .and. &
      ieee_is_finite(settings%absolute_tolerance_k) .and. &
      settings%absolute_tolerance_k > 0.0_dp .and. &
      ieee_is_finite(settings%relative_tolerance) .and. &
      settings%relative_tolerance >= 0.0_dp .and. &
      settings%maximum_iterations > 0
  end function valid_wlf_configuration

  pure function valid_identified_wlf_parameters( &
      fit, pole_distance_k, pole_margin_k) result(valid)
    type(tts_wlf_fit_result_t), intent(in) :: fit
    real(dp), intent(in) :: pole_distance_k
    real(dp), intent(in) :: pole_margin_k
    logical :: valid

    valid = ieee_is_finite(fit%c1) .and. fit%c1 > 0.0_dp .and. &
      ieee_is_finite(fit%c2_k) .and. &
      fit%c2_k > pole_distance_k+pole_margin_k .and. &
      ieee_is_finite(fit%p_c1_over_c2_per_k) .and. &
      fit%p_c1_over_c2_per_k > 0.0_dp .and. &
      ieee_is_finite(fit%q_inverse_c2_per_k) .and. &
      fit%q_inverse_c2_per_k > 0.0_dp
  end function valid_identified_wlf_parameters

  pure function residual_metrics_are_finite( &
      rmse, maximum_absolute_residual, mean_residual) result(valid)
    real(dp), intent(in) :: rmse
    real(dp), intent(in) :: maximum_absolute_residual
    real(dp), intent(in) :: mean_residual
    logical :: valid

    valid = ieee_is_finite(rmse) .and. rmse >= 0.0_dp .and. &
      ieee_is_finite(maximum_absolute_residual) .and. &
      maximum_absolute_residual >= 0.0_dp .and. &
      ieee_is_finite(mean_residual)
  end function residual_metrics_are_finite

end module tms_tts_shift_law_fit
