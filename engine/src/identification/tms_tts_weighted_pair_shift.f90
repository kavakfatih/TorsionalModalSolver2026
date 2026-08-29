module tms_tts_weighted_pair_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_scalar_minimizer, only : scalar_minimizer_result_t, &
    SCALAR_MINIMIZER_SUCCESS, minimize_scalar_brent, &
    is_valid_scalar_minimum_bracket
  use tms_weighted_piecewise_integrals, only : &
    weighted_interval_integral_t, huber_interval_integral_t, &
    integrate_weighted_linear_interval, integrate_huber_linear_interval
  use tms_tts_types, only : tts_isotherm_t, &
    tts_pair_shift_configuration_t, TTS_CHANNEL_STORAGE, &
    TTS_CHANNEL_LOSS, TTS_CHANNEL_JOINT, SHIFT_FROM_JOINT, &
    SHIFT_FROM_STORAGE_ONLY, are_tts_values_machine_equivalent
  use tms_tts_pair_shift, only : get_tts_pair_feasible_shift_domain
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, &
    tts_uncertainty_log_segment_t, &
    tts_uncertainty_objective_evaluation_t, &
    tts_uncertainty_pair_solution_t, TTS_WEIGHTED_L2_OBJECTIVE, &
    TTS_STANDARDIZED_HUBER_OBJECTIVE, WEIGHTED_SHIFT_SUCCESS, &
    WEIGHTED_SHIFT_STORAGE_ONLY, WEIGHTED_SHIFT_INVALID_INPUT, &
    WEIGHTED_SHIFT_NO_OVERLAP, WEIGHTED_SHIFT_INSUFFICIENT_SUPPORT, &
    WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM, &
    WEIGHTED_SHIFT_NUMERICAL_FAILURE
  use tms_tts_uncertainty_validation, only : &
    build_tts_uncertainty_log_segments
  implicit none
  private

  type :: channel_optimization_t
    integer :: status = WEIGHTED_SHIFT_INVALID_INPUT
    logical :: available = .false.
    logical :: any_valid_objective = .false.
    logical :: numerical_failure_observed = .false.
    real(dp) :: shift = 0.0_dp
    real(dp) :: objective = huge(1.0_dp)
    integer :: iteration_count = 0
    integer :: evaluation_count = 0
    type(tts_uncertainty_objective_evaluation_t) :: diagnostics
  end type channel_optimization_t

  public :: evaluate_tts_uncertainty_pair_channel_objective
  public :: evaluate_tts_uncertainty_pair_objective
  public :: identify_tts_uncertainty_pair_shift

contains

  !> Tek storage/loss channel için uncertainty-supported, overlap-width
  !! normalize objective'i değerlendirir. Weighted L2, integral r^2/v;
  !! Huber ise standardized z=r/sqrt(v) üzerinde rho_c(z) kullanır. Moving
  !! curve mevcut V0.8.1 convention'ıyla x_shifted=x+delta_s taşınır.
  pure function evaluate_tts_uncertainty_pair_channel_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, delta_s, &
      channel, objective_mode, huber_c) result(evaluation)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    real(dp), intent(in) :: delta_s
    integer, intent(in) :: channel
    integer, intent(in) :: objective_mode
    real(dp), intent(in) :: huber_c
    type(tts_uncertainty_objective_evaluation_t) :: evaluation

    type(tts_uncertainty_log_segment_t), allocatable :: &
      reference_segments(:)
    type(tts_uncertainty_log_segment_t), allocatable :: moving_segments(:)
    real(dp) :: integral
    real(dp) :: squared_standardized_integral
    real(dp) :: reference_measure
    real(dp) :: moving_measure
    real(dp) :: overlap_measure
    real(dp) :: quadratic_measure
    real(dp) :: tail_measure
    integer :: interval_count
    logical :: numerical_failure

    if (.not. ieee_is_finite(delta_s) .or. &
        .not. ieee_is_finite(huber_c) .or. huber_c <= 0.0_dp) return
    if (channel /= TTS_CHANNEL_STORAGE .and. &
        channel /= TTS_CHANNEL_LOSS) return
    if (objective_mode /= TTS_WEIGHTED_L2_OBJECTIVE .and. &
        objective_mode /= TTS_STANDARDIZED_HUBER_OBJECTIVE) return

    reference_segments = build_tts_uncertainty_log_segments( &
      reference_isotherm, uncertainty_family, channel)
    moving_segments = build_tts_uncertainty_log_segments( &
      moving_isotherm, uncertainty_family, channel)
    if (size(reference_segments) == 0 .or. size(moving_segments) == 0) return

    call integrate_uncertainty_segment_overlaps(reference_segments, &
      moving_segments, delta_s, objective_mode, huber_c, integral, &
      squared_standardized_integral, overlap_measure, quadratic_measure, &
      tail_measure, interval_count, numerical_failure)
    if (numerical_failure) then
      evaluation%numerical_failure = .true.
      return
    end if
    if (interval_count <= 0 .or. overlap_measure <= 0.0_dp .or. &
        .not. ieee_is_finite(integral)) return
    reference_measure = total_segment_measure(reference_segments)
    moving_measure = total_segment_measure(moving_segments)
    if (reference_measure <= 0.0_dp .or. moving_measure <= 0.0_dp) return

    evaluation%objective = integral/overlap_measure
    evaluation%rms_standardized_residual = &
      sqrt(squared_standardized_integral/overlap_measure)
    if (.not. ieee_is_finite(evaluation%objective) .or. &
        .not. ieee_is_finite(evaluation%rms_standardized_residual)) then
      evaluation%numerical_failure = .true.
      evaluation%objective = huge(1.0_dp)
      return
    end if
    evaluation%valid = .true.
    evaluation%overlap_width_decades = overlap_measure
    evaluation%overlap_fraction = overlap_measure / &
      min(reference_measure, moving_measure)
    evaluation%quadratic_width_decades = quadratic_measure
    evaluation%tail_width_decades = tail_measure
    if (overlap_measure > 0.0_dp) &
      evaluation%tail_fraction = tail_measure/overlap_measure
    evaluation%interpolation_interval_count = interval_count
  end function evaluate_tts_uncertainty_pair_channel_objective

  !> Joint objective her channel kendi gerçek support genişliğiyle normalize
  !! edildikten sonra 0.5 storage + 0.5 loss olarak kurulur. Sample point sayısı
  !! channel ağırlığını değiştirmez. Storage-only değerlendirme explicit'tir.
  pure function evaluate_tts_uncertainty_pair_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, delta_s, &
      channel, objective_mode, huber_c) result(evaluation)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    real(dp), intent(in) :: delta_s
    integer, intent(in) :: channel
    integer, intent(in) :: objective_mode
    real(dp), intent(in) :: huber_c
    type(tts_uncertainty_objective_evaluation_t) :: evaluation

    type(tts_uncertainty_objective_evaluation_t) :: storage
    type(tts_uncertainty_objective_evaluation_t) :: loss

    select case (channel)
    case (TTS_CHANNEL_STORAGE, TTS_CHANNEL_LOSS)
      evaluation = evaluate_tts_uncertainty_pair_channel_objective( &
        reference_isotherm, moving_isotherm, uncertainty_family, delta_s, &
        channel, objective_mode, huber_c)
    case (TTS_CHANNEL_JOINT)
      storage = evaluate_tts_uncertainty_pair_channel_objective( &
        reference_isotherm, moving_isotherm, uncertainty_family, delta_s, &
        TTS_CHANNEL_STORAGE, objective_mode, huber_c)
      loss = evaluate_tts_uncertainty_pair_channel_objective( &
        reference_isotherm, moving_isotherm, uncertainty_family, delta_s, &
        TTS_CHANNEL_LOSS, objective_mode, huber_c)
      evaluation%numerical_failure = storage%numerical_failure .or. &
        loss%numerical_failure
      if (.not. storage%valid .or. .not. loss%valid) return
      evaluation%valid = .true.
      evaluation%objective = 0.5_dp*(storage%objective + loss%objective)
      evaluation%overlap_width_decades = storage%overlap_width_decades + &
        loss%overlap_width_decades
      evaluation%overlap_fraction = min(storage%overlap_fraction, &
        loss%overlap_fraction)
      evaluation%quadratic_width_decades = &
        storage%quadratic_width_decades + loss%quadratic_width_decades
      evaluation%tail_width_decades = storage%tail_width_decades + &
        loss%tail_width_decades
      if (evaluation%overlap_width_decades > 0.0_dp) &
        evaluation%tail_fraction = evaluation%tail_width_decades / &
          evaluation%overlap_width_decades
      evaluation%rms_standardized_residual = sqrt(0.5_dp * &
        (storage%rms_standardized_residual**2 + &
        loss%rms_standardized_residual**2))
      evaluation%interpolation_interval_count = &
        storage%interpolation_interval_count + &
        loss%interpolation_interval_count
    end select
  end function evaluate_tts_uncertainty_pair_objective

  !> Bir adjacent pair için weighted-L2 veya standardized-Huber shift'i mevcut
  !! V0.8.1 coarse-scan/strict-interior-bracket/Brent politikasıyla çözer.
  !! Boundary-only ya da flat minimum başarı olarak uydurulmaz.
  pure function identify_tts_uncertainty_pair_shift( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      objective_mode, configuration, huber_c) result(pair_result)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    integer, intent(in) :: objective_mode
    type(tts_pair_shift_configuration_t), intent(in) :: configuration
    real(dp), intent(in) :: huber_c
    type(tts_uncertainty_pair_solution_t) :: pair_result

    type(channel_optimization_t) :: joint
    type(channel_optimization_t) :: loss
    type(channel_optimization_t) :: storage

    pair_result%objective_mode = objective_mode
    if (.not. pair_configuration_is_valid(configuration) .or. &
        .not. ieee_is_finite(huber_c) .or. huber_c <= 0.0_dp .or. &
        (objective_mode /= TTS_WEIGHTED_L2_OBJECTIVE .and. &
        objective_mode /= TTS_STANDARDIZED_HUBER_OBJECTIVE)) return

    storage = optimize_pair_channel(reference_isotherm, moving_isotherm, &
      uncertainty_family, TTS_CHANNEL_STORAGE, objective_mode, &
      configuration, huber_c)
    loss = optimize_pair_channel(reference_isotherm, moving_isotherm, &
      uncertainty_family, TTS_CHANNEL_LOSS, objective_mode, &
      configuration, huber_c)
    joint = optimize_pair_channel(reference_isotherm, moving_isotherm, &
      uncertainty_family, TTS_CHANNEL_JOINT, objective_mode, &
      configuration, huber_c)
    pair_result%storage_shift_available = storage%available
    pair_result%loss_shift_available = loss%available
    pair_result%joint_shift_available = joint%available
    if (storage%available) then
      pair_result%storage_shift = storage%shift
      pair_result%storage_diagnostics = storage%diagnostics
    end if
    if (loss%available) then
      pair_result%loss_shift = loss%shift
      pair_result%loss_diagnostics = loss%diagnostics
    end if

    if (joint%available) then
      pair_result%status = WEIGHTED_SHIFT_SUCCESS
      pair_result%production_channel = SHIFT_FROM_JOINT
      pair_result%shift_available = .true.
      pair_result%shift = joint%shift
      pair_result%joint_shift = joint%shift
      pair_result%objective_minimum = joint%objective
      pair_result%iteration_count = joint%iteration_count
      pair_result%evaluation_count = joint%evaluation_count
      pair_result%production_diagnostics = joint%diagnostics
      return
    end if
    if (storage%available .and. .not. loss%any_valid_objective) then
      pair_result%status = WEIGHTED_SHIFT_STORAGE_ONLY
      pair_result%production_channel = SHIFT_FROM_STORAGE_ONLY
      pair_result%shift_available = .true.
      pair_result%shift = storage%shift
      pair_result%objective_minimum = storage%objective
      pair_result%iteration_count = storage%iteration_count
      pair_result%evaluation_count = storage%evaluation_count
      pair_result%production_diagnostics = storage%diagnostics
      return
    end if
    pair_result%status = merge_failure_status(storage, loss, joint)
    pair_result%iteration_count = joint%iteration_count
    pair_result%evaluation_count = storage%evaluation_count + &
      loss%evaluation_count + joint%evaluation_count
  end function identify_tts_uncertainty_pair_shift

  pure subroutine integrate_uncertainty_segment_overlaps( &
      reference_segments, moving_segments, delta_s, objective_mode, &
      huber_c, integral, squared_standardized_integral, overlap_measure, &
      quadratic_measure, tail_measure, interval_count, numerical_failure)
    type(tts_uncertainty_log_segment_t), intent(in) :: reference_segments(:)
    type(tts_uncertainty_log_segment_t), intent(in) :: moving_segments(:)
    real(dp), intent(in) :: delta_s
    integer, intent(in) :: objective_mode
    real(dp), intent(in) :: huber_c
    real(dp), intent(out) :: integral
    real(dp), intent(out) :: squared_standardized_integral
    real(dp), intent(out) :: overlap_measure
    real(dp), intent(out) :: quadratic_measure
    real(dp), intent(out) :: tail_measure
    integer, intent(out) :: interval_count
    logical, intent(out) :: numerical_failure

    real(dp) :: pair_integral
    real(dp) :: pair_squared_integral
    real(dp) :: pair_overlap
    real(dp) :: pair_quadratic
    real(dp) :: pair_tail
    integer :: i
    integer :: j
    integer :: pair_intervals
    logical :: pair_failure

    integral = 0.0_dp
    squared_standardized_integral = 0.0_dp
    overlap_measure = 0.0_dp
    quadratic_measure = 0.0_dp
    tail_measure = 0.0_dp
    interval_count = 0
    numerical_failure = .false.
    do i = 1, size(reference_segments)
      do j = 1, size(moving_segments)
        call integrate_uncertainty_segment_pair(reference_segments(i), &
          moving_segments(j), delta_s, objective_mode, huber_c, &
          pair_integral, pair_squared_integral, pair_overlap, &
          pair_quadratic, pair_tail, pair_intervals, pair_failure)
        if (pair_failure) then
          numerical_failure = .true.
          return
        end if
        integral = integral + pair_integral
        squared_standardized_integral = &
          squared_standardized_integral + pair_squared_integral
        overlap_measure = overlap_measure + pair_overlap
        quadratic_measure = quadratic_measure + pair_quadratic
        tail_measure = tail_measure + pair_tail
        interval_count = interval_count + pair_intervals
      end do
    end do
  end subroutine integrate_uncertainty_segment_overlaps

  pure subroutine integrate_uncertainty_segment_pair( &
      reference_segment, moving_segment, delta_s, objective_mode, &
      huber_c, integral, squared_standardized_integral, overlap_measure, &
      quadratic_measure, tail_measure, interval_count, numerical_failure)
    type(tts_uncertainty_log_segment_t), intent(in) :: reference_segment
    type(tts_uncertainty_log_segment_t), intent(in) :: moving_segment
    real(dp), intent(in) :: delta_s
    integer, intent(in) :: objective_mode
    real(dp), intent(in) :: huber_c
    real(dp), intent(out) :: integral
    real(dp), intent(out) :: squared_standardized_integral
    real(dp), intent(out) :: overlap_measure
    real(dp), intent(out) :: quadratic_measure
    real(dp), intent(out) :: tail_measure
    integer, intent(out) :: interval_count
    logical, intent(out) :: numerical_failure

    type(weighted_interval_integral_t) :: weighted
    type(huber_interval_integral_t) :: huber
    real(dp), allocatable :: knots(:)
    real(dp) :: h
    real(dp) :: lower_overlap
    real(dp) :: r0
    real(dp) :: r1
    real(dp) :: reference_v
    real(dp) :: moving_v
    real(dp) :: upper_overlap
    real(dp) :: v0
    real(dp) :: v1
    integer :: i
    integer :: knot_count

    integral = 0.0_dp
    squared_standardized_integral = 0.0_dp
    overlap_measure = 0.0_dp
    quadratic_measure = 0.0_dp
    tail_measure = 0.0_dp
    interval_count = 0
    numerical_failure = .false.
    lower_overlap = max(reference_segment%x(1), &
      moving_segment%x(1) + delta_s)
    upper_overlap = min(reference_segment%x(size(reference_segment%x)), &
      moving_segment%x(size(moving_segment%x)) + delta_s)
    if (upper_overlap <= lower_overlap .or. &
        are_tts_values_machine_equivalent(upper_overlap, lower_overlap)) return

    allocate(knots(size(reference_segment%x) + size(moving_segment%x) + 2))
    call build_merged_knots(reference_segment, moving_segment, delta_s, &
      lower_overlap, upper_overlap, knots, knot_count)
    do i = 1, knot_count - 1
      h = knots(i + 1) - knots(i)
      if (h <= 0.0_dp .or. are_tts_values_machine_equivalent( &
          knots(i + 1), knots(i))) cycle
      r0 = interpolate_segment_field(reference_segment, knots(i), .false.) - &
        interpolate_segment_field(moving_segment, &
          knots(i) - delta_s, .false.)
      r1 = interpolate_segment_field(reference_segment, &
        knots(i + 1), .false.) - interpolate_segment_field( &
        moving_segment, knots(i + 1) - delta_s, .false.)
      reference_v = interpolate_segment_field( &
        reference_segment, knots(i), .true.)
      moving_v = interpolate_segment_field( &
        moving_segment, knots(i) - delta_s, .true.)
      v0 = reference_v + moving_v
      reference_v = interpolate_segment_field( &
        reference_segment, knots(i + 1), .true.)
      moving_v = interpolate_segment_field( &
        moving_segment, knots(i + 1) - delta_s, .true.)
      v1 = reference_v + moving_v
      if (.not. ieee_is_finite(v0) .or. .not. ieee_is_finite(v1) .or. &
          v0 <= 0.0_dp .or. v1 <= 0.0_dp) then
        numerical_failure = .true.
        return
      end if
      weighted = integrate_weighted_linear_interval(h, r0, r1, v0, v1)
      if (.not. weighted%valid) then
        numerical_failure = .true.
        return
      end if
      squared_standardized_integral = &
        squared_standardized_integral + weighted%integral
      if (objective_mode == TTS_WEIGHTED_L2_OBJECTIVE) then
        integral = integral + weighted%integral
      else
        huber = integrate_huber_linear_interval( &
          h, r0, r1, v0, v1, huber_c)
        if (.not. huber%valid) then
          numerical_failure = .true.
          return
        end if
        integral = integral + huber%integral
        quadratic_measure = quadratic_measure + &
          huber%quadratic_width_decades
        tail_measure = tail_measure + huber%tail_width_decades
      end if
      overlap_measure = overlap_measure + h
      interval_count = interval_count + 1
    end do
    if (objective_mode == TTS_WEIGHTED_L2_OBJECTIVE) &
      quadratic_measure = overlap_measure
  end subroutine integrate_uncertainty_segment_pair

  pure subroutine build_merged_knots( &
      reference_segment, moving_segment, delta_s, lower_overlap, &
      upper_overlap, knots, knot_count)
    type(tts_uncertainty_log_segment_t), intent(in) :: reference_segment
    type(tts_uncertainty_log_segment_t), intent(in) :: moving_segment
    real(dp), intent(in) :: delta_s
    real(dp), intent(in) :: lower_overlap
    real(dp), intent(in) :: upper_overlap
    real(dp), intent(out) :: knots(:)
    integer, intent(out) :: knot_count
    integer :: i

    knot_count = 2
    knots(1) = lower_overlap
    knots(2) = upper_overlap
    do i = 1, size(reference_segment%x)
      if (reference_segment%x(i) > lower_overlap .and. &
          reference_segment%x(i) < upper_overlap) then
        knot_count = knot_count + 1
        knots(knot_count) = reference_segment%x(i)
      end if
    end do
    do i = 1, size(moving_segment%x)
      if (moving_segment%x(i) + delta_s > lower_overlap .and. &
          moving_segment%x(i) + delta_s < upper_overlap) then
        knot_count = knot_count + 1
        knots(knot_count) = moving_segment%x(i) + delta_s
      end if
    end do
    call sort_and_compact_knots(knots, knot_count)
  end subroutine build_merged_knots

  pure subroutine sort_and_compact_knots(knots, knot_count)
    real(dp), intent(inout) :: knots(:)
    integer, intent(inout) :: knot_count
    real(dp) :: key
    integer :: i
    integer :: j
    integer :: unique_count

    do i = 2, knot_count
      key = knots(i)
      j = i - 1
      do while (j >= 1)
        if (knots(j) <= key) exit
        knots(j + 1) = knots(j)
        j = j - 1
      end do
      knots(j + 1) = key
    end do
    unique_count = 1
    do i = 2, knot_count
      if (.not. are_tts_values_machine_equivalent( &
          knots(i), knots(unique_count))) then
        unique_count = unique_count + 1
        knots(unique_count) = knots(i)
      end if
    end do
    knot_count = unique_count
  end subroutine sort_and_compact_knots

  pure function interpolate_segment_field(segment, x, variance_field) &
      result(value)
    type(tts_uncertainty_log_segment_t), intent(in) :: segment
    real(dp), intent(in) :: x
    logical, intent(in) :: variance_field
    real(dp) :: value
    real(dp) :: alpha
    integer :: i

    if (are_tts_values_machine_equivalent(x, segment%x(1))) then
      if (variance_field) then
        value = segment%variance(1)
      else
        value = segment%y(1)
      end if
      return
    end if
    if (are_tts_values_machine_equivalent( &
        x, segment%x(size(segment%x)))) then
      if (variance_field) then
        value = segment%variance(size(segment%variance))
      else
        value = segment%y(size(segment%y))
      end if
      return
    end if
    do i = 1, size(segment%x) - 1
      if (x >= segment%x(i) .and. x <= segment%x(i + 1)) then
        alpha = (x - segment%x(i))/(segment%x(i + 1) - segment%x(i))
        if (variance_field) then
          value = (1.0_dp - alpha)*segment%variance(i) + &
            alpha*segment%variance(i + 1)
        else
          value = (1.0_dp - alpha)*segment%y(i) + alpha*segment%y(i + 1)
        end if
        return
      end if
    end do
    value = huge(1.0_dp)
  end function interpolate_segment_field

  pure function total_segment_measure(segments) result(measure)
    type(tts_uncertainty_log_segment_t), intent(in) :: segments(:)
    real(dp) :: measure
    integer :: i

    measure = 0.0_dp
    do i = 1, size(segments)
      measure = measure + &
        segments(i)%x(size(segments(i)%x)) - segments(i)%x(1)
    end do
  end function total_segment_measure

  pure function optimize_pair_channel( &
      reference_isotherm, moving_isotherm, uncertainty_family, channel, &
      objective_mode, settings, huber_c) result(solution)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    integer, intent(in) :: channel
    integer, intent(in) :: objective_mode
    type(tts_pair_shift_configuration_t), intent(in) :: settings
    real(dp), intent(in) :: huber_c
    type(channel_optimization_t) :: solution

    real(dp), allocatable :: objective_values(:)
    logical, allocatable :: objective_valid(:)
    real(dp), allocatable :: scan_shifts(:)
    type(tts_uncertainty_objective_evaluation_t) :: evaluation
    type(scalar_minimizer_result_t) :: minimizer
    real(dp) :: interior_lower
    real(dp) :: interior_upper
    real(dp) :: lower_shift
    real(dp) :: margin
    real(dp) :: scale
    real(dp) :: upper_shift
    integer :: bracket_index
    integer :: i
    logical :: domain_valid

    call get_tts_pair_feasible_shift_domain(reference_isotherm, &
      moving_isotherm, lower_shift, upper_shift, domain_valid)
    if (.not. domain_valid) then
      solution%status = WEIGHTED_SHIFT_NO_OVERLAP
      return
    end if
    scale = max(1.0_dp, abs(lower_shift), abs(upper_shift))
    margin = 128.0_dp*epsilon(1.0_dp)*scale
    interior_lower = lower_shift + margin
    interior_upper = upper_shift - margin
    if (interior_lower >= interior_upper) then
      solution%status = WEIGHTED_SHIFT_NO_OVERLAP
      return
    end if

    allocate(scan_shifts(settings%coarse_scan_point_count))
    allocate(objective_values(settings%coarse_scan_point_count))
    allocate(objective_valid(settings%coarse_scan_point_count))
    do i = 1, settings%coarse_scan_point_count
      scan_shifts(i) = interior_lower + real(i - 1, dp) * &
        (interior_upper - interior_lower) / &
        real(settings%coarse_scan_point_count - 1, dp)
      evaluation = evaluate_tts_uncertainty_pair_objective( &
        reference_isotherm, moving_isotherm, uncertainty_family, &
        scan_shifts(i), channel, objective_mode, huber_c)
      objective_valid(i) = evaluation%valid
      objective_values(i) = evaluation%objective
      solution%numerical_failure_observed = &
        solution%numerical_failure_observed .or. &
        evaluation%numerical_failure
    end do
    solution%evaluation_count = settings%coarse_scan_point_count
    solution%any_valid_objective = any(objective_valid)
    if (.not. solution%any_valid_objective) then
      if (solution%numerical_failure_observed) then
        solution%status = WEIGHTED_SHIFT_NUMERICAL_FAILURE
      else
        solution%status = WEIGHTED_SHIFT_INSUFFICIENT_SUPPORT
      end if
      return
    end if

    bracket_index = 0
    do i = 2, settings%coarse_scan_point_count - 1
      if (objective_valid(i - 1) .and. objective_valid(i) .and. &
          objective_valid(i + 1) .and. is_valid_scalar_minimum_bracket( &
            scan_shifts(i - 1), scan_shifts(i), scan_shifts(i + 1), &
            objective_values(i - 1), objective_values(i), &
            objective_values(i + 1))) then
        if (bracket_index == 0) then
          bracket_index = i
        else if (objective_values(i) < objective_values(bracket_index)) then
          bracket_index = i
        end if
      end if
    end do
    if (bracket_index == 0) then
      solution%status = WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM
      return
    end if

    minimizer = minimize_scalar_brent(objective_at_shift, &
      scan_shifts(bracket_index - 1), scan_shifts(bracket_index), &
      scan_shifts(bracket_index + 1), settings%absolute_tolerance, &
      settings%relative_tolerance, settings%maximum_iterations)
    solution%evaluation_count = solution%evaluation_count + &
      minimizer%function_evaluation_count
    solution%iteration_count = minimizer%iteration_count
    if (minimizer%status /= SCALAR_MINIMIZER_SUCCESS) then
      solution%status = WEIGHTED_SHIFT_NUMERICAL_FAILURE
      return
    end if
    evaluation = evaluate_tts_uncertainty_pair_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      minimizer%x_minimum, channel, objective_mode, huber_c)
    solution%evaluation_count = solution%evaluation_count + 1
    if (.not. evaluation%valid) then
      solution%status = WEIGHTED_SHIFT_NUMERICAL_FAILURE
      return
    end if
    solution%available = .true.
    solution%status = WEIGHTED_SHIFT_SUCCESS
    solution%shift = minimizer%x_minimum
    solution%objective = evaluation%objective
    solution%diagnostics = evaluation

  contains

    pure function objective_at_shift(shift) result(value)
      real(dp), intent(in) :: shift
      real(dp) :: value
      type(tts_uncertainty_objective_evaluation_t) :: local_evaluation

      local_evaluation = evaluate_tts_uncertainty_pair_objective( &
        reference_isotherm, moving_isotherm, uncertainty_family, shift, &
        channel, objective_mode, huber_c)
      if (local_evaluation%valid) then
        value = local_evaluation%objective
      else
        value = huge(1.0_dp)
      end if
    end function objective_at_shift

  end function optimize_pair_channel

  pure function pair_configuration_is_valid(settings) result(valid)
    type(tts_pair_shift_configuration_t), intent(in) :: settings
    logical :: valid

    valid = settings%coarse_scan_point_count >= 3 .and. &
      ieee_is_finite(settings%absolute_tolerance) .and. &
      ieee_is_finite(settings%relative_tolerance) .and. &
      settings%absolute_tolerance > 0.0_dp .and. &
      settings%relative_tolerance >= 0.0_dp .and. &
      settings%maximum_iterations > 0
  end function pair_configuration_is_valid

  pure function merge_failure_status(storage, loss, joint) result(status)
    type(channel_optimization_t), intent(in) :: storage
    type(channel_optimization_t), intent(in) :: loss
    type(channel_optimization_t), intent(in) :: joint
    integer :: status

    if (joint%status == WEIGHTED_SHIFT_NUMERICAL_FAILURE .or. &
        storage%status == WEIGHTED_SHIFT_NUMERICAL_FAILURE .or. &
        loss%status == WEIGHTED_SHIFT_NUMERICAL_FAILURE) then
      status = WEIGHTED_SHIFT_NUMERICAL_FAILURE
    else if (joint%status == WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM .or. &
        storage%status == WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM .or. &
        loss%status == WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM) then
      status = WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM
    else if (joint%status == WEIGHTED_SHIFT_NO_OVERLAP .or. &
        storage%status == WEIGHTED_SHIFT_NO_OVERLAP) then
      status = WEIGHTED_SHIFT_NO_OVERLAP
    else
      status = WEIGHTED_SHIFT_INSUFFICIENT_SUPPORT
    end if
  end function merge_failure_status

end module tms_tts_weighted_pair_shift
