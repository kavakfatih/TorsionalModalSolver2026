module tms_tts_pair_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_scalar_minimizer, only : scalar_minimizer_result_t, &
    SCALAR_MINIMIZER_SUCCESS, minimize_scalar_brent, &
    is_valid_scalar_minimum_bracket
  use tms_tts_types, only : tts_isotherm_t, tts_measurement_point_t, &
    tts_log_segment_t, tts_pair_objective_evaluation_t, &
    tts_pair_shift_configuration_t, tts_pair_shift_result_t, &
    TTS_CHANNEL_STORAGE, TTS_CHANNEL_LOSS, TTS_CHANNEL_JOINT, &
    SHIFT_FROM_JOINT, SHIFT_FROM_STORAGE_ONLY, PAIR_SHIFT_SUCCESS, &
    PAIR_SHIFT_STORAGE_ONLY, PAIR_SHIFT_INVALID_INPUT, &
    PAIR_SHIFT_NO_OVERLAP, PAIR_SHIFT_INSUFFICIENT_SUPPORT, &
    PAIR_SHIFT_NO_INTERIOR_MINIMUM, PAIR_SHIFT_OPTIMIZATION_FAILED, &
    is_storage_log_usable, is_loss_log_usable, &
    are_tts_values_machine_equivalent
  implicit none
  private

  type :: channel_optimization_t
    integer :: status = PAIR_SHIFT_INVALID_INPUT
    logical :: available = .false.
    logical :: any_valid_objective = .false.
    real(dp) :: shift = 0.0_dp
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: overlap_width_decades = 0.0_dp
    real(dp) :: overlap_fraction = 0.0_dp
    real(dp) :: curvature = 0.0_dp
    integer :: iteration_count = 0
    integer :: evaluation_count = 0
  end type channel_optimization_t

  public :: build_tts_valid_log_segments
  public :: evaluate_tts_pair_channel_objective
  public :: evaluate_tts_pair_objective
  public :: get_tts_pair_feasible_shift_domain
  public :: identify_tts_pair_shift
  public :: integrate_linear_squared_residual

contains

  !> Bir linear residual'ın [x0,x1] aralığındaki exact L2 integralini verir.
  !! h=x1-x0 [decade], r0/r1 log10-modulus farklarıdır. Model, residual'ın
  !! subinterval içinde lineer olmasıdır: integral=h/3(r0^2+r0*r1+r1^2).
  pure function integrate_linear_squared_residual(h, r0, r1) result(value)
    real(dp), intent(in) :: h
    real(dp), intent(in) :: r0
    real(dp), intent(in) :: r1
    real(dp) :: value

    value = h*(r0*r0 + r0*r1 + r1*r1)/3.0_dp
  end function integrate_linear_squared_residual

  !> Bir isotherm'de quality-invalid noktaların üzerinden interpolation
  !! yapmadan contiguous log-usable segmentler kurar. Storage için VALID G'>0,
  !! loss için VALID G''>0 gerekir. x=log10(f/Hz), y=log10(G/Pa); tek noktalı
  !! koşular gerçek interpolation interval'i olmadığından sonuç dışıdır.
  pure function build_tts_valid_log_segments(isotherm, channel) &
      result(segments)
    type(tts_isotherm_t), intent(in) :: isotherm
    integer, intent(in) :: channel
    type(tts_log_segment_t), allocatable :: segments(:)

    integer :: point_count
    integer :: point_index
    integer :: run_end
    integer :: run_start
    integer :: segment_count
    integer :: segment_index

    if (.not. allocated(isotherm%points)) then
      allocate(segments(0))
      return
    end if
    point_count = size(isotherm%points)
    segment_count = 0
    point_index = 1
    do while (point_index <= point_count)
      if (.not. point_is_log_usable(isotherm%points(point_index), channel)) then
        point_index = point_index + 1
        cycle
      end if
      run_start = point_index
      run_end = point_index
      do while (run_end < point_count)
        if (.not. point_is_log_usable( &
            isotherm%points(run_end + 1), channel)) exit
        run_end = run_end + 1
      end do
      if (run_end - run_start + 1 >= 2) segment_count = segment_count + 1
      point_index = run_end + 1
    end do

    allocate(segments(segment_count))
    point_index = 1
    segment_index = 0
    do while (point_index <= point_count)
      if (.not. point_is_log_usable(isotherm%points(point_index), channel)) then
        point_index = point_index + 1
        cycle
      end if
      run_start = point_index
      run_end = point_index
      do while (run_end < point_count)
        if (.not. point_is_log_usable( &
            isotherm%points(run_end + 1), channel)) exit
        run_end = run_end + 1
      end do
      if (run_end - run_start + 1 >= 2) then
        segment_index = segment_index + 1
        call populate_segment(segments(segment_index), isotherm, channel, &
          run_start, run_end)
      end if
      point_index = run_end + 1
    end do
  end function build_tts_valid_log_segments

  !> Measured log-frequency domain'lerinden relative shift'in açık feasible
  !! aralığını türetir: L_i-U_j < delta_s < U_i-L_j. Arbitrary global bound
  !! kullanılmaz. Bounds s=log10(a_T) boyutsuz koordinatındadır.
  pure subroutine get_tts_pair_feasible_shift_domain( &
      reference_isotherm, moving_isotherm, lower_shift, upper_shift, valid)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    real(dp), intent(out) :: lower_shift
    real(dp), intent(out) :: upper_shift
    logical, intent(out) :: valid

    real(dp) :: reference_lower
    real(dp) :: reference_upper
    real(dp) :: moving_lower
    real(dp) :: moving_upper

    lower_shift = 0.0_dp
    upper_shift = 0.0_dp
    valid = is_isotherm_structurally_valid(reference_isotherm) .and. &
      is_isotherm_structurally_valid(moving_isotherm)
    if (.not. valid) return

    reference_lower = log10(reference_isotherm%points(1)%frequency_hz)
    reference_upper = log10(reference_isotherm%points( &
      size(reference_isotherm%points))%frequency_hz)
    moving_lower = log10(moving_isotherm%points(1)%frequency_hz)
    moving_upper = log10(moving_isotherm%points( &
      size(moving_isotherm%points))%frequency_hz)
    lower_shift = reference_lower - moving_upper
    upper_shift = reference_upper - moving_lower
    valid = ieee_is_finite(lower_shift) .and. ieee_is_finite(upper_shift) &
      .and. lower_shift < upper_shift
  end subroutine get_tts_pair_feasible_shift_domain

  !> Tek storage veya positive-loss channel için exact, normalized,
  !! piecewise-linear log-space residual objective'i hesaplar. Moving curve
  !! x_j_shifted=x_j+delta_s ile taşınır. Yalnız iki curve'ün contiguous valid
  !! segmentlerinin gerçek interval overlap'ları entegre edilir; extrapolation
  !! veya invalid-gap interpolation yapılmaz.
  pure function evaluate_tts_pair_channel_objective( &
      reference_isotherm, moving_isotherm, delta_s, channel) &
      result(evaluation)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    real(dp), intent(in) :: delta_s
    integer, intent(in) :: channel
    type(tts_pair_objective_evaluation_t) :: evaluation

    type(tts_log_segment_t), allocatable :: reference_segments(:)
    type(tts_log_segment_t), allocatable :: moving_segments(:)
    real(dp) :: integral
    real(dp) :: reference_measure
    real(dp) :: moving_measure
    real(dp) :: overlap_measure
    integer :: interval_count

    if (.not. ieee_is_finite(delta_s)) return
    if (channel /= TTS_CHANNEL_STORAGE .and. &
        channel /= TTS_CHANNEL_LOSS) return
    reference_segments = build_tts_valid_log_segments( &
      reference_isotherm, channel)
    moving_segments = build_tts_valid_log_segments(moving_isotherm, channel)
    if (size(reference_segments) == 0 .or. size(moving_segments) == 0) return

    call integrate_segment_overlaps(reference_segments, moving_segments, &
      delta_s, integral, overlap_measure, interval_count)
    if (interval_count <= 0 .or. overlap_measure <= 0.0_dp .or. &
        .not. ieee_is_finite(integral)) return
    reference_measure = total_segment_measure(reference_segments)
    moving_measure = total_segment_measure(moving_segments)
    if (reference_measure <= 0.0_dp .or. moving_measure <= 0.0_dp) return

    evaluation%valid = .true.
    evaluation%objective = integral/overlap_measure
    evaluation%overlap_width_decades = overlap_measure
    evaluation%overlap_fraction = overlap_measure / &
      min(reference_measure, moving_measure)
    evaluation%interpolation_interval_count = interval_count
  end function evaluate_tts_pair_channel_objective

  !> Storage, loss veya joint production objective'ini değerlendirir. Joint
  !! objective yalnız iki channel da mathematical interval support taşıyorsa
  !! J=0.5(J_storage+J_loss) olur; her channel kendi overlap measure'ıyla
  !! önceden normalize edilmiştir.
  pure function evaluate_tts_pair_objective( &
      reference_isotherm, moving_isotherm, delta_s, channel) &
      result(evaluation)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    real(dp), intent(in) :: delta_s
    integer, intent(in) :: channel
    type(tts_pair_objective_evaluation_t) :: evaluation

    type(tts_pair_objective_evaluation_t) :: storage
    type(tts_pair_objective_evaluation_t) :: loss

    select case (channel)
    case (TTS_CHANNEL_STORAGE, TTS_CHANNEL_LOSS)
      evaluation = evaluate_tts_pair_channel_objective( &
        reference_isotherm, moving_isotherm, delta_s, channel)
    case (TTS_CHANNEL_JOINT)
      storage = evaluate_tts_pair_channel_objective(reference_isotherm, &
        moving_isotherm, delta_s, TTS_CHANNEL_STORAGE)
      loss = evaluate_tts_pair_channel_objective(reference_isotherm, &
        moving_isotherm, delta_s, TTS_CHANNEL_LOSS)
      if (.not. storage%valid .or. .not. loss%valid) return
      evaluation%valid = .true.
      evaluation%objective = 0.5_dp*(storage%objective + loss%objective)
      evaluation%overlap_width_decades = min( &
        storage%overlap_width_decades, loss%overlap_width_decades)
      evaluation%overlap_fraction = min( &
        storage%overlap_fraction, loss%overlap_fraction)
      evaluation%interpolation_interval_count = &
        storage%interpolation_interval_count + loss%interpolation_interval_count
    end select
  end function evaluate_tts_pair_objective

  !> İki adjacent measured isotherm arasındaki production horizontal shift'i
  !! deterministic coarse scan ve yalnız valid interior three-point bracket
  !! üzerine kurulan Brent ile belirler. G' ve G'' support varsa joint shift;
  !! positive-valid G'' support yoksa açık storage-only fallback döner.
  !! Başarı full thermorheological simplicity kanıtı değildir.
  pure function identify_tts_pair_shift( &
      reference_isotherm, moving_isotherm, configuration) result(pair_result)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_pair_shift_configuration_t), intent(in), optional :: configuration
    type(tts_pair_shift_result_t) :: pair_result

    type(tts_pair_shift_configuration_t) :: settings
    type(channel_optimization_t) :: joint
    type(channel_optimization_t) :: loss
    type(channel_optimization_t) :: storage

    if (present(configuration)) settings = configuration
    pair_result%reference_isotherm_identifier = &
      reference_isotherm%isotherm_identifier
    pair_result%moving_isotherm_identifier = &
      moving_isotherm%isotherm_identifier
    if (.not. is_pair_configuration_valid(settings) .or. &
        .not. is_isotherm_structurally_valid(reference_isotherm) .or. &
        .not. is_isotherm_structurally_valid(moving_isotherm)) then
      pair_result%status = PAIR_SHIFT_INVALID_INPUT
      return
    end if

    storage = optimize_pair_channel(reference_isotherm, moving_isotherm, &
      TTS_CHANNEL_STORAGE, settings)
    loss = optimize_pair_channel(reference_isotherm, moving_isotherm, &
      TTS_CHANNEL_LOSS, settings)
    joint = optimize_pair_channel(reference_isotherm, moving_isotherm, &
      TTS_CHANNEL_JOINT, settings)
    pair_result%storage_shift_available = storage%available
    pair_result%loss_shift_available = loss%available
    pair_result%joint_shift_available = joint%available
    if (storage%available) then
      pair_result%delta_s_storage = storage%shift
      pair_result%storage_overlap_width_decades = &
        storage%overlap_width_decades
    end if
    if (loss%available) then
      pair_result%delta_s_loss = loss%shift
      pair_result%loss_overlap_width_decades = loss%overlap_width_decades
    end if
    if (storage%available .and. loss%available) then
      pair_result%storage_loss_shift_discrepancy = &
        abs(storage%shift - loss%shift)
    end if

    if (joint%available) then
      pair_result%status = PAIR_SHIFT_SUCCESS
      pair_result%production_channel = SHIFT_FROM_JOINT
      pair_result%shift_available = .true.
      pair_result%delta_s = joint%shift
      pair_result%delta_s_joint = joint%shift
      call copy_production_diagnostics(pair_result, joint)
      return
    end if
    if (storage%available .and. .not. loss%any_valid_objective) then
      pair_result%status = PAIR_SHIFT_STORAGE_ONLY
      pair_result%production_channel = SHIFT_FROM_STORAGE_ONLY
      pair_result%shift_available = .true.
      pair_result%delta_s = storage%shift
      call copy_production_diagnostics(pair_result, storage)
      return
    end if

    pair_result%status = merge_pair_failure_status(storage, loss, joint)
    pair_result%iteration_count = joint%iteration_count
    pair_result%evaluation_count = storage%evaluation_count + &
      loss%evaluation_count + joint%evaluation_count
  end function identify_tts_pair_shift

  pure subroutine populate_segment( &
      segment, isotherm, channel, run_start, run_end)
    type(tts_log_segment_t), intent(out) :: segment
    type(tts_isotherm_t), intent(in) :: isotherm
    integer, intent(in) :: channel
    integer, intent(in) :: run_start
    integer, intent(in) :: run_end
    integer :: local_index
    integer :: source_index

    segment%channel = channel
    allocate(segment%x(run_end - run_start + 1))
    allocate(segment%y(run_end - run_start + 1))
    allocate(segment%source_point_indices(run_end - run_start + 1))
    do source_index = run_start, run_end
      local_index = source_index - run_start + 1
      segment%x(local_index) = log10( &
        isotherm%points(source_index)%frequency_hz)
      if (channel == TTS_CHANNEL_STORAGE) then
        segment%y(local_index) = log10( &
          isotherm%points(source_index)%storage_modulus_pa)
      else
        segment%y(local_index) = log10( &
          isotherm%points(source_index)%loss_modulus_pa)
      end if
      segment%source_point_indices(local_index) = source_index
    end do
  end subroutine populate_segment

  pure elemental function point_is_log_usable(point, channel) result(usable)
    type(tts_measurement_point_t), intent(in) :: point
    integer, intent(in) :: channel
    logical :: usable

    select case (channel)
    case (TTS_CHANNEL_STORAGE)
      usable = is_storage_log_usable(point)
    case (TTS_CHANNEL_LOSS)
      usable = is_loss_log_usable(point)
    case default
      usable = .false.
    end select
  end function point_is_log_usable

  pure subroutine integrate_segment_overlaps( &
      reference_segments, moving_segments, delta_s, integral, &
      overlap_measure, interval_count)
    type(tts_log_segment_t), intent(in) :: reference_segments(:)
    type(tts_log_segment_t), intent(in) :: moving_segments(:)
    real(dp), intent(in) :: delta_s
    real(dp), intent(out) :: integral
    real(dp), intent(out) :: overlap_measure
    integer, intent(out) :: interval_count

    integer :: i
    integer :: j
    real(dp) :: pair_integral
    real(dp) :: pair_measure
    integer :: pair_intervals

    integral = 0.0_dp
    overlap_measure = 0.0_dp
    interval_count = 0
    do i = 1, size(reference_segments)
      do j = 1, size(moving_segments)
        call integrate_segment_pair(reference_segments(i), &
          moving_segments(j), delta_s, pair_integral, pair_measure, &
          pair_intervals)
        integral = integral + pair_integral
        overlap_measure = overlap_measure + pair_measure
        interval_count = interval_count + pair_intervals
      end do
    end do
  end subroutine integrate_segment_overlaps

  pure subroutine integrate_segment_pair( &
      reference_segment, moving_segment, delta_s, integral, &
      overlap_measure, interval_count)
    type(tts_log_segment_t), intent(in) :: reference_segment
    type(tts_log_segment_t), intent(in) :: moving_segment
    real(dp), intent(in) :: delta_s
    real(dp), intent(out) :: integral
    real(dp), intent(out) :: overlap_measure
    integer, intent(out) :: interval_count

    real(dp), allocatable :: knots(:)
    real(dp) :: h
    real(dp) :: lower_overlap
    real(dp) :: r0
    real(dp) :: r1
    real(dp) :: upper_overlap
    integer :: i
    integer :: knot_count

    integral = 0.0_dp
    overlap_measure = 0.0_dp
    interval_count = 0
    lower_overlap = max(reference_segment%x(1), &
      moving_segment%x(1) + delta_s)
    upper_overlap = min(reference_segment%x(size(reference_segment%x)), &
      moving_segment%x(size(moving_segment%x)) + delta_s)
    if (upper_overlap <= lower_overlap .or. &
        are_tts_values_machine_equivalent(upper_overlap, lower_overlap)) return

    allocate(knots(size(reference_segment%x) + size(moving_segment%x) + 2))
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

    do i = 1, knot_count - 1
      h = knots(i + 1) - knots(i)
      if (h <= 0.0_dp .or. &
          are_tts_values_machine_equivalent(knots(i + 1), knots(i))) cycle
      r0 = interpolate_log_segment(reference_segment, knots(i)) - &
        interpolate_log_segment(moving_segment, knots(i) - delta_s)
      r1 = interpolate_log_segment(reference_segment, knots(i + 1)) - &
        interpolate_log_segment(moving_segment, knots(i + 1) - delta_s)
      integral = integral + integrate_linear_squared_residual(h, r0, r1)
      overlap_measure = overlap_measure + h
      interval_count = interval_count + 1
    end do
  end subroutine integrate_segment_pair

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

  pure function interpolate_log_segment(segment, x) result(y)
    type(tts_log_segment_t), intent(in) :: segment
    real(dp), intent(in) :: x
    real(dp) :: y
    real(dp) :: alpha
    integer :: i

    if (are_tts_values_machine_equivalent(x, segment%x(1))) then
      y = segment%y(1)
      return
    end if
    if (are_tts_values_machine_equivalent( &
        x, segment%x(size(segment%x)))) then
      y = segment%y(size(segment%y))
      return
    end if
    do i = 1, size(segment%x) - 1
      if (x >= segment%x(i) .and. x <= segment%x(i + 1)) then
        alpha = (x - segment%x(i))/(segment%x(i + 1) - segment%x(i))
        y = (1.0_dp - alpha)*segment%y(i) + alpha*segment%y(i + 1)
        return
      end if
    end do
    y = huge(1.0_dp)
  end function interpolate_log_segment

  pure function total_segment_measure(segments) result(measure)
    type(tts_log_segment_t), intent(in) :: segments(:)
    real(dp) :: measure
    integer :: i

    measure = 0.0_dp
    do i = 1, size(segments)
      measure = measure + segments(i)%x(size(segments(i)%x)) - segments(i)%x(1)
    end do
  end function total_segment_measure

  pure function optimize_pair_channel( &
      reference_isotherm, moving_isotherm, channel, settings) result(solution)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    integer, intent(in) :: channel
    type(tts_pair_shift_configuration_t), intent(in) :: settings
    type(channel_optimization_t) :: solution

    real(dp), allocatable :: objective_values(:)
    logical, allocatable :: objective_valid(:)
    real(dp), allocatable :: scan_shifts(:)
    type(tts_pair_objective_evaluation_t) :: evaluation
    type(tts_pair_objective_evaluation_t) :: left_curvature
    type(tts_pair_objective_evaluation_t) :: right_curvature
    type(scalar_minimizer_result_t) :: minimizer
    real(dp) :: curvature_step
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
      solution%status = PAIR_SHIFT_NO_OVERLAP
      return
    end if
    scale = max(1.0_dp, abs(lower_shift), abs(upper_shift))
    margin = 128.0_dp*epsilon(1.0_dp)*scale
    interior_lower = lower_shift + margin
    interior_upper = upper_shift - margin
    if (interior_lower >= interior_upper) then
      solution%status = PAIR_SHIFT_NO_OVERLAP
      return
    end if

    allocate(scan_shifts(settings%coarse_scan_point_count))
    allocate(objective_values(settings%coarse_scan_point_count))
    allocate(objective_valid(settings%coarse_scan_point_count))
    do i = 1, settings%coarse_scan_point_count
      scan_shifts(i) = interior_lower + real(i - 1, dp) * &
        (interior_upper - interior_lower) / &
        real(settings%coarse_scan_point_count - 1, dp)
      evaluation = evaluate_tts_pair_objective(reference_isotherm, &
        moving_isotherm, scan_shifts(i), channel)
      objective_valid(i) = evaluation%valid
      objective_values(i) = evaluation%objective
    end do
    solution%evaluation_count = settings%coarse_scan_point_count
    solution%any_valid_objective = any(objective_valid)
    if (.not. solution%any_valid_objective) then
      solution%status = PAIR_SHIFT_INSUFFICIENT_SUPPORT
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
      solution%status = PAIR_SHIFT_NO_INTERIOR_MINIMUM
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
      solution%status = PAIR_SHIFT_OPTIMIZATION_FAILED
      return
    end if
    evaluation = evaluate_tts_pair_objective(reference_isotherm, &
      moving_isotherm, minimizer%x_minimum, channel)
    solution%evaluation_count = solution%evaluation_count + 1
    if (.not. evaluation%valid) then
      solution%status = PAIR_SHIFT_OPTIMIZATION_FAILED
      return
    end if

    solution%available = .true.
    solution%status = PAIR_SHIFT_SUCCESS
    solution%shift = minimizer%x_minimum
    solution%objective = evaluation%objective
    solution%overlap_width_decades = evaluation%overlap_width_decades
    solution%overlap_fraction = evaluation%overlap_fraction
    curvature_step = min(0.25_dp*(interior_upper - interior_lower), &
      max(32.0_dp*sqrt(epsilon(1.0_dp)) * &
        (1.0_dp + abs(solution%shift)), &
        1.0e-4_dp*(interior_upper - interior_lower)))
    if (solution%shift - curvature_step > interior_lower .and. &
        solution%shift + curvature_step < interior_upper) then
      left_curvature = evaluate_tts_pair_objective(reference_isotherm, &
        moving_isotherm, solution%shift - curvature_step, channel)
      right_curvature = evaluate_tts_pair_objective(reference_isotherm, &
        moving_isotherm, solution%shift + curvature_step, channel)
      solution%evaluation_count = solution%evaluation_count + 2
      if (left_curvature%valid .and. right_curvature%valid) then
        solution%curvature = (left_curvature%objective - &
          2.0_dp*solution%objective + right_curvature%objective) / &
          (curvature_step*curvature_step)
      end if
    end if

  contains

    pure function objective_at_shift(shift) result(value)
      real(dp), intent(in) :: shift
      real(dp) :: value
      type(tts_pair_objective_evaluation_t) :: local_evaluation

      local_evaluation = evaluate_tts_pair_objective(reference_isotherm, &
        moving_isotherm, shift, channel)
      if (local_evaluation%valid) then
        value = local_evaluation%objective
      else
        value = huge(1.0_dp)
      end if
    end function objective_at_shift

  end function optimize_pair_channel

  pure function is_pair_configuration_valid(settings) result(valid)
    type(tts_pair_shift_configuration_t), intent(in) :: settings
    logical :: valid

    valid = settings%coarse_scan_point_count >= 3 .and. &
      ieee_is_finite(settings%absolute_tolerance) .and. &
      ieee_is_finite(settings%relative_tolerance) .and. &
      settings%absolute_tolerance > 0.0_dp .and. &
      settings%relative_tolerance >= 0.0_dp .and. &
      settings%maximum_iterations > 0
  end function is_pair_configuration_valid

  pure function is_isotherm_structurally_valid(isotherm) result(valid)
    type(tts_isotherm_t), intent(in) :: isotherm
    logical :: valid
    integer :: i

    valid = allocated(isotherm%points)
    if (.not. valid) return
    valid = size(isotherm%points) >= 2 .and. &
      ieee_is_finite(isotherm%temperature_k) .and. &
      isotherm%temperature_k > 0.0_dp
    if (.not. valid) return
    do i = 1, size(isotherm%points)
      valid = ieee_is_finite(isotherm%points(i)%frequency_hz) .and. &
        isotherm%points(i)%frequency_hz > 0.0_dp
      if (.not. valid) return
      if (i > 1) then
        valid = isotherm%points(i)%frequency_hz > &
          isotherm%points(i - 1)%frequency_hz .and. &
          .not. are_tts_values_machine_equivalent( &
            isotherm%points(i)%frequency_hz, &
            isotherm%points(i - 1)%frequency_hz)
        if (.not. valid) return
      end if
    end do
  end function is_isotherm_structurally_valid

  pure subroutine copy_production_diagnostics(pair_result, solution)
    type(tts_pair_shift_result_t), intent(inout) :: pair_result
    type(channel_optimization_t), intent(in) :: solution

    pair_result%objective_minimum = solution%objective
    pair_result%overlap_width_decades = solution%overlap_width_decades
    pair_result%overlap_fraction = solution%overlap_fraction
    pair_result%objective_curvature = solution%curvature
    pair_result%iteration_count = solution%iteration_count
    pair_result%evaluation_count = solution%evaluation_count
  end subroutine copy_production_diagnostics

  pure function merge_pair_failure_status(storage, loss, joint) result(status)
    type(channel_optimization_t), intent(in) :: storage
    type(channel_optimization_t), intent(in) :: loss
    type(channel_optimization_t), intent(in) :: joint
    integer :: status

    if (joint%status == PAIR_SHIFT_OPTIMIZATION_FAILED .or. &
        storage%status == PAIR_SHIFT_OPTIMIZATION_FAILED .or. &
        loss%status == PAIR_SHIFT_OPTIMIZATION_FAILED) then
      status = PAIR_SHIFT_OPTIMIZATION_FAILED
    else if (joint%status == PAIR_SHIFT_NO_INTERIOR_MINIMUM .or. &
        storage%status == PAIR_SHIFT_NO_INTERIOR_MINIMUM .or. &
        loss%status == PAIR_SHIFT_NO_INTERIOR_MINIMUM) then
      status = PAIR_SHIFT_NO_INTERIOR_MINIMUM
    else if (joint%status == PAIR_SHIFT_NO_OVERLAP .or. &
        storage%status == PAIR_SHIFT_NO_OVERLAP) then
      status = PAIR_SHIFT_NO_OVERLAP
    else
      status = PAIR_SHIFT_INSUFFICIENT_SUPPORT
    end if
  end function merge_pair_failure_status

end module tms_tts_pair_shift
