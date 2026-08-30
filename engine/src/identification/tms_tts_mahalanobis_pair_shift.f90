module tms_tts_mahalanobis_pair_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_scalar_minimizer, only : scalar_minimizer_result_t, &
    SCALAR_MINIMIZER_SUCCESS, minimize_scalar_brent, &
    is_valid_scalar_minimum_bracket
  use tms_bivariate_covariance_integrals, only : &
    bivariate_covariance_interval_integral_t, &
    integrate_bivariate_covariance_linear_interval
  use tms_tts_types, only : tts_isotherm_t, &
    tts_pair_shift_configuration_t, are_tts_values_machine_equivalent
  use tms_tts_pair_shift, only : get_tts_pair_feasible_shift_domain
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t
  use tms_tts_covariance_types, only : &
    tts_dynamic_modulus_covariance_family_t, &
    tts_bivariate_covariance_log_segment_t, &
    tts_covariance_objective_evaluation_t, &
    tts_covariance_shift_solution_t, tts_covariance_pair_solution_t, &
    TTS_MATCHED_DIAGONAL_OBJECTIVE, TTS_MAHALANOBIS_OBJECTIVE, &
    MAHALANOBIS_SHIFT_SUCCESS, MAHALANOBIS_SHIFT_NO_OVERLAP, &
    MAHALANOBIS_SHIFT_INVALID_COVARIANCE, &
    MAHALANOBIS_SHIFT_NO_INTERIOR_MINIMUM, &
    MAHALANOBIS_SHIFT_NUMERICAL_FAILURE
  use tms_tts_covariance_validation, only : &
    build_tts_bivariate_covariance_log_segments
  implicit none
  private

  public :: evaluate_tts_covariance_pair_objective
  public :: identify_tts_covariance_pair_shift

contains

  !> Aynı bivariate common support O_B üzerinde matched-diagonal ve
  !! Mahalanobis objective'lerini değerlendirir. Moving x koordinatı mevcut
  !! V0.8.1 convention'ıyla x_shifted=x+delta_s taşınır. Her iki sonuç
  !! (2*|O_B|) ile normalize edilir; gap genişliği ölçüye katılmaz.
  pure function evaluate_tts_covariance_pair_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      covariance_family, delta_s) result(evaluation)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    real(dp), intent(in) :: delta_s
    type(tts_covariance_objective_evaluation_t) :: evaluation

    type(tts_bivariate_covariance_log_segment_t), allocatable :: &
      reference_segments(:)
    type(tts_bivariate_covariance_log_segment_t), allocatable :: &
      moving_segments(:)
    real(dp) :: diagonal_integral
    real(dp) :: mahalanobis_integral
    real(dp) :: overlap_width
    integer :: interval_count
    logical :: numerical_failure

    if (.not. ieee_is_finite(delta_s)) return
    reference_segments = build_tts_bivariate_covariance_log_segments( &
      reference_isotherm, uncertainty_family, covariance_family)
    moving_segments = build_tts_bivariate_covariance_log_segments( &
      moving_isotherm, uncertainty_family, covariance_family)
    if (size(reference_segments) == 0 .or. size(moving_segments) == 0) return
    call integrate_segment_overlaps(reference_segments, moving_segments, &
      delta_s, diagonal_integral, mahalanobis_integral, overlap_width, &
      interval_count, numerical_failure)
    if (numerical_failure) then
      evaluation%numerical_failure = .true.
      return
    end if
    if (interval_count <= 0 .or. overlap_width <= 0.0_dp) return
    evaluation%matched_diagonal_objective = &
      0.5_dp*diagonal_integral/overlap_width
    evaluation%mahalanobis_objective = &
      0.5_dp*mahalanobis_integral/overlap_width
    if (.not. ieee_is_finite(evaluation%matched_diagonal_objective) .or. &
        .not. ieee_is_finite(evaluation%mahalanobis_objective) .or. &
        evaluation%matched_diagonal_objective < 0.0_dp .or. &
        evaluation%mahalanobis_objective < 0.0_dp) then
      evaluation%numerical_failure = .true.
      evaluation%matched_diagonal_objective = huge(1.0_dp)
      evaluation%mahalanobis_objective = huge(1.0_dp)
      return
    end if
    evaluation%valid = .true.
    evaluation%overlap_width_decades = overlap_width
    evaluation%interpolation_interval_count = interval_count
  end function evaluate_tts_covariance_pair_objective

  !> Matched-diagonal ve Mahalanobis shift'lerini ayrı objective'lerle fakat
  !! aynı O_B support kurallarıyla çözer. Deterministic coarse scan, genuine
  !! strict-interior bracket ve mevcut Brent minimizer reuse edilir. Flat veya
  !! boundary-only minimum başarıya çevrilmez; storage-only fallback yoktur.
  pure function identify_tts_covariance_pair_shift( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      covariance_family, configuration) result(pair_result)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    type(tts_pair_shift_configuration_t), intent(in) :: configuration
    type(tts_covariance_pair_solution_t) :: pair_result

    if (.not. pair_configuration_is_valid(configuration)) return
    pair_result%diagonal_matched = optimize_pair_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      covariance_family, TTS_MATCHED_DIAGONAL_OBJECTIVE, configuration)
    pair_result%mahalanobis = optimize_pair_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      covariance_family, TTS_MAHALANOBIS_OBJECTIVE, configuration)
  end function identify_tts_covariance_pair_shift

  pure subroutine integrate_segment_overlaps( &
      reference_segments, moving_segments, delta_s, diagonal_integral, &
      mahalanobis_integral, overlap_width, interval_count, numerical_failure)
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: &
      reference_segments(:)
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: &
      moving_segments(:)
    real(dp), intent(in) :: delta_s
    real(dp), intent(out) :: diagonal_integral
    real(dp), intent(out) :: mahalanobis_integral
    real(dp), intent(out) :: overlap_width
    integer, intent(out) :: interval_count
    logical, intent(out) :: numerical_failure

    real(dp) :: pair_diagonal
    real(dp) :: pair_mahalanobis
    real(dp) :: pair_width
    integer :: i
    integer :: j
    integer :: pair_intervals
    logical :: pair_failure

    diagonal_integral = 0.0_dp
    mahalanobis_integral = 0.0_dp
    overlap_width = 0.0_dp
    interval_count = 0
    numerical_failure = .false.
    do i = 1, size(reference_segments)
      do j = 1, size(moving_segments)
        call integrate_segment_pair(reference_segments(i), &
          moving_segments(j), delta_s, pair_diagonal, pair_mahalanobis, &
          pair_width, pair_intervals, pair_failure)
        if (pair_failure) then
          numerical_failure = .true.
          return
        end if
        diagonal_integral = diagonal_integral + pair_diagonal
        mahalanobis_integral = mahalanobis_integral + pair_mahalanobis
        overlap_width = overlap_width + pair_width
        interval_count = interval_count + pair_intervals
      end do
    end do
  end subroutine integrate_segment_overlaps

  pure subroutine integrate_segment_pair( &
      reference_segment, moving_segment, delta_s, diagonal_integral, &
      mahalanobis_integral, overlap_width, interval_count, numerical_failure)
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: &
      reference_segment
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: moving_segment
    real(dp), intent(in) :: delta_s
    real(dp), intent(out) :: diagonal_integral
    real(dp), intent(out) :: mahalanobis_integral
    real(dp), intent(out) :: overlap_width
    integer, intent(out) :: interval_count
    logical, intent(out) :: numerical_failure

    type(bivariate_covariance_interval_integral_t) :: interval
    real(dp), allocatable :: knots(:)
    real(dp) :: h
    real(dp) :: lower_overlap
    real(dp) :: moving_x0
    real(dp) :: moving_x1
    real(dp) :: r_l0
    real(dp) :: r_l1
    real(dp) :: r_s0
    real(dp) :: r_s1
    real(dp) :: reference_c0
    real(dp) :: reference_c1
    real(dp) :: reference_lv0
    real(dp) :: reference_lv1
    real(dp) :: reference_sv0
    real(dp) :: reference_sv1
    real(dp) :: upper_overlap
    real(dp) :: covariance0
    real(dp) :: covariance1
    real(dp) :: loss_variance0
    real(dp) :: loss_variance1
    real(dp) :: storage_variance0
    real(dp) :: storage_variance1
    integer :: i
    integer :: knot_count

    diagonal_integral = 0.0_dp
    mahalanobis_integral = 0.0_dp
    overlap_width = 0.0_dp
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
      moving_x0 = knots(i) - delta_s
      moving_x1 = knots(i + 1) - delta_s
      r_s0 = interpolate_field(reference_segment, knots(i), 1) - &
        interpolate_field(moving_segment, moving_x0, 1)
      r_s1 = interpolate_field(reference_segment, knots(i + 1), 1) - &
        interpolate_field(moving_segment, moving_x1, 1)
      r_l0 = interpolate_field(reference_segment, knots(i), 2) - &
        interpolate_field(moving_segment, moving_x0, 2)
      r_l1 = interpolate_field(reference_segment, knots(i + 1), 2) - &
        interpolate_field(moving_segment, moving_x1, 2)
      reference_sv0 = interpolate_field(reference_segment, knots(i), 3)
      reference_sv1 = interpolate_field(reference_segment, knots(i + 1), 3)
      reference_lv0 = interpolate_field(reference_segment, knots(i), 4)
      reference_lv1 = interpolate_field(reference_segment, knots(i + 1), 4)
      reference_c0 = interpolate_field(reference_segment, knots(i), 5)
      reference_c1 = interpolate_field(reference_segment, knots(i + 1), 5)
      storage_variance0 = reference_sv0 + &
        interpolate_field(moving_segment, moving_x0, 3)
      storage_variance1 = reference_sv1 + &
        interpolate_field(moving_segment, moving_x1, 3)
      loss_variance0 = reference_lv0 + &
        interpolate_field(moving_segment, moving_x0, 4)
      loss_variance1 = reference_lv1 + &
        interpolate_field(moving_segment, moving_x1, 4)
      covariance0 = reference_c0 + &
        interpolate_field(moving_segment, moving_x0, 5)
      covariance1 = reference_c1 + &
        interpolate_field(moving_segment, moving_x1, 5)
      interval = integrate_bivariate_covariance_linear_interval(h, &
        r_s0, r_s1, r_l0, r_l1, storage_variance0, storage_variance1, &
        loss_variance0, loss_variance1, covariance0, covariance1)
      if (.not. interval%valid) then
        numerical_failure = .true.
        return
      end if
      diagonal_integral = diagonal_integral + &
        interval%matched_diagonal_integral
      mahalanobis_integral = mahalanobis_integral + &
        interval%mahalanobis_integral
      overlap_width = overlap_width + h
      interval_count = interval_count + 1
    end do
  end subroutine integrate_segment_pair

  pure subroutine build_merged_knots( &
      reference_segment, moving_segment, delta_s, lower_overlap, &
      upper_overlap, knots, knot_count)
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: &
      reference_segment
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: moving_segment
    real(dp), intent(in) :: delta_s
    real(dp), intent(in) :: lower_overlap
    real(dp), intent(in) :: upper_overlap
    real(dp), intent(out) :: knots(:)
    integer, intent(out) :: knot_count
    integer :: i

    knot_count = 0
    call append_knot(knots, knot_count, lower_overlap)
    do i = 1, size(reference_segment%x)
      if (reference_segment%x(i) > lower_overlap .and. &
          reference_segment%x(i) < upper_overlap) &
        call append_knot(knots, knot_count, reference_segment%x(i))
    end do
    do i = 1, size(moving_segment%x)
      if (moving_segment%x(i) + delta_s > lower_overlap .and. &
          moving_segment%x(i) + delta_s < upper_overlap) &
        call append_knot(knots, knot_count, &
          moving_segment%x(i) + delta_s)
    end do
    call append_knot(knots, knot_count, upper_overlap)
    call sort_knots(knots, knot_count)
  end subroutine build_merged_knots

  pure subroutine append_knot(knots, knot_count, value)
    real(dp), intent(inout) :: knots(:)
    integer, intent(inout) :: knot_count
    real(dp), intent(in) :: value
    integer :: i

    do i = 1, knot_count
      if (are_tts_values_machine_equivalent(knots(i), value)) return
    end do
    knot_count = knot_count + 1
    knots(knot_count) = value
  end subroutine append_knot

  pure subroutine sort_knots(knots, knot_count)
    real(dp), intent(inout) :: knots(:)
    integer, intent(in) :: knot_count
    real(dp) :: key
    integer :: i
    integer :: j

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
  end subroutine sort_knots

  pure function interpolate_field(segment, x, field) result(value)
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: segment
    real(dp), intent(in) :: x
    integer, intent(in) :: field
    real(dp) :: value
    real(dp) :: alpha
    real(dp) :: left_value
    real(dp) :: right_value
    integer :: i

    value = huge(1.0_dp)
    if (are_tts_values_machine_equivalent(x, segment%x(1))) then
      value = field_at_index(segment, 1, field)
      return
    end if
    if (are_tts_values_machine_equivalent( &
        x, segment%x(size(segment%x)))) then
      value = field_at_index(segment, size(segment%x), field)
      return
    end if
    do i = 1, size(segment%x) - 1
      if ((x < segment%x(i) .and. .not. &
          are_tts_values_machine_equivalent(x, segment%x(i))) .or. &
          (x > segment%x(i + 1) .and. .not. &
          are_tts_values_machine_equivalent(x, segment%x(i + 1)))) cycle
      alpha = (x - segment%x(i))/(segment%x(i + 1) - segment%x(i))
      alpha = max(0.0_dp, min(1.0_dp, alpha))
      left_value = field_at_index(segment, i, field)
      right_value = field_at_index(segment, i + 1, field)
      if (.not. ieee_is_finite(left_value) .or. &
          .not. ieee_is_finite(right_value)) return
      value = (1.0_dp - alpha)*left_value + alpha*right_value
      return
    end do
  end function interpolate_field

  pure function field_at_index(segment, index, field) result(value)
    type(tts_bivariate_covariance_log_segment_t), intent(in) :: segment
    integer, intent(in) :: index
    integer, intent(in) :: field
    real(dp) :: value

    value = huge(1.0_dp)
    select case (field)
    case (1)
      value = segment%storage_y(index)
    case (2)
      value = segment%loss_y(index)
    case (3)
      value = segment%storage_variance(index)
    case (4)
      value = segment%loss_variance(index)
    case (5)
      value = segment%storage_loss_covariance(index)
    end select
  end function field_at_index

  pure function optimize_pair_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      covariance_family, objective_mode, settings) result(solution)
    type(tts_isotherm_t), intent(in) :: reference_isotherm
    type(tts_isotherm_t), intent(in) :: moving_isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    integer, intent(in) :: objective_mode
    type(tts_pair_shift_configuration_t), intent(in) :: settings
    type(tts_covariance_shift_solution_t) :: solution

    real(dp), allocatable :: objective_values(:)
    logical, allocatable :: objective_valid(:)
    real(dp), allocatable :: scan_shifts(:)
    type(tts_covariance_objective_evaluation_t) :: evaluation
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
    logical :: numerical_failure_observed

    solution%objective_mode = objective_mode
    if (objective_mode /= TTS_MATCHED_DIAGONAL_OBJECTIVE .and. &
        objective_mode /= TTS_MAHALANOBIS_OBJECTIVE) return
    call get_tts_pair_feasible_shift_domain(reference_isotherm, &
      moving_isotherm, lower_shift, upper_shift, domain_valid)
    if (.not. domain_valid) then
      solution%status = MAHALANOBIS_SHIFT_NO_OVERLAP
      return
    end if
    scale = max(1.0_dp, abs(lower_shift), abs(upper_shift))
    margin = 128.0_dp*epsilon(1.0_dp)*scale
    interior_lower = lower_shift + margin
    interior_upper = upper_shift - margin
    if (interior_lower >= interior_upper) then
      solution%status = MAHALANOBIS_SHIFT_NO_OVERLAP
      return
    end if

    allocate(scan_shifts(settings%coarse_scan_point_count))
    allocate(objective_values(settings%coarse_scan_point_count))
    allocate(objective_valid(settings%coarse_scan_point_count))
    numerical_failure_observed = .false.
    do i = 1, settings%coarse_scan_point_count
      scan_shifts(i) = interior_lower + real(i - 1, dp) * &
        (interior_upper - interior_lower) / &
        real(settings%coarse_scan_point_count - 1, dp)
      evaluation = evaluate_tts_covariance_pair_objective( &
        reference_isotherm, moving_isotherm, uncertainty_family, &
        covariance_family, scan_shifts(i))
      objective_valid(i) = evaluation%valid
      objective_values(i) = selected_objective(evaluation, objective_mode)
      numerical_failure_observed = numerical_failure_observed .or. &
        evaluation%numerical_failure
    end do
    solution%evaluation_count = settings%coarse_scan_point_count
    if (.not. any(objective_valid)) then
      if (numerical_failure_observed) then
        solution%status = MAHALANOBIS_SHIFT_NUMERICAL_FAILURE
      else
        solution%status = MAHALANOBIS_SHIFT_INVALID_COVARIANCE
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
        ! Fortran mantıksal ifadelerde kısa devre değerlendirmeyi garanti etmez.
        ! İlk geçerli bracket bulunduğunda sıfır indisli dizi erişimini önlemek
        ! için başlangıç ve karşılaştırma durumları açıkça ayrılır.
        if (bracket_index == 0) then
          bracket_index = i
        else if (objective_values(i) < objective_values(bracket_index)) then
          bracket_index = i
        end if
      end if
    end do
    if (bracket_index == 0) then
      solution%status = MAHALANOBIS_SHIFT_NO_INTERIOR_MINIMUM
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
      solution%status = MAHALANOBIS_SHIFT_NUMERICAL_FAILURE
      return
    end if
    evaluation = evaluate_tts_covariance_pair_objective( &
      reference_isotherm, moving_isotherm, uncertainty_family, &
      covariance_family, minimizer%x_minimum)
    solution%evaluation_count = solution%evaluation_count + 1
    if (.not. evaluation%valid) then
      solution%status = MAHALANOBIS_SHIFT_NUMERICAL_FAILURE
      return
    end if
    solution%status = MAHALANOBIS_SHIFT_SUCCESS
    solution%shift_available = .true.
    solution%shift = minimizer%x_minimum
    solution%objective_minimum = selected_objective(evaluation, objective_mode)
    solution%diagnostics = evaluation

  contains

    pure function objective_at_shift(shift) result(value)
      real(dp), intent(in) :: shift
      real(dp) :: value
      type(tts_covariance_objective_evaluation_t) :: local_evaluation

      local_evaluation = evaluate_tts_covariance_pair_objective( &
        reference_isotherm, moving_isotherm, uncertainty_family, &
        covariance_family, shift)
      if (local_evaluation%valid) then
        value = selected_objective(local_evaluation, objective_mode)
      else
        value = huge(1.0_dp)
      end if
    end function objective_at_shift

  end function optimize_pair_objective

  pure function selected_objective(evaluation, objective_mode) result(value)
    type(tts_covariance_objective_evaluation_t), intent(in) :: evaluation
    integer, intent(in) :: objective_mode
    real(dp) :: value

    if (objective_mode == TTS_MATCHED_DIAGONAL_OBJECTIVE) then
      value = evaluation%matched_diagonal_objective
    else
      value = evaluation%mahalanobis_objective
    end if
  end function selected_objective

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

end module tms_tts_mahalanobis_pair_shift
