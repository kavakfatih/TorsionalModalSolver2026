module tms_tts_master_curve
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_empirical_shift_t, &
    tts_master_cloud_point_t, tts_runtime_master_point_t, &
    tts_master_boundary_diagnostic_t, is_runtime_export_usable, &
    is_storage_log_usable, is_loss_log_usable, &
    are_tts_values_machine_equivalent
  implicit none
  private

  public :: build_tts_master_experimental_cloud
  public :: stitch_tts_runtime_master_table

contains

  !> Bütün original measured point'leri koruyarak shifted experimental master
  !! cloud üretir. f_r=10^(log10(f)+s) [Hz], s=log10(a_T) boyutsuzdur.
  !! G'/G'' ve quality/provenance overwrite edilmez. Extreme log-coordinate
  !! representable değilse success=false döner; clamp uygulanmaz.
  pure subroutine build_tts_master_experimental_cloud( &
      family, empirical_shifts, cloud, success)
    type(tts_material_family_t), intent(in) :: family
    type(tts_empirical_shift_t), intent(in) :: empirical_shifts(:)
    type(tts_master_cloud_point_t), allocatable, intent(out) :: cloud(:)
    logical, intent(out) :: success

    integer :: cloud_index
    integer :: isotherm_index
    integer :: point_index
    integer :: point_total
    logical :: shift_found
    real(dp) :: log10_a_t
    real(dp) :: log10_reduced_frequency

    success = .false.
    point_total = 0
    do isotherm_index = 1, size(family%isotherms)
      if (.not. allocated(family%isotherms(isotherm_index)%points)) return
      point_total = point_total + &
        size(family%isotherms(isotherm_index)%points)
    end do
    allocate(cloud(point_total))
    cloud_index = 0
    do isotherm_index = 1, size(family%isotherms)
      call find_log10_shift(empirical_shifts, isotherm_index, &
        log10_a_t, shift_found)
      if (.not. shift_found) return
      do point_index = 1, size(family%isotherms(isotherm_index)%points)
        cloud_index = cloud_index + 1
        log10_reduced_frequency = log10(family%isotherms( &
          isotherm_index)%points(point_index)%frequency_hz) + log10_a_t
        if (log10_reduced_frequency < log10(tiny(1.0_dp)) .or. &
            log10_reduced_frequency > log10(huge(1.0_dp))) return
        cloud(cloud_index)%source_isotherm_index = isotherm_index
        cloud(cloud_index)%source_point_index = point_index
        cloud(cloud_index)%source_isotherm_identifier = &
          family%isotherms(isotherm_index)%isotherm_identifier
        cloud(cloud_index)%specimen_identifier = &
          family%isotherms(isotherm_index)%specimen_identifier
        cloud(cloud_index)%source_identifier = &
          family%isotherms(isotherm_index)%source_identifier
        cloud(cloud_index)%source_temperature_k = &
          family%isotherms(isotherm_index)%temperature_k
        cloud(cloud_index)%source_frequency_hz = family%isotherms( &
          isotherm_index)%points(point_index)%frequency_hz
        cloud(cloud_index)%log10_a_t = log10_a_t
        cloud(cloud_index)%reduced_frequency_hz = &
          10.0_dp**log10_reduced_frequency
        cloud(cloud_index)%storage_modulus_pa = family%isotherms( &
          isotherm_index)%points(point_index)%storage_modulus_pa
        cloud(cloud_index)%loss_modulus_pa = family%isotherms( &
          isotherm_index)%points(point_index)%loss_modulus_pa
        cloud(cloud_index)%storage_quality = family%isotherms( &
          isotherm_index)%points(point_index)%storage_quality
        cloud(cloud_index)%loss_quality = family%isotherms( &
          isotherm_index)%points(point_index)%loss_quality
        cloud(cloud_index)%contributes_to_validation = &
          is_storage_log_usable(family%isotherms( &
            isotherm_index)%points(point_index)) .or. &
          is_loss_log_usable(family%isotherms( &
            isotherm_index)%points(point_index))
      end do
    end do
    success = .true.
  end subroutine build_tts_master_experimental_cloud

  !> Reference isotherm ile başlayan, temperature-distance priority'li
  !! deterministic stitching uygular. Yeni curve yalnız mevcut reduced
  !! frequency aralığını gerçekten genişleten solver-valid measured noktaları
  !! ekler; iki taraflı extension mümkündür. Overlap duplicate edilmez,
  !! averaging/smoothing yoktur. Machine-equivalent duplicate priority'si
  !! reference, sonra reference'a daha yakın temperature, sonra daha uzaktır.
  pure subroutine stitch_tts_runtime_master_table( &
      family, reference_isotherm_index, empirical_shifts, cloud, &
      runtime_table, boundaries, success)
    type(tts_material_family_t), intent(in) :: family
    integer, intent(in) :: reference_isotherm_index
    type(tts_empirical_shift_t), intent(in) :: empirical_shifts(:)
    type(tts_master_cloud_point_t), intent(inout) :: cloud(:)
    type(tts_runtime_master_point_t), allocatable, intent(out) :: runtime_table(:)
    type(tts_master_boundary_diagnostic_t), allocatable, intent(out) :: boundaries(:)
    logical, intent(out) :: success

    type(tts_runtime_master_point_t), allocatable :: selected(:)
    integer, allocatable :: priority_indices(:)
    integer :: boundary_count
    integer :: cloud_index
    integer :: i
    integer :: isotherm_index
    integer :: max_point_count
    integer :: point_index
    integer :: selected_count
    logical :: shift_found
    real(dp) :: current_maximum
    real(dp) :: current_minimum
    real(dp) :: log10_a_t
    real(dp) :: previous_maximum
    real(dp) :: previous_minimum
    real(dp) :: reduced_frequency

    success = .false.
    if (reference_isotherm_index < 1 .or. &
        reference_isotherm_index > size(family%isotherms)) return
    max_point_count = size(cloud)
    allocate(selected(max_point_count))
    priority_indices = make_stitch_priority_indices( &
      family, reference_isotherm_index)
    selected_count = 0
    current_minimum = huge(1.0_dp)
    current_maximum = 0.0_dp

    do i = 1, size(priority_indices)
      isotherm_index = priority_indices(i)
      previous_minimum = current_minimum
      previous_maximum = current_maximum
      call find_log10_shift(empirical_shifts, isotherm_index, &
        log10_a_t, shift_found)
      if (.not. shift_found) return
      do point_index = 1, size(family%isotherms(isotherm_index)%points)
        if (.not. is_runtime_export_usable( &
            family%isotherms(isotherm_index)%points(point_index))) cycle
        reduced_frequency = calculate_reduced_frequency( &
          family%isotherms(isotherm_index)%points(point_index)%frequency_hz, &
          log10_a_t)
        if (.not. ieee_is_finite(reduced_frequency) .or. &
            reduced_frequency <= 0.0_dp) return
        if (selected_count > 0 .and. &
            isotherm_index /= reference_isotherm_index) then
          if (reduced_frequency > previous_minimum .and. &
              reduced_frequency < previous_maximum) cycle
          if (are_tts_values_machine_equivalent( &
              reduced_frequency, previous_minimum) .or. &
              are_tts_values_machine_equivalent( &
                reduced_frequency, previous_maximum)) cycle
        end if
        if (frequency_is_already_selected( &
            selected, selected_count, reduced_frequency)) cycle
        selected_count = selected_count + 1
        selected(selected_count)%reduced_frequency_hz = reduced_frequency
        selected(selected_count)%storage_modulus_pa = family%isotherms( &
          isotherm_index)%points(point_index)%storage_modulus_pa
        selected(selected_count)%loss_modulus_pa = family%isotherms( &
          isotherm_index)%points(point_index)%loss_modulus_pa
        selected(selected_count)%source_isotherm_index = isotherm_index
        selected(selected_count)%source_point_index = point_index
        selected(selected_count)%source_isotherm_identifier = &
          family%isotherms(isotherm_index)%isotherm_identifier
        current_minimum = min(current_minimum, reduced_frequency)
        current_maximum = max(current_maximum, reduced_frequency)
        cloud_index = find_cloud_point_index(cloud, isotherm_index, point_index)
        if (cloud_index > 0) then
          cloud(cloud_index)%contributes_to_runtime_extension = .true.
        end if
      end do
    end do
    if (selected_count < 2) return

    call sort_runtime_points(selected, selected_count)
    allocate(runtime_table(selected_count))
    runtime_table = selected(1:selected_count)
    do i = 2, selected_count
      if (runtime_table(i)%reduced_frequency_hz <= &
          runtime_table(i - 1)%reduced_frequency_hz .or. &
          are_tts_values_machine_equivalent( &
            runtime_table(i)%reduced_frequency_hz, &
            runtime_table(i - 1)%reduced_frequency_hz)) return
    end do

    boundary_count = 0
    do i = 1, selected_count - 1
      if (runtime_table(i)%source_isotherm_index /= &
          runtime_table(i + 1)%source_isotherm_index) then
        boundary_count = boundary_count + 1
      end if
    end do
    allocate(boundaries(boundary_count))
    boundary_count = 0
    do i = 1, selected_count - 1
      if (runtime_table(i)%source_isotherm_index == &
          runtime_table(i + 1)%source_isotherm_index) cycle
      boundary_count = boundary_count + 1
      boundaries(boundary_count)%left_runtime_point_index = i
      boundaries(boundary_count)%right_runtime_point_index = i + 1
      boundaries(boundary_count)%boundary_gap_decades = &
        log10(runtime_table(i + 1)%reduced_frequency_hz) - &
        log10(runtime_table(i)%reduced_frequency_hz)
      call calculate_boundary_mismatch(family, empirical_shifts, &
        runtime_table(i)%source_isotherm_index, &
        runtime_table(i + 1)%source_isotherm_index, &
        boundaries(boundary_count))
    end do
    success = .true.
  end subroutine stitch_tts_runtime_master_table

  pure subroutine find_log10_shift( &
      empirical_shifts, source_index, log10_a_t, found)
    type(tts_empirical_shift_t), intent(in) :: empirical_shifts(:)
    integer, intent(in) :: source_index
    real(dp), intent(out) :: log10_a_t
    logical, intent(out) :: found
    integer :: i

    found = .false.
    log10_a_t = 0.0_dp
    do i = 1, size(empirical_shifts)
      if (empirical_shifts(i)%source_isotherm_index == source_index) then
        log10_a_t = empirical_shifts(i)%log10_a_t
        found = ieee_is_finite(log10_a_t)
        return
      end if
    end do
  end subroutine find_log10_shift

  pure function calculate_reduced_frequency(frequency_hz, log10_a_t) &
      result(reduced_frequency_hz)
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: log10_a_t
    real(dp) :: reduced_frequency_hz
    real(dp) :: log10_reduced

    log10_reduced = log10(frequency_hz) + log10_a_t
    if (log10_reduced < log10(tiny(1.0_dp)) .or. &
        log10_reduced > log10(huge(1.0_dp))) then
      reduced_frequency_hz = huge(1.0_dp)
    else
      reduced_frequency_hz = 10.0_dp**log10_reduced
    end if
  end function calculate_reduced_frequency

  pure function make_stitch_priority_indices( &
      family, reference_index) result(indices)
    type(tts_material_family_t), intent(in) :: family
    integer, intent(in) :: reference_index
    integer, allocatable :: indices(:)
    integer :: i
    integer :: j
    integer :: key
    real(dp) :: key_distance
    real(dp) :: left_distance

    allocate(indices(size(family%isotherms)))
    indices(1) = reference_index
    j = 1
    do i = 1, size(family%isotherms)
      if (i == reference_index) cycle
      j = j + 1
      indices(j) = i
    end do
    do i = 3, size(indices)
      key = indices(i)
      key_distance = abs(family%isotherms(key)%temperature_k - &
        family%isotherms(reference_index)%temperature_k)
      j = i - 1
      do while (j >= 2)
        left_distance = abs(family%isotherms(indices(j))%temperature_k - &
          family%isotherms(reference_index)%temperature_k)
        if (left_distance < key_distance) exit
        if (are_tts_values_machine_equivalent(left_distance, key_distance) &
            .and. indices(j) < key) exit
        indices(j + 1) = indices(j)
        j = j - 1
      end do
      indices(j + 1) = key
    end do
  end function make_stitch_priority_indices

  pure function frequency_is_already_selected( &
      selected, selected_count, frequency_hz) result(already_selected)
    type(tts_runtime_master_point_t), intent(in) :: selected(:)
    integer, intent(in) :: selected_count
    real(dp), intent(in) :: frequency_hz
    logical :: already_selected
    integer :: i

    already_selected = .false.
    do i = 1, selected_count
      if (are_tts_values_machine_equivalent( &
          selected(i)%reduced_frequency_hz, frequency_hz)) then
        already_selected = .true.
        return
      end if
    end do
  end function frequency_is_already_selected

  pure function find_cloud_point_index( &
      cloud, isotherm_index, point_index) result(index_value)
    type(tts_master_cloud_point_t), intent(in) :: cloud(:)
    integer, intent(in) :: isotherm_index
    integer, intent(in) :: point_index
    integer :: index_value
    integer :: i

    index_value = 0
    do i = 1, size(cloud)
      if (cloud(i)%source_isotherm_index == isotherm_index .and. &
          cloud(i)%source_point_index == point_index) then
        index_value = i
        return
      end if
    end do
  end function find_cloud_point_index

  pure subroutine sort_runtime_points(points, point_count)
    type(tts_runtime_master_point_t), intent(inout) :: points(:)
    integer, intent(in) :: point_count
    type(tts_runtime_master_point_t) :: key
    integer :: i
    integer :: j

    do i = 2, point_count
      key = points(i)
      j = i - 1
      do while (j >= 1)
        if (points(j)%reduced_frequency_hz <= key%reduced_frequency_hz) exit
        points(j + 1) = points(j)
        j = j - 1
      end do
      points(j + 1) = key
    end do
  end subroutine sort_runtime_points

  pure subroutine calculate_boundary_mismatch( &
      family, empirical_shifts, left_source, right_source, diagnostic)
    type(tts_material_family_t), intent(in) :: family
    type(tts_empirical_shift_t), intent(in) :: empirical_shifts(:)
    integer, intent(in) :: left_source
    integer, intent(in) :: right_source
    type(tts_master_boundary_diagnostic_t), intent(inout) :: diagnostic

    real(dp) :: left_loss
    real(dp) :: left_shift
    real(dp) :: left_storage
    real(dp) :: overlap_lower
    real(dp) :: overlap_upper
    real(dp) :: right_loss
    real(dp) :: right_shift
    real(dp) :: right_storage
    real(dp) :: x_query
    logical :: found_left
    logical :: found_right
    logical :: left_evaluated
    logical :: right_evaluated

    call find_log10_shift(empirical_shifts, left_source, left_shift, found_left)
    call find_log10_shift(empirical_shifts, right_source, right_shift, &
      found_right)
    if (.not. found_left .or. .not. found_right) return
    call find_runtime_valid_overlap(family%isotherms(left_source), left_shift, &
      family%isotherms(right_source), right_shift, overlap_lower, &
      overlap_upper, diagnostic%has_overlap)
    if (.not. diagnostic%has_overlap) return
    x_query = 0.5_dp*(overlap_lower + overlap_upper)
    call evaluate_runtime_valid_isotherm(family%isotherms(left_source), &
      left_shift, x_query, left_storage, left_loss, left_evaluated)
    call evaluate_runtime_valid_isotherm(family%isotherms(right_source), &
      right_shift, x_query, right_storage, right_loss, right_evaluated)
    if (.not. left_evaluated .or. .not. right_evaluated) then
      diagnostic%has_overlap = .false.
      return
    end if
    diagnostic%storage_log10_mismatch = abs(left_storage - right_storage)
    diagnostic%loss_log10_mismatch = abs(left_loss - right_loss)
  end subroutine calculate_boundary_mismatch

  pure subroutine find_runtime_valid_overlap( &
      left_isotherm, left_shift, right_isotherm, right_shift, &
      overlap_lower, overlap_upper, found)
    use tms_tts_types, only : tts_isotherm_t
    type(tts_isotherm_t), intent(in) :: left_isotherm
    real(dp), intent(in) :: left_shift
    type(tts_isotherm_t), intent(in) :: right_isotherm
    real(dp), intent(in) :: right_shift
    real(dp), intent(out) :: overlap_lower
    real(dp), intent(out) :: overlap_upper
    logical, intent(out) :: found
    real(dp) :: candidate_lower
    real(dp) :: candidate_upper
    integer :: i
    integer :: j

    found = .false.
    overlap_lower = 0.0_dp
    overlap_upper = 0.0_dp
    do i = 1, size(left_isotherm%points) - 1
      if (.not. is_runtime_export_usable(left_isotherm%points(i)) .or. &
          .not. is_runtime_export_usable(left_isotherm%points(i + 1))) cycle
      do j = 1, size(right_isotherm%points) - 1
        if (.not. is_runtime_export_usable(right_isotherm%points(j)) .or. &
            .not. is_runtime_export_usable( &
              right_isotherm%points(j + 1))) cycle
        candidate_lower = max( &
          log10(left_isotherm%points(i)%frequency_hz) + left_shift, &
          log10(right_isotherm%points(j)%frequency_hz) + right_shift)
        candidate_upper = min( &
          log10(left_isotherm%points(i + 1)%frequency_hz) + left_shift, &
          log10(right_isotherm%points(j + 1)%frequency_hz) + right_shift)
        if (candidate_upper > candidate_lower .and. &
            .not. are_tts_values_machine_equivalent( &
              candidate_upper, candidate_lower)) then
          overlap_lower = candidate_lower
          overlap_upper = candidate_upper
          found = .true.
          return
        end if
      end do
    end do
  end subroutine find_runtime_valid_overlap

  pure subroutine evaluate_runtime_valid_isotherm( &
      isotherm, shift, reduced_log_frequency, storage_log, loss_log, found)
    use tms_tts_types, only : tts_isotherm_t
    type(tts_isotherm_t), intent(in) :: isotherm
    real(dp), intent(in) :: shift
    real(dp), intent(in) :: reduced_log_frequency
    real(dp), intent(out) :: storage_log
    real(dp), intent(out) :: loss_log
    logical, intent(out) :: found
    real(dp) :: alpha
    real(dp) :: x0
    real(dp) :: x1
    integer :: i

    found = .false.
    storage_log = 0.0_dp
    loss_log = 0.0_dp
    do i = 1, size(isotherm%points) - 1
      if (.not. is_runtime_export_usable(isotherm%points(i)) .or. &
          .not. is_runtime_export_usable(isotherm%points(i + 1))) cycle
      x0 = log10(isotherm%points(i)%frequency_hz) + shift
      x1 = log10(isotherm%points(i + 1)%frequency_hz) + shift
      if (reduced_log_frequency < x0 .or. reduced_log_frequency > x1) cycle
      alpha = (reduced_log_frequency - x0)/(x1 - x0)
      storage_log = (1.0_dp - alpha) * &
        log10(isotherm%points(i)%storage_modulus_pa) + alpha * &
        log10(isotherm%points(i + 1)%storage_modulus_pa)
      if (isotherm%points(i)%loss_modulus_pa <= 0.0_dp .or. &
          isotherm%points(i + 1)%loss_modulus_pa <= 0.0_dp) return
      loss_log = (1.0_dp - alpha) * &
        log10(isotherm%points(i)%loss_modulus_pa) + alpha * &
        log10(isotherm%points(i + 1)%loss_modulus_pa)
      found = .true.
      return
    end do
  end subroutine evaluate_runtime_valid_isotherm

end module tms_tts_master_curve
