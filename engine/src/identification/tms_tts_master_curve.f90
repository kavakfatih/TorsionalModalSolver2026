module tms_tts_master_curve
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_empirical_shift_t, &
    tts_master_cloud_point_t, tts_runtime_master_point_t, &
    tts_master_boundary_diagnostic_t, is_runtime_export_usable, &
    is_storage_log_usable, is_loss_log_usable, &
    are_tts_values_machine_equivalent, TTS_IDENTIFICATION_SUCCESS, &
    TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED, &
    TTS_IDENTIFICATION_RUNTIME_DOMAIN_GAP
  implicit none
  private

  !> Solver'ın interpolation yapmasına izin verilen tek bir reduced
  !! log-frequency interval'ını taşır. Alt ve üst sınırlar
  !! x_r=log10(f/Hz)+log10(a_T) boyutsuz koordinatındadır.
  type :: runtime_coverage_interval_t
    real(dp) :: lower_log10_frequency = 0.0_dp
    real(dp) :: upper_log10_frequency = 0.0_dp
  end type runtime_coverage_interval_t

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
    logical :: frequency_success
    real(dp) :: reduced_frequency_hz

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
      if (.not. shift_found) then
        deallocate(cloud)
        return
      end if
      do point_index = 1, size(family%isotherms(isotherm_index)%points)
        cloud_index = cloud_index + 1
        call calculate_reduced_frequency_checked( &
          family%isotherms(isotherm_index)%points(point_index)%frequency_hz, &
          log10_a_t, reduced_frequency_hz, frequency_success)
        if (.not. frequency_success) then
          deallocate(cloud)
          return
        end if
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
        cloud(cloud_index)%reduced_frequency_hz = reduced_frequency_hz
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
      runtime_table, boundaries, success, status)
    type(tts_material_family_t), intent(in) :: family
    integer, intent(in) :: reference_isotherm_index
    type(tts_empirical_shift_t), intent(in) :: empirical_shifts(:)
    type(tts_master_cloud_point_t), intent(inout) :: cloud(:)
    type(tts_runtime_master_point_t), allocatable, intent(out) :: runtime_table(:)
    type(tts_master_boundary_diagnostic_t), allocatable, intent(out) :: boundaries(:)
    logical, intent(out) :: success
    integer, intent(out), optional :: status

    type(runtime_coverage_interval_t), allocatable :: coverage_intervals(:)
    type(runtime_coverage_interval_t), allocatable :: merged_coverage(:)
    type(tts_runtime_master_point_t), allocatable :: selected(:)
    integer, allocatable :: priority_indices(:)
    integer :: boundary_count
    integer :: cloud_index
    integer :: i
    integer :: isotherm_index
    integer :: max_point_count
    integer :: point_index
    integer :: selected_count
    logical :: coverage_success
    logical :: frequency_success
    logical :: shift_found
    real(dp) :: current_maximum
    real(dp) :: current_minimum
    real(dp) :: log10_a_t
    real(dp) :: previous_maximum
    real(dp) :: previous_minimum
    real(dp) :: reduced_frequency

    success = .false.
    if (present(status)) status = TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED
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
        if (.not. point_belongs_to_runtime_valid_interval( &
            family%isotherms(isotherm_index), point_index)) cycle
        call calculate_reduced_frequency_checked( &
          family%isotherms(isotherm_index)%points(point_index)%frequency_hz, &
          log10_a_t, reduced_frequency, frequency_success)
        if (.not. frequency_success) return
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

    ! Runtime provider bütün table domain'ini continuous kabul eder. Bu nedenle
    ! seçilen iki nokta arasındaki her koordinat, en az bir original adjacent
    ! runtime-valid measurement interval'ı tarafından desteklenmelidir.
    call build_runtime_valid_reduced_intervals(family, empirical_shifts, &
      coverage_intervals, coverage_success)
    if (.not. coverage_success) then
      deallocate(runtime_table)
      return
    end if
    call merge_runtime_coverage_intervals(coverage_intervals, merged_coverage)
    if (.not. runtime_domain_is_fully_supported( &
        runtime_table, merged_coverage)) then
      deallocate(runtime_table)
      if (present(status)) status = TTS_IDENTIFICATION_RUNTIME_DOMAIN_GAP
      return
    end if

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
    if (present(status)) status = TTS_IDENTIFICATION_SUCCESS
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

  !> Physical frequency [Hz] ve boyutsuz s=log10(a_T) değerinden
  !! x_r=log10(f)+s, f_r=10^x_r [Hz] hesaplar. Girdi sonlu/pozitif değilse
  !! veya exponentiation Fortran real(dp) normal aralığının dışına çıkarsa
  !! success=false döner. huge/tiny sentinel veya clamp üretilmez.
  pure subroutine calculate_reduced_frequency_checked( &
      frequency_hz, log10_a_t, reduced_frequency_hz, success)
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: log10_a_t
    real(dp), intent(out) :: reduced_frequency_hz
    logical, intent(out) :: success
    real(dp) :: log10_reduced

    reduced_frequency_hz = 0.0_dp
    success = .false.
    if (.not. ieee_is_finite(frequency_hz)) return
    if (frequency_hz <= 0.0_dp) return
    if (.not. ieee_is_finite(log10_a_t)) return
    log10_reduced = log10(frequency_hz) + log10_a_t
    if (.not. ieee_is_finite(log10_reduced)) return
    if (log10_reduced < log10(tiny(1.0_dp)) .or. &
        log10_reduced > log10(huge(1.0_dp))) return
    reduced_frequency_hz = 10.0_dp**log10_reduced
    success = ieee_is_finite(reduced_frequency_hz) .and. &
      reduced_frequency_hz > 0.0_dp
    if (.not. success) reduced_frequency_hz = 0.0_dp
  end subroutine calculate_reduced_frequency_checked

  !> Bir measured point'in en az bir adjacent runtime-valid original interval'a
  !! ait olup olmadığını döndürür. İzole usable point solver interpolation
  !! domain'ini tek başına genişletemez. G''=0, quality VALID ise usable kalır.
  pure function point_belongs_to_runtime_valid_interval( &
      isotherm, point_index) result(belongs)
    use tms_tts_types, only : tts_isotherm_t
    type(tts_isotherm_t), intent(in) :: isotherm
    integer, intent(in) :: point_index
    logical :: belongs

    belongs = .false.
    if (point_index < 1) return
    if (point_index > size(isotherm%points)) return
    if (point_index > 1) then
      if (is_runtime_export_usable(isotherm%points(point_index - 1)) .and. &
          is_runtime_export_usable(isotherm%points(point_index))) then
        belongs = .true.
        return
      end if
    end if
    if (point_index < size(isotherm%points)) then
      if (is_runtime_export_usable(isotherm%points(point_index)) .and. &
          is_runtime_export_usable(isotherm%points(point_index + 1))) then
        belongs = .true.
      end if
    end if
  end function point_belongs_to_runtime_valid_interval

  !> Bütün source isotherm'lerde yalnız iki adjacent endpoint'i de
  !! runtime-export usable olan original interval'ları x_r=x+s koordinatına
  !! taşır. Frekans [Hz], s boyutsuzdur. Numerical range dışı koordinat clean
  !! failure üretir; experimental gap threshold uygulanmaz.
  pure subroutine build_runtime_valid_reduced_intervals( &
      family, empirical_shifts, intervals, success)
    type(tts_material_family_t), intent(in) :: family
    type(tts_empirical_shift_t), intent(in) :: empirical_shifts(:)
    type(runtime_coverage_interval_t), allocatable, intent(out) :: intervals(:)
    logical, intent(out) :: success

    type(runtime_coverage_interval_t), allocatable :: buffer(:)
    integer :: interval_count
    integer :: isotherm_index
    integer :: maximum_interval_count
    integer :: point_index
    logical :: shift_found
    real(dp) :: log10_a_t
    real(dp) :: lower_coordinate
    real(dp) :: upper_coordinate

    success = .false.
    maximum_interval_count = 0
    do isotherm_index = 1, size(family%isotherms)
      if (.not. allocated(family%isotherms(isotherm_index)%points)) return
      maximum_interval_count = maximum_interval_count + &
        max(0, size(family%isotherms(isotherm_index)%points) - 1)
    end do
    allocate(buffer(maximum_interval_count))
    interval_count = 0
    do isotherm_index = 1, size(family%isotherms)
      call find_log10_shift(empirical_shifts, isotherm_index, &
        log10_a_t, shift_found)
      if (.not. shift_found) return
      do point_index = 1, &
          size(family%isotherms(isotherm_index)%points) - 1
        if (.not. is_runtime_export_usable( &
            family%isotherms(isotherm_index)%points(point_index))) cycle
        if (.not. is_runtime_export_usable( &
            family%isotherms(isotherm_index)%points(point_index + 1))) cycle
        lower_coordinate = log10(family%isotherms(isotherm_index) &
          %points(point_index)%frequency_hz) + log10_a_t
        upper_coordinate = log10(family%isotherms(isotherm_index) &
          %points(point_index + 1)%frequency_hz) + log10_a_t
        if (.not. reduced_log_coordinate_is_representable(lower_coordinate)) &
          return
        if (.not. reduced_log_coordinate_is_representable(upper_coordinate)) &
          return
        if (upper_coordinate <= lower_coordinate) return
        interval_count = interval_count + 1
        buffer(interval_count)%lower_log10_frequency = lower_coordinate
        buffer(interval_count)%upper_log10_frequency = upper_coordinate
      end do
    end do
    allocate(intervals(interval_count))
    if (interval_count > 0) intervals = buffer(1:interval_count)
    success = .true.
  end subroutine build_runtime_valid_reduced_intervals

  !> Reduced log-frequency interval'larını alt sınıra göre deterministik
  !! sıralar; overlap veya yalnız machine-equivalent touching endpoint varsa
  !! birleştirir. Bu tolerance physical/experimental gap eşiği değildir.
  pure subroutine merge_runtime_coverage_intervals(intervals, merged)
    type(runtime_coverage_interval_t), intent(in) :: intervals(:)
    type(runtime_coverage_interval_t), allocatable, intent(out) :: merged(:)

    type(runtime_coverage_interval_t), allocatable :: sorted(:)
    type(runtime_coverage_interval_t), allocatable :: work(:)
    type(runtime_coverage_interval_t) :: key
    integer :: i
    integer :: j
    integer :: merged_count

    if (size(intervals) == 0) then
      allocate(merged(0))
      return
    end if
    sorted = intervals
    do i = 2, size(sorted)
      key = sorted(i)
      j = i - 1
      do while (j >= 1)
        if (sorted(j)%lower_log10_frequency <= &
            key%lower_log10_frequency) exit
        sorted(j + 1) = sorted(j)
        j = j - 1
      end do
      sorted(j + 1) = key
    end do

    allocate(work(size(sorted)))
    merged_count = 1
    work(1) = sorted(1)
    do i = 2, size(sorted)
      if (sorted(i)%lower_log10_frequency <= &
          work(merged_count)%upper_log10_frequency .or. &
          runtime_log_coordinates_are_equivalent( &
            sorted(i)%lower_log10_frequency, &
            work(merged_count)%upper_log10_frequency)) then
        work(merged_count)%upper_log10_frequency = max( &
          work(merged_count)%upper_log10_frequency, &
          sorted(i)%upper_log10_frequency)
      else
        merged_count = merged_count + 1
        work(merged_count) = sorted(i)
      end if
    end do
    allocate(merged(merged_count))
    merged = work(1:merged_count)
  end subroutine merge_runtime_coverage_intervals

  !> Strict runtime table'ın [min(log10(f_r)),max(log10(f_r))] domain'inin
  !! tek bir merged measurement-supported coverage interval'ı içinde kalıp
  !! kalmadığını doğrular. Ayrık union parçaları arasındaki hole kabul edilmez.
  pure function runtime_domain_is_fully_supported( &
      runtime_table, merged_coverage) result(supported)
    type(tts_runtime_master_point_t), intent(in) :: runtime_table(:)
    type(runtime_coverage_interval_t), intent(in) :: merged_coverage(:)
    logical :: supported
    integer :: i
    logical :: lower_supported
    logical :: upper_supported
    real(dp) :: runtime_lower
    real(dp) :: runtime_upper

    supported = .false.
    if (size(runtime_table) < 2) return
    runtime_lower = log10(runtime_table(1)%reduced_frequency_hz)
    runtime_upper = log10(runtime_table(size(runtime_table)) &
      %reduced_frequency_hz)
    do i = 1, size(merged_coverage)
      lower_supported = runtime_lower >= &
        merged_coverage(i)%lower_log10_frequency .or. &
        runtime_log_coordinates_are_equivalent(runtime_lower, &
          merged_coverage(i)%lower_log10_frequency)
      upper_supported = runtime_upper <= &
        merged_coverage(i)%upper_log10_frequency .or. &
        runtime_log_coordinates_are_equivalent(runtime_upper, &
          merged_coverage(i)%upper_log10_frequency)
      if (lower_supported .and. upper_supported) then
        supported = .true.
        return
      end if
    end do
  end function runtime_domain_is_fully_supported

  pure elemental function reduced_log_coordinate_is_representable( &
      coordinate) result(representable)
    real(dp), intent(in) :: coordinate
    logical :: representable

    representable = ieee_is_finite(coordinate) .and. &
      coordinate >= log10(tiny(1.0_dp)) .and. &
      coordinate <= log10(huge(1.0_dp))
  end function reduced_log_coordinate_is_representable

  pure elemental function runtime_log_coordinates_are_equivalent(a, b) &
      result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent
    real(dp) :: scale

    scale = max(1.0_dp, abs(a), abs(b))
    equivalent = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
      abs(a - b) <= 64.0_dp*epsilon(1.0_dp)*scale
  end function runtime_log_coordinates_are_equivalent

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
