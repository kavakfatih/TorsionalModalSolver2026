module tms_tts_uncertainty_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_isotherm_t, &
    tts_validation_result_t, &
    TTS_CHANNEL_STORAGE, TTS_CHANNEL_LOSS, &
    TTS_IDENTIFICATION_SUCCESS, is_storage_log_usable, &
    is_loss_log_usable, are_tts_values_machine_equivalent, &
    validate_tts_material_family
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_point_t, &
    tts_dynamic_modulus_uncertainty_family_t, &
    tts_log_uncertainty_result_t, tts_uncertainty_validation_result_t, &
    tts_uncertainty_log_segment_t, TTS_UNCERTAINTY_SUCCESS, &
    TTS_UNCERTAINTY_INVALID_INPUT, TTS_UNCERTAINTY_DATA_MISMATCH, &
    TTS_UNCERTAINTY_NONFINITE_DATA, TTS_UNCERTAINTY_NO_SUPPORT, &
    UNCERTAINTY_SOURCE_UNSPECIFIED, UNCERTAINTY_SOURCE_TYPE_A, &
    UNCERTAINTY_SOURCE_TYPE_B, UNCERTAINTY_SOURCE_COMBINED_STANDARD, &
    UNCERTAINTY_SOURCE_REPEAT_MEASUREMENT
  implicit none
  private

  public :: propagate_log10_standard_uncertainty
  public :: validate_tts_uncertainty_family
  public :: find_tts_uncertainty_point
  public :: build_tts_uncertainty_log_segments

contains

  !> y=log10(G) dönüşümünde first-order GUM sensitivity coefficient ile
  !! u_y=u_G/(G ln(10)) hesaplar. G ve u_G [Pa], u_y boyutsuzdur. G>0,
  !! u_G>0 ve sonluluk gerekir; zero uncertainty için epsilon eklenmez.
  pure function propagate_log10_standard_uncertainty( &
      modulus_pa, standard_uncertainty_pa) result(result)
    real(dp), intent(in) :: modulus_pa
    real(dp), intent(in) :: standard_uncertainty_pa
    type(tts_log_uncertainty_result_t) :: result

    if (.not. ieee_is_finite(modulus_pa) .or. &
        .not. ieee_is_finite(standard_uncertainty_pa)) then
      result%status = TTS_UNCERTAINTY_NONFINITE_DATA
      return
    end if
    if (modulus_pa <= 0.0_dp .or. standard_uncertainty_pa <= 0.0_dp) return
    result%standard_uncertainty = &
      (standard_uncertainty_pa/modulus_pa)/log(10.0_dp)
    result%variance = result%standard_uncertainty**2
    if (.not. ieee_is_finite(result%standard_uncertainty) .or. &
        .not. ieee_is_finite(result%variance) .or. &
        result%standard_uncertainty <= 0.0_dp .or. &
        result%variance <= 0.0_dp) then
      result%status = TTS_UNCERTAINTY_NONFINITE_DATA
      result%standard_uncertainty = 0.0_dp
      result%variance = 0.0_dp
      return
    end if
    result%status = TTS_UNCERTAINTY_SUCCESS
    result%valid = .true.
  end function propagate_log10_standard_uncertainty

  !> Additive uncertainty overlay'i authoritative measurement family'ye
  !! unique physical (T,f) anahtarıyla bağlar. Array index eşlemesi yapılmaz.
  !! Missing measurement uncertainty geçerlidir ve support gap oluşturur;
  !! explicit available fakat nonpositive/nonfinite u_G reddedilir.
  pure function validate_tts_uncertainty_family( &
      measurement_family, uncertainty_family) result(validation)
    type(tts_material_family_t), intent(in) :: measurement_family
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_uncertainty_validation_result_t) :: validation

    type(tts_validation_result_t) :: measurement_validation
    integer :: i
    integer :: j
    integer :: source_isotherm_index
    integer :: source_match_count
    integer :: source_point_index
    integer :: point_match_count
    logical :: uncertainties_valid

    measurement_validation = validate_tts_material_family(measurement_family)
    if (measurement_validation%status /= TTS_IDENTIFICATION_SUCCESS) then
      validation%message = "Authoritative measurement family geçersiz."
      return
    end if
    if (len_trim(uncertainty_family%family_identifier) == 0 .or. &
        len_trim(uncertainty_family%provenance) == 0 .or. &
        .not. allocated(uncertainty_family%isotherms)) then
      validation%message = "Uncertainty family kimliği/provenance/dizisi eksik."
      return
    end if

    do i = 1, size(uncertainty_family%isotherms)
      if (.not. ieee_is_finite( &
          uncertainty_family%isotherms(i)%temperature_k)) then
        validation%status = TTS_UNCERTAINTY_NONFINITE_DATA
        validation%message = "Uncertainty isotherm sıcaklığı sonlu olmalıdır."
        return
      end if
      if (uncertainty_family%isotherms(i)%temperature_k <= 0.0_dp) then
        validation%message = "Uncertainty isotherm sıcaklığı pozitif K olmalı."
        return
      end if
      do j = 1, i - 1
        if (are_tts_values_machine_equivalent( &
            uncertainty_family%isotherms(i)%temperature_k, &
            uncertainty_family%isotherms(j)%temperature_k)) then
          validation%status = TTS_UNCERTAINTY_DATA_MISMATCH
          validation%message = "Uncertainty isotherm sıcaklıkları unique olmalı."
          return
        end if
      end do
      call find_measurement_isotherm(measurement_family, &
        uncertainty_family%isotherms(i)%temperature_k, &
        source_isotherm_index, source_match_count)
      if (source_match_count /= 1) then
        validation%status = TTS_UNCERTAINTY_DATA_MISMATCH
        validation%message = "Uncertainty sıcaklığı unique measurement bulmadı."
        return
      end if
      if (.not. allocated(uncertainty_family%isotherms(i)%points)) cycle

      do j = 1, size(uncertainty_family%isotherms(i)%points)
        if (.not. uncertainty_point_key_is_finite( &
            uncertainty_family%isotherms(i)%points(j))) then
          validation%status = TTS_UNCERTAINTY_NONFINITE_DATA
          validation%message = "Uncertainty (T,f) anahtarı sonlu olmalıdır."
          return
        end if
        if (uncertainty_family%isotherms(i)%points(j)%temperature_k <= &
            0.0_dp .or. &
            uncertainty_family%isotherms(i)%points(j)%frequency_hz <= &
            0.0_dp) then
          validation%message = "Uncertainty (T,f) anahtarı pozitif olmalıdır."
          return
        end if
        if (.not. are_tts_values_machine_equivalent( &
            uncertainty_family%isotherms(i)%points(j)%temperature_k, &
            uncertainty_family%isotherms(i)%temperature_k)) then
          validation%status = TTS_UNCERTAINTY_DATA_MISMATCH
          validation%message = "Point ve isotherm sıcaklık anahtarları farklı."
          return
        end if
        call find_measurement_point( &
          measurement_family%isotherms(source_isotherm_index), &
          uncertainty_family%isotherms(i)%points(j)%frequency_hz, &
          source_point_index, point_match_count)
        if (point_match_count /= 1) then
          validation%status = TTS_UNCERTAINTY_DATA_MISMATCH
          validation%message = "Uncertainty frekansı unique measurement bulmadı."
          return
        end if
        if (count_uncertainty_key_matches(uncertainty_family, &
            uncertainty_family%isotherms(i)%points(j)%temperature_k, &
            uncertainty_family%isotherms(i)%points(j)%frequency_hz) /= 1) then
          validation%status = TTS_UNCERTAINTY_DATA_MISMATCH
          validation%message = "Duplicate uncertainty physical key bulundu."
          return
        end if
        if (.not. uncertainty_source_is_known( &
            uncertainty_family%isotherms(i)%points(j)%uncertainty_source)) then
          validation%message = "Tanımsız uncertainty source enum değeri."
          return
        end if
        call validate_available_uncertainties( &
          uncertainty_family%isotherms(i)%points(j), validation, &
          uncertainties_valid)
        if (.not. uncertainties_valid) return
        validation%matched_point_count = validation%matched_point_count + 1
      end do
    end do

    if (validation%storage_available_point_count == 0 .and. &
        validation%loss_available_point_count == 0) then
      validation%status = TTS_UNCERTAINTY_NO_SUPPORT
      validation%message = "Hiçbir channel için pointwise uncertainty yok."
      return
    end if
    validation%status = TTS_UNCERTAINTY_SUCCESS
    validation%valid = .true.
    validation%message = "Pointwise standard uncertainty overlay doğrulandı."
  end function validate_tts_uncertainty_family

  !> Verilen physical (T,f) anahtarının overlay'deki unique kaydını döndürür.
  !! Machine-equivalence yalnız floating-point representation eşlemesidir;
  !! deney uncertainty tolerance'ı değildir.
  pure subroutine find_tts_uncertainty_point( &
      uncertainty_family, temperature_k, frequency_hz, point, found, unique)
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: frequency_hz
    type(tts_dynamic_modulus_uncertainty_point_t), intent(out) :: point
    logical, intent(out) :: found
    logical, intent(out) :: unique

    integer :: i
    integer :: j
    integer :: match_count

    found = .false.
    unique = .false.
    match_count = 0
    if (.not. allocated(uncertainty_family%isotherms)) return
    do i = 1, size(uncertainty_family%isotherms)
      if (.not. are_tts_values_machine_equivalent( &
          uncertainty_family%isotherms(i)%temperature_k, temperature_k)) cycle
      if (.not. allocated(uncertainty_family%isotherms(i)%points)) cycle
      do j = 1, size(uncertainty_family%isotherms(i)%points)
        if (are_tts_values_machine_equivalent( &
            uncertainty_family%isotherms(i)%points(j)%temperature_k, &
            temperature_k) .and. are_tts_values_machine_equivalent( &
            uncertainty_family%isotherms(i)%points(j)%frequency_hz, &
            frequency_hz)) then
          match_count = match_count + 1
          point = uncertainty_family%isotherms(i)%points(j)
        end if
      end do
    end do
    found = match_count > 0
    unique = match_count == 1
  end subroutine find_tts_uncertainty_point

  !> Measurement-quality ve uncertainty availability kesişiminden contiguous
  !! weighted log segmentleri kurar. Missing/invalid uncertainty run'ı böler;
  !! gap üzerinden interpolation ve extrapolation yapılmaz. Segment içinde
  !! variance u_log10_G^2, x=log10(f) üzerinde endpoint-preserving lineerdir.
  pure function build_tts_uncertainty_log_segments( &
      isotherm, uncertainty_family, channel) result(segments)
    type(tts_isotherm_t), intent(in) :: isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    integer, intent(in) :: channel
    type(tts_uncertainty_log_segment_t), allocatable :: segments(:)

    integer :: point_index
    integer :: run_end
    integer :: run_start
    integer :: segment_count
    integer :: segment_index

    if (.not. allocated(isotherm%points)) then
      allocate(segments(0))
      return
    end if
    segment_count = 0
    point_index = 1
    do while (point_index <= size(isotherm%points))
      if (.not. weighted_point_is_usable(isotherm, point_index, &
          uncertainty_family, channel)) then
        point_index = point_index + 1
        cycle
      end if
      run_start = point_index
      run_end = point_index
      do while (run_end < size(isotherm%points))
        if (.not. weighted_point_is_usable(isotherm, run_end + 1, &
            uncertainty_family, channel)) exit
        run_end = run_end + 1
      end do
      if (run_end - run_start + 1 >= 2) segment_count = segment_count + 1
      point_index = run_end + 1
    end do

    allocate(segments(segment_count))
    point_index = 1
    segment_index = 0
    do while (point_index <= size(isotherm%points))
      if (.not. weighted_point_is_usable(isotherm, point_index, &
          uncertainty_family, channel)) then
        point_index = point_index + 1
        cycle
      end if
      run_start = point_index
      run_end = point_index
      do while (run_end < size(isotherm%points))
        if (.not. weighted_point_is_usable(isotherm, run_end + 1, &
            uncertainty_family, channel)) exit
        run_end = run_end + 1
      end do
      if (run_end - run_start + 1 >= 2) then
        segment_index = segment_index + 1
        call populate_uncertainty_segment(segments(segment_index), &
          isotherm, uncertainty_family, channel, run_start, run_end)
      end if
      point_index = run_end + 1
    end do
  end function build_tts_uncertainty_log_segments

  pure function weighted_point_is_usable( &
      isotherm, point_index, uncertainty_family, channel) result(usable)
    type(tts_isotherm_t), intent(in) :: isotherm
    integer, intent(in) :: point_index
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    integer, intent(in) :: channel
    logical :: usable

    type(tts_dynamic_modulus_uncertainty_point_t) :: uncertainty_point
    type(tts_log_uncertainty_result_t) :: propagated
    logical :: found
    logical :: unique
    real(dp) :: modulus_pa
    real(dp) :: uncertainty_pa

    usable = .false.
    if (channel == TTS_CHANNEL_STORAGE) then
      if (.not. is_storage_log_usable(isotherm%points(point_index))) return
      modulus_pa = isotherm%points(point_index)%storage_modulus_pa
    else if (channel == TTS_CHANNEL_LOSS) then
      if (.not. is_loss_log_usable(isotherm%points(point_index))) return
      modulus_pa = isotherm%points(point_index)%loss_modulus_pa
    else
      return
    end if
    call find_tts_uncertainty_point(uncertainty_family, &
      isotherm%temperature_k, isotherm%points(point_index)%frequency_hz, &
      uncertainty_point, found, unique)
    if (.not. found .or. .not. unique) return
    if (channel == TTS_CHANNEL_STORAGE) then
      if (.not. uncertainty_point%storage_uncertainty_available) return
      uncertainty_pa = uncertainty_point%storage_standard_uncertainty_pa
    else
      if (.not. uncertainty_point%loss_uncertainty_available) return
      uncertainty_pa = uncertainty_point%loss_standard_uncertainty_pa
    end if
    propagated = propagate_log10_standard_uncertainty( &
      modulus_pa, uncertainty_pa)
    usable = propagated%valid
  end function weighted_point_is_usable

  pure subroutine populate_uncertainty_segment( &
      segment, isotherm, uncertainty_family, channel, run_start, run_end)
    type(tts_uncertainty_log_segment_t), intent(out) :: segment
    type(tts_isotherm_t), intent(in) :: isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    integer, intent(in) :: channel
    integer, intent(in) :: run_start
    integer, intent(in) :: run_end

    type(tts_dynamic_modulus_uncertainty_point_t) :: uncertainty_point
    type(tts_log_uncertainty_result_t) :: propagated
    logical :: found
    logical :: unique
    real(dp) :: modulus_pa
    real(dp) :: uncertainty_pa
    integer :: local_index
    integer :: source_index

    segment%channel = channel
    allocate(segment%x(run_end - run_start + 1))
    allocate(segment%y(run_end - run_start + 1))
    allocate(segment%variance(run_end - run_start + 1))
    allocate(segment%source_point_indices(run_end - run_start + 1))
    allocate(segment%uncertainty_point_indices(run_end - run_start + 1))
    do source_index = run_start, run_end
      local_index = source_index - run_start + 1
      call find_tts_uncertainty_point(uncertainty_family, &
        isotherm%temperature_k, isotherm%points(source_index)%frequency_hz, &
        uncertainty_point, found, unique)
      if (channel == TTS_CHANNEL_STORAGE) then
        modulus_pa = isotherm%points(source_index)%storage_modulus_pa
        uncertainty_pa = uncertainty_point%storage_standard_uncertainty_pa
      else
        modulus_pa = isotherm%points(source_index)%loss_modulus_pa
        uncertainty_pa = uncertainty_point%loss_standard_uncertainty_pa
      end if
      propagated = propagate_log10_standard_uncertainty( &
        modulus_pa, uncertainty_pa)
      segment%x(local_index) = &
        log10(isotherm%points(source_index)%frequency_hz)
      segment%y(local_index) = log10(modulus_pa)
      segment%variance(local_index) = propagated%variance
      segment%source_point_indices(local_index) = source_index
      segment%uncertainty_point_indices(local_index) = &
        find_uncertainty_source_index(uncertainty_family, &
          isotherm%temperature_k, &
          isotherm%points(source_index)%frequency_hz)
    end do
  end subroutine populate_uncertainty_segment

  pure subroutine find_measurement_isotherm( &
      family, temperature_k, source_index, match_count)
    type(tts_material_family_t), intent(in) :: family
    real(dp), intent(in) :: temperature_k
    integer, intent(out) :: source_index
    integer, intent(out) :: match_count
    integer :: i

    source_index = 0
    match_count = 0
    do i = 1, size(family%isotherms)
      if (are_tts_values_machine_equivalent( &
          family%isotherms(i)%temperature_k, temperature_k)) then
        source_index = i
        match_count = match_count + 1
      end if
    end do
  end subroutine find_measurement_isotherm

  pure subroutine find_measurement_point( &
      isotherm, frequency_hz, source_index, match_count)
    type(tts_isotherm_t), intent(in) :: isotherm
    real(dp), intent(in) :: frequency_hz
    integer, intent(out) :: source_index
    integer, intent(out) :: match_count
    integer :: i

    source_index = 0
    match_count = 0
    do i = 1, size(isotherm%points)
      if (are_tts_values_machine_equivalent( &
          isotherm%points(i)%frequency_hz, frequency_hz)) then
        source_index = i
        match_count = match_count + 1
      end if
    end do
  end subroutine find_measurement_point

  pure function count_uncertainty_key_matches( &
      family, temperature_k, frequency_hz) result(match_count)
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: family
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: frequency_hz
    integer :: match_count
    integer :: i
    integer :: j

    match_count = 0
    do i = 1, size(family%isotherms)
      if (.not. allocated(family%isotherms(i)%points)) cycle
      do j = 1, size(family%isotherms(i)%points)
        if (are_tts_values_machine_equivalent( &
            family%isotherms(i)%points(j)%temperature_k, temperature_k) .and. &
            are_tts_values_machine_equivalent( &
              family%isotherms(i)%points(j)%frequency_hz, frequency_hz)) &
          match_count = match_count + 1
      end do
    end do
  end function count_uncertainty_key_matches

  pure function find_uncertainty_source_index( &
      family, temperature_k, frequency_hz) result(source_index)
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: family
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: frequency_hz
    integer :: source_index
    integer :: i
    integer :: j
    integer :: flat_index

    source_index = 0
    flat_index = 0
    do i = 1, size(family%isotherms)
      if (.not. allocated(family%isotherms(i)%points)) cycle
      do j = 1, size(family%isotherms(i)%points)
        flat_index = flat_index + 1
        if (are_tts_values_machine_equivalent( &
            family%isotherms(i)%points(j)%temperature_k, temperature_k) .and. &
            are_tts_values_machine_equivalent( &
              family%isotherms(i)%points(j)%frequency_hz, frequency_hz)) then
          source_index = flat_index
          return
        end if
      end do
    end do
  end function find_uncertainty_source_index

  pure function uncertainty_point_key_is_finite(point) result(valid)
    type(tts_dynamic_modulus_uncertainty_point_t), intent(in) :: point
    logical :: valid

    valid = ieee_is_finite(point%temperature_k) .and. &
      ieee_is_finite(point%frequency_hz)
  end function uncertainty_point_key_is_finite

  pure function uncertainty_source_is_known(source) result(known)
    integer, intent(in) :: source
    logical :: known

    known = source == UNCERTAINTY_SOURCE_UNSPECIFIED .or. &
      source == UNCERTAINTY_SOURCE_TYPE_A .or. &
      source == UNCERTAINTY_SOURCE_TYPE_B .or. &
      source == UNCERTAINTY_SOURCE_COMBINED_STANDARD .or. &
      source == UNCERTAINTY_SOURCE_REPEAT_MEASUREMENT
  end function uncertainty_source_is_known

  pure subroutine validate_available_uncertainties( &
      point, validation, valid)
    type(tts_dynamic_modulus_uncertainty_point_t), intent(in) :: point
    type(tts_uncertainty_validation_result_t), intent(inout) :: validation
    logical, intent(out) :: valid

    valid = .false.
    if (point%storage_uncertainty_available) then
      if (.not. ieee_is_finite(point%storage_standard_uncertainty_pa)) then
        validation%status = TTS_UNCERTAINTY_NONFINITE_DATA
        validation%message = "Storage standard uncertainty sonlu olmalıdır."
        return
      end if
      if (point%storage_standard_uncertainty_pa <= 0.0_dp) then
        validation%status = TTS_UNCERTAINTY_INVALID_INPUT
        validation%message = "Storage standard uncertainty pozitif olmalıdır."
        return
      end if
      validation%storage_available_point_count = &
        validation%storage_available_point_count + 1
    end if
    if (point%loss_uncertainty_available) then
      if (.not. ieee_is_finite(point%loss_standard_uncertainty_pa)) then
        validation%status = TTS_UNCERTAINTY_NONFINITE_DATA
        validation%message = "Loss standard uncertainty sonlu olmalıdır."
        return
      end if
      if (point%loss_standard_uncertainty_pa <= 0.0_dp) then
        validation%status = TTS_UNCERTAINTY_INVALID_INPUT
        validation%message = "Loss standard uncertainty pozitif olmalıdır."
        return
      end if
      validation%loss_available_point_count = &
        validation%loss_available_point_count + 1
    end if
    valid = .true.
  end subroutine validate_available_uncertainties

end module tms_tts_uncertainty_validation
