module tms_tts_covariance_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_isotherm_t, &
    tts_validation_result_t, TTS_IDENTIFICATION_SUCCESS, &
    is_storage_log_usable, is_loss_log_usable, &
    are_tts_values_machine_equivalent, validate_tts_material_family
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, &
    tts_dynamic_modulus_uncertainty_point_t, &
    tts_uncertainty_validation_result_t
  use tms_tts_uncertainty_validation, only : &
    validate_tts_uncertainty_family, find_tts_uncertainty_point
  use tms_tts_covariance_types, only : &
    tts_covariance_matrix_2x2_t, tts_dynamic_modulus_covariance_point_t, &
    tts_dynamic_modulus_covariance_family_t, &
    tts_covariance_matrix_validation_t, &
    tts_covariance_validation_result_t, &
    tts_bivariate_covariance_log_segment_t, &
    COVARIANCE_SOURCE_UNSPECIFIED, COVARIANCE_SOURCE_DIRECT, &
    COVARIANCE_SOURCE_REPEAT_MEASUREMENT, &
    COVARIANCE_SOURCE_MAGNITUDE_PHASE_PROPAGATION, &
    COVARIANCE_SOURCE_CALIBRATION_MODEL, TTS_COVARIANCE_SUCCESS, &
    TTS_COVARIANCE_INVALID_INPUT, TTS_COVARIANCE_DATA_MISMATCH, &
    TTS_COVARIANCE_INVALID_MATRIX, TTS_COVARIANCE_SINGULAR_MATRIX, &
    TTS_COVARIANCE_ILL_CONDITIONED, &
    TTS_COVARIANCE_NO_BIVARIATE_SUPPORT
  use tms_tts_covariance_propagation, only : &
    validate_covariance_matrix_2x2, propagate_log10_modulus_covariance
  use tms_tts_covariance_types, only : &
    tts_log_covariance_propagation_result_t
  implicit none
  private

  public :: validate_tts_covariance_family
  public :: find_tts_covariance_point
  public :: build_tts_bivariate_covariance_log_segments
  public :: interpolate_covariance_matrix_entries

contains

  !> Covariance overlay'i authoritative measurement ve V0.8.4 uncertainty
  !! overlay'lerine unique physical (T,f) anahtarıyla bağlar. Available tam
  !! matris SPD olmalı ve diyagonalleri mevcut u_G^2 ile machine-equivalent
  !! olmalıdır; iki çelişkili diyagonal kaynak sessizce seçilmez.
  pure function validate_tts_covariance_family( &
      measurement_family, uncertainty_family, covariance_family) &
      result(validation)
    type(tts_material_family_t), intent(in) :: measurement_family
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    type(tts_covariance_validation_result_t) :: validation

    type(tts_validation_result_t) :: measurement_validation
    type(tts_uncertainty_validation_result_t) :: uncertainty_validation
    type(tts_covariance_matrix_validation_t) :: matrix_validation
    type(tts_dynamic_modulus_uncertainty_point_t) :: uncertainty_point
    integer :: i
    integer :: j
    integer :: measurement_isotherm_index
    integer :: measurement_isotherm_matches
    integer :: measurement_point_index
    integer :: measurement_point_matches
    logical :: uncertainty_found
    logical :: uncertainty_unique

    measurement_validation = validate_tts_material_family(measurement_family)
    if (measurement_validation%status /= TTS_IDENTIFICATION_SUCCESS) then
      validation%message = "Authoritative measurement family geçersiz."
      return
    end if
    uncertainty_validation = validate_tts_uncertainty_family( &
      measurement_family, uncertainty_family)
    if (.not. uncertainty_validation%valid) then
      validation%message = "V0.8.4 uncertainty overlay geçersiz."
      return
    end if
    if (len_trim(covariance_family%family_identifier) == 0 .or. &
        len_trim(covariance_family%provenance) == 0 .or. &
        .not. allocated(covariance_family%isotherms)) then
      validation%message = "Covariance family kimliği/provenance/dizisi eksik."
      return
    end if

    do i = 1, size(covariance_family%isotherms)
      if (.not. ieee_is_finite( &
          covariance_family%isotherms(i)%temperature_k) .or. &
          covariance_family%isotherms(i)%temperature_k <= 0.0_dp) then
        validation%message = "Covariance isotherm sıcaklığı pozitif K olmalı."
        return
      end if
      do j = 1, i - 1
        if (are_tts_values_machine_equivalent( &
            covariance_family%isotherms(i)%temperature_k, &
            covariance_family%isotherms(j)%temperature_k)) then
          validation%status = TTS_COVARIANCE_DATA_MISMATCH
          validation%message = "Covariance isotherm sıcaklıkları unique olmalı."
          return
        end if
      end do
      call find_measurement_isotherm(measurement_family, &
        covariance_family%isotherms(i)%temperature_k, &
        measurement_isotherm_index, measurement_isotherm_matches)
      if (measurement_isotherm_matches /= 1) then
        validation%status = TTS_COVARIANCE_DATA_MISMATCH
        validation%message = "Covariance sıcaklığı unique measurement bulmadı."
        return
      end if
      if (.not. allocated(covariance_family%isotherms(i)%points)) cycle

      do j = 1, size(covariance_family%isotherms(i)%points)
        if (.not. covariance_key_is_valid( &
            covariance_family%isotherms(i)%points(j))) then
          validation%message = "Covariance (T,f) anahtarı sonlu/pozitif olmalı."
          return
        end if
        if (.not. are_tts_values_machine_equivalent( &
            covariance_family%isotherms(i)%points(j)%temperature_k, &
            covariance_family%isotherms(i)%temperature_k)) then
          validation%status = TTS_COVARIANCE_DATA_MISMATCH
          validation%message = "Covariance point/isotherm sıcaklıkları farklı."
          return
        end if
        call find_measurement_point( &
          measurement_family%isotherms(measurement_isotherm_index), &
          covariance_family%isotherms(i)%points(j)%frequency_hz, &
          measurement_point_index, measurement_point_matches)
        if (measurement_point_matches /= 1) then
          validation%status = TTS_COVARIANCE_DATA_MISMATCH
          validation%message = "Covariance frekansı unique measurement bulmadı."
          return
        end if
        if (count_covariance_key_matches(covariance_family, &
            covariance_family%isotherms(i)%points(j)%temperature_k, &
            covariance_family%isotherms(i)%points(j)%frequency_hz) /= 1) then
          validation%status = TTS_COVARIANCE_DATA_MISMATCH
          validation%message = "Duplicate covariance physical key bulundu."
          return
        end if
        if (.not. covariance_source_is_known( &
            covariance_family%isotherms(i)%points(j) &
              %covariance_source_kind)) then
          validation%message = "Tanımsız covariance source enum değeri."
          return
        end if
        validation%matched_point_count = validation%matched_point_count + 1
        if (.not. covariance_family%isotherms(i)%points(j) &
            %covariance_available) then
          validation%covariance_gap_count = &
            validation%covariance_gap_count + 1
          cycle
        end if

        matrix_validation = validate_covariance_matrix_2x2( &
          covariance_family%isotherms(i)%points(j)%covariance)
        if (.not. matrix_validation%covariance_valid) then
          validation%status = matrix_validation%status
          validation%message = "Available covariance matrisi SPD değildir."
          return
        end if
        if (.not. matrix_validation &
            %covariance_numerically_well_conditioned) then
          validation%status = TTS_COVARIANCE_ILL_CONDITIONED
          validation%ill_conditioned_covariance_count = &
            validation%ill_conditioned_covariance_count + 1
          validation%message = &
            "Available covariance matrisi numerical olarak near-singular."
          return
        end if
        call find_tts_uncertainty_point(uncertainty_family, &
          covariance_family%isotherms(i)%points(j)%temperature_k, &
          covariance_family%isotherms(i)%points(j)%frequency_hz, &
          uncertainty_point, uncertainty_found, uncertainty_unique)
        if (.not. uncertainty_found .or. .not. uncertainty_unique .or. &
            .not. uncertainty_point%storage_uncertainty_available .or. &
            .not. uncertainty_point%loss_uncertainty_available) then
          validation%covariance_gap_count = &
            validation%covariance_gap_count + 1
          cycle
        end if
        if (.not. covariance_diagonal_is_consistent( &
            covariance_family%isotherms(i)%points(j)%covariance, &
            uncertainty_point)) then
          validation%status = TTS_COVARIANCE_DATA_MISMATCH
          validation%message = &
            "Covariance diyagonali V0.8.4 standard uncertainty ile uyuşmuyor."
          return
        end if
        validation%covariance_available_point_count = &
          validation%covariance_available_point_count + 1
      end do
    end do

    if (validation%covariance_available_point_count == 0) then
      validation%status = TTS_COVARIANCE_NO_BIVARIATE_SUPPORT
      validation%message = "Bivariate covariance support noktası bulunamadı."
      return
    end if
    validation%status = TTS_COVARIANCE_SUCCESS
    validation%valid = .true.
    validation%message = &
      "Point-local physical covariance overlay ve diyagonaller doğrulandı."
  end function validate_tts_covariance_family

  !> Physical (T [K], f [Hz]) anahtarındaki covariance kaydını dizi sırasından
  !! bağımsız bulur. `unique` duplicate kayıtların sessiz seçilmesini önler.
  pure subroutine find_tts_covariance_point( &
      covariance_family, temperature_k, frequency_hz, point, found, unique)
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: frequency_hz
    type(tts_dynamic_modulus_covariance_point_t), intent(out) :: point
    logical, intent(out) :: found
    logical, intent(out) :: unique

    integer :: i
    integer :: j
    integer :: match_count

    match_count = 0
    if (allocated(covariance_family%isotherms)) then
      do i = 1, size(covariance_family%isotherms)
        if (.not. allocated(covariance_family%isotherms(i)%points)) cycle
        do j = 1, size(covariance_family%isotherms(i)%points)
          if (are_tts_values_machine_equivalent( &
              covariance_family%isotherms(i)%points(j)%temperature_k, &
              temperature_k) .and. are_tts_values_machine_equivalent( &
              covariance_family%isotherms(i)%points(j)%frequency_hz, &
              frequency_hz)) then
            point = covariance_family%isotherms(i)%points(j)
            match_count = match_count + 1
          end if
        end do
      end do
    end if
    found = match_count > 0
    unique = match_count == 1
  end subroutine find_tts_covariance_point

  !> Measurement, iki V0.8.4 uncertainty channel'ı ve SPD covariance
  !! kesişiminden contiguous O_B segmentleri kurar. Missing covariance run'ı
  !! böler; gap üzerinden interpolation/extrapolation yapılmaz. Matris
  !! elemanları x=log10(f) koordinatında endpoint-preserving lineerdir.
  pure function build_tts_bivariate_covariance_log_segments( &
      isotherm, uncertainty_family, covariance_family) result(segments)
    type(tts_isotherm_t), intent(in) :: isotherm
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    type(tts_bivariate_covariance_log_segment_t), allocatable :: segments(:)

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
      if (.not. bivariate_point_is_usable(isotherm, point_index, &
          uncertainty_family, covariance_family)) then
        point_index = point_index + 1
        cycle
      end if
      run_start = point_index
      run_end = point_index
      do while (run_end < size(isotherm%points))
        if (.not. bivariate_point_is_usable(isotherm, run_end + 1, &
            uncertainty_family, covariance_family)) exit
        run_end = run_end + 1
      end do
      if (run_end > run_start) segment_count = segment_count + 1
      point_index = run_end + 1
    end do

    allocate(segments(segment_count))
    point_index = 1
    segment_index = 0
    do while (point_index <= size(isotherm%points))
      if (.not. bivariate_point_is_usable(isotherm, point_index, &
          uncertainty_family, covariance_family)) then
        point_index = point_index + 1
        cycle
      end if
      run_start = point_index
      run_end = point_index
      do while (run_end < size(isotherm%points))
        if (.not. bivariate_point_is_usable(isotherm, run_end + 1, &
            uncertainty_family, covariance_family)) exit
        run_end = run_end + 1
      end do
      if (run_end > run_start) then
        segment_index = segment_index + 1
        call populate_bivariate_segment(segments(segment_index), isotherm, &
          covariance_family, run_start, run_end)
      end if
      point_index = run_end + 1
    end do
  end function build_tts_bivariate_covariance_log_segments

  !> İki SPD endpoint covariance matrisinin convex entrywise interpolation'ını
  !! Sigma(alpha)=(1-alpha)Sigma0+alpha Sigma1 ile hesaplar. alpha [0,1]
  !! boyutsuzdur; bu deterministic TMS26 politikası fizik yasası değildir.
  pure subroutine interpolate_covariance_matrix_entries( &
      matrix0, matrix1, alpha, matrix, valid)
    type(tts_covariance_matrix_2x2_t), intent(in) :: matrix0
    type(tts_covariance_matrix_2x2_t), intent(in) :: matrix1
    real(dp), intent(in) :: alpha
    type(tts_covariance_matrix_2x2_t), intent(out) :: matrix
    logical, intent(out) :: valid
    type(tts_covariance_matrix_validation_t) :: validation

    valid = .false.
    if (.not. ieee_is_finite(alpha) .or. alpha < 0.0_dp .or. &
        alpha > 1.0_dp) return
    matrix%storage_variance_pa2 = (1.0_dp - alpha) * &
      matrix0%storage_variance_pa2 + alpha*matrix1%storage_variance_pa2
    matrix%loss_variance_pa2 = (1.0_dp - alpha) * &
      matrix0%loss_variance_pa2 + alpha*matrix1%loss_variance_pa2
    matrix%storage_loss_covariance_pa2 = (1.0_dp - alpha) * &
      matrix0%storage_loss_covariance_pa2 + &
      alpha*matrix1%storage_loss_covariance_pa2
    validation = validate_covariance_matrix_2x2(matrix)
    valid = validation%covariance_numerically_well_conditioned
  end subroutine interpolate_covariance_matrix_entries

  pure function bivariate_point_is_usable( &
      isotherm, point_index, uncertainty_family, covariance_family) &
      result(usable)
    type(tts_isotherm_t), intent(in) :: isotherm
    integer, intent(in) :: point_index
    type(tts_dynamic_modulus_uncertainty_family_t), intent(in) :: &
      uncertainty_family
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    logical :: usable

    type(tts_dynamic_modulus_uncertainty_point_t) :: uncertainty_point
    type(tts_dynamic_modulus_covariance_point_t) :: covariance_point
    type(tts_covariance_matrix_validation_t) :: matrix_validation
    type(tts_log_covariance_propagation_result_t) :: propagated
    logical :: found
    logical :: unique

    usable = .false.
    if (.not. is_storage_log_usable(isotherm%points(point_index)) .or. &
        .not. is_loss_log_usable(isotherm%points(point_index))) return
    call find_tts_uncertainty_point(uncertainty_family, &
      isotherm%temperature_k, isotherm%points(point_index)%frequency_hz, &
      uncertainty_point, found, unique)
    if (.not. found .or. .not. unique .or. &
        .not. uncertainty_point%storage_uncertainty_available .or. &
        .not. uncertainty_point%loss_uncertainty_available) return
    call find_tts_covariance_point(covariance_family, &
      isotherm%temperature_k, isotherm%points(point_index)%frequency_hz, &
      covariance_point, found, unique)
    if (.not. found .or. .not. unique .or. &
        .not. covariance_point%covariance_available) return
    if (.not. covariance_diagonal_is_consistent( &
        covariance_point%covariance, uncertainty_point)) return
    matrix_validation = validate_covariance_matrix_2x2( &
      covariance_point%covariance)
    if (.not. matrix_validation%covariance_numerically_well_conditioned) &
      return
    propagated = propagate_log10_modulus_covariance( &
      isotherm%points(point_index)%storage_modulus_pa, &
      isotherm%points(point_index)%loss_modulus_pa, &
      covariance_point%covariance)
    usable = propagated%valid
  end function bivariate_point_is_usable

  pure subroutine populate_bivariate_segment( &
      segment, isotherm, covariance_family, &
      run_start, run_end)
    type(tts_bivariate_covariance_log_segment_t), intent(out) :: segment
    type(tts_isotherm_t), intent(in) :: isotherm
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: &
      covariance_family
    integer, intent(in) :: run_start
    integer, intent(in) :: run_end

    type(tts_dynamic_modulus_covariance_point_t) :: covariance_point
    type(tts_log_covariance_propagation_result_t) :: propagated
    logical :: found
    logical :: unique
    integer :: local_index
    integer :: source_index
    integer :: n

    n = run_end - run_start + 1
    allocate(segment%x(n), segment%storage_y(n), segment%loss_y(n))
    allocate(segment%storage_variance(n), segment%loss_variance(n))
    allocate(segment%storage_loss_covariance(n))
    allocate(segment%source_point_indices(n))
    allocate(segment%covariance_point_indices(n))
    do source_index = run_start, run_end
      local_index = source_index - run_start + 1
      call find_tts_covariance_point(covariance_family, &
        isotherm%temperature_k, isotherm%points(source_index)%frequency_hz, &
        covariance_point, found, unique)
      if (.not. found .or. .not. unique) return
      propagated = propagate_log10_modulus_covariance( &
        isotherm%points(source_index)%storage_modulus_pa, &
        isotherm%points(source_index)%loss_modulus_pa, &
        covariance_point%covariance)
      segment%x(local_index) = &
        log10(isotherm%points(source_index)%frequency_hz)
      segment%storage_y(local_index) = &
        log10(isotherm%points(source_index)%storage_modulus_pa)
      segment%loss_y(local_index) = &
        log10(isotherm%points(source_index)%loss_modulus_pa)
      segment%storage_variance(local_index) = &
        propagated%storage_variance
      segment%loss_variance(local_index) = propagated%loss_variance
      segment%storage_loss_covariance(local_index) = &
        propagated%storage_loss_covariance
      segment%source_point_indices(local_index) = source_index
      segment%covariance_point_indices(local_index) = &
        find_covariance_source_index(covariance_family, &
          isotherm%temperature_k, &
          isotherm%points(source_index)%frequency_hz)
    end do
  end subroutine populate_bivariate_segment

  pure function covariance_diagonal_is_consistent( &
      covariance, uncertainty_point) result(consistent)
    type(tts_covariance_matrix_2x2_t), intent(in) :: covariance
    type(tts_dynamic_modulus_uncertainty_point_t), intent(in) :: &
      uncertainty_point
    logical :: consistent
    real(dp) :: storage_variance
    real(dp) :: loss_variance

    consistent = .false.
    if (.not. uncertainty_point%storage_uncertainty_available .or. &
        .not. uncertainty_point%loss_uncertainty_available) return
    storage_variance = &
      uncertainty_point%storage_standard_uncertainty_pa**2
    loss_variance = uncertainty_point%loss_standard_uncertainty_pa**2
    if (.not. ieee_is_finite(storage_variance) .or. &
        .not. ieee_is_finite(loss_variance)) return
    consistent = are_tts_values_machine_equivalent( &
      covariance%storage_variance_pa2, storage_variance) .and. &
      are_tts_values_machine_equivalent( &
        covariance%loss_variance_pa2, loss_variance)
  end function covariance_diagonal_is_consistent

  pure function covariance_key_is_valid(point) result(valid)
    type(tts_dynamic_modulus_covariance_point_t), intent(in) :: point
    logical :: valid

    valid = ieee_is_finite(point%temperature_k) .and. &
      ieee_is_finite(point%frequency_hz) .and. &
      point%temperature_k > 0.0_dp .and. point%frequency_hz > 0.0_dp
  end function covariance_key_is_valid

  pure function covariance_source_is_known(source) result(known)
    integer, intent(in) :: source
    logical :: known

    known = source == COVARIANCE_SOURCE_UNSPECIFIED .or. &
      source == COVARIANCE_SOURCE_DIRECT .or. &
      source == COVARIANCE_SOURCE_REPEAT_MEASUREMENT .or. &
      source == COVARIANCE_SOURCE_MAGNITUDE_PHASE_PROPAGATION .or. &
      source == COVARIANCE_SOURCE_CALIBRATION_MODEL
  end function covariance_source_is_known

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

  pure function count_covariance_key_matches( &
      family, temperature_k, frequency_hz) result(match_count)
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: family
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
  end function count_covariance_key_matches

  pure function find_covariance_source_index( &
      family, temperature_k, frequency_hz) result(source_index)
    type(tts_dynamic_modulus_covariance_family_t), intent(in) :: family
    real(dp), intent(in) :: temperature_k
    real(dp), intent(in) :: frequency_hz
    integer :: source_index
    integer :: flat_index
    integer :: i
    integer :: j

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
  end function find_covariance_source_index

end module tms_tts_covariance_validation
