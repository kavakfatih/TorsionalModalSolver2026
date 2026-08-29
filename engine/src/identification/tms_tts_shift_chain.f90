module tms_tts_shift_chain
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_pair_shift_configuration_t, tts_pair_shift_result_t, &
    tts_shift_chain_result_t, tts_validation_result_t, &
    TTS_IDENTIFICATION_SUCCESS, TTS_IDENTIFICATION_REFERENCE_NOT_FOUND, &
    TTS_IDENTIFICATION_CHAIN_BROKEN, validate_tts_material_family
  use tms_tts_pair_shift, only : identify_tts_pair_shift
  implicit none
  private

  public :: build_tts_shift_chain

contains

  !> Measured isotherm'leri temperature [K] sırasına göre yalnız indeks
  !! düzeyinde düzenler ve explicit reference'tan iki yönde adjacent-pair
  !! shifting uygular. Relative delta_s değerleri kümülatif olarak
  !! s_i=log10(a_T(T_i)) tablosuna çevrilir; reference s=0,a_T=1'dir. Bir
  !! zorunlu link çözülemezse non-adjacent bridge yapılmaz ve CHAIN_BROKEN
  !! döner. Original family ve point sırası değiştirilmez.
  pure function build_tts_shift_chain( &
      family, reference_isotherm_identifier, configuration) result(chain)
    type(tts_material_family_t), intent(in) :: family
    character(len=*), intent(in) :: reference_isotherm_identifier
    type(tts_pair_shift_configuration_t), intent(in), optional :: configuration
    type(tts_shift_chain_result_t) :: chain

    type(tts_pair_shift_configuration_t) :: settings
    type(tts_pair_shift_result_t) :: pair
    type(tts_validation_result_t) :: validation
    integer, allocatable :: sorted_indices(:)
    real(dp), allocatable :: absolute_shifts(:)
    integer :: current_index
    integer :: i
    integer :: next_index
    integer :: pair_slot
    integer :: reference_position
    integer :: reference_source_index
    logical :: shift_factor_success

    validation = validate_tts_material_family(family)
    if (.not. validation%valid) then
      chain%status = validation%status
      return
    end if
    if (present(configuration)) settings = configuration
    reference_source_index = find_reference_index( &
      family, reference_isotherm_identifier)
    if (reference_source_index == 0) then
      chain%status = TTS_IDENTIFICATION_REFERENCE_NOT_FOUND
      return
    end if

    sorted_indices = make_temperature_sorted_indices(family)
    reference_position = 0
    do i = 1, size(sorted_indices)
      if (sorted_indices(i) == reference_source_index) reference_position = i
    end do
    if (reference_position == 0) then
      chain%status = TTS_IDENTIFICATION_REFERENCE_NOT_FOUND
      return
    end if

    allocate(absolute_shifts(size(family%isotherms)), source=0.0_dp)
    allocate(chain%pair_shift_results(size(family%isotherms) - 1))
    chain%reference_isotherm_index = reference_source_index
    pair_slot = 0

    ! Colder direction: reference'tan azalan temperature'a adjacent ilerler.
    current_index = reference_source_index
    do i = reference_position - 1, 1, -1
      next_index = sorted_indices(i)
      pair = identify_tts_pair_shift(family%isotherms(current_index), &
        family%isotherms(next_index), settings)
      pair%reference_isotherm_index = current_index
      pair%moving_isotherm_index = next_index
      pair_slot = pair_slot + 1
      chain%pair_shift_results(pair_slot) = pair
      if (.not. pair%shift_available) then
        chain%status = TTS_IDENTIFICATION_CHAIN_BROKEN
        return
      end if
      absolute_shifts(next_index) = absolute_shifts(current_index) + pair%delta_s
      current_index = next_index
    end do

    ! Hotter direction: reference'tan artan temperature'a adjacent ilerler.
    current_index = reference_source_index
    do i = reference_position + 1, size(sorted_indices)
      next_index = sorted_indices(i)
      pair = identify_tts_pair_shift(family%isotherms(current_index), &
        family%isotherms(next_index), settings)
      pair%reference_isotherm_index = current_index
      pair%moving_isotherm_index = next_index
      pair_slot = pair_slot + 1
      chain%pair_shift_results(pair_slot) = pair
      if (.not. pair%shift_available) then
        chain%status = TTS_IDENTIFICATION_CHAIN_BROKEN
        return
      end if
      absolute_shifts(next_index) = absolute_shifts(current_index) + pair%delta_s
      current_index = next_index
    end do

    allocate(chain%empirical_shifts(size(family%isotherms)))
    do i = 1, size(sorted_indices)
      current_index = sorted_indices(i)
      chain%empirical_shifts(i)%source_isotherm_index = current_index
      chain%empirical_shifts(i)%source_isotherm_identifier = &
        family%isotherms(current_index)%isotherm_identifier
      chain%empirical_shifts(i)%temperature_k = &
        family%isotherms(current_index)%temperature_k
      if (current_index == reference_source_index) then
        chain%empirical_shifts(i)%log10_a_t = 0.0_dp
        chain%empirical_shifts(i)%a_t = 1.0_dp
      else
        chain%empirical_shifts(i)%log10_a_t = absolute_shifts(current_index)
        call calculate_shift_factor_checked(absolute_shifts(current_index), &
          chain%empirical_shifts(i)%a_t, shift_factor_success)
        if (.not. shift_factor_success) then
          chain%status = TTS_IDENTIFICATION_CHAIN_BROKEN
          return
        end if
      end if
      if (.not. ieee_is_finite(chain%empirical_shifts(i)%a_t) .or. &
          chain%empirical_shifts(i)%a_t <= 0.0_dp) then
        chain%status = TTS_IDENTIFICATION_CHAIN_BROKEN
        return
      end if
    end do

    chain%status = TTS_IDENTIFICATION_SUCCESS
    chain%available = .true.
  end function build_tts_shift_chain

  !> Primary boyutsuz s=log10(a_T) değerinden derived a_T=10^s üretir.
  !! s sonlu değilse veya exponentiation real(dp) normal aralığını aşarsa
  !! success=false döner; zero/huge sentinel ile geçerli shift taklit edilmez.
  pure subroutine calculate_shift_factor_checked(log10_a_t, a_t, success)
    real(dp), intent(in) :: log10_a_t
    real(dp), intent(out) :: a_t
    logical, intent(out) :: success

    a_t = 0.0_dp
    success = .false.
    if (.not. ieee_is_finite(log10_a_t)) return
    if (log10_a_t < log10(tiny(1.0_dp))) return
    if (log10_a_t > log10(huge(1.0_dp))) return
    a_t = 10.0_dp**log10_a_t
    success = ieee_is_finite(a_t) .and. a_t > 0.0_dp
    if (.not. success) a_t = 0.0_dp
  end subroutine calculate_shift_factor_checked

  pure function find_reference_index(family, identifier) result(index_value)
    type(tts_material_family_t), intent(in) :: family
    character(len=*), intent(in) :: identifier
    integer :: index_value
    integer :: i

    index_value = 0
    do i = 1, size(family%isotherms)
      if (trim(family%isotherms(i)%isotherm_identifier) == &
          trim(identifier)) then
        index_value = i
        return
      end if
    end do
  end function find_reference_index

  pure function make_temperature_sorted_indices(family) result(indices)
    type(tts_material_family_t), intent(in) :: family
    integer, allocatable :: indices(:)
    integer :: i
    integer :: j
    integer :: key

    allocate(indices(size(family%isotherms)))
    do i = 1, size(indices)
      indices(i) = i
    end do
    do i = 2, size(indices)
      key = indices(i)
      j = i - 1
      do while (j >= 1)
        if (family%isotherms(indices(j))%temperature_k <= &
            family%isotherms(key)%temperature_k) exit
        indices(j + 1) = indices(j)
        j = j - 1
      end do
      indices(j + 1) = key
    end do
  end function make_temperature_sorted_indices

end module tms_tts_shift_chain
