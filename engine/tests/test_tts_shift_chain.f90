program test_tts_shift_chain
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_shift_chain_result_t, &
    TTS_IDENTIFICATION_SUCCESS, TTS_IDENTIFICATION_CHAIN_BROKEN, &
    TTS_IDENTIFICATION_REFERENCE_NOT_FOUND, MEASUREMENT_VALID, &
    MEASUREMENT_UNAVAILABLE
  use tms_tts_shift_chain, only : build_tts_shift_chain
  use tms_tts_test_support, only : make_exact_trs_family, &
    find_empirical_shift, assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures(5) = &
    [253.15_dp, 273.15_dp, 293.15_dp, 313.15_dp, 333.15_dp]
  real(dp), parameter :: known_shifts(5) = &
    [2.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, -2.0_dp]
  type(tts_material_family_t) :: family
  type(tts_shift_chain_result_t) :: chain
  integer :: i

  family = make_exact_trs_family(temperatures, known_shifts)
  chain = build_tts_shift_chain(family, "ISO-3")
  call assert_true(chain%status == TTS_IDENTIFICATION_SUCCESS .and. &
    chain%available, "Reference-anchored shift chain kurulamadı.")
  call assert_true(chain%reference_isotherm_index == 3 .and. &
    size(chain%pair_shift_results) == 4 .and. &
    size(chain%empirical_shifts) == 5, &
    "Shift chain result boyutları/reference index hatalı.")
  do i = 1, 5
    call assert_close(find_empirical_shift(chain%empirical_shifts, i), &
      known_shifts(i), 1.2e-6_dp, &
      "Cold/hot cumulative absolute shift hatalı.")
  end do
  call assert_close(find_empirical_shift(chain%empirical_shifts, 3), &
    0.0_dp, 0.0_dp, "Reference shift exact sıfır değil.")
  call assert_close(chain%empirical_shifts(3)%a_t, 1.0_dp, 0.0_dp, &
    "Reference a_T exact bir değil.")
  do i = 2, size(chain%empirical_shifts)
    call assert_true(chain%empirical_shifts(i)%temperature_k > &
      chain%empirical_shifts(i - 1)%temperature_k, &
      "Empirical shift table temperature-sorted değil.")
  end do

  ! Reference'tan outward ilerleyen zorunlu son link support taşımıyorsa chain
  ! kırılmalı; non-adjacent bridge veya partial runtime table yapılmamalıdır.
  do i = 1, size(family%isotherms(5)%points)
    family%isotherms(5)%points(i)%storage_quality = MEASUREMENT_UNAVAILABLE
    family%isotherms(5)%points(i)%loss_quality = MEASUREMENT_UNAVAILABLE
  end do
  chain = build_tts_shift_chain(family, "ISO-3")
  call assert_true(chain%status == TTS_IDENTIFICATION_CHAIN_BROKEN .and. &
    .not. chain%available .and. .not. allocated(chain%empirical_shifts), &
    "Broken adjacent chain complete shift table gibi raporlandı.")

  chain = build_tts_shift_chain(family, "SYNTHETIC-REFERENCE")
  call assert_true(chain%status == TTS_IDENTIFICATION_REFERENCE_NOT_FOUND &
    .and. .not. chain%available, &
    "Olmayan explicit measured reference sessizce seçildi.")

  ! Pair objective matematiksel olarak s=+310 bulabilse de derived
  ! a_T=10^s real(dp) aralığını aşar. Shift chain fake huge a_T üretmeden clean
  ! failure vermelidir; primary log10(a_T) architecture'ı değişmez.
  family = make_exact_trs_family( &
    [293.15_dp, 313.15_dp], [0.0_dp, 0.0_dp])
  call populate_extreme_shift_pair(family, 300.0_dp, -10.0_dp)
  chain = build_tts_shift_chain(family, "ISO-1")
  call assert_true(chain%status == TTS_IDENTIFICATION_CHAIN_BROKEN .and. &
    .not. chain%available, &
    "Overflow a_T shift chain içinde geçerli sayıldı.")

  ! s=-310 için 10^s normal representable domain'in altındadır. Underflow
  ! zero/sentinel a_T ile gizlenmeden aynı clean failure yolunu kullanmalıdır.
  call populate_extreme_shift_pair(family, -10.0_dp, 300.0_dp)
  chain = build_tts_shift_chain(family, "ISO-1")
  call assert_true(chain%status == TTS_IDENTIFICATION_CHAIN_BROKEN .and. &
    .not. chain%available, &
    "Underflow a_T shift chain içinde geçerli sayıldı.")

  print *, "V0.8.1 reference-anchored hot/cold shift chain doğrulandı."

contains

  !> İki üç-noktalı power-law curve'ü log-frequency ekseninde verilen
  !! merkezlere yerleştirir. G'/G'' [Pa] pointwise aynı kaldığından analytical
  !! relative shift reference_center-moving_center olur; amaç yalnız derived
  !! a_T exponentiation range kontrolünü bağımsız biçimde zorlamaktır.
  subroutine populate_extreme_shift_pair( &
      target_family, reference_center, moving_center)
    type(tts_material_family_t), intent(inout) :: target_family
    real(dp), intent(in) :: reference_center
    real(dp), intent(in) :: moving_center

    real(dp) :: centers(2)
    real(dp) :: local_coordinate
    integer :: isotherm_index
    integer :: point_index

    centers = [reference_center, moving_center]
    do isotherm_index = 1, 2
      if (allocated(target_family%isotherms(isotherm_index)%points)) then
        deallocate(target_family%isotherms(isotherm_index)%points)
      end if
      allocate(target_family%isotherms(isotherm_index)%points(3))
      do point_index = 1, 3
        local_coordinate = real(point_index - 2, dp)
        target_family%isotherms(isotherm_index)%points(point_index) &
          %frequency_hz = 10.0_dp**(centers(isotherm_index) + local_coordinate)
        target_family%isotherms(isotherm_index)%points(point_index) &
          %storage_modulus_pa = 10.0_dp**(6.0_dp + 0.1_dp*local_coordinate)
        target_family%isotherms(isotherm_index)%points(point_index) &
          %loss_modulus_pa = 10.0_dp**(5.0_dp + 0.05_dp*local_coordinate)
        target_family%isotherms(isotherm_index)%points(point_index) &
          %storage_quality = MEASUREMENT_VALID
        target_family%isotherms(isotherm_index)%points(point_index) &
          %loss_quality = MEASUREMENT_VALID
      end do
    end do
  end subroutine populate_extreme_shift_pair
end program test_tts_shift_chain
