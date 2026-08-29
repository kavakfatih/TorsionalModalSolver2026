module tms_sample_statistics
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: SAMPLE_STATISTICS_SUCCESS = 0
  integer, parameter, public :: SAMPLE_STATISTICS_EMPTY = 1
  integer, parameter, public :: SAMPLE_STATISTICS_NONFINITE_INPUT = 2
  integer, parameter, public :: SAMPLE_STATISTICS_NONFINITE_RESULT = 3

  !> Normal dağılımla tutarlı robust scale için kullanılan sabittir.
  !! scaled_MAD=MAD/0.6744897501960817 bağıntısındaki paydadır. Bu yalnız
  !! matematiksel bir ölçek dönüşümüdür; standard deviation değildir.
  real(dp), parameter, public :: NORMAL_MAD_CONSISTENCY_DENOMINATOR = &
    0.6744897501960817_dp

  !> Sonlu scalar örneklemin descriptive statistics sonucudur. Bütün değerler
  !! girdiyle aynı fiziksel birimdedir. n=1 için location statistics vardır;
  !! sample SD ve SE yoktur. Sayısal alanların ancak ilgili availability flag
  !! true olduğunda anlamlı olduğu açıkça taşınır. Coefficient of variation
  !! bilinçli olarak veri modeline dahil edilmez.
  type, public :: sample_statistics_t
    integer :: status = SAMPLE_STATISTICS_EMPTY
    integer :: sample_count = 0
    logical :: location_statistics_available = .false.
    logical :: spread_statistics_available = .false.
    real(dp) :: mean = 0.0_dp
    real(dp) :: sample_standard_deviation = 0.0_dp
    real(dp) :: standard_error = 0.0_dp
    real(dp) :: median = 0.0_dp
    real(dp) :: median_absolute_deviation = 0.0_dp
    real(dp) :: scaled_median_absolute_deviation = 0.0_dp
    real(dp) :: minimum = 0.0_dp
    real(dp) :: maximum = 0.0_dp
    real(dp) :: absolute_mean_median_difference = 0.0_dp
  end type sample_statistics_t

  public :: calculate_sample_statistics

contains

  !> Sonlu x_i örneklerinden arithmetic mean, n-1 paydalı sample SD,
  !! SE=SD/sqrt(n), median, MAD ve normal-consistent scaled MAD hesaplar.
  !! Girdiler herhangi bir scalar fiziksel birimde olabilir; bütün çıktılar
  !! aynı birimdedir. Welford güncellemesi cancellation riskini azaltır.
  !! n=0 unavailable, n=1 spread unavailable olur; sıfır spread yalnız n>=2
  !! sabit örneklemde gerçek descriptive sonuçtur. Ayrıntılı model için
  !! docs/mathematics/tts_repeatability_bootstrap.md belgesine bakınız.
  pure function calculate_sample_statistics(values) result(statistics)
    real(dp), intent(in) :: values(:)
    type(sample_statistics_t) :: statistics

    real(dp), allocatable :: absolute_deviations(:)
    real(dp), allocatable :: sorted_values(:)
    real(dp) :: delta
    real(dp) :: delta_after_update
    real(dp) :: m2
    integer :: i

    statistics%sample_count = size(values)
    if (size(values) == 0) then
      statistics%status = SAMPLE_STATISTICS_EMPTY
      return
    end if
    if (.not. all(ieee_is_finite(values))) then
      statistics%status = SAMPLE_STATISTICS_NONFINITE_INPUT
      return
    end if

    statistics%mean = 0.0_dp
    m2 = 0.0_dp
    do i = 1, size(values)
      delta = values(i) - statistics%mean
      if (.not. ieee_is_finite(delta)) then
        statistics%status = SAMPLE_STATISTICS_NONFINITE_RESULT
        return
      end if
      statistics%mean = statistics%mean + delta/real(i, dp)
      delta_after_update = values(i) - statistics%mean
      m2 = m2 + delta*delta_after_update
      if (.not. ieee_is_finite(statistics%mean) .or. &
          .not. ieee_is_finite(m2)) then
        statistics%status = SAMPLE_STATISTICS_NONFINITE_RESULT
        return
      end if
    end do

    sorted_values = values
    call sort_ascending(sorted_values)
    statistics%median = median_of_sorted_values(sorted_values)
    statistics%minimum = sorted_values(1)
    statistics%maximum = sorted_values(size(sorted_values))
    allocate(absolute_deviations(size(values)))
    absolute_deviations = abs(values - statistics%median)
    call sort_ascending(absolute_deviations)
    statistics%median_absolute_deviation = &
      median_of_sorted_values(absolute_deviations)
    statistics%scaled_median_absolute_deviation = &
      statistics%median_absolute_deviation / &
      NORMAL_MAD_CONSISTENCY_DENOMINATOR
    statistics%absolute_mean_median_difference = &
      abs(statistics%mean - statistics%median)

    if (.not. ieee_is_finite(statistics%median) .or. &
        .not. ieee_is_finite( &
          statistics%scaled_median_absolute_deviation) .or. &
        .not. ieee_is_finite( &
          statistics%absolute_mean_median_difference)) then
      statistics%status = SAMPLE_STATISTICS_NONFINITE_RESULT
      return
    end if

    statistics%location_statistics_available = .true.
    if (size(values) >= 2) then
      statistics%sample_standard_deviation = &
        sqrt(max(0.0_dp, m2/real(size(values) - 1, dp)))
      statistics%standard_error = statistics%sample_standard_deviation / &
        sqrt(real(size(values), dp))
      if (.not. ieee_is_finite(statistics%sample_standard_deviation) .or. &
          .not. ieee_is_finite(statistics%standard_error)) then
        statistics%status = SAMPLE_STATISTICS_NONFINITE_RESULT
        statistics%location_statistics_available = .false.
        return
      end if
      statistics%spread_statistics_available = .true.
    end if

    statistics%status = SAMPLE_STATISTICS_SUCCESS
  end function calculate_sample_statistics

  !> Sonlu bir diziyi artan sıraya koyan deterministic insertion sort'tur.
  !! Matematiksel değerleri değiştirmez; sample quantile/median hazırlığıdır.
  pure subroutine sort_ascending(values)
    real(dp), intent(inout) :: values(:)

    real(dp) :: key
    integer :: i
    integer :: j

    do i = 2, size(values)
      key = values(i)
      j = i - 1
      do while (j >= 1)
        if (values(j) <= key) exit
        values(j + 1) = values(j)
        j = j - 1
      end do
      values(j + 1) = key
    end do
  end subroutine sort_ascending

  !> Artan sıralı sonlu örneklemin median değerini döndürür. Tek n için orta
  !! sıra istatistiği, çift n için iki orta değerin arithmetic mean'i alınır.
  pure function median_of_sorted_values(values) result(median)
    real(dp), intent(in) :: values(:)
    real(dp) :: median
    integer :: upper_middle

    upper_middle = size(values)/2 + 1
    if (mod(size(values), 2) == 1) then
      median = values(upper_middle)
    else
      median = values(upper_middle - 1) + &
        0.5_dp*(values(upper_middle) - values(upper_middle - 1))
    end if
  end function median_of_sorted_values

end module tms_sample_statistics
