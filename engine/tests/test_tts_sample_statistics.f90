program test_tts_sample_statistics
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tms_kinds, only : dp
  use tms_sample_statistics, only : sample_statistics_t, &
    SAMPLE_STATISTICS_SUCCESS, SAMPLE_STATISTICS_EMPTY, &
    SAMPLE_STATISTICS_NONFINITE_INPUT, &
    NORMAL_MAD_CONSISTENCY_DENOMINATOR, calculate_sample_statistics
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  type(sample_statistics_t) :: statistics
  real(dp), allocatable :: empty_values(:)
  real(dp) :: nonfinite_values(2)

  statistics = calculate_sample_statistics([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp])
  call assert_true(statistics%status == SAMPLE_STATISTICS_SUCCESS .and. &
    statistics%sample_count == 4 .and. &
    statistics%location_statistics_available .and. &
    statistics%spread_statistics_available, &
    "Dört-sample descriptive availability semantics hatalı.")
  call assert_close(statistics%mean, 2.5_dp, 1.0e-15_dp, &
    "Arithmetic sample mean hatalı.")
  call assert_close(statistics%sample_standard_deviation, &
    sqrt(5.0_dp/3.0_dp), 1.0e-15_dp, &
    "Sample SD n-1 paydasıyla hesaplanmadı.")
  call assert_close(statistics%standard_error, &
    sqrt(5.0_dp/3.0_dp)/2.0_dp, 1.0e-15_dp, &
    "Standard error SD/sqrt(n) bağıntısını sağlamıyor.")
  call assert_close(statistics%median, 2.5_dp, 1.0e-15_dp, &
    "Çift örneklem median hesabı hatalı.")
  call assert_close(statistics%median_absolute_deviation, 1.0_dp, &
    1.0e-15_dp, "MAD hesabı hatalı.")
  call assert_close(statistics%scaled_median_absolute_deviation, &
    1.0_dp/NORMAL_MAD_CONSISTENCY_DENOMINATOR, 1.0e-15_dp, &
    "Normal-consistent scaled MAD hesabı hatalı.")
  call assert_close(statistics%minimum, 1.0_dp, 0.0_dp, &
    "Sample minimum hatalı.")
  call assert_close(statistics%maximum, 4.0_dp, 0.0_dp, &
    "Sample maximum hatalı.")
  call assert_close(statistics%absolute_mean_median_difference, 0.0_dp, &
    0.0_dp, "Absolute mean-median difference hatalı.")

  statistics = calculate_sample_statistics([3.0_dp, 3.0_dp, 3.0_dp])
  call assert_true(statistics%spread_statistics_available, &
    "n>=2 sabit örneklem spread availability taşımıyor.")
  call assert_close(statistics%sample_standard_deviation, 0.0_dp, 0.0_dp, &
    "Sabit örneklemin sample spread'i sıfır değil.")

  statistics = calculate_sample_statistics([-2.0_dp, -1.0_dp, &
    1.0_dp, 2.0_dp])
  call assert_close(statistics%mean, 0.0_dp, 1.0e-15_dp, &
    "Signed near-zero sample mean hatalı.")
  call assert_close(statistics%median, 0.0_dp, 1.0e-15_dp, &
    "Signed sample median hatalı.")

  statistics = calculate_sample_statistics([5.0_dp])
  call assert_true(statistics%status == SAMPLE_STATISTICS_SUCCESS .and. &
    statistics%location_statistics_available .and. &
    .not. statistics%spread_statistics_available, &
    "n=1 sample için fake SD/SE availability üretildi.")
  call assert_close(statistics%mean, 5.0_dp, 0.0_dp, &
    "Tek sample descriptive value korunmadı.")

  allocate(empty_values(0))
  statistics = calculate_sample_statistics(empty_values)
  call assert_true(statistics%status == SAMPLE_STATISTICS_EMPTY .and. &
    .not. statistics%location_statistics_available .and. &
    .not. statistics%spread_statistics_available, &
    "Empty sample unavailable dönmedi.")

  nonfinite_values = [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)]
  statistics = calculate_sample_statistics(nonfinite_values)
  call assert_true(statistics%status == SAMPLE_STATISTICS_NONFINITE_INPUT .and. &
    .not. statistics%location_statistics_available, &
    "Nonfinite sample clean status ile reddedilmedi.")

  print *, "V0.8.3 sample mean/SD/SE/median/MAD semantics doğrulandı."
end program test_tts_sample_statistics
