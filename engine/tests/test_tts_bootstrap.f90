program test_tts_bootstrap
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  use tms_deterministic_rng, only : deterministic_rng_t
  use tms_bootstrap, only : bootstrap_configuration_t, &
    cluster_bootstrap_plan_t, bootstrap_mean_result_t, BOOTSTRAP_SUCCESS, &
    BOOTSTRAP_INVALID_INPUT, create_cluster_bootstrap_plan, &
    calculate_cluster_bootstrap_mean, calculate_type7_quantile
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  type(deterministic_rng_t) :: rng
  type(bootstrap_configuration_t) :: configuration
  type(cluster_bootstrap_plan_t) :: plan_a
  type(cluster_bootstrap_plan_t) :: plan_b
  type(cluster_bootstrap_plan_t) :: plan_c
  type(bootstrap_mean_result_t) :: result_a
  type(bootstrap_mean_result_t) :: result_b
  type(bootstrap_mean_result_t) :: partial_result
  integer(int64) :: rng_value
  real(dp) :: quantile
  integer :: draw_index
  integer :: quantile_status
  logical :: has_replacement_duplicate
  logical :: usable_three(3)
  logical :: usable_partial(5)

  call rng%initialize(1_int64)
  call rng%next_int31(rng_value)
  call assert_true(rng_value == 16807_int64, &
    "Park-Miller portable sequence ilk değeri hatalı.")
  call rng%next_int31(rng_value)
  call assert_true(rng_value == 282475249_int64, &
    "Park-Miller portable sequence ikinci değeri hatalı.")
  call rng%next_int31(rng_value)
  call assert_true(rng_value == 1622650073_int64, &
    "Park-Miller portable sequence üçüncü değeri hatalı.")

  configuration%draw_count = 200
  configuration%confidence_level = 0.95_dp
  configuration%seed = 73_int64
  plan_a = create_cluster_bootstrap_plan([1, 2, 3], configuration)
  plan_b = create_cluster_bootstrap_plan([1, 2, 3], configuration)
  call assert_true(plan_a%status == BOOTSTRAP_SUCCESS .and. &
    plan_a%available .and. all(plan_a%campaign_indices == &
      plan_b%campaign_indices), &
    "Aynı seed aynı cluster draw planını üretmedi.")

  configuration%seed = 74_int64
  plan_c = create_cluster_bootstrap_plan([1, 2, 3], configuration)
  call assert_true(any(plan_a%campaign_indices /= plan_c%campaign_indices), &
    "Farklı deterministic seed farklı sequence üretmedi.")
  call assert_true(all(plan_a%campaign_indices >= 1) .and. &
    all(plan_a%campaign_indices <= 3), &
    "Cluster draw campaign population dışına çıktı.")
  has_replacement_duplicate = .false.
  do draw_index = 1, size(plan_a%campaign_indices, 2)
    if (count(plan_a%campaign_indices(:, draw_index) == &
        plan_a%campaign_indices(1, draw_index)) > 1) then
      has_replacement_duplicate = .true.
      exit
    end if
  end do
  call assert_true(has_replacement_duplicate, &
    "Bootstrap replacement sampling davranışı gözlenmedi.")

  call calculate_type7_quantile([0.0_dp, 10.0_dp], 0.25_dp, quantile, &
    quantile_status)
  call assert_true(quantile_status == BOOTSTRAP_SUCCESS, &
    "Valid Type-7 quantile unavailable döndü.")
  call assert_close(quantile, 2.5_dp, 1.0e-15_dp, &
    "Hyndman-Fan Type-7 interpolation hatalı.")
  call calculate_type7_quantile([0.0_dp, 10.0_dp], 0.5_dp, quantile, &
    quantile_status)
  call assert_close(quantile, 5.0_dp, 1.0e-15_dp, &
    "Type-7 median quantile hatalı.")
  call calculate_type7_quantile([0.0_dp, 10.0_dp], -0.1_dp, quantile, &
    quantile_status)
  call assert_true(quantile_status == BOOTSTRAP_INVALID_INPUT, &
    "Geçersiz quantile probability reddedilmedi.")

  ! A ve B aynı complete-campaign draw planını kullanır. B=10*A coupling'i
  ! her draw mean'inde korunmalıdır; quantity'ler bağımsız resample edilirse
  ! bu mandatory cluster-dependence regresyonu başarısız olur.
  usable_three = .true.
  result_a = calculate_cluster_bootstrap_mean( &
    [1.0_dp, 2.0_dp, 3.0_dp], usable_three, plan_a)
  result_b = calculate_cluster_bootstrap_mean( &
    [10.0_dp, 20.0_dp, 30.0_dp], usable_three, plan_a)
  call assert_true(result_a%interval%available .and. &
    result_b%interval%available .and. all(result_a%draw_available) .and. &
    all(result_b%draw_available), &
    "Complete coupled bootstrap draw'ları unavailable oldu.")
  do draw_index = 1, size(result_a%draw_means)
    call assert_close(result_b%draw_means(draw_index), &
      10.0_dp*result_a%draw_means(draw_index), 1.0e-15_dp, &
      "Whole-campaign draw coupling A/B quantities arasında kırıldı.")
  end do

  configuration%draw_count = 300
  configuration%seed = 91_int64
  plan_a = create_cluster_bootstrap_plan([1, 2, 3, 4, 5], configuration)
  usable_partial = [.true., .true., .false., .false., .false.]
  partial_result = calculate_cluster_bootstrap_mean( &
    [1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], usable_partial, plan_a)
  call assert_true(partial_result%interval%available .and. &
    partial_result%interval%valid_bootstrap_draw_count > 0 .and. &
    partial_result%interval%unavailable_bootstrap_draw_count > 0 .and. &
    partial_result%interval%valid_bootstrap_draw_count + &
      partial_result%interval%unavailable_bootstrap_draw_count == 300, &
    "Partial availability draw accounting hatalı.")

  configuration%confidence_level = 0.0_dp
  plan_a = create_cluster_bootstrap_plan([1, 2, 3], configuration)
  call assert_true(plan_a%status == BOOTSTRAP_INVALID_INPUT, &
    "confidence<=0 bootstrap configuration reddedilmedi.")
  configuration%confidence_level = 1.0_dp
  plan_a = create_cluster_bootstrap_plan([1, 2, 3], configuration)
  call assert_true(plan_a%status == BOOTSTRAP_INVALID_INPUT, &
    "confidence>=1 bootstrap configuration reddedilmedi.")
  configuration%confidence_level = 0.95_dp
  configuration%draw_count = 0
  plan_a = create_cluster_bootstrap_plan([1, 2, 3], configuration)
  call assert_true(plan_a%status == BOOTSTRAP_INVALID_INPUT, &
    "Nonpositive bootstrap draw count reddedilmedi.")

  print *, "V0.8.3 deterministic campaign-cluster bootstrap doğrulandı."
end program test_tts_bootstrap
