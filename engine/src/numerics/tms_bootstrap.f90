module tms_bootstrap
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  use tms_deterministic_rng, only : deterministic_rng_t
  implicit none
  private

  integer, parameter, public :: BOOTSTRAP_SUCCESS = 0
  integer, parameter, public :: BOOTSTRAP_INVALID_INPUT = 1
  integer, parameter, public :: BOOTSTRAP_INSUFFICIENT_CAMPAIGNS = 2
  integer, parameter, public :: BOOTSTRAP_NO_VALID_DRAWS = 3

  integer, parameter, public :: DEFAULT_BOOTSTRAP_DRAW_COUNT = 1000
  real(dp), parameter, public :: DEFAULT_BOOTSTRAP_CONFIDENCE_LEVEL = 0.95_dp
  integer(int64), parameter, public :: DEFAULT_BOOTSTRAP_SEED = &
    20260803_int64

  !> Nonparametric cluster bootstrap numerical ayarlarıdır. Draw count ve
  !! confidence level engineering kabul eşiği değildir. Seed explicit ve
  !! portable deterministic RNG state'inin başlangıcıdır.
  type, public :: bootstrap_configuration_t
    integer :: draw_count = DEFAULT_BOOTSTRAP_DRAW_COUNT
    real(dp) :: confidence_level = DEFAULT_BOOTSTRAP_CONFIDENCE_LEVEL
    integer(int64) :: seed = DEFAULT_BOOTSTRAP_SEED
  end type bootstrap_configuration_t

  !> Complete campaign population'ından replacement ile üretilmiş ortak draw
  !! planıdır. campaign_indices(:,draw), bütün derived quantities için aynen
  !! kullanılmalıdır; pair/isotherm/frequency düzeyinde bağımsız resampling
  !! yapılmaz.
  type, public :: cluster_bootstrap_plan_t
    integer :: status = BOOTSTRAP_INVALID_INPUT
    logical :: available = .false.
    integer :: population_count = 0
    integer :: requested_draw_count = 0
    real(dp) :: confidence_level = 0.0_dp
    integer(int64) :: seed = 0_int64
    integer, allocatable :: population_campaign_indices(:)
    integer, allocatable :: campaign_indices(:,:)
  end type cluster_bootstrap_plan_t

  !> Cohort mean için percentile bootstrap confidence interval sonucudur.
  !! lower/upper girdi quantity ile aynı birimdedir. Availability, requested,
  !! valid ve unavailable draw sayılarıyla birlikte taşınır; unavailable draw
  !! hiçbir zaman sıfır/NaN placeholder olarak quantile'a girmez.
  type, public :: bootstrap_interval_t
    integer :: status = BOOTSTRAP_INVALID_INPUT
    logical :: available = .false.
    integer :: requested_bootstrap_draw_count = 0
    integer :: valid_bootstrap_draw_count = 0
    integer :: unavailable_bootstrap_draw_count = 0
    real(dp) :: confidence_level = 0.0_dp
    integer(int64) :: seed = 0_int64
    real(dp) :: lower = 0.0_dp
    real(dp) :: upper = 0.0_dp
  end type bootstrap_interval_t

  !> Test ve coupled-quantity doğrulaması için draw-level mean izini taşır.
  !! draw_available=false olan sayısal hücreler kullanılmamalıdır.
  type, public :: bootstrap_mean_result_t
    type(bootstrap_interval_t) :: interval
    real(dp), allocatable :: draw_means(:)
    logical, allocatable :: draw_available(:)
  end type bootstrap_mean_result_t

  public :: is_valid_bootstrap_configuration
  public :: create_cluster_bootstrap_plan
  public :: calculate_cluster_bootstrap_mean
  public :: calculate_type7_quantile

contains

  !> Bootstrap ayarlarının draw_count>0 ve 0<confidence<1 koşullarını sınar.
  !! Confidence boyutsuzdur; bu validation material acceptance belirtmez.
  pure function is_valid_bootstrap_configuration(configuration) result(valid)
    type(bootstrap_configuration_t), intent(in) :: configuration
    logical :: valid

    valid = configuration%draw_count > 0 .and. &
      ieee_is_finite(configuration%confidence_level) .and. &
      configuration%confidence_level > 0.0_dp .and. &
      configuration%confidence_level < 1.0_dp
  end function is_valid_bootstrap_configuration

  !> Eligible independent complete campaign kimliklerinden R x B ortak
  !! cluster draw planı üretir. Her draw R kez replacement sampling yapar.
  !! Aynı plan bütün quantities'e verilerek campaign içi dependence korunur.
  !! En az iki bağımsız campaign yoksa uncertainty planı unavailable olur.
  pure function create_cluster_bootstrap_plan( &
      campaign_indices, configuration) result(plan)
    integer, intent(in) :: campaign_indices(:)
    type(bootstrap_configuration_t), intent(in) :: configuration
    type(cluster_bootstrap_plan_t) :: plan

    type(deterministic_rng_t) :: rng
    integer :: draw_index
    integer :: population_position
    integer :: sampled_position
    integer :: sample_index
    logical :: valid_index

    plan%population_count = size(campaign_indices)
    plan%requested_draw_count = configuration%draw_count
    plan%confidence_level = configuration%confidence_level
    plan%seed = configuration%seed
    if (.not. is_valid_bootstrap_configuration(configuration) .or. &
        any(campaign_indices <= 0) .or. &
        campaign_indices_have_duplicates(campaign_indices)) then
      plan%status = BOOTSTRAP_INVALID_INPUT
      return
    end if
    if (size(campaign_indices) < 2) then
      plan%status = BOOTSTRAP_INSUFFICIENT_CAMPAIGNS
      return
    end if

    plan%population_campaign_indices = campaign_indices
    allocate(plan%campaign_indices(size(campaign_indices), &
      configuration%draw_count))
    call rng%initialize(configuration%seed)
    do draw_index = 1, configuration%draw_count
      do sample_index = 1, size(campaign_indices)
        call rng%next_index(size(campaign_indices), sampled_position, &
          valid_index)
        if (.not. valid_index) then
          plan%status = BOOTSTRAP_INVALID_INPUT
          deallocate(plan%campaign_indices)
          return
        end if
        population_position = sampled_position
        plan%campaign_indices(sample_index, draw_index) = &
          campaign_indices(population_position)
      end do
    end do
    plan%status = BOOTSTRAP_SUCCESS
    plan%available = .true.
  end function create_cluster_bootstrap_plan

  !> Aynı cluster draw planı altında quantity-specific usable campaign
  !! değerlerinden cohort mean dağılımını ve percentile CI'yı üretir. Global
  !! population'da en az iki usable independent campaign aranır; her draw'da
  !! en az iki usable sampled value gerekir. Missing derived results draw'dan
  !! sonra elenir ve hiçbir placeholder istatistiğe girmez.
  pure function calculate_cluster_bootstrap_mean( &
      values, usable, plan) result(result)
    real(dp), intent(in) :: values(:)
    logical, intent(in) :: usable(:)
    type(cluster_bootstrap_plan_t), intent(in) :: plan
    type(bootstrap_mean_result_t) :: result

    real(dp), allocatable :: valid_draw_means(:)
    real(dp) :: alpha
    real(dp) :: delta
    real(dp) :: draw_mean
    integer :: campaign_index
    integer :: draw_index
    integer :: quantile_status
    integer :: sample_index
    integer :: usable_count
    integer :: usable_population_count
    logical :: finite_draw

    result%interval%requested_bootstrap_draw_count = &
      plan%requested_draw_count
    result%interval%confidence_level = plan%confidence_level
    result%interval%seed = plan%seed
    if (.not. plan%available .or. plan%status /= BOOTSTRAP_SUCCESS .or. &
        size(values) /= size(usable) .or. &
        .not. allocated(plan%campaign_indices) .or. &
        .not. allocated(plan%population_campaign_indices)) then
      result%interval%status = BOOTSTRAP_INVALID_INPUT
      return
    end if
    if (any(plan%population_campaign_indices > size(values))) then
      result%interval%status = BOOTSTRAP_INVALID_INPUT
      return
    end if
    if (any(usable .and. .not. ieee_is_finite(values))) then
      result%interval%status = BOOTSTRAP_INVALID_INPUT
      return
    end if

    usable_population_count = 0
    do sample_index = 1, size(plan%population_campaign_indices)
      campaign_index = plan%population_campaign_indices(sample_index)
      if (usable(campaign_index)) usable_population_count = &
        usable_population_count + 1
    end do
    if (usable_population_count < 2) then
      result%interval%status = BOOTSTRAP_INSUFFICIENT_CAMPAIGNS
      result%interval%unavailable_bootstrap_draw_count = &
        plan%requested_draw_count
      return
    end if

    allocate(result%draw_means(plan%requested_draw_count), &
      result%draw_available(plan%requested_draw_count))
    result%draw_means = 0.0_dp
    result%draw_available = .false.
    do draw_index = 1, plan%requested_draw_count
      draw_mean = 0.0_dp
      usable_count = 0
      finite_draw = .true.
      do sample_index = 1, size(plan%campaign_indices, 1)
        campaign_index = plan%campaign_indices(sample_index, draw_index)
        if (campaign_index < 1 .or. campaign_index > size(values)) then
          finite_draw = .false.
          exit
        end if
        if (.not. usable(campaign_index)) cycle
        usable_count = usable_count + 1
        delta = values(campaign_index) - draw_mean
        if (.not. ieee_is_finite(delta)) then
          finite_draw = .false.
          exit
        end if
        draw_mean = draw_mean + delta/real(usable_count, dp)
        if (.not. ieee_is_finite(draw_mean)) then
          finite_draw = .false.
          exit
        end if
      end do
      if (finite_draw .and. usable_count >= 2) then
        result%draw_means(draw_index) = draw_mean
        result%draw_available(draw_index) = .true.
      end if
    end do

    result%interval%valid_bootstrap_draw_count = &
      count(result%draw_available)
    result%interval%unavailable_bootstrap_draw_count = &
      plan%requested_draw_count - &
      result%interval%valid_bootstrap_draw_count
    if (result%interval%valid_bootstrap_draw_count == 0) then
      result%interval%status = BOOTSTRAP_NO_VALID_DRAWS
      return
    end if

    valid_draw_means = pack(result%draw_means, result%draw_available)
    alpha = 1.0_dp - plan%confidence_level
    call calculate_type7_quantile(valid_draw_means, 0.5_dp*alpha, &
      result%interval%lower, quantile_status)
    if (quantile_status /= BOOTSTRAP_SUCCESS) then
      result%interval%status = quantile_status
      return
    end if
    call calculate_type7_quantile(valid_draw_means, &
      1.0_dp - 0.5_dp*alpha, result%interval%upper, quantile_status)
    if (quantile_status /= BOOTSTRAP_SUCCESS) then
      result%interval%status = quantile_status
      return
    end if
    result%interval%status = BOOTSTRAP_SUCCESS
    result%interval%available = .true.
  end function calculate_cluster_bootstrap_mean

  !> Hyndman-Fan Type-7 quantile hesaplar: h=1+(n-1)p, j=floor(h),
  !! gamma=h-j ve Q=(1-gamma)x_j+gamma*x_(j+1). p boyutsuz [0,1],
  !! values/quantile aynı birimdedir. Input sort edilmeden kopyalanır.
  pure subroutine calculate_type7_quantile(values, probability, quantile, &
      status)
    real(dp), intent(in) :: values(:)
    real(dp), intent(in) :: probability
    real(dp), intent(out) :: quantile
    integer, intent(out) :: status

    real(dp), allocatable :: sorted_values(:)
    real(dp) :: gamma
    real(dp) :: h
    integer :: lower_index

    quantile = 0.0_dp
    status = BOOTSTRAP_INVALID_INPUT
    if (size(values) == 0 .or. .not. all(ieee_is_finite(values)) .or. &
        .not. ieee_is_finite(probability) .or. probability < 0.0_dp .or. &
        probability > 1.0_dp) return

    sorted_values = values
    call sort_ascending(sorted_values)
    h = 1.0_dp + real(size(values) - 1, dp)*probability
    lower_index = int(floor(h))
    gamma = h - real(lower_index, dp)
    if (lower_index >= size(values)) then
      quantile = sorted_values(size(values))
    else
      quantile = sorted_values(lower_index) + gamma * &
        (sorted_values(lower_index + 1) - sorted_values(lower_index))
    end if
    if (.not. ieee_is_finite(quantile)) return
    status = BOOTSTRAP_SUCCESS
  end subroutine calculate_type7_quantile

  pure function campaign_indices_have_duplicates(indices) result(duplicate)
    integer, intent(in) :: indices(:)
    logical :: duplicate
    integer :: i
    integer :: j

    duplicate = .false.
    do i = 2, size(indices)
      do j = 1, i - 1
        if (indices(i) == indices(j)) then
          duplicate = .true.
          return
        end if
      end do
    end do
  end function campaign_indices_have_duplicates

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

end module tms_bootstrap
