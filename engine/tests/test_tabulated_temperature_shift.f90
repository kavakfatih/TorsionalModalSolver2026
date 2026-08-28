module tms_test_inconsistent_temperature_shift
  use tms_kinds, only : dp
  use tms_temperature_shift_types, only : temperature_shift_domain_t, &
    temperature_shift_evaluation_t, WLF_TEMPERATURE_SHIFT
  use tms_temperature_shift_provider, only : temperature_shift_provider_t, &
    validate_temperature_shift_domain
  implicit none
  private

  !> Ortak provider-boundary doğrulamasını sınamak için kasıtlı olarak
  !! s=0 ile tutarsız a_T=2 döndüren test double'ıdır; fizik modeli değildir.
  type, extends(temperature_shift_provider_t), public :: &
      inconsistent_temperature_shift_provider_t
    private
    type(temperature_shift_domain_t) :: domain
  contains
    procedure, public :: evaluate => evaluate_inconsistent_shift
    procedure, public :: validate => validate_inconsistent_shift
    procedure, public :: get_domain => get_inconsistent_shift_domain
  end type inconsistent_temperature_shift_provider_t

  public :: create_inconsistent_temperature_shift_provider

contains

  pure function create_inconsistent_temperature_shift_provider() &
      result(provider)
    type(inconsistent_temperature_shift_provider_t) :: provider

    provider%domain = temperature_shift_domain_t( &
      minimum_temperature_k=273.15_dp, &
      maximum_temperature_k=313.15_dp, &
      reference_temperature_k=293.15_dp)
  end function create_inconsistent_temperature_shift_provider

  pure function evaluate_inconsistent_shift(self, temperature_k) &
      result(evaluation)
    class(inconsistent_temperature_shift_provider_t), intent(in) :: self
    real(dp), intent(in) :: temperature_k
    type(temperature_shift_evaluation_t) :: evaluation

    call self%validate()
    evaluation%shift_model_kind = WLF_TEMPERATURE_SHIFT
    evaluation%operating_temperature_k = temperature_k
    evaluation%reference_temperature_k = &
      self%domain%reference_temperature_k
    evaluation%log10_a_t = 0.0_dp
    evaluation%a_t = 2.0_dp
  end function evaluate_inconsistent_shift

  pure subroutine validate_inconsistent_shift(self)
    class(inconsistent_temperature_shift_provider_t), intent(in) :: self

    call validate_temperature_shift_domain(self%domain)
  end subroutine validate_inconsistent_shift

  pure function get_inconsistent_shift_domain(self) result(domain)
    class(inconsistent_temperature_shift_provider_t), intent(in) :: self
    type(temperature_shift_domain_t) :: domain

    call self%validate()
    domain = self%domain
  end function get_inconsistent_shift_domain

end module tms_test_inconsistent_temperature_shift

program test_tabulated_temperature_shift
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_quiet_nan, ieee_positive_inf
  use tms_kinds, only : dp
  use tms_temperature_shift_types, only : temperature_shift_evaluation_t, &
    TABULATED_LOG10_SHIFT
  use tms_temperature_shift_provider, only : evaluate_temperature_shift, &
    validate_temperature_shift_evaluation
  use tms_tabulated_temperature_shift, only : &
    tabulated_log10_shift_provider_t, &
    create_tabulated_temperature_shift_provider
  use tms_test_inconsistent_temperature_shift, only : &
    inconsistent_temperature_shift_provider_t, &
    create_inconsistent_temperature_shift_provider
  implicit none

  real(dp), parameter :: tolerance = 2.0e-13_dp
  character(len=80) :: validation_case

  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_exact_and_reference_points()
  call test_log10_shift_interpolation()
  call test_nonmonotonic_shift_data()

  print *, "V0.8 tabulated log10 temperature-shift provider doğrulandı."

contains

  !> Stored sıcaklık noktalarının exact olarak bulunduğunu ve T_ref için
  !! s=log10(a_T)=0, a_T=1 olduğunu doğrular. T [K], s ve a_T boyutsuzdur.
  subroutine test_exact_and_reference_points()
    type(tabulated_log10_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: evaluation
    real(dp) :: temperature_k(3)
    real(dp) :: log10_a_t(3)

    provider = make_provider()
    evaluation = evaluate_temperature_shift(provider, 293.15_dp)

    if (evaluation%shift_model_kind /= TABULATED_LOG10_SHIFT .or. &
        .not. evaluation%exact_temperature_point) then
      error stop "Reference sıcaklığı exact tabulated shift noktası değil."
    end if
    call assert_close(evaluation%operating_temperature_k, &
      293.15_dp, 0.0_dp, "Shift trace operating sıcaklığı hatalı.")
    call assert_close(evaluation%reference_temperature_k, &
      293.15_dp, 0.0_dp, "Shift trace reference sıcaklığı hatalı.")
    call assert_close(evaluation%log10_a_t, 0.0_dp, 0.0_dp, &
      "Reference sıcaklığında log10(a_T) sıfır değil.")
    call assert_close(evaluation%a_t, 1.0_dp, 0.0_dp, &
      "Reference sıcaklığında a_T bir değil.")

    ! Bir ULP fark yalnız representation seviyesindedir ve aynı stored
    ! sıcaklık noktasını seçmelidir; bu bir fiziksel tolerance değildir.
    evaluation = evaluate_temperature_shift( &
      provider, nearest(293.15_dp, 1.0_dp))
    if (.not. evaluation%exact_temperature_point) then
      error stop "Bir ULP sıcaklık farkı exact point kabul edilmedi."
    end if
    call assert_close(evaluation%log10_a_t, 0.0_dp, 0.0_dp, &
      "Numerical exact reference shift değeri değişti.")

    ! Ölçüm girdisindeki machine-zero reference shift kabul edilir; runtime
    ! identity ise drift yaratmamak için tam s=0 ve a_T=1 olarak döndürülür.
    temperature_k = [273.15_dp, 293.15_dp, 313.15_dp]
    log10_a_t = [1.0_dp, 32.0_dp*epsilon(1.0_dp), -2.0_dp]
    provider = create_tabulated_temperature_shift_provider( &
      temperature_k, log10_a_t, 293.15_dp)
    evaluation = evaluate_temperature_shift(provider, 293.15_dp)
    call assert_close(evaluation%log10_a_t, 0.0_dp, 0.0_dp, &
      "Machine-zero reference shift canonical sıfıra çevrilmedi.")
    call assert_close(evaluation%a_t, 1.0_dp, 0.0_dp, &
      "Machine-zero reference shift canonical a_T=1 üretmedi.")
  end subroutine test_exact_and_reference_points

  !> Tabulated modelin a_T değerini değil authoritative s=log10(a_T)
  !! değerini lineer interpolate ettiğini analitik midpoint ile doğrular.
  !! 293.15 K'de s=0 ve 313.15 K'de s=-2 için 303.15 K'de s=-1,
  !! dolayısıyla a_T=0.1 olmalıdır.
  subroutine test_log10_shift_interpolation()
    type(tabulated_log10_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: evaluation

    provider = make_provider()
    evaluation = evaluate_temperature_shift(provider, 303.15_dp)

    if (evaluation%exact_temperature_point) then
      error stop "Ara sıcaklık exact tabulated shift noktası işaretlendi."
    end if
    call assert_close(evaluation%lower_temperature_k, &
      293.15_dp, tolerance, "Shift lower temperature bracket hatalı.")
    call assert_close(evaluation%upper_temperature_k, &
      313.15_dp, tolerance, "Shift upper temperature bracket hatalı.")
    call assert_close(evaluation%interpolation_alpha, &
      0.5_dp, tolerance, "Shift temperature interpolation alpha hatalı.")
    call assert_close(evaluation%log10_a_t, -1.0_dp, tolerance, &
      "log10(a_T) midpoint interpolation sonucu hatalı.")
    call assert_close(evaluation%a_t, 0.1_dp, tolerance, &
      "Derived a_T=10**s sonucu hatalı.")
  end subroutine test_log10_shift_interpolation

  !> Ölçülmüş shift data için monotonicity'nin data-quality kuralı olarak
  !! dayatılmadığını doğrular. Sıcaklıklar strictly increasing kalırken
  !! s değerleri yön değiştirebilir; provider yalnız lineer s interpolation yapar.
  subroutine test_nonmonotonic_shift_data()
    real(dp) :: temperature_k(3)
    real(dp) :: log10_a_t(3)
    type(tabulated_log10_shift_provider_t) :: provider
    type(temperature_shift_evaluation_t) :: evaluation

    temperature_k = [273.15_dp, 293.15_dp, 313.15_dp]
    log10_a_t = [0.5_dp, 0.0_dp, 0.25_dp]
    provider = create_tabulated_temperature_shift_provider( &
      temperature_k, log10_a_t, 293.15_dp)
    evaluation = evaluate_temperature_shift(provider, 303.15_dp)

    call assert_close(evaluation%log10_a_t, 0.125_dp, tolerance, &
      "Nonmonotonic measured shift data interpolation'ı hatalı.")
  end subroutine test_nonmonotonic_shift_data

  !> Constructor ve query error-stop yollarını bağımsız CTest süreçleri için
  !! tetikler. Extrapolation, invalid SI sıcaklıkları ve bozuk tablo reddedilir.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    real(dp) :: temperature_k(3)
    real(dp) :: log10_a_t(3)
    real(dp) :: nan_value
    real(dp) :: infinity
    type(tabulated_log10_shift_provider_t) :: provider
    type(inconsistent_temperature_shift_provider_t) :: &
      inconsistent_provider
    type(temperature_shift_evaluation_t) :: evaluation

    temperature_k = [273.15_dp, 293.15_dp, 313.15_dp]
    log10_a_t = [1.0_dp, 0.0_dp, -2.0_dp]
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)

    select case (case_name)
    case ("less_than_two_points")
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k(1:1), log10_a_t(1:1), 273.15_dp)
    case ("size_mismatch")
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t(1:2), 293.15_dp)
    case ("zero_temperature")
      temperature_k(1) = 0.0_dp
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("negative_temperature")
      temperature_k(1) = -1.0_dp
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("nan_temperature")
      temperature_k(1) = nan_value
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("infinite_temperature")
      temperature_k(3) = infinity
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("duplicate_temperature")
      temperature_k(2) = temperature_k(1)
      log10_a_t(1) = 0.0_dp
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, temperature_k(1))
    case ("numerical_duplicate_temperature")
      temperature_k(2) = nearest(temperature_k(1), 1.0_dp)
      log10_a_t(1) = 0.0_dp
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, temperature_k(1))
    case ("unordered_temperature")
      temperature_k(2) = 270.0_dp
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("nan_shift")
      log10_a_t(1) = nan_value
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("infinite_shift")
      log10_a_t(3) = infinity
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("invalid_reference_temperature")
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 0.0_dp)
    case ("missing_reference_temperature")
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 300.0_dp)
    case ("nonzero_reference_shift")
      log10_a_t(2) = 1.0e-6_dp
      provider = create_tabulated_temperature_shift_provider( &
        temperature_k, log10_a_t, 293.15_dp)
    case ("below_domain")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 273.0_dp)
    case ("above_domain")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 313.2_dp)
    case ("inconsistent_shift_factor")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 293.15_dp)
      evaluation%a_t = 2.0_dp
      call validate_temperature_shift_evaluation( &
        provider, 293.15_dp, evaluation)
    case ("invalid_shift_model_kind")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 293.15_dp)
      evaluation%shift_model_kind = 99
      call validate_temperature_shift_evaluation( &
        provider, 293.15_dp, evaluation)
    case ("invalid_temperature_bracket")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 303.15_dp)
      evaluation%lower_temperature_k = 310.0_dp
      call validate_temperature_shift_evaluation( &
        provider, 303.15_dp, evaluation)
    case ("nonfinite_temperature_bracket")
      provider = make_provider()
      evaluation = evaluate_temperature_shift(provider, 303.15_dp)
      evaluation%lower_temperature_k = nan_value
      call validate_temperature_shift_evaluation( &
        provider, 303.15_dp, evaluation)
    case ("inconsistent_provider_evaluation")
      inconsistent_provider = &
        create_inconsistent_temperature_shift_provider()
      evaluation = evaluate_temperature_shift( &
        inconsistent_provider, 293.15_dp)
    case default
      error stop "Bilinmeyen tabulated temperature-shift selector."
    end select
  end subroutine exercise_invalid_case

  pure function make_provider() result(provider)
    type(tabulated_log10_shift_provider_t) :: provider
    real(dp) :: temperature_k(3)
    real(dp) :: log10_a_t(3)

    temperature_k = [273.15_dp, 293.15_dp, 313.15_dp]
    log10_a_t = [1.0_dp, 0.0_dp, -2.0_dp]
    provider = create_tabulated_temperature_shift_provider( &
      temperature_k, log10_a_t, 293.15_dp)
  end function make_provider

  subroutine assert_close(actual, expected, relative_tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: relative_tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(relative_tolerance) .or. &
        relative_tolerance < 0.0_dp .or. abs(actual-expected) > &
        relative_tolerance*max(1.0_dp, abs(expected))) then
      error stop message
    end if
  end subroutine assert_close

end program test_tabulated_temperature_shift
