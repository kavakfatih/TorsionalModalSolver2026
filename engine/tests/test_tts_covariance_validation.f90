program test_tts_covariance_validation
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
    ieee_positive_inf
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t
  use tms_tts_covariance_types, only : tts_covariance_matrix_2x2_t, &
    tts_covariance_matrix_validation_t, &
    tts_dynamic_modulus_covariance_family_t, &
    tts_dynamic_modulus_covariance_point_t, &
    tts_covariance_validation_result_t, &
    tts_bivariate_covariance_log_segment_t, TTS_COVARIANCE_SUCCESS, &
    TTS_COVARIANCE_DATA_MISMATCH, TTS_COVARIANCE_SINGULAR_MATRIX, &
    TTS_COVARIANCE_ILL_CONDITIONED
  use tms_tts_covariance_propagation, only : &
    validate_covariance_matrix_2x2
  use tms_tts_covariance_validation, only : &
    validate_tts_covariance_family, &
    build_tts_bivariate_covariance_log_segments
  use tms_tts_test_support, only : make_exact_trs_family, &
    populate_isotherm_from_log_grid, assert_true, assert_close
  use tms_tts_uncertainty_test_support, only : &
    make_relative_uncertainty_family
  use tms_tts_covariance_test_support, only : make_covariance_family
  implicit none

  type(tts_material_family_t) :: family
  type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty
  type(tts_dynamic_modulus_covariance_family_t) :: covariance_family
  type(tts_dynamic_modulus_covariance_family_t) :: changed_family
  type(tts_dynamic_modulus_covariance_point_t) :: temporary_point
  type(tts_covariance_matrix_2x2_t) :: matrix
  type(tts_covariance_matrix_validation_t) :: matrix_validation
  type(tts_covariance_validation_result_t) :: validation
  type(tts_bivariate_covariance_log_segment_t), allocatable :: segments(:)
  real(dp), parameter :: temperatures(2) = [293.15_dp, 313.15_dp]
  real(dp) :: five_point_grid(5)
  integer :: i
  integer :: n

  ! Symmetric 2x2 storage formunda pozitif diyagonaller, det>0 ve machine-safe
  ! reciprocal condition birlikte SPD kabulü üretmelidir.
  matrix%storage_variance_pa2 = 4.0_dp
  matrix%loss_variance_pa2 = 9.0_dp
  matrix%storage_loss_covariance_pa2 = 3.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(matrix_validation%status == TTS_COVARIANCE_SUCCESS .and. &
    matrix_validation%covariance_valid .and. &
    matrix_validation%covariance_numerically_well_conditioned, &
    "Geçerli symmetric SPD covariance reddedildi.")
  call assert_close(matrix_validation%determinant, 27.0_dp, 2.0e-15_dp, &
    "SPD determinant tanısı hatalı.")
  call assert_close(matrix_validation%correlation, 0.5_dp, 2.0e-15_dp, &
    "SPD correlation tanısı hatalı.")

  ! Geçersiz diyagonal, nonfinite, determinant<=0, perfect correlation ve
  ! materially near-singular covariance hiçbir jitter/clipping ile kurtarılmaz.
  matrix%storage_variance_pa2 = 0.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "Sıfır storage variance reddedilmedi.")
  matrix%storage_variance_pa2 = -1.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "Negatif storage variance reddedilmedi.")
  matrix%storage_variance_pa2 = 4.0_dp
  matrix%loss_variance_pa2 = 0.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "Sıfır loss variance reddedilmedi.")
  matrix%loss_variance_pa2 = -9.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "Negatif loss variance reddedilmedi.")
  matrix%storage_variance_pa2 = ieee_value(1.0_dp, ieee_quiet_nan)
  matrix%loss_variance_pa2 = 9.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "NaN variance reddedilmedi.")
  matrix%storage_variance_pa2 = 4.0_dp
  matrix%loss_variance_pa2 = ieee_value(1.0_dp, ieee_positive_inf)
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "Inf variance reddedilmedi.")
  matrix%loss_variance_pa2 = 9.0_dp
  matrix%storage_loss_covariance_pa2 = &
    ieee_value(1.0_dp, ieee_quiet_nan)
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "NaN covariance reddedilmedi.")
  matrix%storage_loss_covariance_pa2 = &
    ieee_value(1.0_dp, ieee_positive_inf)
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "Inf covariance reddedilmedi.")
  matrix%storage_loss_covariance_pa2 = 7.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(.not. matrix_validation%covariance_valid, &
    "Negatif determinant covariance reddedilmedi.")
  matrix%storage_loss_covariance_pa2 = 6.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(matrix_validation%status == TTS_COVARIANCE_SINGULAR_MATRIX, &
    "rho=+1 singular covariance açık status vermedi.")
  matrix%storage_loss_covariance_pa2 = -6.0_dp
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(matrix_validation%status == TTS_COVARIANCE_SINGULAR_MATRIX, &
    "rho=-1 singular covariance açık status vermedi.")
  matrix%storage_loss_covariance_pa2 = 6.0_dp*(1.0_dp - 1.0e-10_dp)
  matrix_validation = validate_covariance_matrix_2x2(matrix)
  call assert_true(matrix_validation%status == &
    TTS_COVARIANCE_ILL_CONDITIONED .and. &
    matrix_validation%covariance_valid .and. &
    .not. matrix_validation%covariance_numerically_well_conditioned, &
    "Near-singular covariance machine conditioning ile reddedilmedi.")

  ! Physical-key eşlemesi array sırasından bağımsızdır; full covariance
  ! diyagonalleri V0.8.4 u_G^2 ile machine-equivalent olmak zorundadır.
  family = make_exact_trs_family(temperatures, [0.0_dp, 0.5_dp])
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  covariance_family = make_covariance_family(family, uncertainty, 0.25_dp)
  validation = validate_tts_covariance_family( &
    family, uncertainty, covariance_family)
  call assert_true(validation%valid .and. &
    validation%status == TTS_COVARIANCE_SUCCESS, &
    "Physical-key covariance family doğrulanamadı.")
  do i = 1, size(covariance_family%isotherms)
    n = size(covariance_family%isotherms(i)%points)
    temporary_point = covariance_family%isotherms(i)%points(1)
    covariance_family%isotherms(i)%points(1) = &
      covariance_family%isotherms(i)%points(n)
    covariance_family%isotherms(i)%points(n) = temporary_point
  end do
  validation = validate_tts_covariance_family( &
    family, uncertainty, covariance_family)
  call assert_true(validation%valid, &
    "Covariance physical matching array-order'a bağımlı kaldı.")

  changed_family = covariance_family
  changed_family%isotherms(1)%points(1)%covariance &
    %storage_variance_pa2 = 1.1_dp * &
    changed_family%isotherms(1)%points(1)%covariance &
      %storage_variance_pa2
  validation = validate_tts_covariance_family( &
    family, uncertainty, changed_family)
  call assert_true(validation%status == TTS_COVARIANCE_DATA_MISMATCH, &
    "Çelişkili V0.8.4/V0.8.5 diagonal kaynakları reddedilmedi.")

  changed_family = covariance_family
  changed_family%isotherms(1)%points = [ &
    changed_family%isotherms(1)%points, &
    changed_family%isotherms(1)%points(1)]
  validation = validate_tts_covariance_family( &
    family, uncertainty, changed_family)
  call assert_true(validation%status == TTS_COVARIANCE_DATA_MISMATCH, &
    "Duplicate covariance physical key reddedilmedi.")

  ! 1,2,5,10,20 Hz fixture'ında 5 Hz covariance gap'i support'u [1,2] ve
  ! [10,20] olarak böler; 2->10 arasında interpolation bridge oluşmaz.
  five_point_grid = [0.0_dp, log10(2.0_dp), log10(5.0_dp), &
    1.0_dp, log10(20.0_dp)]
  family = make_exact_trs_family(temperatures, [0.0_dp, 0.5_dp])
  do i = 1, size(family%isotherms)
    call populate_isotherm_from_log_grid(family%isotherms(i), &
      five_point_grid, merge(0.0_dp, 0.5_dp, i == 1), 0.25_dp, 0.15_dp)
  end do
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.03_dp)
  covariance_family = make_covariance_family(family, uncertainty, 0.2_dp)
  covariance_family%isotherms(1)%points(3)%covariance_available = .false.
  segments = build_tts_bivariate_covariance_log_segments( &
    family%isotherms(1), uncertainty, covariance_family)
  call assert_true(size(segments) == 2, &
    "Covariance gap contiguous bivariate support'u ikiye bölmedi.")
  call assert_close(segments(1)%x(1), 0.0_dp, 0.0_dp, &
    "İlk covariance segment alt sınırı hatalı.")
  call assert_close(segments(1)%x(size(segments(1)%x)), log10(2.0_dp), &
    2.0e-15_dp, "İlk covariance segment üst sınırı hatalı.")
  call assert_close(segments(2)%x(1), 1.0_dp, 2.0e-15_dp, &
    "İkinci covariance segment alt sınırı gap üzerinden taşındı.")
  call assert_close(segments(2)%x(size(segments(2)%x)), log10(20.0_dp), &
    2.0e-15_dp, "İkinci covariance segment üst sınırı hatalı.")

  print *, "V0.8.5 covariance validation, matching ve gaps doğrulandı."
end program test_tts_covariance_validation
