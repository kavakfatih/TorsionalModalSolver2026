program test_tts_mahalanobis_pair_shift
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_pair_shift_configuration_t
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t
  use tms_tts_covariance_types, only : tts_covariance_matrix_2x2_t, &
    tts_dynamic_modulus_covariance_family_t, &
    tts_covariance_objective_evaluation_t, &
    tts_covariance_pair_solution_t, MAHALANOBIS_SHIFT_SUCCESS, &
    MAHALANOBIS_SHIFT_INVALID_COVARIANCE, &
    MAHALANOBIS_SHIFT_NO_INTERIOR_MINIMUM
  use tms_tts_covariance_propagation, only : squared_mahalanobis_2x2
  use tms_tts_mahalanobis_pair_shift, only : &
    evaluate_tts_covariance_pair_objective, &
    identify_tts_covariance_pair_shift
  use tms_tts_test_support, only : make_exact_trs_family, &
    make_generalized_maxwell_trs_family, replace_loss_truth_shift, &
    assert_true, assert_close
  use tms_tts_uncertainty_test_support, only : &
    make_relative_uncertainty_family
  use tms_tts_covariance_test_support, only : &
    make_covariance_family, scale_covariance_inputs
  implicit none

  type(tts_material_family_t) :: family
  type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty
  type(tts_dynamic_modulus_uncertainty_family_t) :: scaled_uncertainty
  type(tts_dynamic_modulus_covariance_family_t) :: covariance_zero
  type(tts_dynamic_modulus_covariance_family_t) :: covariance_positive
  type(tts_dynamic_modulus_covariance_family_t) :: covariance_negative
  type(tts_dynamic_modulus_covariance_family_t) :: scaled_covariance
  type(tts_pair_shift_configuration_t) :: configuration
  type(tts_covariance_pair_solution_t) :: zero_solution
  type(tts_covariance_pair_solution_t) :: positive_solution
  type(tts_covariance_pair_solution_t) :: negative_solution
  type(tts_covariance_pair_solution_t) :: scaled_solution
  type(tts_covariance_objective_evaluation_t) :: evaluation
  type(tts_covariance_objective_evaluation_t) :: scaled_evaluation
  type(tts_covariance_matrix_2x2_t) :: point_covariance
  real(dp), parameter :: temperatures(2) = [293.15_dp, 313.15_dp]
  real(dp), parameter :: truth_shift = 0.65_dp
  real(dp) :: same_direction
  real(dp) :: opposite_direction
  logical :: valid
  integer :: i

  ! Low-level covariance ellipse physics: pozitif rho same-sign residual
  ! yönünü, negatif rho opposite-sign yönünü daha düşük d_M^2 ile ağırlıklar.
  point_covariance%storage_variance_pa2 = 1.0_dp
  point_covariance%loss_variance_pa2 = 1.0_dp
  point_covariance%storage_loss_covariance_pa2 = 0.6_dp
  call squared_mahalanobis_2x2( &
    1.0_dp, 1.0_dp, point_covariance, same_direction, valid)
  call assert_true(valid, "Pozitif rho same-direction point solve başarısız.")
  call squared_mahalanobis_2x2( &
    1.0_dp, -1.0_dp, point_covariance, opposite_direction, valid)
  call assert_true(valid .and. same_direction < opposite_direction, &
    "Pozitif rho covariance ellipse yönü ters.")
  point_covariance%storage_loss_covariance_pa2 = -0.6_dp
  call squared_mahalanobis_2x2( &
    1.0_dp, 1.0_dp, point_covariance, same_direction, valid)
  call squared_mahalanobis_2x2( &
    1.0_dp, -1.0_dp, point_covariance, opposite_direction, valid)
  call assert_true(valid .and. same_direction > opposite_direction, &
    "Negatif rho covariance ellipse yönü tersine dönmedi.")

  ! Zero off-diagonal covariance'da Mahalanobis ve matched-support diagonal
  ! aynı O_B objective/minimizer'ına indirgenmelidir.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  covariance_zero = make_covariance_family(family, uncertainty, 0.0_dp)
  zero_solution = identify_tts_covariance_pair_shift( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_zero, configuration)
  call assert_true(zero_solution%diagonal_matched%status == &
    MAHALANOBIS_SHIFT_SUCCESS .and. &
    zero_solution%mahalanobis%status == MAHALANOBIS_SHIFT_SUCCESS, &
    "Zero-covariance pair shift çözülemedi.")
  call assert_close(zero_solution%mahalanobis%shift, &
    zero_solution%diagonal_matched%shift, 4.0e-8_dp, &
    "Zero covariance diagonal-limit minimizer invariant bozuldu.")
  call assert_close(zero_solution%mahalanobis%objective_minimum, &
    zero_solution%diagonal_matched%objective_minimum, 2.0e-10_dp, &
    "Zero covariance diagonal-limit objective invariant bozuldu.")

  ! Storage ve loss kanallarının tercih ettiği shift'ler ayrıldığında kontrollü
  ! positive/negative rho Mahalanobis residual geometry'sini öngörülen zıt
  ! yönlerde değiştirir; yalnız "shift changed" kontrolü yapılmaz.
  family = make_exact_trs_family(temperatures, [0.0_dp, 0.60_dp])
  call replace_loss_truth_shift(family%isotherms(2), 0.20_dp)
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  covariance_zero = make_covariance_family(family, uncertainty, 0.0_dp)
  covariance_positive = make_covariance_family(family, uncertainty, 0.5_dp)
  covariance_negative = make_covariance_family(family, uncertainty, -0.5_dp)
  zero_solution = identify_tts_covariance_pair_shift( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_zero, configuration)
  positive_solution = identify_tts_covariance_pair_shift( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_positive, configuration)
  negative_solution = identify_tts_covariance_pair_shift( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_negative, configuration)
  call assert_true(zero_solution%mahalanobis%shift_available .and. &
    positive_solution%mahalanobis%shift_available .and. &
    negative_solution%mahalanobis%shift_available, &
    "Correlation-direction pair fixture çözülemedi.")
  call assert_true(positive_solution%mahalanobis%shift > &
    zero_solution%mahalanobis%shift .and. &
    zero_solution%mahalanobis%shift > negative_solution%mahalanobis%shift, &
    "Positive/negative covariance pair shift yönü analytical beklentiye ters.")

  ! Sigma->lambda^2 Sigma objective'i 1/lambda^2 ölçekler; matched ve
  ! Mahalanobis minimizer'ları değişmez.
  scaled_uncertainty = uncertainty
  scaled_covariance = covariance_positive
  call scale_covariance_inputs( &
    scaled_uncertainty, scaled_covariance, 7.0_dp)
  scaled_solution = identify_tts_covariance_pair_shift( &
    family%isotherms(1), family%isotherms(2), scaled_uncertainty, &
    scaled_covariance, configuration)
  call assert_close(scaled_solution%mahalanobis%shift, &
    positive_solution%mahalanobis%shift, 4.0e-8_dp, &
    "Covariance scale Mahalanobis minimizer'ını değiştirdi.")
  call assert_close(scaled_solution%diagonal_matched%shift, &
    positive_solution%diagonal_matched%shift, 4.0e-8_dp, &
    "Covariance scale matched-diagonal minimizer'ını değiştirdi.")
  evaluation = evaluate_tts_covariance_pair_objective( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_positive, 0.4_dp)
  scaled_evaluation = evaluate_tts_covariance_pair_objective( &
    family%isotherms(1), family%isotherms(2), scaled_uncertainty, &
    scaled_covariance, 0.4_dp)
  call assert_close(scaled_evaluation%mahalanobis_objective, &
    evaluation%mahalanobis_objective/49.0_dp, 3.0e-10_dp, &
    "Mahalanobis objective covariance scale ile 1/lambda^2 değişmedi.")
  call assert_close(scaled_evaluation%matched_diagonal_objective, &
    evaluation%matched_diagonal_objective/49.0_dp, 3.0e-10_dp, &
    "Matched objective covariance scale ile 1/lambda^2 değişmedi.")

  ! Missing covariance gerçek gap'tir; kalan iki contiguous parça objective
  ! kurabilir fakat gap üzerinden interpolation yapılmadığı için interval
  ! sayısı ve overlap full support'tan küçülür.
  covariance_positive%isotherms(1)%points(4)%covariance_available = .false.
  covariance_positive%isotherms(2)%points(4)%covariance_available = .false.
  evaluation = evaluate_tts_covariance_pair_objective( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_positive, 0.4_dp)
  call assert_true(evaluation%valid .and. &
    evaluation%overlap_width_decades < 4.0_dp, &
    "Covariance gap bivariate common support'u daraltmadı.")

  evaluation = evaluate_tts_covariance_pair_objective( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_zero, 100.0_dp)
  call assert_true(.not. evaluation%valid .and. &
    evaluation%overlap_width_decades <= 0.0_dp, &
    "No-overlap covariance objective false support üretti.")

  ! Flat objective genuine strict interior bracket taşımaz; boundary/flat
  ! çözüm başarıya çevrilmez.
  family = make_exact_trs_family( &
    temperatures, [0.0_dp, truth_shift], 0.0_dp, 0.0_dp)
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  covariance_zero = make_covariance_family(family, uncertainty, 0.0_dp)
  zero_solution = identify_tts_covariance_pair_shift( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_zero, configuration)
  call assert_true(zero_solution%mahalanobis%status == &
    MAHALANOBIS_SHIFT_NO_INTERIOR_MINIMUM .and. &
    .not. zero_solution%mahalanobis%shift_available, &
    "Flat Mahalanobis objective false interior success üretti.")

  ! Loss uncertainty/covariance tamamen unavailable ise bivariate O_B yoktur.
  ! V0.8.5 storage-only Mahalanobis fallback üretmemelidir.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  covariance_zero = make_covariance_family(family, uncertainty, 0.0_dp)
  do i = 1, size(uncertainty%isotherms(1)%points)
    uncertainty%isotherms(1)%points(i)%loss_uncertainty_available = .false.
    uncertainty%isotherms(2)%points(i)%loss_uncertainty_available = .false.
    covariance_zero%isotherms(1)%points(i)%covariance_available = .false.
    covariance_zero%isotherms(2)%points(i)%covariance_available = .false.
  end do
  zero_solution = identify_tts_covariance_pair_shift( &
    family%isotherms(1), family%isotherms(2), uncertainty, &
    covariance_zero, configuration)
  call assert_true(zero_solution%mahalanobis%status == &
    MAHALANOBIS_SHIFT_INVALID_COVARIANCE .and. &
    .not. zero_solution%mahalanobis%shift_available .and. &
    .not. zero_solution%diagonal_matched%shift_available, &
    "Loss unavailable iken storage-only Mahalanobis fallback üretildi.")

  print *, "V0.8.5 Mahalanobis pair shift ve failure semantics doğrulandı."
end program test_tts_mahalanobis_pair_shift
