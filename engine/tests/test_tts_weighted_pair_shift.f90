program test_tts_weighted_pair_shift
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, &
    tts_pair_shift_configuration_t, tts_pair_shift_result_t, &
    TTS_CHANNEL_STORAGE, TTS_CHANNEL_LOSS, TTS_CHANNEL_JOINT
  use tms_tts_pair_shift, only : identify_tts_pair_shift
  use tms_tts_uncertainty_types, only : &
    tts_dynamic_modulus_uncertainty_family_t, &
    tts_uncertainty_pair_solution_t, &
    tts_uncertainty_objective_evaluation_t, TTS_WEIGHTED_L2_OBJECTIVE, &
    TTS_STANDARDIZED_HUBER_OBJECTIVE, WEIGHTED_SHIFT_SUCCESS, &
    WEIGHTED_SHIFT_STORAGE_ONLY, WEIGHTED_SHIFT_INSUFFICIENT_SUPPORT, &
    WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM
  use tms_tts_weighted_pair_shift, only : &
    evaluate_tts_uncertainty_pair_objective, &
    identify_tts_uncertainty_pair_shift
  use tms_tts_test_support, only : make_exact_trs_family, &
    make_generalized_maxwell_trs_family, assert_true, assert_close
  use tms_tts_uncertainty_test_support, only : &
    make_relative_uncertainty_family, scale_uncertainty_family
  implicit none

  type(tts_material_family_t) :: family
  type(tts_dynamic_modulus_uncertainty_family_t) :: uncertainty
  type(tts_dynamic_modulus_uncertainty_family_t) :: scaled_uncertainty
  type(tts_pair_shift_configuration_t) :: configuration
  type(tts_pair_shift_result_t) :: baseline
  type(tts_uncertainty_pair_solution_t) :: weighted
  type(tts_uncertainty_pair_solution_t) :: scaled_weighted
  type(tts_uncertainty_pair_solution_t) :: huber
  type(tts_uncertainty_objective_evaluation_t) :: storage_objective
  type(tts_uncertainty_objective_evaluation_t) :: loss_objective
  type(tts_uncertainty_objective_evaluation_t) :: joint_objective
  real(dp), parameter :: temperatures(2) = [293.15_dp, 313.15_dp]
  real(dp), parameter :: truth_shift = 0.65_dp
  integer :: i

  ! Sabit relative uncertainty, bütün noktalarda aynı u_log10 ve residual
  ! variance üretir. Bu durumda weighted-L2 objective yalnız sabit ölçeklenir;
  ! minimizer authoritative unweighted baseline ile aynı kalmalıdır.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  baseline = identify_tts_pair_shift( &
    family%isotherms(1), family%isotherms(2), configuration)
  weighted = identify_tts_uncertainty_pair_shift(family%isotherms(1), &
    family%isotherms(2), uncertainty, TTS_WEIGHTED_L2_OBJECTIVE, &
    configuration, 1.345_dp)
  call assert_true(weighted%status == WEIGHTED_SHIFT_SUCCESS .and. &
    weighted%shift_available, "Homoscedastic weighted pair çözülemedi.")
  call assert_close(weighted%shift, baseline%delta_s, 2.0e-6_dp, &
    "Homoscedastic weighted minimizer baseline'dan ayrıldı.")

  ! Bütün u_G değerlerinin ortak lambda ile ölçeklenmesi weighted objective'i
  ! 1/lambda^2 değiştirir fakat argmin shift'i değiştirmez.
  scaled_uncertainty = uncertainty
  call scale_uncertainty_family(scaled_uncertainty, 7.0_dp)
  scaled_weighted = identify_tts_uncertainty_pair_shift( &
    family%isotherms(1), family%isotherms(2), scaled_uncertainty, &
    TTS_WEIGHTED_L2_OBJECTIVE, configuration, 1.345_dp)
  call assert_close(scaled_weighted%shift, weighted%shift, 3.0e-7_dp, &
    "Global uncertainty ölçeği weighted shift minimizer'ını değiştirdi.")
  call assert_close(scaled_weighted%objective_minimum, &
    weighted%objective_minimum/49.0_dp, 2.0e-8_dp, &
    "Weighted objective ortak uncertainty ölçeğinde 1/lambda^2 değişmedi.")

  ! Heteroscedastic fixture: high-frequency moving bölgesine deterministic
  ! bias eklenir ve aynı bölgenin uncertainty'si büyütülür. Inverse-variance
  ! semantics weighted shift'i known truth'e unweighted L2'den yaklaştırmalıdır.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  do i = 27, 35
    family%isotherms(2)%points(i)%storage_modulus_pa = &
      family%isotherms(2)%points(i)%storage_modulus_pa*10.0_dp**0.40_dp
    family%isotherms(2)%points(i)%loss_modulus_pa = &
      family%isotherms(2)%points(i)%loss_modulus_pa*10.0_dp**0.40_dp
  end do
  baseline = identify_tts_pair_shift( &
    family%isotherms(1), family%isotherms(2), configuration)
  uncertainty = make_relative_uncertainty_family(family, 0.01_dp, 0.01_dp)
  do i = 27, 35
    uncertainty%isotherms(2)%points(i) &
      %storage_standard_uncertainty_pa = &
      0.50_dp*family%isotherms(2)%points(i)%storage_modulus_pa
    uncertainty%isotherms(2)%points(i)%loss_standard_uncertainty_pa = &
      0.50_dp*family%isotherms(2)%points(i)%loss_modulus_pa
  end do
  weighted = identify_tts_uncertainty_pair_shift(family%isotherms(1), &
    family%isotherms(2), uncertainty, TTS_WEIGHTED_L2_OBJECTIVE, &
    configuration, 1.345_dp)
  call assert_true(baseline%shift_available .and. weighted%shift_available, &
    "Heteroscedastic fixture shift sonuçlarını üretmedi.")
  call assert_true(abs(weighted%shift - truth_shift) < &
    abs(baseline%delta_s - truth_shift), &
    "Weighted shift biased region'i inverse-variance ile azaltmadı.")

  ! Storage channel çok daha fazla sample point taşısa da joint objective ham
  ! point count ile değil 0.5*storage+0.5*loss olarak kurulmalıdır.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  do i = 1, size(family%isotherms(2)%points)
    family%isotherms(2)%points(i)%loss_modulus_pa = &
      family%isotherms(2)%points(i)%loss_modulus_pa*10.0_dp**0.08_dp
  end do
  do i = 13, size(uncertainty%isotherms(1)%points)
    uncertainty%isotherms(1)%points(i)%loss_uncertainty_available = .false.
    uncertainty%isotherms(2)%points(i)%loss_uncertainty_available = .false.
  end do
  storage_objective = evaluate_tts_uncertainty_pair_objective( &
    family%isotherms(1), family%isotherms(2), uncertainty, truth_shift, &
    TTS_CHANNEL_STORAGE, TTS_WEIGHTED_L2_OBJECTIVE, 1.345_dp)
  loss_objective = evaluate_tts_uncertainty_pair_objective( &
    family%isotherms(1), family%isotherms(2), uncertainty, truth_shift, &
    TTS_CHANNEL_LOSS, TTS_WEIGHTED_L2_OBJECTIVE, 1.345_dp)
  joint_objective = evaluate_tts_uncertainty_pair_objective( &
    family%isotherms(1), family%isotherms(2), uncertainty, truth_shift, &
    TTS_CHANNEL_JOINT, TTS_WEIGHTED_L2_OBJECTIVE, 1.345_dp)
  call assert_true(storage_objective%valid .and. loss_objective%valid .and. &
    joint_objective%valid, "Channel-balance objective support'u eksik.")
  call assert_close(joint_objective%objective, 0.5_dp * &
    (storage_objective%objective + loss_objective%objective), 2.0e-13_dp, &
    "Joint weighted objective equal channel weighting'i korumadı.")

  ! Loss uncertainty tamamen unavailable ise G'' fabricate edilmez; storage
  ! weighted ve Huber çözümleri açık STORAGE_ONLY status ile kalır.
  family = make_generalized_maxwell_trs_family( &
    temperatures, [0.0_dp, truth_shift])
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  do i = 1, size(uncertainty%isotherms(1)%points)
    uncertainty%isotherms(1)%points(i)%loss_uncertainty_available = .false.
    uncertainty%isotherms(2)%points(i)%loss_uncertainty_available = .false.
  end do
  weighted = identify_tts_uncertainty_pair_shift(family%isotherms(1), &
    family%isotherms(2), uncertainty, TTS_WEIGHTED_L2_OBJECTIVE, &
    configuration, 1.345_dp)
  huber = identify_tts_uncertainty_pair_shift(family%isotherms(1), &
    family%isotherms(2), uncertainty, TTS_STANDARDIZED_HUBER_OBJECTIVE, &
    configuration, 1.345_dp)
  call assert_true(weighted%status == WEIGHTED_SHIFT_STORAGE_ONLY .and. &
    weighted%shift_available .and. .not. weighted%joint_shift_available, &
    "Weighted storage-only fallback explicit değil.")
  call assert_true(huber%status == WEIGHTED_SHIFT_STORAGE_ONLY .and. &
    huber%shift_available .and. .not. huber%joint_shift_available, &
    "Huber storage-only fallback explicit değil.")

  ! Tek uncertainty noktası interpolation interval'i kuramaz; clean
  ! insufficient-support status gerekir.
  do i = 1, size(uncertainty%isotherms(1)%points)
    uncertainty%isotherms(1)%points(i)%storage_uncertainty_available = &
      i == 1
    uncertainty%isotherms(2)%points(i)%storage_uncertainty_available = &
      i == 1
  end do
  weighted = identify_tts_uncertainty_pair_shift(family%isotherms(1), &
    family%isotherms(2), uncertainty, TTS_WEIGHTED_L2_OBJECTIVE, &
    configuration, 1.345_dp)
  call assert_true(weighted%status == WEIGHTED_SHIFT_INSUFFICIENT_SUPPORT .and. &
    .not. weighted%shift_available, &
    "Tek-nokta uncertainty support'u clean failure vermedi.")

  ! Flat ve boundary-only objectives gerçek strict interior bracket taşımaz;
  ! Brent başarısı gibi raporlanmamalıdır.
  family = make_exact_trs_family( &
    temperatures, [0.0_dp, truth_shift], 0.0_dp, 0.0_dp)
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  weighted = identify_tts_uncertainty_pair_shift(family%isotherms(1), &
    family%isotherms(2), uncertainty, TTS_WEIGHTED_L2_OBJECTIVE, &
    configuration, 1.345_dp)
  call assert_true(weighted%status == WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM, &
    "Flat weighted objective false interior success üretti.")
  family = make_exact_trs_family(temperatures, [0.0_dp, 0.0_dp])
  do i = 1, size(family%isotherms(2)%points)
    family%isotherms(2)%points(i)%storage_modulus_pa = &
      100.0_dp*family%isotherms(2)%points(i)%storage_modulus_pa
    family%isotherms(2)%points(i)%loss_modulus_pa = &
      100.0_dp*family%isotherms(2)%points(i)%loss_modulus_pa
  end do
  uncertainty = make_relative_uncertainty_family(family, 0.02_dp, 0.02_dp)
  weighted = identify_tts_uncertainty_pair_shift(family%isotherms(1), &
    family%isotherms(2), uncertainty, TTS_WEIGHTED_L2_OBJECTIVE, &
    configuration, 1.345_dp)
  call assert_true(weighted%status == WEIGHTED_SHIFT_NO_INTERIOR_MINIMUM, &
    "Boundary-only weighted optimum false success üretti.")
  joint_objective = evaluate_tts_uncertainty_pair_objective( &
    family%isotherms(1), family%isotherms(2), uncertainty, 100.0_dp, &
    TTS_CHANNEL_JOINT, TTS_WEIGHTED_L2_OBJECTIVE, 1.345_dp)
  call assert_true(.not. joint_objective%valid .and. &
    joint_objective%overlap_width_decades <= 0.0_dp, &
    "No-overlap weighted objective false support üretti.")

  print *, "V0.8.4 weighted pair shift invariants ve failures doğrulandı."
end program test_tts_weighted_pair_shift
