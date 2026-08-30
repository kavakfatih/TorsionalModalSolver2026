module tms_tts_covariance_propagation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_tts_covariance_types, only : &
    tts_covariance_matrix_2x2_t, tts_covariance_matrix_validation_t, &
    tts_polar_covariance_propagation_result_t, &
    tts_log_covariance_propagation_result_t, TTS_COVARIANCE_SUCCESS, &
    TTS_COVARIANCE_INVALID_INPUT, TTS_COVARIANCE_INVALID_MATRIX, &
    TTS_COVARIANCE_SINGULAR_MATRIX, TTS_COVARIANCE_ILL_CONDITIONED, &
    TTS_COVARIANCE_NUMERICAL_FAILURE
  implicit none
  private

  public :: validate_covariance_matrix_2x2
  public :: propagate_magnitude_phase_covariance
  public :: propagate_log10_modulus_covariance
  public :: squared_mahalanobis_2x2

contains

  !> Sigma=[[v_s,c],[c,v_l]] matrisinin SPD ve numerical conditioning
  !! durumunu değerlendirir. Girdiler aynı koordinat sisteminde covariance
  !! birimindedir; test ölçekten bağımsız normalize matris üzerinde yapılır.
  !! Varsayım: matris simetrisi üç bağımsız alanla yapısal olarak sağlanır.
  pure function validate_covariance_matrix_2x2(matrix) result(validation)
    type(tts_covariance_matrix_2x2_t), intent(in) :: matrix
    type(tts_covariance_matrix_validation_t) :: validation

    real(dp) :: a
    real(dp) :: b
    real(dp) :: c
    real(dp) :: determinant_normalized
    real(dp) :: eigen_discriminant
    real(dp) :: eigenvalue_maximum
    real(dp) :: eigenvalue_minimum
    real(dp) :: scale
    real(dp) :: singular_tolerance

    if (.not. ieee_is_finite(matrix%storage_variance_pa2) .or. &
        .not. ieee_is_finite(matrix%loss_variance_pa2) .or. &
        .not. ieee_is_finite(matrix%storage_loss_covariance_pa2)) return
    if (matrix%storage_variance_pa2 <= 0.0_dp .or. &
        matrix%loss_variance_pa2 <= 0.0_dp) return

    scale = max(matrix%storage_variance_pa2, &
      matrix%loss_variance_pa2, &
      abs(matrix%storage_loss_covariance_pa2))
    if (.not. ieee_is_finite(scale) .or. scale <= 0.0_dp) return
    a = matrix%storage_variance_pa2/scale
    b = matrix%loss_variance_pa2/scale
    c = matrix%storage_loss_covariance_pa2/scale
    determinant_normalized = a*b - c*c
    singular_tolerance = 32.0_dp*epsilon(1.0_dp) * &
      max(a*b, c*c, tiny(1.0_dp))
    validation%correlation = c/sqrt(a*b)

    if (determinant_normalized < -singular_tolerance .or. &
        abs(validation%correlation) > 1.0_dp + &
          32.0_dp*epsilon(1.0_dp)) then
      validation%status = TTS_COVARIANCE_INVALID_MATRIX
      return
    end if
    if (determinant_normalized <= singular_tolerance .or. &
        abs(validation%correlation) >= 1.0_dp) then
      validation%status = TTS_COVARIANCE_SINGULAR_MATRIX
      return
    end if

    eigen_discriminant = sqrt((a - b)*(a - b) + 4.0_dp*c*c)
    eigenvalue_maximum = 0.5_dp*(a + b + eigen_discriminant)
    eigenvalue_minimum = determinant_normalized/eigenvalue_maximum
    validation%reciprocal_condition_estimate = &
      eigenvalue_minimum/eigenvalue_maximum
    validation%covariance_valid = .true.
    if (scale <= sqrt(huge(1.0_dp)) .and. &
        determinant_normalized <= huge(1.0_dp)/(scale*scale)) then
      validation%determinant = determinant_normalized*scale*scale
    else
      validation%determinant = huge(1.0_dp)
    end if
    if (validation%reciprocal_condition_estimate <= sqrt(epsilon(1.0_dp))) then
      validation%status = TTS_COVARIANCE_ILL_CONDITIONED
      return
    end if
    validation%status = TTS_COVARIANCE_SUCCESS
    validation%covariance_numerically_well_conditioned = .true.
  end function validate_covariance_matrix_2x2

  !> G*=M*exp(i*delta) modelinin Jacobian'ı ile Sigma_G=J*Sigma_Mdelta*J^T
  !! hesaplar. M [Pa], delta [rad], Var(M) [Pa^2], Var(delta) [rad^2] ve
  !! Cov(M,delta) [Pa rad]'dır. Model first-order'dır; probability dönüşümü
  !! exact kabul edilmez. M>0 ve input covariance SPD olmalıdır.
  pure function propagate_magnitude_phase_covariance( &
      magnitude_pa, phase_angle_rad, magnitude_variance_pa2, &
      phase_variance_rad2, magnitude_phase_covariance_pa_rad) result(result)
    real(dp), intent(in) :: magnitude_pa
    real(dp), intent(in) :: phase_angle_rad
    real(dp), intent(in) :: magnitude_variance_pa2
    real(dp), intent(in) :: phase_variance_rad2
    real(dp), intent(in) :: magnitude_phase_covariance_pa_rad
    type(tts_polar_covariance_propagation_result_t) :: result

    real(dp) :: cosine
    real(dp) :: correlation
    real(dp) :: sine
    type(tts_covariance_matrix_validation_t) :: validation

    if (.not. ieee_is_finite(magnitude_pa) .or. &
        .not. ieee_is_finite(phase_angle_rad) .or. &
        .not. ieee_is_finite(magnitude_variance_pa2) .or. &
        .not. ieee_is_finite(phase_variance_rad2) .or. &
        .not. ieee_is_finite(magnitude_phase_covariance_pa_rad)) return
    if (magnitude_pa <= 0.0_dp .or. magnitude_variance_pa2 <= 0.0_dp .or. &
        phase_variance_rad2 <= 0.0_dp) return
    correlation = (magnitude_phase_covariance_pa_rad / &
      sqrt(magnitude_variance_pa2))/sqrt(phase_variance_rad2)
    if (.not. ieee_is_finite(correlation) .or. &
        1.0_dp - correlation*correlation <= &
          sqrt(epsilon(1.0_dp))) then
      result%status = TTS_COVARIANCE_ILL_CONDITIONED
      return
    end if

    cosine = cos(phase_angle_rad)
    sine = sin(phase_angle_rad)
    result%storage_modulus_pa = magnitude_pa*cosine
    result%loss_modulus_pa = magnitude_pa*sine
    result%covariance%storage_variance_pa2 = &
      cosine*cosine*magnitude_variance_pa2 + &
      magnitude_pa*magnitude_pa*sine*sine*phase_variance_rad2 - &
      2.0_dp*magnitude_pa*cosine*sine* &
        magnitude_phase_covariance_pa_rad
    result%covariance%loss_variance_pa2 = &
      sine*sine*magnitude_variance_pa2 + &
      magnitude_pa*magnitude_pa*cosine*cosine*phase_variance_rad2 + &
      2.0_dp*magnitude_pa*sine*cosine* &
        magnitude_phase_covariance_pa_rad
    result%covariance%storage_loss_covariance_pa2 = &
      sine*cosine*magnitude_variance_pa2 + &
      magnitude_pa*(cosine*cosine - sine*sine)* &
        magnitude_phase_covariance_pa_rad - &
      magnitude_pa*magnitude_pa*sine*cosine*phase_variance_rad2
    if (.not. ieee_is_finite(result%storage_modulus_pa) .or. &
        .not. ieee_is_finite(result%loss_modulus_pa) .or. &
        .not. ieee_is_finite(result%covariance%storage_variance_pa2) .or. &
        .not. ieee_is_finite(result%covariance%loss_variance_pa2) .or. &
        .not. ieee_is_finite( &
          result%covariance%storage_loss_covariance_pa2)) then
      result%status = TTS_COVARIANCE_NUMERICAL_FAILURE
      return
    end if
    validation = validate_covariance_matrix_2x2(result%covariance)
    if (.not. validation%covariance_numerically_well_conditioned) then
      result%status = validation%status
      return
    end if
    result%status = TTS_COVARIANCE_SUCCESS
    result%valid = .true.
  end function propagate_magnitude_phase_covariance

  !> y=[log10(G'),log10(G'')] için D=diag(1/(G ln(10))) ve
  !! Sigma_y=D*Sigma_G*D^T first-order propagation'ını hesaplar. G değerleri
  !! [Pa], Sigma_G [Pa^2], output covariance boyutsuzdur. G',G''>0 şarttır;
  !! zero loss modulus epsilon ile değiştirilmez.
  pure function propagate_log10_modulus_covariance( &
      storage_modulus_pa, loss_modulus_pa, covariance) result(result)
    real(dp), intent(in) :: storage_modulus_pa
    real(dp), intent(in) :: loss_modulus_pa
    type(tts_covariance_matrix_2x2_t), intent(in) :: covariance
    type(tts_log_covariance_propagation_result_t) :: result

    real(dp) :: storage_sensitivity
    real(dp) :: loss_sensitivity
    type(tts_covariance_matrix_2x2_t) :: log_matrix
    type(tts_covariance_matrix_validation_t) :: input_validation
    type(tts_covariance_matrix_validation_t) :: output_validation

    if (.not. ieee_is_finite(storage_modulus_pa) .or. &
        .not. ieee_is_finite(loss_modulus_pa) .or. &
        storage_modulus_pa <= 0.0_dp .or. loss_modulus_pa <= 0.0_dp) return
    input_validation = validate_covariance_matrix_2x2(covariance)
    if (.not. input_validation%covariance_numerically_well_conditioned) then
      result%status = input_validation%status
      return
    end if

    storage_sensitivity = 1.0_dp/(storage_modulus_pa*log(10.0_dp))
    loss_sensitivity = 1.0_dp/(loss_modulus_pa*log(10.0_dp))
    log_matrix%storage_variance_pa2 = &
      covariance%storage_variance_pa2*storage_sensitivity**2
    log_matrix%loss_variance_pa2 = &
      covariance%loss_variance_pa2*loss_sensitivity**2
    log_matrix%storage_loss_covariance_pa2 = &
      covariance%storage_loss_covariance_pa2 * &
      storage_sensitivity*loss_sensitivity
    if (.not. ieee_is_finite(log_matrix%storage_variance_pa2) .or. &
        .not. ieee_is_finite(log_matrix%loss_variance_pa2) .or. &
        .not. ieee_is_finite( &
          log_matrix%storage_loss_covariance_pa2)) then
      result%status = TTS_COVARIANCE_NUMERICAL_FAILURE
      return
    end if
    output_validation = validate_covariance_matrix_2x2(log_matrix)
    if (.not. output_validation%covariance_numerically_well_conditioned) then
      result%status = output_validation%status
      return
    end if
    result%storage_variance = log_matrix%storage_variance_pa2
    result%loss_variance = log_matrix%loss_variance_pa2
    result%storage_loss_covariance = &
      log_matrix%storage_loss_covariance_pa2
    result%correlation = output_validation%correlation
    result%status = TTS_COVARIANCE_SUCCESS
    result%valid = .true.
  end function propagate_log10_modulus_covariance

  !> Point diagnostic için 2x2 Cholesky çözümüyle d_M^2=r^T Sigma^-1 r
  !! hesaplar. Residual bileşenleri covariance koordinatıyla aynı birimdedir;
  !! sonuç boyutsuzdur. SPD/conditioning geçmezse valid=.false. döner.
  pure subroutine squared_mahalanobis_2x2( &
      storage_residual, loss_residual, covariance, value, valid)
    real(dp), intent(in) :: storage_residual
    real(dp), intent(in) :: loss_residual
    type(tts_covariance_matrix_2x2_t), intent(in) :: covariance
    real(dp), intent(out) :: value
    logical, intent(out) :: valid

    real(dp) :: l11
    real(dp) :: l21
    real(dp) :: l22
    real(dp) :: z1
    real(dp) :: z2
    type(tts_covariance_matrix_validation_t) :: validation

    value = 0.0_dp
    valid = .false.
    if (.not. ieee_is_finite(storage_residual) .or. &
        .not. ieee_is_finite(loss_residual)) return
    validation = validate_covariance_matrix_2x2(covariance)
    if (.not. validation%covariance_numerically_well_conditioned) return
    l11 = sqrt(covariance%storage_variance_pa2)
    l21 = covariance%storage_loss_covariance_pa2/l11
    l22 = sqrt(covariance%loss_variance_pa2 - l21*l21)
    if (.not. ieee_is_finite(l22) .or. l22 <= 0.0_dp) return
    z1 = storage_residual/l11
    z2 = (loss_residual - l21*z1)/l22
    value = z1*z1 + z2*z2
    valid = ieee_is_finite(value) .and. value >= 0.0_dp
  end subroutine squared_mahalanobis_2x2

end module tms_tts_covariance_propagation
