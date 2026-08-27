module tms_frf
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : pi
  implicit none
  private

  public :: calculate_rotational_receptance
  public :: calculate_rotational_mobility
  public :: calculate_rotational_accelerance

contains

  !> Tanımlı tek torque input/output kanalı için rotational receptance hesaplar.
  !!
  !! Fiziksel açıklama: Receptance, output angular-displacement genliğinin aynı
  !! load case içindeki tanımlı input torque genliğine transfer oranıdır.
  !! Matematiksel model: H_thetaT=theta_hat/T_hat.
  !! Girdiler: Complex peak theta_hat [rad] ve sıfır olmayan complex peak
  !! T_hat [N*m]. Çıktı: H_thetaT [rad/(N*m)].
  !! Varsayımlar ve geçerlilik: Tek input channel açıkça tanımlanmıştır;
  !! eşzamanlı bağımsız yüklerin toplamına bölünen genel response FRF sayılmaz.
  pure function calculate_rotational_receptance( &
      angular_response, input_torque) result(receptance)
    complex(dp), intent(in) :: angular_response
    complex(dp), intent(in) :: input_torque
    complex(dp) :: receptance

    call validate_complex_value(angular_response, "Angular response")
    call validate_nonzero_input_torque(input_torque)
    receptance = angular_response/input_torque
    call validate_complex_value(receptance, "Rotational receptance")
  end function calculate_rotational_receptance

  !> Tanımlı torque input/output kanalı için rotational mobility hesaplar.
  !!
  !! Fiziksel açıklama: Mobility, complex angular velocity'nin input torque'a
  !! transfer oranıdır. exp(+i*omega*t) convention'ında türev +i*omega'dır.
  !! Matematiksel model: H_omegaT=i*omega*theta_hat/T_hat.
  !! Girdiler: theta_hat [rad], T_hat [N*m], sonlu pozitif f [Hz]. Çıktı:
  !! H_omegaT [rad/(s*N*m)]. Küçük genlikli steady-state hareket kabul edilir.
  pure function calculate_rotational_mobility( &
      angular_response, input_torque, frequency_hz) result(mobility)
    complex(dp), intent(in) :: angular_response
    complex(dp), intent(in) :: input_torque
    real(dp), intent(in) :: frequency_hz
    complex(dp) :: mobility

    real(dp) :: angular_frequency

    call validate_frequency(frequency_hz)
    angular_frequency = 2.0_dp*pi*frequency_hz
    if (.not. ieee_is_finite(angular_frequency)) then
      error stop "FRF açısal frekansı sonlu olmalıdır."
    end if
    mobility = cmplx(0.0_dp, angular_frequency, kind=dp) * &
      calculate_rotational_receptance(angular_response, input_torque)
    call validate_complex_value(mobility, "Rotational mobility")
  end function calculate_rotational_mobility

  !> Tanımlı torque input/output kanalı için rotational accelerance hesaplar.
  !!
  !! Fiziksel açıklama: Accelerance, complex angular acceleration'ın input
  !! torque'a transfer oranıdır. exp(+i*omega*t) convention'ında ikinci türev
  !! -omega^2 katsayısı üretir.
  !! Matematiksel model: H_alphaT=-omega^2*theta_hat/T_hat.
  !! Girdiler theta_hat [rad], T_hat [N*m], sonlu pozitif f [Hz]; çıktı
  !! H_alphaT [rad/(s^2*N*m)] değeridir.
  pure function calculate_rotational_accelerance( &
      angular_response, input_torque, frequency_hz) result(accelerance)
    complex(dp), intent(in) :: angular_response
    complex(dp), intent(in) :: input_torque
    real(dp), intent(in) :: frequency_hz
    complex(dp) :: accelerance

    real(dp) :: angular_frequency
    real(dp) :: angular_frequency_squared

    call validate_frequency(frequency_hz)
    angular_frequency = 2.0_dp*pi*frequency_hz
    angular_frequency_squared = angular_frequency*angular_frequency
    if (.not. ieee_is_finite(angular_frequency_squared)) then
      error stop "FRF açısal frekans karesi sonlu olmalıdır."
    end if
    accelerance = -angular_frequency_squared * &
      calculate_rotational_receptance(angular_response, input_torque)
    call validate_complex_value(accelerance, "Rotational accelerance")
  end function calculate_rotational_accelerance

  !> FRF frekansının [Hz] sonlu ve pozitif olduğunu doğrular.
  pure subroutine validate_frequency(frequency_hz)
    real(dp), intent(in) :: frequency_hz

    if (.not. ieee_is_finite(frequency_hz) .or. frequency_hz <= 0.0_dp) then
      error stop "FRF frekansı sonlu ve pozitif olmalıdır."
    end if
  end subroutine validate_frequency

  !> Tek FRF input torque kanalının sonlu ve sıfırdan farklı olduğunu sınar.
  pure subroutine validate_nonzero_input_torque(input_torque)
    complex(dp), intent(in) :: input_torque

    call validate_complex_value(input_torque, "FRF input torque")
    if (.not. abs(input_torque) > 0.0_dp) then
      error stop "FRF input torque genliği sıfır olamaz."
    end if
  end subroutine validate_nonzero_input_torque

  !> Bir complex büyüklüğün gerçek ve sanal bileşenlerini sonluluk açısından
  !! doğrular. Fiziksel birim çağıran yordamın sözleşmesiyle belirlenir.
  pure subroutine validate_complex_value(value, quantity_name)
    complex(dp), intent(in) :: value
    character(len=*), intent(in) :: quantity_name

    if (.not. ieee_is_finite(real(value, dp)) .or. &
        .not. ieee_is_finite(aimag(value))) then
      error stop quantity_name//" sonlu complex değer olmalıdır."
    end if
  end subroutine validate_complex_value

end module tms_frf
