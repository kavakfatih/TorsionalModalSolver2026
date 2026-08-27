module tms_dynamic_stiffness
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    validate_stiffness_matrix, get_stiffness_matrix_values, &
    get_stiffness_matrix_size
  use tms_loss_stiffness_matrix, only : loss_stiffness_matrix_t, &
    validate_loss_stiffness_matrix, get_loss_stiffness_matrix_values, &
    get_loss_stiffness_matrix_size
  use tms_damping_matrix, only : damping_matrix_t, &
    validate_damping_matrix, get_damping_matrix_values, &
    get_damping_matrix_size
  use tms_mass_matrix, only : mass_matrix_t, validate_mass_matrix, &
    get_mass_matrix_values, get_mass_matrix_size
  implicit none
  private

  !> Tek bir pozitif uyarım frekansındaki reduced complex dynamic stiffness
  !! matrisini taşır. Katsayılar [N*m/rad] birimindedir ve private depolama,
  !! solver backend'inden bağımsız tam mantıksal Z matrisini korur.
  type, public :: dynamic_stiffness_matrix_t
    private
    complex(dp), allocatable :: values(:, :)
    real(dp) :: frequency_hz = 0.0_dp
  end type dynamic_stiffness_matrix_t

  public :: build_dynamic_stiffness
  public :: validate_dynamic_stiffness
  public :: get_dynamic_stiffness_values
  public :: get_dynamic_stiffness_size
  public :: get_dynamic_stiffness_frequency_hz

contains

  !> Reduced torsional matrislerden frequency-domain dynamic stiffness üretir.
  !!
  !! Fiziksel açıklama: Harmonic convention
  !! theta(t)=Re{theta_hat*exp(i*omega*t)} altında atalet, depolama rijitliği,
  !! yapısal kayıp rijitliği ve viskoz sönüm aynı denge denkleminde birleşir.
  !! Matematiksel model: Z=K'-omega^2*M+i*(K''+omega*C), omega=2*pi*f.
  !! Girdiler: K' ve K'' [N*m/rad], C [N*m*s/rad], M [kg*m^2] ve sonlu,
  !! pozitif f [Hz]. Çıktı: [N*m/rad] birimli complex symmetric Z matrisi.
  !! Varsayımlar ve geçerlilik: Küçük genlikli lineer, doğrudan full-order ve
  !! frozen-property çözüm kullanılır. Girdi matrisleri değiştirilmez. 0 Hz,
  !! boyut uyuşmazlığı veya sonlu olmayan ara/çıktı değerleri reddedilir.
  pure function build_dynamic_stiffness( &
      storage_stiffness, loss_stiffness, damping, mass, frequency_hz) &
      result(dynamic_stiffness)
    type(stiffness_matrix_t), intent(in) :: storage_stiffness
    type(loss_stiffness_matrix_t), intent(in) :: loss_stiffness
    type(damping_matrix_t), intent(in) :: damping
    type(mass_matrix_t), intent(in) :: mass
    real(dp), intent(in) :: frequency_hz
    type(dynamic_stiffness_matrix_t) :: dynamic_stiffness

    real(dp), allocatable :: damping_values(:, :)
    real(dp), allocatable :: loss_values(:, :)
    real(dp), allocatable :: mass_values(:, :)
    real(dp), allocatable :: storage_values(:, :)
    real(dp) :: angular_frequency
    real(dp) :: angular_frequency_squared
    integer :: matrix_size

    call validate_stiffness_matrix(storage_stiffness)
    call validate_loss_stiffness_matrix(loss_stiffness)
    call validate_damping_matrix(damping)
    call validate_mass_matrix(mass)

    matrix_size = get_stiffness_matrix_size(storage_stiffness)
    if (get_loss_stiffness_matrix_size(loss_stiffness) /= matrix_size .or. &
        get_damping_matrix_size(damping) /= matrix_size .or. &
        get_mass_matrix_size(mass) /= matrix_size) then
      error stop "Dynamic stiffness matrislerinin boyutları aynı olmalıdır."
    end if

    if (.not. ieee_is_finite(frequency_hz) .or. frequency_hz <= 0.0_dp) then
      error stop "Harmonic uyarım frekansı sonlu ve pozitif olmalıdır."
    end if

    angular_frequency = 2.0_dp*pi*frequency_hz
    if (.not. ieee_is_finite(angular_frequency)) then
      error stop "Harmonic açısal frekans sonlu sayı aralığında olmalıdır."
    end if
    angular_frequency_squared = angular_frequency*angular_frequency
    if (.not. ieee_is_finite(angular_frequency_squared)) then
      error stop "Harmonic açısal frekansın karesi sonlu olmalıdır."
    end if

    storage_values = get_stiffness_matrix_values(storage_stiffness)
    loss_values = get_loss_stiffness_matrix_values(loss_stiffness)
    damping_values = get_damping_matrix_values(damping)
    mass_values = get_mass_matrix_values(mass)

    allocate(dynamic_stiffness%values(matrix_size, matrix_size))
    dynamic_stiffness%values = cmplx( &
      storage_values-angular_frequency_squared*mass_values, &
      loss_values+angular_frequency*damping_values, kind=dp)
    dynamic_stiffness%frequency_hz = frequency_hz

    call validate_dynamic_stiffness(dynamic_stiffness)
  end function build_dynamic_stiffness

  !> Dynamic stiffness veri sözleşmesini doğrular.
  !! Matematiksel açıklama: Z kare, sonlu ve complex symmetric olmalıdır;
  !! Z^T=Z kontrolü conjugation kullanılmadan yapılır. Frekans [Hz] sonlu ve
  !! pozitiftir. Hermitian olma koşulu aranmaz. Geçersiz veri reddedilir.
  pure subroutine validate_dynamic_stiffness(dynamic_stiffness)
    type(dynamic_stiffness_matrix_t), intent(in) :: dynamic_stiffness

    integer :: column
    integer :: row

    if (.not. allocated(dynamic_stiffness%values)) then
      error stop "Dynamic stiffness kullanılmadan önce oluşturulmalıdır."
    end if
    if (size(dynamic_stiffness%values, 1) /= &
        size(dynamic_stiffness%values, 2)) then
      error stop "Dynamic stiffness matrisi kare olmalıdır."
    end if
    if (.not. ieee_is_finite(dynamic_stiffness%frequency_hz) .or. &
        dynamic_stiffness%frequency_hz <= 0.0_dp) then
      error stop "Dynamic stiffness frekansı sonlu ve pozitif olmalıdır."
    end if

    do column = 1, size(dynamic_stiffness%values, 2)
      do row = 1, size(dynamic_stiffness%values, 1)
        if (.not. ieee_is_finite(real(dynamic_stiffness%values(row, column), dp)) .or. &
            .not. ieee_is_finite(aimag(dynamic_stiffness%values(row, column)))) then
          error stop "Dynamic stiffness yalnız sonlu complex katsayılar içermelidir."
        end if
        if (abs(dynamic_stiffness%values(row, column) - &
            dynamic_stiffness%values(column, row)) > &
            128.0_dp*epsilon(1.0_dp)*max( &
              1.0_dp, abs(dynamic_stiffness%values(row, column)), &
              abs(dynamic_stiffness%values(column, row)))) then
          error stop "Dynamic stiffness complex symmetric olmalıdır."
        end if
      end do
    end do
  end subroutine validate_dynamic_stiffness

  !> Tam mantıksal dynamic stiffness katsayılarının bağımsız kopyasını verir.
  !! Çıktı [N*m/rad] birimli complex(dp) kare matristir; private storage bu
  !! kopya üzerinden değiştirilemez.
  pure function get_dynamic_stiffness_values(dynamic_stiffness) result(values)
    type(dynamic_stiffness_matrix_t), intent(in) :: dynamic_stiffness
    complex(dp), allocatable :: values(:, :)

    call validate_dynamic_stiffness(dynamic_stiffness)
    values = dynamic_stiffness%values
  end function get_dynamic_stiffness_values

  !> Dynamic stiffness matrisinin boyutsuz denklem sayısını döndürür.
  pure function get_dynamic_stiffness_size(dynamic_stiffness) result(matrix_size)
    type(dynamic_stiffness_matrix_t), intent(in) :: dynamic_stiffness
    integer :: matrix_size

    call validate_dynamic_stiffness(dynamic_stiffness)
    matrix_size = size(dynamic_stiffness%values, 1)
  end function get_dynamic_stiffness_size

  !> Dynamic stiffness'in oluşturulduğu uyarım frekansını [Hz] döndürür.
  pure function get_dynamic_stiffness_frequency_hz(dynamic_stiffness) &
      result(frequency_hz)
    type(dynamic_stiffness_matrix_t), intent(in) :: dynamic_stiffness
    real(dp) :: frequency_hz

    call validate_dynamic_stiffness(dynamic_stiffness)
    frequency_hz = dynamic_stiffness%frequency_hz
  end function get_dynamic_stiffness_frequency_hz

end module tms_dynamic_stiffness
