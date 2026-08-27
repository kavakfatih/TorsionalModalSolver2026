module tms_harmonic_response
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_dof_types, only : TORSIONAL_ROTATION, is_supported_dof_type
  use tms_torsional_element, only : torsional_element_t, &
    validate_torsional_element
  use tms_constraint_manager, only : active_dof_map_t, &
    active_dof_map_entry_t, validate_active_dof_map, &
    get_active_dof_map_entries
  use tms_complex_linear_solution, only : COMPLEX_SOLVE_SOLVED, &
    COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED, COMPLEX_SOLVE_SINGULAR
  implicit none
  private

  !> Bir explicit frequency sweep'in backend-neutral complex torsional
  !! response ve sayısal tanılarını taşır. Complex genlikler exp(+i*omega*t)
  !! convention'ında PEAK değerlerdir; storage private tutulur.
  type, public :: harmonic_response_t
    private
    real(dp), allocatable :: frequencies_hz(:)
    integer, allocatable :: solution_statuses(:)
    logical, allocatable :: response_available(:)
    real(dp), allocatable :: reciprocal_condition_numbers(:)
    real(dp), allocatable :: relative_residuals(:)
    real(dp), allocatable :: forward_error_bounds(:)
    real(dp), allocatable :: backward_errors(:)
    complex(dp), allocatable :: reduced_complex_response(:, :)
    complex(dp), allocatable :: physical_complex_response(:, :)
    character(len=:), allocatable :: backend_identity
  end type harmonic_response_t

  public :: create_harmonic_response
  public :: validate_harmonic_response
  public :: get_harmonic_frequency_count
  public :: get_harmonic_frequencies_hz
  public :: get_harmonic_solution_statuses
  public :: get_harmonic_response_availability
  public :: get_harmonic_reciprocal_condition_numbers
  public :: get_harmonic_relative_residuals
  public :: get_harmonic_forward_error_bounds
  public :: get_harmonic_backward_errors
  public :: get_reduced_complex_response
  public :: get_physical_complex_response
  public :: get_reduced_response_at_frequency
  public :: get_physical_response_at_frequency
  public :: get_harmonic_backend_identity
  public :: get_physical_response_magnitudes
  public :: get_physical_response_phases_rad
  public :: get_physical_angular_velocities
  public :: get_physical_angular_accelerations
  public :: get_physical_node_response
  public :: calculate_element_relative_angle
  public :: calculate_element_dynamic_torque
  public :: calculate_element_average_dissipated_power
  public :: calculate_element_dissipated_energy_per_cycle
  public :: get_element_relative_angle_response
  public :: get_element_dynamic_torque_response
  public :: get_element_dynamic_torque_magnitudes
  public :: get_element_average_dissipated_power_response
  public :: get_element_dissipated_energy_response
  public :: validate_harmonic_frequency_sweep

contains

  !> Harmonic sweep sonuçlarını doğrulanmış private result nesnesine kopyalar.
  !!
  !! Fiziksel açıklama: Her sütun tek bir pozitif frequency point'in reduced ve
  !! physical angular response genliğini [rad] taşır. Complex değerler peak'tir.
  !! Matematiksel açıklama: Status, RCOND, residual, FERR ve BERR aynı frequency
  !! indeksiyle eşlenir. SINGULAR noktada unique response olmadığından response
  !! ve solution-only tanılar complex/real NaN unavailable sentinel taşır;
  !! bunlar çözüm olarak yorumlanamaz.
  !! Girdiler bağımsız kopyalanır. Çıktı harmonic_response_t değeridir.
  pure function create_harmonic_response( &
      frequencies_hz, solution_statuses, reciprocal_condition_numbers, &
      relative_residuals, forward_error_bounds, backward_errors, &
      reduced_complex_response, physical_complex_response, backend_identity) &
      result(response)
    real(dp), intent(in) :: frequencies_hz(:)
    integer, intent(in) :: solution_statuses(:)
    real(dp), intent(in) :: reciprocal_condition_numbers(:)
    real(dp), intent(in) :: relative_residuals(:)
    real(dp), intent(in) :: forward_error_bounds(:)
    real(dp), intent(in) :: backward_errors(:)
    complex(dp), intent(in) :: reduced_complex_response(:, :)
    complex(dp), intent(in) :: physical_complex_response(:, :)
    character(len=*), intent(in) :: backend_identity
    type(harmonic_response_t) :: response

    integer :: frequency_index

    response%frequencies_hz = frequencies_hz
    response%solution_statuses = solution_statuses
    response%reciprocal_condition_numbers = reciprocal_condition_numbers
    response%relative_residuals = relative_residuals
    response%forward_error_bounds = forward_error_bounds
    response%backward_errors = backward_errors
    response%reduced_complex_response = reduced_complex_response
    response%physical_complex_response = physical_complex_response
    response%backend_identity = trim(backend_identity)
    allocate(response%response_available(size(solution_statuses)))
    do frequency_index = 1, size(solution_statuses)
      response%response_available(frequency_index) = &
        is_solved_status(solution_statuses(frequency_index))
    end do

    call validate_harmonic_response(response)
  end function create_harmonic_response

  !> Harmonic result boyut, status, sonluluk ve unavailable invariantlarını
  !! doğrular. Solved noktalarda tüm response/tanılar sonlu ve negatif olmayan,
  !! singular noktalarda solution-only alanlar NaN olmalıdır. RCOND her durumda
  !! sonlu, boyutsuz ve negatif olmayandır. LAPACK tahmini roundoff nedeniyle
  !! bir değerini çok az aşabileceğinden üst sınır fiziksel invariant sayılmaz.
  !! Geçersiz sonuç reddedilir.
  pure subroutine validate_harmonic_response(response)
    type(harmonic_response_t), intent(in) :: response

    integer :: frequency_count
    integer :: frequency_index

    if (.not. allocated(response%frequencies_hz) .or. &
        .not. allocated(response%solution_statuses) .or. &
        .not. allocated(response%response_available) .or. &
        .not. allocated(response%reciprocal_condition_numbers) .or. &
        .not. allocated(response%relative_residuals) .or. &
        .not. allocated(response%forward_error_bounds) .or. &
        .not. allocated(response%backward_errors) .or. &
        .not. allocated(response%reduced_complex_response) .or. &
        .not. allocated(response%physical_complex_response) .or. &
        .not. allocated(response%backend_identity)) then
      error stop "Harmonic response kullanılmadan önce oluşturulmalıdır."
    end if

    frequency_count = size(response%frequencies_hz)
    if (frequency_count == 0) then
      error stop "Harmonic response en az bir frequency point içermelidir."
    end if
    if (size(response%solution_statuses) /= frequency_count .or. &
        size(response%response_available) /= frequency_count .or. &
        size(response%reciprocal_condition_numbers) /= frequency_count .or. &
        size(response%relative_residuals) /= frequency_count .or. &
        size(response%forward_error_bounds) /= frequency_count .or. &
        size(response%backward_errors) /= frequency_count .or. &
        size(response%reduced_complex_response, 2) /= frequency_count .or. &
        size(response%physical_complex_response, 2) /= frequency_count) then
      error stop "Harmonic response frequency boyutları tutarsız."
    end if
    if (size(response%reduced_complex_response, 1) == 0 .or. &
        size(response%physical_complex_response, 1) == 0) then
      error stop "Harmonic response en az bir aktif ve fiziksel DOF gerektirir."
    end if
    if (len_trim(response%backend_identity) == 0) then
      error stop "Harmonic response solver backend kimliği boş olamaz."
    end if

    call validate_harmonic_frequency_sweep(response%frequencies_hz)
    do frequency_index = 1, frequency_count
      if (.not. is_known_status(response%solution_statuses(frequency_index))) then
        error stop "Harmonic response bilinmeyen solution status içeriyor."
      end if
      if (response%response_available(frequency_index) .neqv. &
          is_solved_status(response%solution_statuses(frequency_index))) then
        error stop "Harmonic response availability ile status tutarsız."
      end if
      if (.not. ieee_is_finite( &
          response%reciprocal_condition_numbers(frequency_index)) .or. &
          response%reciprocal_condition_numbers(frequency_index) < 0.0_dp) then
        error stop "Harmonic response RCOND sonlu ve negatif olmayan olmalıdır."
      end if

      if (response%response_available(frequency_index)) then
        call require_finite_complex_vector( &
          response%reduced_complex_response(:, frequency_index))
        call require_finite_complex_vector( &
          response%physical_complex_response(:, frequency_index))
        call require_finite_nonnegative( &
          response%relative_residuals(frequency_index), "Relative residual")
        call require_finite_nonnegative( &
          response%forward_error_bounds(frequency_index), "FERR")
        call require_finite_nonnegative( &
          response%backward_errors(frequency_index), "BERR")
      else
        call require_unavailable_complex_vector( &
          response%reduced_complex_response(:, frequency_index))
        call require_unavailable_complex_vector( &
          response%physical_complex_response(:, frequency_index))
        if (.not. ieee_is_nan(response%relative_residuals(frequency_index)) .or. &
            .not. ieee_is_nan(response%forward_error_bounds(frequency_index)) .or. &
            .not. ieee_is_nan(response%backward_errors(frequency_index))) then
          error stop "Singular harmonic noktada solution tanıları unavailable olmalıdır."
        end if
      end if
    end do
  end subroutine validate_harmonic_response

  !> Harmonic sweep'teki boyutsuz frequency point sayısını döndürür.
  pure function get_harmonic_frequency_count(response) result(frequency_count)
    type(harmonic_response_t), intent(in) :: response
    integer :: frequency_count

    call validate_harmonic_response(response)
    frequency_count = size(response%frequencies_hz)
  end function get_harmonic_frequency_count

  !> Explicit harmonic frequency grid'inin [Hz] bağımsız kopyasını döndürür.
  pure function get_harmonic_frequencies_hz(response) result(frequencies_hz)
    type(harmonic_response_t), intent(in) :: response
    real(dp), allocatable :: frequencies_hz(:)

    call validate_harmonic_response(response)
    frequencies_hz = response%frequencies_hz
  end function get_harmonic_frequencies_hz

  !> Frequency point başına backend-neutral solution status dizisini döndürür.
  pure function get_harmonic_solution_statuses(response) result(statuses)
    type(harmonic_response_t), intent(in) :: response
    integer, allocatable :: statuses(:)

    call validate_harmonic_response(response)
    statuses = response%solution_statuses
  end function get_harmonic_solution_statuses

  !> Her frequency point'te unique response bulunup bulunmadığını döndürür.
  pure function get_harmonic_response_availability(response) result(availability)
    type(harmonic_response_t), intent(in) :: response
    logical, allocatable :: availability(:)

    call validate_harmonic_response(response)
    availability = response%response_available
  end function get_harmonic_response_availability

  !> Frequency point başına boyutsuz reciprocal condition estimate döndürür.
  pure function get_harmonic_reciprocal_condition_numbers(response) &
      result(values)
    type(harmonic_response_t), intent(in) :: response
    real(dp), allocatable :: values(:)

    call validate_harmonic_response(response)
    values = response%reciprocal_condition_numbers
  end function get_harmonic_reciprocal_condition_numbers

  !> Frequency point başına backend-independent boyutsuz residual döndürür.
  !! Singular noktadaki değer NaN unavailable sentinel'dır.
  pure function get_harmonic_relative_residuals(response) result(values)
    type(harmonic_response_t), intent(in) :: response
    real(dp), allocatable :: values(:)

    call validate_harmonic_response(response)
    values = response%relative_residuals
  end function get_harmonic_relative_residuals

  !> Tek RHS için frequency point başına ZSYSVX forward error bound döndürür.
  !! Singular noktadaki değer NaN unavailable sentinel'dır.
  pure function get_harmonic_forward_error_bounds(response) result(values)
    type(harmonic_response_t), intent(in) :: response
    real(dp), allocatable :: values(:)

    call validate_harmonic_response(response)
    values = response%forward_error_bounds
  end function get_harmonic_forward_error_bounds

  !> Tek RHS için frequency point başına ZSYSVX componentwise backward error
  !! değerini döndürür. Singular noktadaki değer NaN unavailable sentinel'dır.
  pure function get_harmonic_backward_errors(response) result(values)
    type(harmonic_response_t), intent(in) :: response
    real(dp), allocatable :: values(:)

    call validate_harmonic_response(response)
    values = response%backward_errors
  end function get_harmonic_backward_errors

  !> Active-equation sıralı reduced complex response matrisinin [rad]
  !! bağımsız kopyasını döndürür. Sütunlar frequency sırasındadır. Singular
  !! sütunlar çözüm değildir ve NaN unavailable sentinel taşır; status birlikte
  !! incelenmelidir.
  pure function get_reduced_complex_response(response) result(values)
    type(harmonic_response_t), intent(in) :: response
    complex(dp), allocatable :: values(:, :)

    call validate_harmonic_response(response)
    values = response%reduced_complex_response
  end function get_reduced_complex_response

  !> Physical-DOF sıralı complex response matrisinin [rad] bağımsız kopyasını
  !! döndürür. Singular sütunlar NaN unavailable sentinel'dır.
  pure function get_physical_complex_response(response) result(values)
    type(harmonic_response_t), intent(in) :: response
    complex(dp), allocatable :: values(:, :)

    call validate_harmonic_response(response)
    values = response%physical_complex_response
  end function get_physical_complex_response

  !> Tek frequency point'in reduced complex response vektörünü [rad] döndürür.
  !! Singular noktada unique solution bulunmadığından clean diagnostic üretir.
  pure function get_reduced_response_at_frequency(response, frequency_index) &
      result(values)
    type(harmonic_response_t), intent(in) :: response
    integer, intent(in) :: frequency_index
    complex(dp), allocatable :: values(:)

    call require_available_frequency(response, frequency_index)
    values = response%reduced_complex_response(:, frequency_index)
  end function get_reduced_response_at_frequency

  !> Tek frequency point'in physical complex response vektörünü [rad] verir.
  !! Constraint uygulanmış bileşenler sıfır complex phasor'dur.
  pure function get_physical_response_at_frequency(response, frequency_index) &
      result(values)
    type(harmonic_response_t), intent(in) :: response
    integer, intent(in) :: frequency_index
    complex(dp), allocatable :: values(:)

    call require_available_frequency(response, frequency_index)
    values = response%physical_complex_response(:, frequency_index)
  end function get_physical_response_at_frequency

  !> Harmonic çözümün backend kimliğinin bağımsız kopyasını döndürür.
  pure function get_harmonic_backend_identity(response) result(identity)
    type(harmonic_response_t), intent(in) :: response
    character(len=:), allocatable :: identity

    call validate_harmonic_response(response)
    identity = response%backend_identity
  end function get_harmonic_backend_identity

  !> Physical complex angular response genliklerini abs(theta_hat) [rad]
  !! olarak döndürür. Singular sütunlar NaN unavailable kalır.
  pure function get_physical_response_magnitudes(response) result(magnitudes)
    type(harmonic_response_t), intent(in) :: response
    real(dp), allocatable :: magnitudes(:, :)

    integer :: frequency_index

    call validate_harmonic_response(response)
    allocate(magnitudes( &
      size(response%physical_complex_response, 1), &
      size(response%physical_complex_response, 2)))
    do frequency_index = 1, size(magnitudes, 2)
      if (response%response_available(frequency_index)) then
        magnitudes(:, frequency_index) = &
          abs(response%physical_complex_response(:, frequency_index))
      else
        magnitudes(:, frequency_index) = &
          real(response%physical_complex_response(:, frequency_index), dp)
      end if
    end do
  end function get_physical_response_magnitudes

  !> Physical complex angular response fazını atan2(imag,real) ile [rad]
  !! döndürür. Phase unwrapping uygulanmaz; singular sütunlar NaN'dır.
  pure function get_physical_response_phases_rad(response) result(phases_rad)
    type(harmonic_response_t), intent(in) :: response
    real(dp), allocatable :: phases_rad(:, :)

    integer :: column
    integer :: row

    call validate_harmonic_response(response)
    allocate(phases_rad( &
      size(response%physical_complex_response, 1), &
      size(response%physical_complex_response, 2)))
    do column = 1, size(phases_rad, 2)
      do row = 1, size(phases_rad, 1)
        phases_rad(row, column) = atan2( &
          aimag(response%physical_complex_response(row, column)), &
          real(response%physical_complex_response(row, column), dp))
      end do
    end do
  end function get_physical_response_phases_rad

  !> Physical angular velocity phasor'larını [rad/s] döndürür.
  !! exp(+i*omega*t) convention'ında omega_hat=i*omega*theta_hat kullanılır;
  !! yeni solver çağrısı yapılmaz ve singular sütunlar unavailable kalır.
  pure function get_physical_angular_velocities(response) result(velocities)
    type(harmonic_response_t), intent(in) :: response
    complex(dp), allocatable :: velocities(:, :)

    real(dp) :: angular_frequency
    integer :: frequency_index

    call validate_harmonic_response(response)
    allocate(velocities( &
      size(response%physical_complex_response, 1), &
      size(response%physical_complex_response, 2)))
    do frequency_index = 1, size(velocities, 2)
      if (response%response_available(frequency_index)) then
        angular_frequency = 2.0_dp*pi*response%frequencies_hz(frequency_index)
        velocities(:, frequency_index) = &
          cmplx(0.0_dp, angular_frequency, kind=dp) * &
          response%physical_complex_response(:, frequency_index)
        call require_finite_complex_vector(velocities(:, frequency_index))
      else
        velocities(:, frequency_index) = &
          response%physical_complex_response(:, frequency_index)
      end if
    end do
  end function get_physical_angular_velocities

  !> Physical angular acceleration phasor'larını [rad/s^2] döndürür.
  !! exp(+i*omega*t) convention'ında alpha_hat=-omega^2*theta_hat kullanılır;
  !! singular sütunlar unavailable kalır.
  pure function get_physical_angular_accelerations(response) &
      result(accelerations)
    type(harmonic_response_t), intent(in) :: response
    complex(dp), allocatable :: accelerations(:, :)

    real(dp) :: angular_frequency
    real(dp) :: angular_frequency_squared
    integer :: frequency_index

    call validate_harmonic_response(response)
    allocate(accelerations( &
      size(response%physical_complex_response, 1), &
      size(response%physical_complex_response, 2)))
    do frequency_index = 1, size(accelerations, 2)
      if (response%response_available(frequency_index)) then
        angular_frequency = 2.0_dp*pi*response%frequencies_hz(frequency_index)
        angular_frequency_squared = angular_frequency*angular_frequency
        accelerations(:, frequency_index) = &
          -angular_frequency_squared * &
          response%physical_complex_response(:, frequency_index)
        call require_finite_complex_vector(accelerations(:, frequency_index))
      else
        accelerations(:, frequency_index) = &
          response%physical_complex_response(:, frequency_index)
      end if
    end do
  end function get_physical_angular_accelerations

  !> Bir fiziksel (node_id,dof_type) kanalının bütün sweep response'unu [rad]
  !! döndürür. Physical DOF haritası ile result satır sırası doğrulanır.
  pure function get_physical_node_response( &
      response, mapping, node_id, dof_type) result(values)
    type(harmonic_response_t), intent(in) :: response
    type(active_dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in) :: dof_type
    complex(dp), allocatable :: values(:)

    integer :: physical_dof_id

    call validate_harmonic_response(response)
    physical_dof_id = find_physical_dof_id(response, mapping, node_id, dof_type)
    values = response%physical_complex_response(physical_dof_id, :)
  end function get_physical_node_response

  !> Oriented element relative complex angle değerini hesaplar.
  !! Matematiksel model: Delta_theta_hat=theta_i_hat-theta_j_hat; orientation
  !! node_i -> node_j'dir. Girdiler/çıktı [rad] complex peak genliklerdir.
  pure function calculate_element_relative_angle( &
      node_i_response, node_j_response) result(relative_angle)
    complex(dp), intent(in) :: node_i_response
    complex(dp), intent(in) :: node_j_response
    complex(dp) :: relative_angle

    call require_finite_complex(node_i_response, "Node i response")
    call require_finite_complex(node_j_response, "Node j response")
    relative_angle = node_i_response-node_j_response
    call require_finite_complex(relative_angle, "Element relative angle")
  end function calculate_element_relative_angle

  !> Elementin oriented complex internal generalized torque'unu hesaplar.
  !!
  !! Fiziksel açıklama: Storage, structural-loss ve viscous kanalları aynı
  !! relative harmonic rotation üzerinde etki eder fakat K'' ile c ayrı kalır.
  !! Matematiksel model: T_e=[K'+i(K''+omega*c)]*Delta_theta_hat.
  !! Girdiler: Geçerli element, f>0 [Hz], relative angle [rad]. Çıktı: Complex
  !! peak T_e [N*m]. Matrix-equilibrium işareti node_i için +T_e, node_j için
  !! -T_e'dir. Küçük genlikli frozen-property model kabul edilir.
  pure function calculate_element_dynamic_torque( &
      element, frequency_hz, relative_angle) result(dynamic_torque)
    type(torsional_element_t), intent(in) :: element
    real(dp), intent(in) :: frequency_hz
    complex(dp), intent(in) :: relative_angle
    complex(dp) :: dynamic_torque

    real(dp) :: angular_frequency
    real(dp) :: dissipative_coefficient

    call validate_element_response_input(element, frequency_hz, relative_angle)
    angular_frequency = 2.0_dp*pi*frequency_hz
    dissipative_coefficient = element%loss_stiffness_nm_per_rad + &
      angular_frequency*element%damping_nms_per_rad
    if (.not. ieee_is_finite(dissipative_coefficient)) then
      error stop "Element dissipative coefficient sonlu olmalıdır."
    end if
    dynamic_torque = cmplx( &
      element%stiffness_nm_per_rad, dissipative_coefficient, kind=dp) * &
      relative_angle
    call require_finite_complex(dynamic_torque, "Element dynamic torque")
  end function calculate_element_dynamic_torque

  !> Passive elementin ortalama dissipated power değerini hesaplar.
  !! Matematiksel model: P_avg=(omega/2)*(K''+omega*c)*|Delta theta|^2.
  !! Girdiler element, f>0 [Hz], relative complex peak angle [rad]; çıktı [W].
  !! PEAK-amplitude convention kullanılır. Negatif fiziksel enerji kabul
  !! edilmez; roundoff düzeyindeki negatiflik scale-aware toleransla sıfırlanır.
  pure function calculate_element_average_dissipated_power( &
      element, frequency_hz, relative_angle) result(average_power_w)
    type(torsional_element_t), intent(in) :: element
    real(dp), intent(in) :: frequency_hz
    complex(dp), intent(in) :: relative_angle
    real(dp) :: average_power_w

    real(dp) :: angular_frequency
    real(dp) :: dissipative_coefficient
    real(dp) :: magnitude_squared

    call validate_element_response_input(element, frequency_hz, relative_angle)
    angular_frequency = 2.0_dp*pi*frequency_hz
    dissipative_coefficient = element%loss_stiffness_nm_per_rad + &
      angular_frequency*element%damping_nms_per_rad
    magnitude_squared = stable_complex_magnitude_squared(relative_angle)
    average_power_w = 0.5_dp*angular_frequency*dissipative_coefficient * &
      magnitude_squared
    call enforce_passive_energy(average_power_w, "Average dissipated power")
  end function calculate_element_average_dissipated_power

  !> Passive elementin bir harmonic çevrimde kaybettiği enerjiyi hesaplar.
  !! Matematiksel model: E_cycle=pi*(K''+omega*c)*|Delta theta|^2.
  !! Girdiler element, f>0 [Hz], relative complex peak angle [rad]; çıktı
  !! [J/cycle]. PEAK-amplitude convention ve lineer passive model kullanılır.
  pure function calculate_element_dissipated_energy_per_cycle( &
      element, frequency_hz, relative_angle) result(energy_j_per_cycle)
    type(torsional_element_t), intent(in) :: element
    real(dp), intent(in) :: frequency_hz
    complex(dp), intent(in) :: relative_angle
    real(dp) :: energy_j_per_cycle

    real(dp) :: angular_frequency
    real(dp) :: dissipative_coefficient
    real(dp) :: magnitude_squared

    call validate_element_response_input(element, frequency_hz, relative_angle)
    angular_frequency = 2.0_dp*pi*frequency_hz
    dissipative_coefficient = element%loss_stiffness_nm_per_rad + &
      angular_frequency*element%damping_nms_per_rad
    magnitude_squared = stable_complex_magnitude_squared(relative_angle)
    energy_j_per_cycle = pi*dissipative_coefficient*magnitude_squared
    call enforce_passive_energy(energy_j_per_cycle, "Dissipated energy")
  end function calculate_element_dissipated_energy_per_cycle

  !> Bir elementin oriented relative angular response'unu bütün sweep için
  !! [rad] döndürür. Singular frequency point'ler NaN unavailable kalır.
  pure function get_element_relative_angle_response( &
      response, mapping, element) result(relative_angles)
    type(harmonic_response_t), intent(in) :: response
    type(active_dof_map_t), intent(in) :: mapping
    type(torsional_element_t), intent(in) :: element
    complex(dp), allocatable :: relative_angles(:)

    complex(dp), allocatable :: node_i_response(:)
    complex(dp), allocatable :: node_j_response(:)
    integer :: frequency_index

    call validate_torsional_element(element)
    node_i_response = get_physical_node_response( &
      response, mapping, element%node_i_id, TORSIONAL_ROTATION)
    node_j_response = get_physical_node_response( &
      response, mapping, element%node_j_id, TORSIONAL_ROTATION)
    allocate(relative_angles(size(node_i_response)))
    do frequency_index = 1, size(relative_angles)
      if (response%response_available(frequency_index)) then
        relative_angles(frequency_index) = calculate_element_relative_angle( &
          node_i_response(frequency_index), node_j_response(frequency_index))
      else
        relative_angles(frequency_index) = node_i_response(frequency_index)
      end if
    end do
  end function get_element_relative_angle_response

  !> Element complex internal torque response'unu [N*m] bütün sweep için
  !! döndürür. node_i için +T_e, node_j için -T_e işareti kullanılır.
  pure function get_element_dynamic_torque_response( &
      response, mapping, element) result(dynamic_torques)
    type(harmonic_response_t), intent(in) :: response
    type(active_dof_map_t), intent(in) :: mapping
    type(torsional_element_t), intent(in) :: element
    complex(dp), allocatable :: dynamic_torques(:)

    complex(dp), allocatable :: relative_angles(:)
    integer :: frequency_index

    relative_angles = get_element_relative_angle_response( &
      response, mapping, element)
    allocate(dynamic_torques(size(relative_angles)))
    do frequency_index = 1, size(dynamic_torques)
      if (response%response_available(frequency_index)) then
        dynamic_torques(frequency_index) = calculate_element_dynamic_torque( &
          element, response%frequencies_hz(frequency_index), &
          relative_angles(frequency_index))
      else
        dynamic_torques(frequency_index) = relative_angles(frequency_index)
      end if
    end do
  end function get_element_dynamic_torque_response

  !> Element internal torque response'unun transmitted magnitude değerlerini
  !! [N*m] olarak döndürür.
  !! Matematiksel model her solved noktada |T_hat_e|=abs(T_hat_e)'dir; complex
  !! torque node_i->node_j orientation'ını taşırken magnitude işaretsizdir.
  !! Girdiler harmonic result, aynı physical mapping ve geçerli elemandır.
  !! Singular frequency point'ler çözüm uydurulmadan NaN unavailable kalır.
  pure function get_element_dynamic_torque_magnitudes( &
      response, mapping, element) result(torque_magnitudes)
    type(harmonic_response_t), intent(in) :: response
    type(active_dof_map_t), intent(in) :: mapping
    type(torsional_element_t), intent(in) :: element
    real(dp), allocatable :: torque_magnitudes(:)

    complex(dp), allocatable :: dynamic_torques(:)
    integer :: frequency_index

    dynamic_torques = get_element_dynamic_torque_response( &
      response, mapping, element)
    allocate(torque_magnitudes(size(dynamic_torques)))
    do frequency_index = 1, size(torque_magnitudes)
      if (response%response_available(frequency_index)) then
        torque_magnitudes(frequency_index) = &
          abs(dynamic_torques(frequency_index))
      else
        torque_magnitudes(frequency_index) = &
          real(dynamic_torques(frequency_index), dp)
      end if
    end do
  end function get_element_dynamic_torque_magnitudes

  !> Element ortalama dissipated power sweep'ini [W] döndürür. Singular
  !! frequency point'ler NaN unavailable kalır.
  pure function get_element_average_dissipated_power_response( &
      response, mapping, element) result(average_power_w)
    type(harmonic_response_t), intent(in) :: response
    type(active_dof_map_t), intent(in) :: mapping
    type(torsional_element_t), intent(in) :: element
    real(dp), allocatable :: average_power_w(:)

    complex(dp), allocatable :: relative_angles(:)
    integer :: frequency_index

    relative_angles = get_element_relative_angle_response( &
      response, mapping, element)
    allocate(average_power_w(size(relative_angles)))
    do frequency_index = 1, size(average_power_w)
      if (response%response_available(frequency_index)) then
        average_power_w(frequency_index) = &
          calculate_element_average_dissipated_power( &
            element, response%frequencies_hz(frequency_index), &
            relative_angles(frequency_index))
      else
        average_power_w(frequency_index) = &
          real(relative_angles(frequency_index), dp)
      end if
    end do
  end function get_element_average_dissipated_power_response

  !> Element dissipated energy-per-cycle sweep'ini [J/cycle] döndürür.
  !! Singular frequency point'ler NaN unavailable kalır.
  pure function get_element_dissipated_energy_response( &
      response, mapping, element) result(energy_j_per_cycle)
    type(harmonic_response_t), intent(in) :: response
    type(active_dof_map_t), intent(in) :: mapping
    type(torsional_element_t), intent(in) :: element
    real(dp), allocatable :: energy_j_per_cycle(:)

    complex(dp), allocatable :: relative_angles(:)
    integer :: frequency_index

    relative_angles = get_element_relative_angle_response( &
      response, mapping, element)
    allocate(energy_j_per_cycle(size(relative_angles)))
    do frequency_index = 1, size(energy_j_per_cycle)
      if (response%response_available(frequency_index)) then
        energy_j_per_cycle(frequency_index) = &
          calculate_element_dissipated_energy_per_cycle( &
            element, response%frequencies_hz(frequency_index), &
            relative_angles(frequency_index))
      else
        energy_j_per_cycle(frequency_index) = &
          real(relative_angles(frequency_index), dp)
      end if
    end do
  end function get_element_dissipated_energy_response

  !> Explicit frequency sweep sözleşmesini doğrular: dizi boş olmayan, sonlu,
  !! pozitif ve kesin artan [Hz] değerlerden oluşmalıdır.
  pure subroutine validate_harmonic_frequency_sweep(frequencies_hz)
    real(dp), intent(in) :: frequencies_hz(:)

    integer :: frequency_index

    if (size(frequencies_hz) == 0) then
      error stop "Harmonic frequency sweep boş olamaz."
    end if
    do frequency_index = 1, size(frequencies_hz)
      if (.not. ieee_is_finite(frequencies_hz(frequency_index)) .or. &
          frequencies_hz(frequency_index) <= 0.0_dp) then
        error stop "Harmonic frekanslar sonlu ve pozitif olmalıdır."
      end if
    end do
    do frequency_index = 2, size(frequencies_hz)
      if (frequencies_hz(frequency_index) <= &
          frequencies_hz(frequency_index-1)) then
        error stop "Harmonic frequency sweep kesin artan olmalıdır."
      end if
    end do
  end subroutine validate_harmonic_frequency_sweep

  !> Status değerinin unique solution içerip içermediğini sınar.
  pure function is_solved_status(status) result(is_solved)
    integer, intent(in) :: status
    logical :: is_solved

    is_solved = status == COMPLEX_SOLVE_SOLVED .or. &
      status == COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED
  end function is_solved_status

  !> Status değerinin V0.6 solver durum kümesinde olduğunu sınar.
  pure function is_known_status(status) result(is_known)
    integer, intent(in) :: status
    logical :: is_known

    is_known = is_solved_status(status) .or. status == COMPLEX_SOLVE_SINGULAR
  end function is_known_status

  !> İstenen frequency indeksinin geçerli ve response-available olduğunu sınar.
  pure subroutine require_available_frequency(response, frequency_index)
    type(harmonic_response_t), intent(in) :: response
    integer, intent(in) :: frequency_index

    call validate_harmonic_response(response)
    if (frequency_index < 1 .or. &
        frequency_index > size(response%frequencies_hz)) then
      error stop "Harmonic response frequency indeksi geçersiz."
    end if
    if (.not. response%response_available(frequency_index)) then
      error stop "Singular harmonic noktada unique response mevcut değildir."
    end if
  end subroutine require_available_frequency

  !> Physical DOF anahtarının result satır kimliğini active haritada bulur.
  pure function find_physical_dof_id( &
      response, mapping, node_id, dof_type) result(physical_dof_id)
    type(harmonic_response_t), intent(in) :: response
    type(active_dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in) :: dof_type
    integer :: physical_dof_id

    type(active_dof_map_entry_t), allocatable :: entries(:)
    integer :: entry_index

    call validate_active_dof_map(mapping)
    if (node_id <= 0 .or. .not. is_supported_dof_type(dof_type)) then
      error stop "Physical harmonic response DOF anahtarı geçersiz."
    end if
    entries = get_active_dof_map_entries(mapping)
    if (size(entries) /= size(response%physical_complex_response, 1)) then
      error stop "Harmonic response ile active DOF haritası boyutu uyumsuz."
    end if

    physical_dof_id = 0
    do entry_index = 1, size(entries)
      if (entries(entry_index)%physical_dof%node_id == node_id .and. &
          entries(entry_index)%physical_dof%dof_type == dof_type) then
        physical_dof_id = entries(entry_index)%physical_dof_id
        exit
      end if
    end do
    if (physical_dof_id == 0) then
      error stop "Physical harmonic response DOF haritada bulunamadı."
    end if
  end function find_physical_dof_id

  !> Element response hesabının ortak element/frequency/angle girdilerini
  !! doğrular. Birimler sırasıyla element sözleşmesi, [Hz] ve [rad]'dır.
  pure subroutine validate_element_response_input( &
      element, frequency_hz, relative_angle)
    type(torsional_element_t), intent(in) :: element
    real(dp), intent(in) :: frequency_hz
    complex(dp), intent(in) :: relative_angle

    call validate_torsional_element(element)
    if (.not. ieee_is_finite(frequency_hz) .or. frequency_hz <= 0.0_dp) then
      error stop "Element harmonic frekansı sonlu ve pozitif olmalıdır."
    end if
    call require_finite_complex(relative_angle, "Element relative angle")
  end subroutine validate_element_response_input

  !> Complex değerin gerçek ve sanal bileşenlerinin sonluluğunu doğrular.
  pure subroutine require_finite_complex(value, quantity_name)
    complex(dp), intent(in) :: value
    character(len=*), intent(in) :: quantity_name

    if (.not. ieee_is_finite(real(value, dp)) .or. &
        .not. ieee_is_finite(aimag(value))) then
      error stop quantity_name//" sonlu complex değer olmalıdır."
    end if
  end subroutine require_finite_complex

  !> Complex vektörün bütün bileşenlerinin sonlu olduğunu doğrular.
  pure subroutine require_finite_complex_vector(values)
    complex(dp), intent(in) :: values(:)

    integer :: value_index

    do value_index = 1, size(values)
      call require_finite_complex(values(value_index), "Complex response")
    end do
  end subroutine require_finite_complex_vector

  !> Singular response vektörünün çözüm yerine explicit complex NaN sentinel
  !! taşıdığını doğrular.
  pure subroutine require_unavailable_complex_vector(values)
    complex(dp), intent(in) :: values(:)

    integer :: value_index

    do value_index = 1, size(values)
      if (.not. ieee_is_nan(real(values(value_index), dp)) .or. &
          .not. ieee_is_nan(aimag(values(value_index)))) then
        error stop "Singular harmonic response unavailable sentinel taşımalıdır."
      end if
    end do
  end subroutine require_unavailable_complex_vector

  !> Sayısal tanının sonlu ve negatif olmayan boyutsuz değerini doğrular.
  pure subroutine require_finite_nonnegative(value, quantity_name)
    real(dp), intent(in) :: value
    character(len=*), intent(in) :: quantity_name

    if (.not. ieee_is_finite(value) .or. value < 0.0_dp) then
      error stop quantity_name//" sonlu ve negatif olmayan değer olmalıdır."
    end if
  end subroutine require_finite_nonnegative

  !> Complex genliğin karesini taşma riskini azaltan ölçekli normla hesaplar.
  !! Matematiksel model: |z|^2=Re(z)^2+Im(z)^2. Çıktı, girdinin biriminin
  !! karesindedir; sonlu temsil edilemiyorsa clean diagnostic üretilir.
  pure function stable_complex_magnitude_squared(value) result(magnitude_squared)
    complex(dp), intent(in) :: value
    real(dp) :: magnitude_squared

    real(dp) :: scale
    real(dp) :: scaled_imaginary
    real(dp) :: scaled_real

    call require_finite_complex(value, "Complex magnitude girdisi")
    scale = max(abs(real(value, dp)), abs(aimag(value)))
    if (.not. scale > 0.0_dp) then
      magnitude_squared = 0.0_dp
      return
    end if
    scaled_real = real(value, dp)/scale
    scaled_imaginary = aimag(value)/scale
    magnitude_squared = scale*scale * &
      (scaled_real*scaled_real+scaled_imaginary*scaled_imaginary)
    if (.not. ieee_is_finite(magnitude_squared)) then
      error stop "Complex genlik karesi sonlu sayı aralığında olmalıdır."
    end if
  end function stable_complex_magnitude_squared

  !> Passive enerji/güç sonucunu scale-aware roundoff toleransıyla doğrular.
  !! Negatif fiziksel değer reddedilir; yalnız sıfır çevresindeki yuvarlama
  !! sapması sıfıra kırpılır. Birim çağıran yordamın [W] veya [J/cycle]
  !! sözleşmesiyle belirlenir.
  pure subroutine enforce_passive_energy(value, quantity_name)
    real(dp), intent(inout) :: value
    character(len=*), intent(in) :: quantity_name

    real(dp) :: tolerance

    if (.not. ieee_is_finite(value)) then
      error stop quantity_name//" sonlu sayı aralığında olmalıdır."
    end if
    tolerance = 64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(value))
    if (value < -tolerance) then
      error stop quantity_name//" passive element için negatif olamaz."
    end if
    if (value < 0.0_dp) value = 0.0_dp
  end subroutine enforce_passive_energy

end module tms_harmonic_response
