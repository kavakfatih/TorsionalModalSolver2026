module tms_complex_linear_solver
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_complex_linear_problem, only : complex_linear_problem_t, &
    validate_complex_linear_problem, get_complex_linear_problem_order, &
    get_complex_linear_rhs_count, get_complex_linear_coefficient_matrix, &
    get_complex_linear_right_hand_sides
  use tms_complex_linear_solution, only : complex_linear_solution_t, &
    COMPLEX_SOLVE_SINGULAR, initialize_available_complex_linear_solution, &
    initialize_singular_complex_linear_solution, &
    validate_complex_linear_solution
  use tms_lapack_zsysvx_backend, only : solve_with_lapack_zsysvx
  implicit none
  private

  public :: solve_complex_linear_problem
  public :: calculate_complex_relative_residuals

contains

  !> Backend-neutral complex symmetric doğrusal solver facade'ıdır.
  !!
  !! Matematiksel açıklama: A*X=B problemi reference LAPACK ZSYSVX backend'ine
  !! gönderilir. Backend status'u SINGULAR ise çözüm uydurulmaz; diğer
  !! durumlarda her RHS için
  !! rho=||A*x-b||_2/(||A||_inf*||x||_2+||b||_2)
  !! backend'den bağımsız olarak hesaplanır.
  !! Fiziksel açıklama: Harmonic kullanımda A=Z [N*m/rad], B complex peak
  !! torque [N*m], X complex peak angular response [rad] olur. Facade harmonic
  !! frekans, assembly veya physical recovery bilmez.
  !! Girdi immutable private problem, çıktı private backend-neutral solution
  !! değeridir. External numerical backend çağrısı nedeniyle PURE değildir.
  function solve_complex_linear_problem(problem) result(solution)
    type(complex_linear_problem_t), intent(in) :: problem
    type(complex_linear_solution_t) :: solution

    complex(dp), allocatable :: response(:, :)
    real(dp), allocatable :: backward_errors(:)
    real(dp), allocatable :: forward_error_bounds(:)
    real(dp), allocatable :: relative_residuals(:)
    real(dp) :: reciprocal_condition_number
    character(len=:), allocatable :: backend_identity
    integer :: status

    call validate_complex_linear_problem(problem)
    call solve_with_lapack_zsysvx( &
      problem, status, response, reciprocal_condition_number, &
      forward_error_bounds, backward_errors, backend_identity)

    if (status == COMPLEX_SOLVE_SINGULAR) then
      call initialize_singular_complex_linear_solution( &
        solution, get_complex_linear_problem_order(problem), &
        get_complex_linear_rhs_count(problem), reciprocal_condition_number, &
        backend_identity)
    else
      relative_residuals = calculate_complex_relative_residuals( &
        problem, response)
      call initialize_available_complex_linear_solution( &
        solution, status, response, reciprocal_condition_number, &
        forward_error_bounds, backward_errors, relative_residuals, &
        backend_identity)
    end if

    call validate_complex_linear_solution(solution)
  end function solve_complex_linear_problem

  !> Her right-hand side için backend-independent complex relative residual
  !! hesaplar.
  !!
  !! Matematiksel model:
  !!   r_j=A*x_j-b_j
  !!   rho_j=||r_j||_2/(||A||_inf*||x_j||_2+||b_j||_2).
  !! Normlar overflow/underflow riskini azaltan ölçekli hesaplarla bulunur.
  !! Pay ve payda birlikte sıfırsa rho=0 kabul edilir; sıfır payda ile pozitif
  !! residual geçersizdir.
  !! Fiziksel açıklama: Harmonic kullanımda pay ve payda [N*m], rho ise
  !! boyutsuzdur. Girdiler private A/B problem ve aynı n x nrhs boyutunda sonlu
  !! X'tir. Çıktı nrhs uzunluğunda sonlu, negatif olmayan residual dizisidir.
  pure function calculate_complex_relative_residuals(problem, response) &
      result(residuals)
    type(complex_linear_problem_t), intent(in) :: problem
    complex(dp), intent(in) :: response(:, :)
    real(dp), allocatable :: residuals(:)

    complex(dp), allocatable :: coefficient_matrix(:, :)
    complex(dp), allocatable :: residual_vector(:)
    complex(dp), allocatable :: right_hand_sides(:, :)
    real(dp) :: coefficient_norm
    real(dp) :: denominator
    real(dp) :: response_norm
    real(dp) :: residual_norm
    real(dp) :: rhs_norm
    integer :: rhs_index

    call validate_complex_linear_problem(problem)
    if (size(response, 1) /= get_complex_linear_problem_order(problem) .or. &
        size(response, 2) /= get_complex_linear_rhs_count(problem)) then
      error stop "Complex residual response boyutu problem A/B boyutuyla uyumsuz."
    end if
    if (.not. all_complex_finite(response)) then
      error stop "Complex residual yalnız sonlu response değerleri kabul eder."
    end if

    coefficient_matrix = get_complex_linear_coefficient_matrix(problem)
    right_hand_sides = get_complex_linear_right_hand_sides(problem)
    coefficient_norm = complex_matrix_infinity_norm(coefficient_matrix)
    allocate(residuals(size(response, 2)))

    do rhs_index = 1, size(response, 2)
      residual_vector = matmul( &
        coefficient_matrix, response(:, rhs_index)) - &
        right_hand_sides(:, rhs_index)
      if (.not. all_complex_vector_finite(residual_vector)) then
        error stop "Complex doğrusal residual sonlu sayı aralığında kalmalıdır."
      end if

      residual_norm = stable_complex_vector_two_norm(residual_vector)
      response_norm = stable_complex_vector_two_norm(response(:, rhs_index))
      rhs_norm = stable_complex_vector_two_norm( &
        right_hand_sides(:, rhs_index))
      denominator = safe_nonnegative_sum( &
        safe_nonnegative_product(coefficient_norm, response_norm), rhs_norm)

      if (denominator <= 0.0_dp) then
        if (residual_norm <= 0.0_dp) then
          residuals(rhs_index) = 0.0_dp
        else
          error stop "Complex relative residual paydası sıfırken payı pozitiftir."
        end if
      else
        residuals(rhs_index) = residual_norm / denominator
      end if

      if (.not. ieee_is_finite(residuals(rhs_index)) .or. &
          residuals(rhs_index) < 0.0_dp) then
        error stop "Complex relative residual sonlu ve negatif olmayan olmalıdır."
      end if
    end do
  end function calculate_complex_relative_residuals

  !> Sonlu complex matrisin infinity normunu ölçekli satır toplamlarıyla
  !! hesaplar. Çıktı matris katsayısıyla aynı fiziksel birimdedir.
  pure function complex_matrix_infinity_norm(matrix) result(norm_value)
    complex(dp), intent(in) :: matrix(:, :)
    real(dp) :: norm_value

    real(dp) :: row_norm
    real(dp) :: row_scale
    integer :: row

    if (.not. all_complex_finite(matrix)) then
      error stop "Complex matris normu yalnız sonlu katsayılarla hesaplanabilir."
    end if
    norm_value = 0.0_dp
    do row = 1, size(matrix, 1)
      row_scale = maxval(abs(matrix(row, :)))
      if (row_scale <= 0.0_dp) cycle
      row_norm = safe_nonnegative_product( &
        row_scale, sum(abs(matrix(row, :)) / row_scale))
      norm_value = max(norm_value, row_norm)
    end do
  end function complex_matrix_infinity_norm

  !> Sonlu complex vektörün iki normunu ölçekli biçimde hesaplar.
  !! Çıktı vektör bileşenleriyle aynı fiziksel birimdedir.
  pure function stable_complex_vector_two_norm(vector) result(norm_value)
    complex(dp), intent(in) :: vector(:)
    real(dp) :: norm_value

    real(dp) :: scale

    if (.not. all_complex_vector_finite(vector)) then
      error stop "Complex vektör normu yalnız sonlu bileşenlerle hesaplanabilir."
    end if
    if (size(vector) == 0) then
      norm_value = 0.0_dp
      return
    end if
    scale = maxval(abs(vector))
    if (scale <= 0.0_dp) then
      norm_value = 0.0_dp
    else
      norm_value = safe_nonnegative_product( &
        scale, sqrt(sum(abs(vector / scale)**2)))
    end if
  end function stable_complex_vector_two_norm

  !> Negatif olmayan sonlu iki büyüklüğün overflow'a dirençli çarpımını
  !! döndürür. Girdi/çıktı birimleri çağıran norm ifadesine bağlıdır.
  pure function safe_nonnegative_product(first, second) result(product_value)
    real(dp), intent(in) :: first
    real(dp), intent(in) :: second
    real(dp) :: product_value

    if (.not. ieee_is_finite(first) .or. .not. ieee_is_finite(second) .or. &
        first < 0.0_dp .or. second < 0.0_dp) then
      error stop "Güvenli complex norm çarpımı sonlu negatif olmayan girdi ister."
    end if
    if (first <= 0.0_dp .or. second <= 0.0_dp) then
      product_value = 0.0_dp
    else if (first > huge(1.0_dp) / second) then
      product_value = huge(1.0_dp)
    else
      product_value = first * second
    end if
  end function safe_nonnegative_product

  !> Negatif olmayan sonlu iki büyüklüğün overflow'a dirençli toplamını
  !! döndürür. Girdi/çıktı birimleri çağıran norm ifadesine bağlıdır.
  pure function safe_nonnegative_sum(first, second) result(sum_value)
    real(dp), intent(in) :: first
    real(dp), intent(in) :: second
    real(dp) :: sum_value

    if (.not. ieee_is_finite(first) .or. .not. ieee_is_finite(second) .or. &
        first < 0.0_dp .or. second < 0.0_dp) then
      error stop "Güvenli complex norm toplamı sonlu negatif olmayan girdi ister."
    end if
    if (first > huge(1.0_dp) - second) then
      sum_value = huge(1.0_dp)
    else
      sum_value = first + second
    end if
  end function safe_nonnegative_sum

  !> Complex matrisin bütün reel/sanal bileşenlerinin sonluluğunu değerlendirir.
  pure function all_complex_finite(values) result(is_finite)
    complex(dp), intent(in) :: values(:, :)
    logical :: is_finite

    is_finite = all(ieee_is_finite(real(values, dp))) .and. &
      all(ieee_is_finite(aimag(values)))
  end function all_complex_finite

  !> Complex vektörün bütün reel/sanal bileşenlerinin sonluluğunu değerlendirir.
  pure function all_complex_vector_finite(values) result(is_finite)
    complex(dp), intent(in) :: values(:)
    logical :: is_finite

    is_finite = all(ieee_is_finite(real(values, dp))) .and. &
      all(ieee_is_finite(aimag(values)))
  end function all_complex_vector_finite

end module tms_complex_linear_solver
