program test_complex_linear_solver
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
    ieee_positive_inf, ieee_is_finite
  use tms_kinds, only : dp
  use tms_complex_linear_problem, only : complex_linear_problem_t, &
    create_complex_linear_problem, get_complex_linear_problem_order, &
    get_complex_linear_rhs_count, get_complex_linear_coefficient_matrix, &
    get_complex_linear_right_hand_sides
  use tms_complex_linear_solution, only : complex_linear_solution_t, &
    COMPLEX_SOLVE_SOLVED, COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED, &
    COMPLEX_SOLVE_SINGULAR, get_complex_linear_solution_status, &
    has_complex_linear_response, get_complex_linear_solution_equation_count, &
    get_complex_linear_solution_rhs_count, get_complex_linear_response, &
    get_complex_linear_reciprocal_condition_number, &
    get_complex_linear_forward_error_bounds, &
    get_complex_linear_backward_errors, &
    get_complex_linear_relative_residuals, &
    get_complex_linear_backend_identity
  use tms_complex_linear_solver, only : solve_complex_linear_problem
  implicit none

  real(dp), parameter :: solution_tolerance = 1.0e-11_dp
  real(dp), parameter :: residual_tolerance = 1.0e-11_dp
  character(len=64) :: validation_case

  ! Geçersiz problem veya status-gated erişim vakaları ayrı CTest süreçlerinde
  ! error stop üretmelidir. Bilinmeyen selector normal döner ve WILL_FAIL
  ! kaydını başarısız kılar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_known_nonhermitian_solution_and_immutability()
  call test_multiple_right_hand_sides()
  call test_exact_singular_status()
  call test_deterministic_ill_conditioned_status()

  print *, "V0.6 backend-neutral complex linear solver doğrulamaları başarılı."

contains

  !> Complex symmetric fakat Hermitian olmayan 2x2 problemin bilinen çözümünü
  !! doğrular.
  !!
  !! Matematiksel model A*X=B'dir; A^T=A ve A^H/=A. Bilinen
  !! x=[1+2i,-0.5+0.25i]^T için B bağımsız sabitlerle verilir. Test çözüm,
  !! status, RCOND, FERR/BERR, backend-independent residual ve backend kimliğini
  !! sınar. Original A/B hem çağıran dizilerinde hem private problemde exact
  !! değişmeden kalmalıdır. Harmonic bağlamda A [N*m/rad], B [N*m], X [rad]
  !! olarak yorumlanabilir.
  subroutine test_known_nonhermitian_solution_and_immutability()
    type(complex_linear_problem_t) :: problem
    type(complex_linear_solution_t) :: solution
    complex(dp) :: coefficient_matrix(2, 2)
    complex(dp) :: right_hand_side(2, 1)
    complex(dp) :: expected_response(2, 1)
    complex(dp), allocatable :: coefficient_after(:, :)
    complex(dp), allocatable :: coefficient_before(:, :)
    complex(dp), allocatable :: response(:, :)
    complex(dp), allocatable :: rhs_after(:, :)
    complex(dp), allocatable :: rhs_before(:, :)
    character(len=:), allocatable :: backend_identity

    coefficient_matrix = reshape([ &
      cmplx(2.0_dp, 1.0_dp, kind=dp), &
      cmplx(1.0_dp, 2.0_dp, kind=dp), &
      cmplx(1.0_dp, 2.0_dp, kind=dp), &
      cmplx(3.0_dp, -0.5_dp, kind=dp)], [2, 2])
    right_hand_side(:, 1) = [ &
      cmplx(-1.0_dp, 4.25_dp, kind=dp), &
      cmplx(-4.375_dp, 5.0_dp, kind=dp)]
    expected_response(:, 1) = [ &
      cmplx(1.0_dp, 2.0_dp, kind=dp), &
      cmplx(-0.5_dp, 0.25_dp, kind=dp)]

    if (any(abs(coefficient_matrix - transpose(coefficient_matrix)) > &
        0.0_dp)) then
      error stop "Sentetik Z matrisi exact complex symmetric değil."
    end if
    if (maxval(abs(coefficient_matrix - &
        conjg(transpose(coefficient_matrix)))) <= 1.0_dp) then
      error stop "Sentetik Z matrisi Hermitian olmamalıdır."
    end if

    problem = create_complex_linear_problem( &
      coefficient_matrix, right_hand_side)
    coefficient_before = get_complex_linear_coefficient_matrix(problem)
    rhs_before = get_complex_linear_right_hand_sides(problem)
    solution = solve_complex_linear_problem(problem)
    coefficient_after = get_complex_linear_coefficient_matrix(problem)
    rhs_after = get_complex_linear_right_hand_sides(problem)

    if (get_complex_linear_solution_status(solution) /= &
        COMPLEX_SOLVE_SOLVED) then
      error stop "İyi koşullu complex problem SOLVED status üretmedi."
    end if
    if (.not. has_complex_linear_response(solution)) then
      error stop "SOLVED complex problem response taşımıyor."
    end if
    if (get_complex_linear_problem_order(problem) /= 2 .or. &
        get_complex_linear_rhs_count(problem) /= 1 .or. &
        get_complex_linear_solution_equation_count(solution) /= 2 .or. &
        get_complex_linear_solution_rhs_count(solution) /= 1) then
      error stop "Complex problem/solution boyut metadata'sı yanlış."
    end if

    response = get_complex_linear_response(solution)
    call assert_complex_matrix_close( &
      response, expected_response, solution_tolerance, &
      "Complex symmetric non-Hermitian bilinen çözüm yanlış.")
    call assert_available_diagnostics(solution, 1)
    backend_identity = get_complex_linear_backend_identity(solution)
    if (index(backend_identity, "ZSYSVX") == 0) then
      error stop "Complex solution ZSYSVX backend kimliğini taşımıyor."
    end if

    if (any(abs(coefficient_after - coefficient_before) > 0.0_dp) .or. &
        any(abs(rhs_after - rhs_before) > 0.0_dp) .or. &
        any(abs(coefficient_matrix - coefficient_before) > 0.0_dp) .or. &
        any(abs(right_hand_side - rhs_before) > 0.0_dp)) then
      error stop "ZSYSVX original A/B girdilerini değiştirdi."
    end if
  end subroutine test_known_nonhermitian_solution_and_immutability

  !> Tek coefficient matrix ile iki bağımsız RHS'nin aynı çağrıda çözüldüğünü
  !! doğrular. Test ZSYSVX NRHS>1 yolunu, X sütun eşleşmesini ve FERR/BERR ile
  !! backend-independent residual dizilerinin RHS başına korunmasını sınar.
  !! Bu regression harmonic public API tek load case kullansa bile low-level
  !! abstraction'ın gelecekte unit-torque FRF sütunlarını engellememesini sağlar.
  subroutine test_multiple_right_hand_sides()
    type(complex_linear_problem_t) :: problem
    type(complex_linear_solution_t) :: solution
    complex(dp) :: coefficient_matrix(2, 2)
    complex(dp) :: expected_response(2, 2)
    complex(dp) :: right_hand_sides(2, 2)
    complex(dp), allocatable :: response(:, :)

    coefficient_matrix = reshape([ &
      cmplx(3.0_dp, 0.5_dp, kind=dp), &
      cmplx(-0.5_dp, 0.75_dp, kind=dp), &
      cmplx(-0.5_dp, 0.75_dp, kind=dp), &
      cmplx(2.0_dp, 1.5_dp, kind=dp)], [2, 2])
    expected_response = reshape([ &
      cmplx(1.0_dp, -0.5_dp, kind=dp), &
      cmplx(0.25_dp, 0.75_dp, kind=dp), &
      cmplx(-2.0_dp, 1.0_dp, kind=dp), &
      cmplx(0.5_dp, -1.5_dp, kind=dp)], [2, 2])
    right_hand_sides = matmul(coefficient_matrix, expected_response)

    problem = create_complex_linear_problem( &
      coefficient_matrix, right_hand_sides)
    solution = solve_complex_linear_problem(problem)
    if (get_complex_linear_solution_status(solution) /= &
        COMPLEX_SOLVE_SOLVED) then
      error stop "Multiple-RHS complex problem SOLVED status üretmedi."
    end if

    response = get_complex_linear_response(solution)
    call assert_complex_matrix_close( &
      response, expected_response, solution_tolerance, &
      "Multiple-RHS complex çözümler bilinen sütunlarla uyuşmuyor.")
    call assert_available_diagnostics(solution, 2)
  end subroutine test_multiple_right_hand_sides

  !> Exact singular complex symmetric factorization'ın analysis-state olarak
  !! döndüğünü doğrular. A=diag(1+i,0) için unique çözüm yoktur; ZSYSVX
  !! 1<=INFO<=N yolunu SINGULAR'a çevirmeli, RCOND=0 saklamalı ve undefined
  !! X/FERR/BERR değerlerini response gibi sunmamalıdır. Test crash veya error
  !! stop beklemez; sweep katmanı bu status'tan sonra devam edebilmelidir.
  subroutine test_exact_singular_status()
    type(complex_linear_problem_t) :: problem
    type(complex_linear_solution_t) :: solution
    complex(dp) :: coefficient_matrix(2, 2)
    complex(dp) :: right_hand_side(2, 1)

    coefficient_matrix = cmplx(0.0_dp, 0.0_dp, kind=dp)
    coefficient_matrix(1, 1) = cmplx(1.0_dp, 1.0_dp, kind=dp)
    right_hand_side(:, 1) = [ &
      cmplx(1.0_dp, 0.0_dp, kind=dp), &
      cmplx(1.0_dp, 0.0_dp, kind=dp)]

    problem = create_complex_linear_problem( &
      coefficient_matrix, right_hand_side)
    solution = solve_complex_linear_problem(problem)

    if (get_complex_linear_solution_status(solution) /= &
        COMPLEX_SOLVE_SINGULAR) then
      error stop "Exact singular ZSYSVX problemi SINGULAR status üretmedi."
    end if
    if (has_complex_linear_response(solution)) then
      error stop "Exact singular complex solution fabricated response taşıyor."
    end if
    if (get_complex_linear_reciprocal_condition_number(solution) > 0.0_dp) then
      error stop "Exact singular complex solution RCOND=0 taşımıyor."
    end if
  end subroutine test_exact_singular_status

  !> Çalışma hassasiyetinde deterministic kötü koşullu diagonal problemi
  !! doğrular. A=diag(1,epsilon^2) nonsingular, fakat RCOND machine epsilon'un
  !! çok altındadır; ZSYSVX INFO=N+1 üretirken X ve error bounds geçerlidir.
  !! B=[1,epsilon^2]^T için exact X=[1,1]^T beklenir. Exact RCOND sayısı vendor
  !! bağımlı zorlanmaz; yalnız finite, pozitif ve epsilon altı olması sınanır.
  subroutine test_deterministic_ill_conditioned_status()
    type(complex_linear_problem_t) :: problem
    type(complex_linear_solution_t) :: solution
    complex(dp) :: coefficient_matrix(2, 2)
    complex(dp) :: expected_response(2, 1)
    complex(dp) :: right_hand_side(2, 1)
    complex(dp), allocatable :: response(:, :)
    real(dp) :: small_diagonal
    real(dp) :: rcond

    small_diagonal = epsilon(1.0_dp)**2
    coefficient_matrix = cmplx(0.0_dp, 0.0_dp, kind=dp)
    coefficient_matrix(1, 1) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    coefficient_matrix(2, 2) = cmplx(small_diagonal, 0.0_dp, kind=dp)
    expected_response(:, 1) = [ &
      cmplx(1.0_dp, 0.0_dp, kind=dp), &
      cmplx(1.0_dp, 0.0_dp, kind=dp)]
    right_hand_side(:, 1) = [ &
      cmplx(1.0_dp, 0.0_dp, kind=dp), &
      cmplx(small_diagonal, 0.0_dp, kind=dp)]

    problem = create_complex_linear_problem( &
      coefficient_matrix, right_hand_side)
    solution = solve_complex_linear_problem(problem)
    if (get_complex_linear_solution_status(solution) /= &
        COMPLEX_SOLVE_SOLVED_ILL_CONDITIONED) then
      error stop "epsilon^2 diagonal problem ILL_CONDITIONED status üretmedi."
    end if
    if (.not. has_complex_linear_response(solution)) then
      error stop "ILL_CONDITIONED sonuç hesaplanmış response'u kaybetti."
    end if

    rcond = get_complex_linear_reciprocal_condition_number(solution)
    if (.not. ieee_is_finite(rcond) .or. rcond <= 0.0_dp .or. &
        rcond >= epsilon(1.0_dp)) then
      error stop "ILL_CONDITIONED RCOND finite, pozitif ve epsilon altı değil."
    end if
    response = get_complex_linear_response(solution)
    call assert_complex_matrix_close( &
      response, expected_response, solution_tolerance, &
      "ILL_CONDITIONED problem hesaplanmış çözümü korumadı.")
    call assert_available_diagnostics(solution, 1)
  end subroutine test_deterministic_ill_conditioned_status

  !> Complex problem validation ve status-gated erişim hata yollarını üretim
  !! API'si üzerinden çalıştırır. NaN/Inf değerleri reel veya sanal bileşende
  !! sınanır; bütün vakalar solver backend'ine geçmeden clean diagnostic veya
  !! singular response getter'ında bilinçli API misuse tanısı üretmelidir.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(complex_linear_problem_t) :: problem
    type(complex_linear_solution_t) :: solution
    complex(dp), allocatable :: unused_response(:, :)
    complex(dp) :: matrix(2, 2)
    complex(dp) :: rhs(2, 1)
    real(dp) :: nan_value
    real(dp) :: infinity_value

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity_value = ieee_value(0.0_dp, ieee_positive_inf)
    matrix = complex_identity_matrix(2)
    rhs = cmplx(1.0_dp, 0.0_dp, kind=dp)

    select case (case_name)
      case ("zero_order")
        problem = create_complex_linear_problem( &
          reshape([complex(dp) ::], [0, 0]), &
          reshape([complex(dp) ::], [0, 1]))
      case ("zero_rhs")
        problem = create_complex_linear_problem( &
          complex_identity_matrix(2), reshape([complex(dp) ::], [2, 0]))
      case ("nonsquare_coefficient")
        problem = create_complex_linear_problem( &
          reshape([ &
            cmplx(1.0_dp, 0.0_dp, kind=dp), &
            cmplx(0.0_dp, 0.0_dp, kind=dp), &
            cmplx(0.0_dp, 0.0_dp, kind=dp), &
            cmplx(1.0_dp, 0.0_dp, kind=dp), &
            cmplx(0.0_dp, 0.0_dp, kind=dp), &
            cmplx(0.0_dp, 0.0_dp, kind=dp)], [2, 3]), rhs)
      case ("rhs_dimension_mismatch")
        problem = create_complex_linear_problem( &
          complex_identity_matrix(2), &
          reshape([ &
            cmplx(1.0_dp, 0.0_dp, kind=dp), &
            cmplx(1.0_dp, 0.0_dp, kind=dp), &
            cmplx(1.0_dp, 0.0_dp, kind=dp)], [3, 1]))
      case ("nonsymmetric_coefficient")
        matrix(1, 2) = cmplx(2.0_dp, 1.0_dp, kind=dp)
        problem = create_complex_linear_problem(matrix, rhs)
      case ("nan_coefficient_real")
        matrix(1, 1) = cmplx(nan_value, 0.0_dp, kind=dp)
        problem = create_complex_linear_problem(matrix, rhs)
      case ("infinite_coefficient_imaginary")
        matrix(1, 1) = cmplx(1.0_dp, infinity_value, kind=dp)
        problem = create_complex_linear_problem(matrix, rhs)
      case ("nan_rhs_imaginary")
        rhs(1, 1) = cmplx(1.0_dp, nan_value, kind=dp)
        problem = create_complex_linear_problem(matrix, rhs)
      case ("infinite_rhs_real")
        rhs(1, 1) = cmplx(infinity_value, 0.0_dp, kind=dp)
        problem = create_complex_linear_problem(matrix, rhs)
      case ("singular_response_access")
        matrix = cmplx(0.0_dp, 0.0_dp, kind=dp)
        matrix(1, 1) = cmplx(1.0_dp, 1.0_dp, kind=dp)
        problem = create_complex_linear_problem(matrix, rhs)
        solution = solve_complex_linear_problem(problem)
        unused_response = get_complex_linear_response(solution)
      case default
        print *, "Bilinmeyen complex solver doğrulama selector'ı: ", case_name
        return
    end select
  end subroutine exercise_invalid_case

  !> Çözümü bulunan status için RCOND, FERR, BERR ve independent residual
  !! tanılarının boyut, sonluluk ve işaret koşullarını doğrular. Exact FERR/BERR
  !! sayıları LAPACK vendor'ları arasında zorlanmaz.
  subroutine assert_available_diagnostics(solution, expected_rhs_count)
    type(complex_linear_solution_t), intent(in) :: solution
    integer, intent(in) :: expected_rhs_count

    real(dp), allocatable :: backward_errors(:)
    real(dp), allocatable :: forward_error_bounds(:)
    real(dp), allocatable :: relative_residuals(:)
    real(dp) :: rcond

    rcond = get_complex_linear_reciprocal_condition_number(solution)
    forward_error_bounds = get_complex_linear_forward_error_bounds(solution)
    backward_errors = get_complex_linear_backward_errors(solution)
    relative_residuals = get_complex_linear_relative_residuals(solution)

    if (.not. ieee_is_finite(rcond) .or. rcond < 0.0_dp) then
      error stop "Complex solution RCOND sonlu ve negatif olmayan değil."
    end if
    if (size(forward_error_bounds) /= expected_rhs_count .or. &
        size(backward_errors) /= expected_rhs_count .or. &
        size(relative_residuals) /= expected_rhs_count) then
      error stop "Complex per-RHS diagnostic dizilerinin boyutu yanlış."
    end if
    if (.not. all(ieee_is_finite(forward_error_bounds)) .or. &
        .not. all(ieee_is_finite(backward_errors)) .or. &
        .not. all(ieee_is_finite(relative_residuals)) .or. &
        any(forward_error_bounds < 0.0_dp) .or. &
        any(backward_errors < 0.0_dp) .or. &
        any(relative_residuals < 0.0_dp)) then
      error stop "Complex solution diagnostics finite/nonnegative değil."
    end if
    if (any(relative_residuals > residual_tolerance)) then
      error stop "Backend-independent complex relative residual toleransı aştı."
    end if
  end subroutine assert_available_diagnostics

  !> İki complex matrisin ölçeğe göre göreli farkını sınar. Girdi matrisler
  !! aynı fiziksel birimde, tolerance boyutsuzdur.
  subroutine assert_complex_matrix_close( &
      actual, expected, tolerance, message)
    complex(dp), intent(in) :: actual(:, :)
    complex(dp), intent(in) :: expected(:, :)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    real(dp) :: difference_scale
    real(dp) :: reference_scale

    if (any(shape(actual) /= shape(expected))) then
      error stop message
    end if
    difference_scale = maxval(abs(actual - expected))
    reference_scale = max(1.0_dp, maxval(abs(expected)))
    if (difference_scale > tolerance * reference_scale) then
      error stop message
    end if
  end subroutine assert_complex_matrix_close

  !> n x n complex birim matrisi üretir. Katsayılar boyutsuzdur.
  pure function complex_identity_matrix(order) result(matrix)
    integer, intent(in) :: order
    complex(dp) :: matrix(order, order)

    integer :: index

    matrix = cmplx(0.0_dp, 0.0_dp, kind=dp)
    do index = 1, order
      matrix(index, index) = cmplx(1.0_dp, 0.0_dp, kind=dp)
    end do
  end function complex_identity_matrix

end program test_complex_linear_solver
