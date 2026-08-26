program test_generalized_eigen_solver
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
    ieee_positive_inf, ieee_negative_inf, ieee_is_finite
  use tms_kinds, only : dp
  use tms_generalized_eigen_problem, only : generalized_eigen_problem_t, &
    create_generalized_eigen_problem, get_generalized_eigen_stiffness, &
    get_generalized_eigen_mass
  use tms_eigen_solution, only : eigen_solution_t, initialize_eigen_solution, &
    get_eigen_mode_count, get_eigenvalues, get_eigenvectors, &
    get_eigen_backend_identity
  use tms_generalized_eigen_solver, only : solve_generalized_eigen_problem
  implicit none

  real(dp), parameter :: relative_tolerance = 1.0e-10_dp
  real(dp), parameter :: residual_tolerance = 1.0e-11_dp
  character(len=64) :: validation_case

  ! Geçersiz generalized problem vakaları ayrı CTest süreçlerinde error stop
  ! üretmelidir. Bilinmeyen selector normal döner ve WILL_FAIL testini bozar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call test_singular_stiffness_and_ordering()
  call test_repeated_eigenvalue_subspace()
  call test_input_matrix_immutability()
  call test_partial_eigen_solution_contract()

  print *, "V0.5 generalized eigen solver doğrulamaları başarılı."

contains

  !> Singular symmetric positive-semidefinite K matrisinin fiziksel bir
  !! rigid-body mode olarak kabul edildiğini ve elastic eigenpair'in doğru
  !! sıralandığını doğrular. Model K=k[[1,-1],[-1,1]], M=diag(J1,J2) olup
  !! beklenen lambda=[0,k(1/J1+1/J2)] [1/s^2]'dir. Test üretim solver'ını
  !! çağırır; eigenvalue sırasını, mass normalization ve residual'ı sınar.
  subroutine test_singular_stiffness_and_ordering()
    real(dp), parameter :: stiffness_value = 10.0_dp
    real(dp), parameter :: inertia_one = 2.0_dp
    real(dp), parameter :: inertia_two = 3.0_dp
    type(generalized_eigen_problem_t) :: problem
    type(eigen_solution_t) :: solution
    real(dp) :: stiffness(2, 2)
    real(dp) :: mass(2, 2)
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: eigenvectors(:, :)
    character(len=:), allocatable :: backend
    integer :: mode_index

    stiffness = reshape([ &
      stiffness_value, -stiffness_value, &
      -stiffness_value, stiffness_value], [2, 2])
    mass = reshape([inertia_one, 0.0_dp, 0.0_dp, inertia_two], [2, 2])

    problem = create_generalized_eigen_problem(stiffness, mass)
    solution = solve_generalized_eigen_problem(problem)
    eigenvalues = get_eigenvalues(solution)
    eigenvectors = get_eigenvectors(solution)
    backend = get_eigen_backend_identity(solution)

    if (get_eigen_mode_count(solution) /= 2) then
      error stop "İki DOF generalized problem iki eigenpair üretmedi."
    end if
    if (index(backend, "DSYGV") == 0) then
      error stop "Eigen solution DSYGV reference backend kimliğini taşımıyor."
    end if

    call assert_absolute_close( &
      eigenvalues(1), 0.0_dp, 1.0e-12_dp, &
      "Singular K rigid-body eigenvalue üretmedi.")
    call assert_relative_close( &
      eigenvalues(2), stiffness_value * &
        (1.0_dp / inertia_one + 1.0_dp / inertia_two), &
      relative_tolerance, "İki ataletli elastic eigenvalue yanlış.")

    do mode_index = 1, 2
      call assert_mass_normalized( &
        eigenvectors(:, mode_index), mass, relative_tolerance)
      call assert_relative_residual( &
        stiffness, mass, eigenvalues(mode_index), &
        eigenvectors(:, mode_index), residual_tolerance)
    end do
    call assert_modal_stiffness_diagonalization( &
      stiffness, eigenvalues, eigenvectors, relative_tolerance)
  end subroutine test_singular_stiffness_and_ordering

  !> İki repeated ve bir ayrı eigenvalue içeren sentetik SPD-M problemini
  !! doğrular. Repeated lambda için individual eigenvector unique değildir;
  !! test exact vektör kıyaslamaz. Bunun yerine multiplicity, residual,
  !! M-orthogonality ve beklenen ilk-iki-koordinat eigenspace overlap ölçütünü
  !! kullanır. K/M birimleri sırasıyla [N*m/rad] ve [kg*m^2]'dir.
  subroutine test_repeated_eigenvalue_subspace()
    type(generalized_eigen_problem_t) :: problem
    type(eigen_solution_t) :: solution
    real(dp) :: stiffness(3, 3)
    real(dp) :: mass(3, 3)
    real(dp) :: expected_basis(3, 2)
    real(dp) :: overlap(2, 2)
    real(dp) :: gram(3, 3)
    real(dp), allocatable :: eigenvalues(:)
    real(dp), allocatable :: eigenvectors(:, :)
    integer :: mode_index

    stiffness = diagonal_matrix([2.0_dp, 2.0_dp, 9.0_dp])
    mass = diagonal_matrix([1.0_dp, 1.0_dp, 3.0_dp])
    problem = create_generalized_eigen_problem(stiffness, mass)
    solution = solve_generalized_eigen_problem(problem)
    eigenvalues = get_eigenvalues(solution)
    eigenvectors = get_eigenvectors(solution)

    call assert_relative_close( &
      eigenvalues(1), 2.0_dp, relative_tolerance, &
      "Repeated ilk eigenvalue yanlış.")
    call assert_relative_close( &
      eigenvalues(2), 2.0_dp, relative_tolerance, &
      "Repeated ikinci eigenvalue yanlış.")
    call assert_relative_close( &
      eigenvalues(3), 3.0_dp, relative_tolerance, &
      "Repeated problem ayrı eigenvalue yanlış.")

    gram = matmul(transpose(eigenvectors), matmul(mass, eigenvectors))
    call assert_matrix_close( &
      gram, identity_matrix(3), relative_tolerance, &
      "Repeated problem mode'ları M-orthonormal değil.")

    expected_basis = 0.0_dp
    expected_basis(1, 1) = 1.0_dp
    expected_basis(2, 2) = 1.0_dp
    overlap = matmul( &
      transpose(expected_basis), matmul(mass, eigenvectors(:, 1:2)))
    call assert_matrix_close( &
      matmul(transpose(overlap), overlap), identity_matrix(2), &
      relative_tolerance, &
      "Repeated eigenvalue için hesaplanan eigenspace eşdeğer değil.")

    do mode_index = 1, 3
      call assert_relative_residual( &
        stiffness, mass, eigenvalues(mode_index), &
        eigenvectors(:, mode_index), residual_tolerance)
    end do
    call assert_modal_stiffness_diagonalization( &
      stiffness, eigenvalues, eigenvectors, relative_tolerance)
  end subroutine test_repeated_eigenvalue_subspace

  !> LAPACK DSYGV'nin A/B matrislerini overwrite etmesine rağmen problem
  !! contract'ındaki original K/M kopyalarının değişmediğini doğrular.
  !! Karşılaştırma toleranssız exact yapılır; solver yalnız backend çalışma
  !! kopyalarını değiştirebilir.
  subroutine test_input_matrix_immutability()
    type(generalized_eigen_problem_t) :: problem
    type(eigen_solution_t) :: unused_solution
    real(dp) :: stiffness(2, 2)
    real(dp) :: mass(2, 2)
    real(dp), allocatable :: stiffness_before(:, :)
    real(dp), allocatable :: mass_before(:, :)
    real(dp), allocatable :: stiffness_after(:, :)
    real(dp), allocatable :: mass_after(:, :)

    stiffness = reshape([4.0_dp, -1.0_dp, -1.0_dp, 2.0_dp], [2, 2])
    mass = diagonal_matrix([0.5_dp, 1.5_dp])
    problem = create_generalized_eigen_problem(stiffness, mass)
    stiffness_before = get_generalized_eigen_stiffness(problem)
    mass_before = get_generalized_eigen_mass(problem)

    unused_solution = solve_generalized_eigen_problem(problem)

    stiffness_after = get_generalized_eigen_stiffness(problem)
    mass_after = get_generalized_eigen_mass(problem)
    if (any(abs(stiffness_after - stiffness_before) > 0.0_dp) .or. &
        any(abs(mass_after - mass_before) > 0.0_dp) .or. &
        any(abs(stiffness - stiffness_before) > 0.0_dp) .or. &
        any(abs(mass - mass_before) > 0.0_dp)) then
      error stop "DSYGV original K/M matrislerini değiştirdi."
    end if
  end subroutine test_input_matrix_immutability

  !> Backend-neutral solution contract'ının gelecekteki partial-spectrum
  !! backend'ler için 1<=m<=n mode kabul ettiğini doğrular. Bu test DSYGV'nin
  !! V0.5'te bütün spectrum'u çözme davranışını değiştirmez.
  subroutine test_partial_eigen_solution_contract()
    type(eigen_solution_t) :: solution
    real(dp) :: eigenvectors(3, 2)

    eigenvectors = reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp], [3, 2])
    call initialize_eigen_solution( &
      solution, [1.0_dp, 2.0_dp], eigenvectors, "synthetic partial backend")
    if (get_eigen_mode_count(solution) /= 2) then
      error stop "Partial-spectrum eigen solution mode sayısını korumadı."
    end if
  end subroutine test_partial_eigen_solution_contract

  !> Generalized problem ve backend önkoşullarının geçersiz sentetik
  !! matrisleri reddettiğini çalıştırır. K [N*m/rad], M [kg*m^2] ve bütün
  !! scalar özel değerler IEEE sınıflarındadır. Non-SPD M vakaları DSYGV INFO
  !! çevirisine ulaşmalı; structural vakalar backend çağrısından önce durmalıdır.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(generalized_eigen_problem_t) :: problem
    type(eigen_solution_t) :: solution
    real(dp) :: nan_value
    real(dp) :: positive_infinity
    real(dp) :: negative_infinity

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    positive_infinity = ieee_value(0.0_dp, ieee_positive_inf)
    negative_infinity = ieee_value(0.0_dp, ieee_negative_inf)

    select case (case_name)
      case ("nonsquare_stiffness")
        problem = create_generalized_eigen_problem( &
          reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
            [2, 3]), identity_matrix(2))
      case ("nonsquare_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), &
          reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
            [2, 3]))
      case ("dimension_mismatch")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), identity_matrix(3))
      case ("nonsymmetric_stiffness")
        problem = create_generalized_eigen_problem( &
          reshape([1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp], [2, 2]), &
          identity_matrix(2))
      case ("nonsymmetric_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), &
          reshape([1.0_dp, 0.0_dp, 2.0_dp, 1.0_dp], [2, 2]))
      case ("nan_stiffness")
        problem = create_generalized_eigen_problem( &
          diagonal_matrix([nan_value, 1.0_dp]), identity_matrix(2))
      case ("positive_infinity_stiffness")
        problem = create_generalized_eigen_problem( &
          diagonal_matrix([positive_infinity, 1.0_dp]), identity_matrix(2))
      case ("negative_infinity_stiffness")
        problem = create_generalized_eigen_problem( &
          diagonal_matrix([negative_infinity, 1.0_dp]), identity_matrix(2))
      case ("nan_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), diagonal_matrix([nan_value, 1.0_dp]))
      case ("positive_infinity_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), &
          diagonal_matrix([positive_infinity, 1.0_dp]))
      case ("negative_infinity_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), &
          diagonal_matrix([negative_infinity, 1.0_dp]))
      case ("zero_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
            [2, 2]))
        solution = solve_generalized_eigen_problem(problem)
      case ("negative_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), diagonal_matrix([-1.0_dp, 1.0_dp]))
        solution = solve_generalized_eigen_problem(problem)
      case ("indefinite_mass")
        problem = create_generalized_eigen_problem( &
          identity_matrix(2), &
          reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2]))
        solution = solve_generalized_eigen_problem(problem)
      case ("zero_order")
        problem = create_generalized_eigen_problem( &
          reshape([real(dp) ::], [0, 0]), reshape([real(dp) ::], [0, 0]))
        solution = solve_generalized_eigen_problem(problem)
      case ("excessive_mode_count")
        call initialize_eigen_solution( &
          solution, [1.0_dp, 2.0_dp, 3.0_dp], &
          reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
            [2, 3]), "invalid synthetic backend")
      case default
        print *, "Bilinmeyen eigen doğrulama selector'ı: ", case_name
        return
    end select

    print *, "Geçersiz generalized eigenproblem beklenmedik biçimde kabul edildi."
  end subroutine exercise_invalid_case

  !> Pozitif skaler sonuçları boyutsuz bağıl toleransla karşılaştırır.
  pure subroutine assert_relative_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp .or. &
        abs(expected) <= tiny(1.0_dp)) then
      error stop "Bağıl eigen assertion girdileri geçersiz."
    end if
    if (abs(actual - expected) / abs(expected) > tolerance) error stop message
  end subroutine assert_relative_close

  !> Aynı birimdeki iki skaler sonucu mutlak toleransla karşılaştırır.
  pure subroutine assert_absolute_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (.not. ieee_is_finite(actual) .or. &
        .not. ieee_is_finite(expected) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) then
      error stop "Mutlak eigen assertion girdileri geçersiz."
    end if
    if (abs(actual - expected) > tolerance) error stop message
  end subroutine assert_absolute_close

  !> İki sonlu matrisin maksimum katsayı farkını mutlak toleransla sınar.
  pure subroutine assert_matrix_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:, :)
    real(dp), intent(in) :: expected(:, :)
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: message

    if (any(shape(actual) /= shape(expected))) error stop message
    if (.not. all(ieee_is_finite(actual)) .or. &
        .not. all(ieee_is_finite(expected)) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) then
      error stop "Matris eigen assertion girdileri geçersiz."
    end if
    if (maxval(abs(actual - expected)) > tolerance) error stop message
  end subroutine assert_matrix_close

  !> Bir mode shape'in phi^T M phi=1 kütle normalizasyonunu doğrular.
  pure subroutine assert_mass_normalized(mode, mass, tolerance)
    real(dp), intent(in) :: mode(:)
    real(dp), intent(in) :: mass(:, :)
    real(dp), intent(in) :: tolerance
    real(dp) :: mass_norm

    mass_norm = dot_product(mode, matmul(mass, mode))
    call assert_absolute_close( &
      mass_norm, 1.0_dp, tolerance, "DSYGV mode shape M-normalized değil.")
  end subroutine assert_mass_normalized

  !> M-orthonormal modal bazın rijitlik metriğinde özdeğerleri köşegenleştirdiğini
  !! doğrular. Matematiksel model Phi^T*K*Phi=diag(lambda) olup K [N*m/rad],
  !! lambda [1/s^2] ve Phi mass-normalized modal genliklerdir. Bu bağımsız test,
  !! üretim residual yordamını tekrar kullanmadan mode eşleşmesi, orthogonality
  !! ve repeated eigenspace davranışını birlikte sınar.
  pure subroutine assert_modal_stiffness_diagonalization( &
      stiffness, eigenvalues, eigenvectors, tolerance)
    real(dp), intent(in) :: stiffness(:, :)
    real(dp), intent(in) :: eigenvalues(:)
    real(dp), intent(in) :: eigenvectors(:, :)
    real(dp), intent(in) :: tolerance

    real(dp) :: expected(size(eigenvalues), size(eigenvalues))
    real(dp) :: projected(size(eigenvalues), size(eigenvalues))
    real(dp) :: reference_scale
    integer :: mode_index

    if (size(stiffness, 1) /= size(stiffness, 2) .or. &
        size(eigenvectors, 1) /= size(stiffness, 1) .or. &
        size(eigenvectors, 2) /= size(eigenvalues) .or. &
        .not. all(ieee_is_finite(stiffness)) .or. &
        .not. all(ieee_is_finite(eigenvalues)) .or. &
        .not. all(ieee_is_finite(eigenvectors)) .or. &
        .not. ieee_is_finite(tolerance) .or. tolerance < 0.0_dp) then
      error stop "Modal rijitlik köşegenleştirme assertion girdileri geçersiz."
    end if

    expected = 0.0_dp
    do mode_index = 1, size(eigenvalues)
      expected(mode_index, mode_index) = eigenvalues(mode_index)
    end do
    projected = matmul( &
      transpose(eigenvectors), matmul(stiffness, eigenvectors))
    reference_scale = max(1.0_dp, maxval(abs(eigenvalues)))
    if (maxval(abs(projected - expected)) / reference_scale > tolerance) then
      error stop "Modal baz Phi^T*K*Phi=diag(lambda) koşulunu sağlamıyor."
    end if
  end subroutine assert_modal_stiffness_diagonalization

  !> Dimensionless scaled eigenpair residual değerini bağımsız hesaplar.
  !! rho=||K phi-lambda M phi||2 /
  !! ((||K||inf+|lambda| ||M||inf)||phi||2). Test girdileri finite, M/K
  !! boyutları uyumlu ve payda pozitiftir.
  pure subroutine assert_relative_residual( &
      stiffness, mass, eigenvalue, mode, tolerance)
    real(dp), intent(in) :: stiffness(:, :)
    real(dp), intent(in) :: mass(:, :)
    real(dp), intent(in) :: eigenvalue
    real(dp), intent(in) :: mode(:)
    real(dp), intent(in) :: tolerance
    real(dp) :: denominator
    real(dp) :: residual

    denominator = (matrix_infinity_norm(stiffness) + &
      abs(eigenvalue) * matrix_infinity_norm(mass)) * vector_two_norm(mode)
    if (denominator <= 0.0_dp) then
      error stop "Eigen residual testi pozitif payda gerektirir."
    end if
    residual = vector_two_norm( &
      matmul(stiffness, mode) - eigenvalue * matmul(mass, mode)) / denominator
    if (.not. ieee_is_finite(residual) .or. residual > tolerance) then
      error stop "Generalized eigenpair relative residual sınırı aşıldı."
    end if
  end subroutine assert_relative_residual

  !> Kare matris için maksimum mutlak satır toplamı normunu döndürür.
  pure function matrix_infinity_norm(matrix) result(norm_value)
    real(dp), intent(in) :: matrix(:, :)
    real(dp) :: norm_value

    norm_value = maxval(sum(abs(matrix), dim=2))
  end function matrix_infinity_norm

  !> Vektörün Euclidean iki normunu ölçeklenmiş intrinsic NORM2 ile döndürür.
  pure function vector_two_norm(vector) result(norm_value)
    real(dp), intent(in) :: vector(:)
    real(dp) :: norm_value

    norm_value = norm2(vector)
  end function vector_two_norm

  !> Verilen köşegen katsayılarından n x n diagonal matris oluşturur.
  pure function diagonal_matrix(diagonal) result(matrix)
    real(dp), intent(in) :: diagonal(:)
    real(dp) :: matrix(size(diagonal), size(diagonal))
    integer :: index

    matrix = 0.0_dp
    do index = 1, size(diagonal)
      matrix(index, index) = diagonal(index)
    end do
  end function diagonal_matrix

  !> Verilen boyutta boyutsuz birim matris oluşturur.
  pure function identity_matrix(order) result(matrix)
    integer, intent(in) :: order
    real(dp) :: matrix(order, order)
    integer :: index

    matrix = 0.0_dp
    do index = 1, order
      matrix(index, index) = 1.0_dp
    end do
  end function identity_matrix

end program test_generalized_eigen_solver
