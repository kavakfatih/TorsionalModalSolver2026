module tms_modal_analysis
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_constants, only : pi
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    get_stiffness_matrix_values
  use tms_mass_matrix, only : mass_matrix_t, get_mass_matrix_values
  use tms_reduced_system, only : reduced_torsional_system_t, &
    get_reduced_stiffness, get_reduced_mass, &
    get_reduced_active_dof_count, recover_mode_shape
  use tms_generalized_eigen_problem, only : generalized_eigen_problem_t, &
    create_generalized_eigen_problem, get_generalized_eigen_problem_order, &
    get_generalized_eigen_stiffness, get_generalized_eigen_mass
  use tms_eigen_solution, only : eigen_solution_t, get_eigen_mode_count, &
    get_eigenvalues, get_eigenvectors, get_eigen_backend_identity
  use tms_generalized_eigen_solver, only : solve_generalized_eigen_problem
  use tms_modal_validation, only : AUTO_RIGID_TOLERANCE_MULTIPLIER, &
    calculate_auto_rigid_eigenvalue_tolerance, &
    is_rigid_eigenvalue, &
    normalize_and_validate_mass_modes, &
    calculate_relative_eigenpair_residuals, &
    calculate_mass_orthogonality_error
  use tms_modal_result, only : modal_result_t, create_modal_result, &
    RIGID_MODE, ELASTIC_MODE
  implicit none
  private

  !> Modal analizde backend'den bağımsız ortak seçimleri taşır.
  !!
  !! requested_mode_count boyutsuzdur; `0` bütün özçiftlerin döndürülmesini,
  !! pozitif değer ise artan özdeğer sırasındaki ilk N mode'u ister. DSYGV V0.5
  !! backend'i yine bütün özçiftleri çözer; seçim ortak post-processing katmanında
  !! yapılır. auto_rigid_tolerance_multiplier, sabit bir Hz eşiği değil,
  !! machine epsilon ve problem özdeğer ölçeğiyle kullanılan boyutsuz katsayıdır.
  type, public :: modal_analysis_options_t
    integer :: requested_mode_count = 0
    real(dp) :: auto_rigid_tolerance_multiplier = &
      AUTO_RIGID_TOLERANCE_MULTIPLIER
  end type modal_analysis_options_t

  public :: analyze_reduced_torsional_system
  public :: AUTO_RIGID_TOLERANCE_MULTIPLIER

contains

  !> Constraint-aware reduced torsional sistemi genel reel simetrik özdeğer
  !! problemi olarak çözer ve fiziksel mode shape'leri geri kazanır.
  !!
  !! Fiziksel açıklama: Lineer, sönümsüz ve frozen-property sistem için
  !! K_r*phi=lambda*M_r*phi çözülür. lambda=omega^2 ve f=omega/(2*pi)
  !! bağıntıları kullanılır. Rijit-cisim modları fiziksel hata değildir.
  !! Matematiksel açıklama: K_r simetrik ve positive-semidefinite olabilir;
  !! M_r simetrik positive-definite olmalıdır. Özvektörler M-normalize edilir,
  !! boyutsuz göreli residual ve ||Phi^T*M*Phi-I||_F hesaplanır. Reduced her
  !! sütun `recover_mode_shape` ile phi=P*phi_r biçiminde physical DOF uzayına
  !! açılır; prescribed q_p modal şekle eklenmez.
  !! Girdi: K_r [N*m/rad], M_r [kg*m^2] ve recovery eşlemesini taşıyan reduced
  !! sistem; opsiyonel boyutsuz mode sayısı ve AUTO tolerans katsayısı.
  !! Çıktı: lambda [1/s^2], omega [rad/s], f [Hz], sınıf, residual, reduced ve
  !! physical şekiller ile backend metadata'sını taşıyan modal_result_t.
  !! Varsayımlar ve geçerlilik: En az bir active DOF gerekir. Significantly
  !! negative lambda kararsız/geçersiz model olarak reddedilir; AUTO tolerans
  !! içindeki küçük negatif lambda rijit mode olarak sıfır frekansa eşlenir.
  function analyze_reduced_torsional_system(reduced_system, options) &
      result(modal_result)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    type(modal_analysis_options_t), intent(in), optional :: options
    type(modal_result_t) :: modal_result

    type(modal_analysis_options_t) :: effective_options
    type(stiffness_matrix_t) :: reduced_stiffness
    type(mass_matrix_t) :: reduced_mass
    type(generalized_eigen_problem_t) :: eigenproblem
    type(eigen_solution_t) :: eigen_solution
    character(len=:), allocatable :: backend_identity
    real(dp), allocatable :: all_eigenvalues(:)
    real(dp), allocatable :: all_mode_shapes(:, :)
    real(dp), allocatable :: angular_frequencies_rad_s(:)
    real(dp), allocatable :: eigenproblem_mass_after(:, :)
    real(dp), allocatable :: eigenproblem_stiffness_after(:, :)
    real(dp), allocatable :: frequencies_hz(:)
    real(dp), allocatable :: mass_values(:, :)
    real(dp), allocatable :: physical_mode_shapes(:, :)
    real(dp), allocatable :: recovered_mode(:)
    real(dp), allocatable :: reduced_mode_shapes(:, :)
    real(dp), allocatable :: relative_residuals(:)
    real(dp), allocatable :: selected_eigenvalues(:)
    real(dp), allocatable :: stiffness_values(:, :)
    integer, allocatable :: mode_classifications(:)
    integer :: active_dof_count
    integer :: available_mode_count
    integer :: mode_count
    integer :: mode_index
    real(dp) :: mass_orthogonality_error
    real(dp) :: rigid_eigenvalue_tolerance

    effective_options = modal_analysis_options_t()
    if (present(options)) effective_options = options
    call validate_modal_analysis_options(effective_options)

    active_dof_count = get_reduced_active_dof_count(reduced_system)
    if (active_dof_count == 0) then
      error stop &
        "Modal analiz için aktif DOF yoktur; 0x0 sistem DSYGV'ye gönderilemez."
    end if

    reduced_stiffness = get_reduced_stiffness(reduced_system)
    reduced_mass = get_reduced_mass(reduced_system)
    stiffness_values = get_stiffness_matrix_values(reduced_stiffness)
    mass_values = get_mass_matrix_values(reduced_mass)

    ! LAPACK K ve M çalışma dizilerini overwrite eder. Eigen problem factory'si
    ! bağımsız kopya oluşturur; aşağıdaki ilk diziler residual ve doğrulama için
    ! authoritative original Kr/Mr olarak korunur.
    eigenproblem = create_generalized_eigen_problem( &
      stiffness_values, mass_values)
    if (get_generalized_eigen_problem_order(eigenproblem) /= &
        active_dof_count) then
      error stop "Genelleştirilmiş özdeğer problem boyutu active DOF ile uyumsuz."
    end if

    eigen_solution = solve_generalized_eigen_problem(eigenproblem)

    ! Backend'in çalışma kopyası sınırını ihlal ederek problem girdilerini
    ! değiştirmesi sessizce kabul edilmez. Karşılaştırma tamdır; çözüm işlemi
    ! original Kr/Mr üzerinde hiçbir aritmetik güncelleme yapmamalıdır.
    eigenproblem_stiffness_after = &
      get_generalized_eigen_stiffness(eigenproblem)
    eigenproblem_mass_after = get_generalized_eigen_mass(eigenproblem)
    if (any(abs(eigenproblem_stiffness_after - stiffness_values) > 0.0_dp) .or. &
        any(abs(eigenproblem_mass_after - mass_values) > 0.0_dp)) then
      error stop "Özdeğer backend'i original K/M problem girdilerini değiştirdi."
    end if

    available_mode_count = get_eigen_mode_count(eigen_solution)
    all_eigenvalues = get_eigenvalues(eigen_solution)
    all_mode_shapes = get_eigenvectors(eigen_solution)
    backend_identity = get_eigen_backend_identity(eigen_solution)

    if (available_mode_count <= 0 .or. &
        available_mode_count > active_dof_count .or. &
        size(all_eigenvalues) /= available_mode_count .or. &
        size(all_mode_shapes, 1) /= active_dof_count .or. &
        size(all_mode_shapes, 2) /= available_mode_count) then
      error stop "Özdeğer çözümü active DOF sınırında tutarlı özçiftler vermedi."
    end if
    if (.not. all(ieee_is_finite(all_eigenvalues)) .or. &
        .not. all(ieee_is_finite(all_mode_shapes))) then
      error stop "Özdeğer çözümü yalnız sonlu lambda ve mode shape üretmelidir."
    end if

    ! ITYPE=1 DSYGV normal olarak B=M-normalize vektör verir. Önce bu özellik
    ! doğrulanır; yalnız makine-hassasiyeti üstündeki diagonal norm sapmalarında
    ! belgeli sayısal cleanup uygulanır ve global M-ortogonallik yeniden sınanır.
    call normalize_and_validate_mass_modes(mass_values, all_mode_shapes)

    rigid_eigenvalue_tolerance = &
      calculate_auto_rigid_eigenvalue_tolerance( &
      stiffness_values, mass_values, all_eigenvalues, &
      effective_options%auto_rigid_tolerance_multiplier)
    mode_count = available_mode_count
    if (effective_options%requested_mode_count > 0) then
      if (effective_options%requested_mode_count > available_mode_count) then
        error stop "İstenen mode sayısı backend'in sağladığı mode sayısını aşıyor."
      end if
      mode_count = effective_options%requested_mode_count
    end if

    selected_eigenvalues = all_eigenvalues(1:mode_count)
    reduced_mode_shapes = all_mode_shapes(:, 1:mode_count)
    allocate(angular_frequencies_rad_s(mode_count))
    allocate(frequencies_hz(mode_count))
    allocate(mode_classifications(mode_count))

    do mode_index = 1, mode_count
      if (is_rigid_eigenvalue( &
          selected_eigenvalues(mode_index), &
          rigid_eigenvalue_tolerance)) then
        mode_classifications(mode_index) = RIGID_MODE
        angular_frequencies_rad_s(mode_index) = 0.0_dp
        frequencies_hz(mode_index) = 0.0_dp
      else
        mode_classifications(mode_index) = ELASTIC_MODE
        angular_frequencies_rad_s(mode_index) = &
          sqrt(selected_eigenvalues(mode_index))
        frequencies_hz(mode_index) = &
          angular_frequencies_rad_s(mode_index) / (2.0_dp * pi)
      end if
    end do

    relative_residuals = calculate_relative_eigenpair_residuals( &
      stiffness_values, mass_values, selected_eigenvalues, &
      reduced_mode_shapes)
    mass_orthogonality_error = calculate_mass_orthogonality_error( &
      mass_values, reduced_mode_shapes)

    do mode_index = 1, mode_count
      recovered_mode = recover_mode_shape( &
        reduced_system, reduced_mode_shapes(:, mode_index))
      if (mode_index == 1) then
        allocate(physical_mode_shapes(size(recovered_mode), mode_count))
      end if
      physical_mode_shapes(:, mode_index) = recovered_mode
    end do

    modal_result = create_modal_result( &
      selected_eigenvalues, angular_frequencies_rad_s, frequencies_hz, &
      mode_classifications, relative_residuals, reduced_mode_shapes, &
      physical_mode_shapes, mass_orthogonality_error, &
      rigid_eigenvalue_tolerance, backend_identity)
  end function analyze_reduced_torsional_system

  !> Modal analiz seçeneklerinin backend'den bağımsız giriş sözleşmesini sınar.
  !! requested_mode_count boyutsuz ve negatif olmayan; AUTO rijit tolerans
  !! katsayısı sonlu, boyutsuz ve pozitif olmalıdır. Çıktı üretmez.
  pure subroutine validate_modal_analysis_options(options)
    type(modal_analysis_options_t), intent(in) :: options

    if (options%requested_mode_count < 0) then
      error stop "İstenen modal mode sayısı negatif olamaz; 0 bütün mode'lardır."
    end if
    if (.not. ieee_is_finite(options%auto_rigid_tolerance_multiplier) .or. &
        options%auto_rigid_tolerance_multiplier <= 0.0_dp) then
      error stop "AUTO rijit mode tolerans katsayısı sonlu ve pozitif olmalıdır."
    end if
  end subroutine validate_modal_analysis_options

end module tms_modal_analysis
