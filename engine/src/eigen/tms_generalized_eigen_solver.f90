module tms_generalized_eigen_solver
  use tms_generalized_eigen_problem, only : generalized_eigen_problem_t, &
    validate_generalized_eigen_problem, get_generalized_eigen_problem_order
  use tms_eigen_solution, only : eigen_solution_t, validate_eigen_solution
  use tms_lapack_dsygv_backend, only : solve_with_lapack_dsygv
  implicit none
  private

  public :: solve_generalized_eigen_problem

contains

  !> Backend-neutral generalized symmetric eigen solver facade'ıdır.
  !!
  !! Fiziksel açıklama: Reduced torsional K [N*m/rad] ve M [kg*m^2]
  !! matrislerinin K*phi=lambda*M*phi modal problemini çözer. Çıktı lambda
  !! [1/s^2] ile reduced mode shape'leri taşır; frequency conversion, rigid-mode
  !! classification, residual ve physical recovery modal katmanın görevidir.
  !! Matematiksel açıklama: K reel simetrik ve singular/PSD olabilir; M reel
  !! simetrik positive definite olmalıdır. V0.5 reference backend DSYGV,
  !! ITYPE=1 ile bütün eigenpair'leri artan sırada hesaplar.
  !! Girdi private storage taşıyan generalized_eigen_problem_t, çıktı
  !! backend-neutral eigen_solution_t değeridir. Sıfır aktif DOF clean diagnostic
  !! ile reddedilir ve LAPACK çağrılmaz. Bu facade external backend çağırdığı
  !! için PURE değildir; ileride sparse/Lanczos backend eklenmesi public modal
  !! analysis contract'ını değiştirmemelidir.
  function solve_generalized_eigen_problem(problem) result(solution)
    type(generalized_eigen_problem_t), intent(in) :: problem
    type(eigen_solution_t) :: solution

    call validate_generalized_eigen_problem(problem)
    if (get_generalized_eigen_problem_order(problem) == 0) then
      error stop &
        "Genelleştirilmiş modal çözüm için aktif serbestlik derecesi bulunmuyor."
    end if

    solution = solve_with_lapack_dsygv(problem)
    call validate_eigen_solution(solution)
  end function solve_generalized_eigen_problem

end module tms_generalized_eigen_solver
