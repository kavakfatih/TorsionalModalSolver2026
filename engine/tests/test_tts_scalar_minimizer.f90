program test_tts_scalar_minimizer
  use tms_kinds, only : dp
  use tms_scalar_minimizer, only : scalar_minimizer_result_t, &
    SCALAR_MINIMIZER_SUCCESS, SCALAR_MINIMIZER_INVALID_BRACKET, &
    minimize_scalar_brent, is_valid_scalar_minimum_bracket
  use tms_tts_test_support, only : assert_true, assert_close
  implicit none

  type(scalar_minimizer_result_t) :: result

  call assert_true(is_valid_scalar_minimum_bracket( &
    -1.0_dp, 1.0_dp, 3.0_dp, quadratic(-1.0_dp), quadratic(1.0_dp), &
    quadratic(3.0_dp)), "Geçerli Brent bracket tanınmadı.")
  result = minimize_scalar_brent(quadratic, -1.0_dp, 1.0_dp, 3.0_dp)
  call assert_true(result%status == SCALAR_MINIMIZER_SUCCESS .and. &
    result%converged, "Brent quadratic truth üzerinde yakınsamadı.")
  call assert_close(result%x_minimum, 1.25_dp, 5.0e-8_dp, &
    "Brent quadratic minimum koordinatı hatalı.")
  call assert_close(result%f_minimum, 2.0_dp, 5.0e-12_dp, &
    "Brent quadratic minimum objective değeri hatalı.")
  call assert_true(result%iteration_count > 0 .and. &
    result%function_evaluation_count >= 4, &
    "Brent iteration/evaluation tanısı tutulmadı.")

  result = minimize_scalar_brent(monotonic, -1.0_dp, 0.0_dp, 1.0_dp)
  call assert_true(result%status == SCALAR_MINIMIZER_INVALID_BRACKET .and. &
    .not. result%converged, &
    "Interior minimum içermeyen bracket başarı sayıldı.")

  print *, "V0.8.1 generic Brent scalar minimizer doğrulandı."

contains

  !> Bilinen minimumu x=1.25 olan boyutsuz quadratic numerical referanstır.
  pure function quadratic(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    value = (x - 1.25_dp)**2 + 2.0_dp
  end function quadratic

  !> Boundary minimum üreten monoton objective, invalid bracket regresyonudur.
  pure function monotonic(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value

    value = x
  end function monotonic

end program test_tts_scalar_minimizer
