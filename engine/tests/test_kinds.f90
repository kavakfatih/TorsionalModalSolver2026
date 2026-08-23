program test_kinds
  use tms_kinds, only : dp
  implicit none

  real(dp) :: sample_value

  ! Basit atama ve tür sorgusu, Fortran derleyicisi ile modül bağını doğrular.
  sample_value = 1.0_dp

  if (precision(sample_value) < 15) then
    error stop "dp türü beklenen çift hassasiyeti sağlamıyor."
  end if

  print *, "Fortran derleyicisi ve dp türü doğrulandı."
end program test_kinds
