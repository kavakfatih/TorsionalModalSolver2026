module tms_kinds
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  private

  ! Gerçel sayılar için proje genelinde kullanılacak çift hassasiyetli türdür.
  integer, parameter, public :: dp = real64

end module tms_kinds
