module tms_deterministic_rng
  use, intrinsic :: iso_fortran_env, only : int64
  use tms_kinds, only : dp
  implicit none
  private

  integer(int64), parameter :: park_miller_modulus = 2147483647_int64
  integer(int64), parameter :: park_miller_multiplier = 16807_int64
  integer(int64), parameter :: park_miller_quotient = 127773_int64
  integer(int64), parameter :: park_miller_remainder = 2836_int64

  !> Yerel ve explicit seed'li Park-Miller pseudo-random state'idir. Schrage
  !! güncellemesi signed integer overflow oluşturmadan tanımlı int64 aritmetiği
  !! kullanır. Cryptographic amaç taşımaz; bootstrap testlerinin platformlar
  !! arasında deterministic olmasını sağlar.
  type, public :: deterministic_rng_t
    private
    integer(int64) :: state = 1_int64
  contains
    procedure, public :: initialize => initialize_rng
    procedure, public :: next_int31 => next_rng_int31
    procedure, public :: next_uniform => next_rng_uniform
    procedure, public :: next_index => next_rng_index
    procedure, public :: current_state => current_rng_state
  end type deterministic_rng_t

contains

  !> Explicit integer seed'i 1..m-1 Park-Miller state aralığına taşır.
  !! Seed fiziksel bir nicelik değildir; aynı seed aynı bootstrap dizisini
  !! üretir. Negatif ve sıfır seed güvenli biçimde normalize edilir.
  pure subroutine initialize_rng(self, seed)
    class(deterministic_rng_t), intent(inout) :: self
    integer(int64), intent(in) :: seed

    self%state = modulo(seed, park_miller_modulus - 1_int64)
    if (self%state == 0_int64) self%state = 1_int64
  end subroutine initialize_rng

  !> Park-Miller minimal-standard recurrence'in bir sonraki 31-bit pozitif
  !! durumunu üretir. Schrage ayrıştırması a*state mod m hesabındaki taşmayı
  !! engeller; çıktı boyutsuz integer pseudo-random değerdir.
  pure subroutine next_rng_int31(self, value)
    class(deterministic_rng_t), intent(inout) :: self
    integer(int64), intent(out) :: value

    integer(int64) :: candidate
    integer(int64) :: high_part
    integer(int64) :: low_part

    high_part = self%state/park_miller_quotient
    low_part = modulo(self%state, park_miller_quotient)
    candidate = park_miller_multiplier*low_part - &
      park_miller_remainder*high_part
    if (candidate <= 0_int64) candidate = candidate + park_miller_modulus
    self%state = candidate
    value = self%state
  end subroutine next_rng_int31

  !> Bir sonraki Park-Miller durumunu açık (0,1) aralığındaki boyutsuz real
  !! değere dönüştürür. Bu değer bootstrap sampling içindir; cryptographic
  !! veya physical randomness iddiası taşımaz.
  pure subroutine next_rng_uniform(self, value)
    class(deterministic_rng_t), intent(inout) :: self
    real(dp), intent(out) :: value

    integer(int64) :: integer_value

    call self%next_int31(integer_value)
    value = real(integer_value, dp)/real(park_miller_modulus, dp)
  end subroutine next_rng_uniform

  !> 1..upper_bound aralığından replacement sampling için deterministic bir
  !! index üretir. upper_bound bir campaign population boyutudur; geçersiz
  !! sınırda index=0 ve valid=false döner, programı durdurmaz.
  pure subroutine next_rng_index(self, upper_bound, index, valid)
    class(deterministic_rng_t), intent(inout) :: self
    integer, intent(in) :: upper_bound
    integer, intent(out) :: index
    logical, intent(out) :: valid

    integer(int64) :: integer_value

    index = 0
    valid = upper_bound > 0
    if (.not. valid) return
    call self%next_int31(integer_value)
    index = 1 + int(modulo(integer_value, int(upper_bound, int64)))
  end subroutine next_rng_index

  !> Test/provenance amacıyla güncel deterministic RNG state'ini döndürür.
  pure function current_rng_state(self) result(state)
    class(deterministic_rng_t), intent(in) :: self
    integer(int64) :: state

    state = self%state
  end function current_rng_state

end module tms_deterministic_rng
