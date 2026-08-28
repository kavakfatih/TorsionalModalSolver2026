module tms_temperature_shift_types
  use tms_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: TEMPERATURE_SHIFT_NONE = 0
  integer, parameter, public :: TABULATED_LOG10_SHIFT = 1
  integer, parameter, public :: WLF_TEMPERATURE_SHIFT = 2
  integer, parameter, public :: ARRHENIUS_TEMPERATURE_SHIFT = 3

  !> Bir sıcaklık kaydırma modelinin doğrulanmış kullanım aralığını taşır.
  !! Bütün sıcaklıklar mutlak sıcaklık olarak SI birimi K cinsindedir.
  type, public :: temperature_shift_domain_t
    real(dp) :: minimum_temperature_k = 0.0_dp
    real(dp) :: maximum_temperature_k = 0.0_dp
    real(dp) :: reference_temperature_k = 0.0_dp
  end type temperature_shift_domain_t

  !> Bir operating sıcaklığındaki yatay kaydırma sonucunu ve izini taşır.
  !! TMS26 convention'ı a_T=tau(T)/tau(T_ref), f_r=a_T*f ve
  !! t_r=t/a_T biçimindedir. Authoritative değer log10(a_T), a_T ise bundan
  !! türetilen boyutsuz değerdir. Sıcaklık bracket alanları yalnız tabulated
  !! model için anlamlıdır ve K birimindedir; alpha boyutsuzdur.
  type, public :: temperature_shift_evaluation_t
    integer :: shift_model_kind = TEMPERATURE_SHIFT_NONE
    real(dp) :: operating_temperature_k = 0.0_dp
    real(dp) :: reference_temperature_k = 0.0_dp
    real(dp) :: log10_a_t = 0.0_dp
    real(dp) :: a_t = 1.0_dp
    logical :: has_temperature_bracket = .false.
    logical :: exact_temperature_point = .false.
    real(dp) :: lower_temperature_k = 0.0_dp
    real(dp) :: upper_temperature_k = 0.0_dp
    real(dp) :: interpolation_alpha = 0.0_dp
  end type temperature_shift_evaluation_t

end module tms_temperature_shift_types
