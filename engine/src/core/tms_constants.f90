module tms_constants
  use tms_kinds, only : dp
  implicit none
  private

  !> Birimsiz pi sabiti, dairesel geometri ve açısal hesaplarda kullanılır.
  real(dp), parameter, public :: pi = &
    3.141592653589793238462643383279502884197_dp

  !> Milimetre cinsinden uzunluğu metreye çevirmek için kullanılan çarpandır.
  !! Fiziksel eşitlik: 1 mm = 10^-3 m.
  real(dp), parameter, public :: mm_to_m_factor = 1.0e-3_dp

  !> Megapaskal cinsinden basınç veya gerilmeyi paskala çeviren çarpandır.
  !! Fiziksel eşitlik: 1 MPa = 10^6 Pa.
  real(dp), parameter, public :: mpa_to_pa_factor = 1.0e6_dp

  !> Derece cinsinden düzlem açısını radyana çevirmek için kullanılan çarpandır.
  !! Matematiksel eşitlik: 1 derece = pi / 180 radyan.
  real(dp), parameter, public :: degree_to_radian_factor = pi / 180.0_dp

  !> Evrensel gaz sabiti; bir mol ideal sistem için sıcaklık ile enerji
  !! ölçeğini ilişkilendirir. Arrhenius sıcaklık kaydırma modelinde aktivasyon
  !! enerjisini boyutsuz üsse dönüştürür. SI birimi J/(mol K)'dir.
  real(dp), parameter, public :: universal_gas_constant_j_per_mol_k = &
    8.31446261815324_dp

end module tms_constants
