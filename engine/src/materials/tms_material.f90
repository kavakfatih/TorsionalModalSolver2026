module tms_material
  use tms_kinds, only : dp
  implicit none
  private

  ! Malzeme adları için ayrılan sabit karakter uzunluğudur.
  integer, parameter, public :: material_name_length = 128

  !> Dinamik elastomer malzemenin bir çalışma noktasındaki verilerini taşır.
  !! G' enerji depolayan elastik bileşeni, G'' ise kayıp bileşenini temsil eder.
  !! Modüller Pa, yoğunluk kg/m^3, sıcaklık K ve frekans Hz cinsindendir.
  type, public :: dynamic_rubber_material_t
    character(len=material_name_length) :: name = ""
    real(dp) :: density_kg_m3 = 0.0_dp
    real(dp) :: storage_shear_modulus_pa = 0.0_dp
    real(dp) :: loss_shear_modulus_pa = 0.0_dp
    real(dp) :: temperature_k = 0.0_dp
    real(dp) :: frequency_hz = 0.0_dp
  end type dynamic_rubber_material_t

end module tms_material
