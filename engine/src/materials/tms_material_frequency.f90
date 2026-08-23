module tms_material_frequency
  use tms_dynamic_modulus, only : dynamic_shear_modulus
  implicit none
  private

  !> Bir elastomer malzemenin belirli frekans ve sıcaklıktaki dinamik
  !! kayma modülü veri noktasını temsil eder.
  !!
  !! Kalıtılan değişkenler storage_modulus [Pa], loss_modulus [Pa],
  !! frequency [Hz] ve temperature [K] değerleridir. Bu yapı ileride DMA
  !! deneyleri ile Dewesoft veya Brüel & Kjær modal test verilerinin aynı
  !! fiziksel alanlarda eşleştirilmesine temel sağlar. Bu sürümde veri
  !! interpolasyonu, eğri uydurma ve frekans dönüşümü yapılmaz.
  type, extends(dynamic_shear_modulus), public :: material_frequency_point
  end type material_frequency_point

end module tms_material_frequency
