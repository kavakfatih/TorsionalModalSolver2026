module tms_local_matrix
  use tms_kinds, only : dp
  implicit none
  private

  !> İki yerel torsional serbestlik derecesine ait küçük matris katkısını
  !! taşır.
  !!
  !! Fiziksel anlam: Satır ve sütunlar elemanın sırasıyla i ve j uçlarındaki
  !! açısal koordinatlarla ilişkilidir. Bu sürümde tür, torsional rijitlik
  !! matrisinin dört katsayısını saklamak için kullanılır.
  !! Matematiksel anlam: value, yerel koordinat sırası [theta_i, theta_j]
  !! olan 2x2 katsayı matrisidir.
  !! Birimler: Tür birimi kendi başına kodlamaz; birim, matrisi üreten yordamın
  !! sözleşmesinden gelir. Lokal torsional rijitlik için tüm değerler
  !! [N*m/rad] SI birimindedir.
  type, public :: local_matrix_2x2
    real(dp) :: value(2, 2) = 0.0_dp
  end type local_matrix_2x2

end module tms_local_matrix
