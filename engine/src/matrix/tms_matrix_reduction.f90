module tms_matrix_reduction
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    extract_stiffness_principal_submatrix
  use tms_mass_matrix, only : mass_matrix_t, &
    extract_mass_principal_submatrix
  implicit none
  private

  public :: reduce_matrix

  !> Full torsional matrisi retained full equation indekslerine göre indirger.
  !! Generic arayüz fiziksel matrix türüne göre rijitlik veya atalet
  !! gerçekleştirimini seçer; somut dense/sparse storage türünü dışarı açmaz.
  interface reduce_matrix
    module procedure reduce_stiffness_matrix
    module procedure reduce_mass_matrix
  end interface reduce_matrix

contains

  !> Full torsional rijitlik matrisini active equation uzayına indirger.
  !!
  !! Fiziksel açıklama: Kinematik olarak elenen torsional DOF'ların satır ve
  !! sütunları full rijitlik denklem takımından çıkarılır.
  !! Matematiksel açıklama: K_active=P^T*K_full*P. P, her active equation için
  !! karşılık gelen full equation logical indeksini seçen Boolean prolongation
  !! operatörüdür; uygulama P matrisini açıkça oluşturmaz.
  !! Girdiler: K_full [N*m/rad] ve active sıra ile verilmiş, benzersiz bir
  !! tabanlı retained full equation indeksleri [-]. Çıktı: K_active
  !! [N*m/rad]. Boş retained dizi 0x0 matris üretir.
  !! Varsayımlar ve geçerlilik: Yalnız DOF seçimine dayalı homojen kinematik
  !! eliminasyon uygulanır; çok noktalı constraint dönüşümü bu kapsamda yoktur.
  pure function reduce_stiffness_matrix( &
      full_matrix, retained_full_equation_indices) result(active_matrix)
    type(stiffness_matrix_t), intent(in) :: full_matrix
    integer, intent(in) :: retained_full_equation_indices(:)
    type(stiffness_matrix_t) :: active_matrix

    active_matrix = extract_stiffness_principal_submatrix( &
      full_matrix, retained_full_equation_indices)
  end function reduce_stiffness_matrix

  !> Full torsional atalet matrisini active equation uzayına indirger.
  !!
  !! Fiziksel açıklama: Kısıtlanan torsional DOF'ların polar atalet denklem
  !! katkıları active hareket denklem takımının dışında bırakılır.
  !! Matematiksel açıklama: M_active=P^T*M_full*P. P, retained full equation
  !! indeksleriyle örtük tanımlanan Boolean prolongation operatörüdür.
  !! Girdiler: M_full [kg*m^2] ve active sıra ile verilmiş, benzersiz bir
  !! tabanlı retained full equation indeksleri [-]. Çıktı: M_active
  !! [kg*m^2]. Boş retained dizi 0x0 matris üretir.
  !! Varsayımlar ve geçerlilik: Matrix katsayıları değiştirilmez; yalnız
  !! principal alt matris seçilir.
  pure function reduce_mass_matrix( &
      full_matrix, retained_full_equation_indices) result(active_matrix)
    type(mass_matrix_t), intent(in) :: full_matrix
    integer, intent(in) :: retained_full_equation_indices(:)
    type(mass_matrix_t) :: active_matrix

    active_matrix = extract_mass_principal_submatrix( &
      full_matrix, retained_full_equation_indices)
  end function reduce_mass_matrix

end module tms_matrix_reduction
