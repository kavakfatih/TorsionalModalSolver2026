module tms_torsional_element
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> İki torsional düğüm arasındaki lineer bağlantının topolojik ve fiziksel
  !! özelliklerini taşır.
  type, public :: torsional_element_t
    !> Elemanın sistem içindeki benzersiz pozitif kimliği [-]. Matematiksel
    !! olarak bağlantı topolojisindeki eleman etiketidir.
    integer :: id = 0

    !> Bağlantının ilk ucundaki düğüm kimliği [-]. Bu değer DOF sıra numarası
    !! değil, torsional_node_t kimliğine başvurudur.
    integer :: node_i_id = 0

    !> Bağlantının ikinci ucundaki düğüm kimliği [-]. İlk uçtan farklı ve
    !! sistemde tanımlı bir torsional_node_t kimliğine başvurmalıdır.
    integer :: node_j_id = 0

    !> İki düğüm arasındaki lineer burulma rijitliği K [N*m/rad]. Fiziksel
    !! olarak bağıl açıya karşı geri çağırıcı momenti, matematiksel olarak
    !! T_K = K*(theta_i-theta_j) bağıntısındaki katsayıyı temsil eder.
    real(dp) :: stiffness_nm_per_rad = 0.0_dp

    !> İki düğüm arasındaki eşdeğer viskoz burulma sönüm katsayısı c
    !! [N*m*s/rad]. Matematiksel olarak bağıl açısal hıza karşı
    !! T_c = c*(theta_i'-theta_j') moment katsayısıdır.
    real(dp) :: damping_nms_per_rad = 0.0_dp
  end type torsional_element_t

  public :: validate_torsional_element

contains

  !> Torsional elemanın temel fiziksel ve topolojik önkoşullarını doğrular.
  !!
  !! Fiziksel açıklama: Lineer bağlantı pozitif rijitlik ve pasif model için
  !! negatif olmayan viskoz sönüm taşır; elemanın iki ucu farklı düğümlerdir.
  !! Matematiksel açıklama: id > 0, node_i > 0, node_j > 0,
  !! node_i /= node_j, sonlu K > 0 ve sonlu c >= 0 koşulları uygulanır.
  !! Girdi: Kimlikleri boyutsuz, K değeri [N*m/rad] ve c değeri
  !! [N*m*s/rad] olan torsional_element_t değeridir.
  !! Çıktı üretmez; geçersiz eleman error stop ile reddedilir.
  !! Varsayımlar ve geçerlilik: Rijitlik ve sönüm lineerdir. Düğüm
  !! başvurularının sistemde bulunması ana sistem doğrulamasında denetlenir.
  pure subroutine validate_torsional_element(element)
    type(torsional_element_t), intent(in) :: element

    if (element%id <= 0) then
      error stop "Torsional eleman kimliği pozitif olmalıdır."
    end if

    if (element%node_i_id <= 0 .or. element%node_j_id <= 0) then
      error stop "Torsional eleman düğüm kimlikleri pozitif olmalıdır."
    end if

    if (element%node_i_id == element%node_j_id) then
      error stop "Torsional eleman iki farklı düğümü bağlamalıdır."
    end if

    if (.not. ieee_is_finite(element%stiffness_nm_per_rad) .or. &
        element%stiffness_nm_per_rad <= 0.0_dp) then
      error stop "Torsional eleman rijitliği sonlu ve pozitif olmalıdır."
    end if

    if (.not. ieee_is_finite(element%damping_nms_per_rad) .or. &
        element%damping_nms_per_rad < 0.0_dp) then
      error stop "Torsional eleman viskoz sönümü sonlu olmalı ve negatif olmamalıdır."
    end if
  end subroutine validate_torsional_element

end module tms_torsional_element
