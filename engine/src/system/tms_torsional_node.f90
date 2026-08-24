module tms_torsional_node
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> Ayrık-parametreli torsional sistemde tek bir fiziksel dönme düğümünü
  !! temsil eder. Düğüm kimliği topolojik etikettir; doğrudan DOF sıra numarası
  !! olmak zorunda değildir.
  type, public :: torsional_node_t
    !> Düğümün sistem içindeki benzersiz pozitif kimliği [-]. Matematiksel
    !! olarak bağlantı topolojisindeki düğüm etiketidir.
    integer :: id = 0

    !> Düğüme yığılmış polar kütle ataleti J [kg*m^2]. Fiziksel olarak açısal
    !! ivmeye karşı direnci, matematiksel olarak gelecekteki M matrisinin ilgili
    !! atalet katkısını temsil eder.
    real(dp) :: polar_inertia_kg_m2 = 0.0_dp

    !> Analiz başlangıcındaki açısal konum theta_0 [rad]. Matematiksel olarak
    !! başlangıç koşulu bileşenidir; V0.2.3 statik veri modeli bunu çözmez.
    real(dp) :: initial_angle_rad = 0.0_dp

    !> Düğüm dönmesinin kinematik olarak sabitlenip sabitlenmediğini gösteren
    !! boyutsuz mantıksal bilgidir. Doğruysa düğüm aktif DOF sayısına katılmaz.
    logical :: constrained = .false.
  end type torsional_node_t

  public :: validate_torsional_node

contains

  !> Torsional düğümün temel fiziksel ve topolojik önkoşullarını doğrular.
  !!
  !! Fiziksel açıklama: Rijit bir dönel gövdeyi temsil eden düğüm pozitif polar
  !! kütle ataletine sahip olmalıdır.
  !! Matematiksel açıklama: id > 0, sonlu J > 0 ve sonlu theta_0 koşulları
  !! uygulanır.
  !! Girdi: Boyutsuz düğüm kimliği, J [kg*m^2], başlangıç açısı [rad] ve
  !! boyutsuz sabitlenmişlik bilgisini taşıyan torsional_node_t değeridir.
  !! Çıktı üretmez; geçersiz düğüm error stop ile reddedilir.
  !! Varsayımlar ve geçerlilik: Atalet düğümde yığılmıştır. Başlangıç açısı
  !! sonlu olmalı, fakat bu veri katmanında ayrıca büyüklük sınırı uygulanmaz.
  pure subroutine validate_torsional_node(node)
    type(torsional_node_t), intent(in) :: node

    if (node%id <= 0) then
      error stop "Torsional düğüm kimliği pozitif olmalıdır."
    end if

    if (.not. ieee_is_finite(node%polar_inertia_kg_m2) .or. &
        node%polar_inertia_kg_m2 <= 0.0_dp) then
      error stop "Torsional düğüm polar ataleti sonlu ve pozitif olmalıdır."
    end if

    if (.not. ieee_is_finite(node%initial_angle_rad)) then
      error stop "Torsional düğüm başlangıç açısı sonlu olmalıdır."
    end if
  end subroutine validate_torsional_node

end module tms_torsional_node
