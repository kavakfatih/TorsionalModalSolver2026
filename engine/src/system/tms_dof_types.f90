module tms_dof_types
  implicit none
  private

  !> Ortak dönme ekseni çevresindeki torsional açısal serbestlik derecesi.
  !! Değer boyutsuz bir tür kimliğidir; açısal sonuçların SI birimi radyandır.
  integer, parameter, public :: TORSIONAL_ROTATION = 1

  !> Topolojideki fiziksel serbestlik derecesini numaralandırmadan bağımsız
  !! olarak tanımlar.
  !!
  !! Fiziksel anlam: node_id fiziksel düğümü, dof_type ise o düğümdeki hareket
  !! bileşenini belirtir. TMS26'nın mevcut ayrık torsional modelinde her düğüm
  !! yalnız TORSIONAL_ROTATION türünde tek bir açısal DOF taşır.
  !! Matematiksel anlam: (node_id,dof_type) ikilisi fiziksel koordinatın
  !! benzersiz anahtarıdır; equation veya matrix indeksi değildir.
  !! Birimler: Her iki alan da boyutsuz tam sayı kimliktir.
  type, public :: physical_dof_t
    integer :: node_id = 0
    integer :: dof_type = 0
  end type physical_dof_t

  public :: validate_physical_dof
  public :: is_supported_dof_type
  public :: are_physical_dofs_equal
  public :: operator(==)

  interface operator(==)
    module procedure are_physical_dofs_equal
  end interface operator(==)

contains

  !> Fiziksel torsional DOF anahtarının geçerli olduğunu doğrular.
  !!
  !! Fiziksel açıklama: Her DOF pozitif kimlikli bir fiziksel düğüme ve bu
  !! sürümün desteklediği torsional dönme bileşenine bağlı olmalıdır.
  !! Matematiksel açıklama: node_id > 0 ve dof_type = TORSIONAL_ROTATION
  !! koşulları uygulanır.
  !! Girdi: Boyutsuz kimlikler taşıyan physical_dof_t. Çıktı üretmez;
  !! geçersiz anahtar error stop ile reddedilir.
  !! Varsayımlar ve geçerlilik: Düğümün gerçekten bir sistemde bulunması üst
  !! seviye sistem veya DOF haritası doğrulamasının sorumluluğundadır.
  pure subroutine validate_physical_dof(dof)
    type(physical_dof_t), intent(in) :: dof

    if (dof%node_id <= 0) then
      error stop "Fiziksel DOF düğüm kimliği pozitif olmalıdır."
    end if

    if (.not. is_supported_dof_type(dof%dof_type)) then
      error stop "Fiziksel DOF türü desteklenen torsional dönme türü olmalıdır."
    end if
  end subroutine validate_physical_dof

  !> Verilen boyutsuz DOF tür kimliğinin bu sürümde desteklenip
  !! desteklenmediğini döndürür. V0.4.0 yalnız TORSIONAL_ROTATION kabul eder;
  !! gelecekteki DOF türleri bu merkezi seçim noktasına eklenir.
  pure elemental function is_supported_dof_type(dof_type) &
      result(is_supported)
    integer, intent(in) :: dof_type
    logical :: is_supported

    is_supported = dof_type == TORSIONAL_ROTATION
  end function is_supported_dof_type

  !> İki fiziksel DOF anahtarının aynı koordinatı gösterip göstermediğini sınar.
  !!
  !! Matematiksel model: Eşitlik ancak node_id ve dof_type alanlarının ikisi de
  !! eşitse doğrudur. Girdiler ve mantıksal çıktı boyutsuzdur.
  !! Varsayımlar ve geçerlilik: Bu yapısal karşılaştırma girdileri doğrulamaz;
  !! geçersiz fakat özdeş iki anahtar da yapısal olarak eşit kabul edilir.
  pure elemental function are_physical_dofs_equal(first, second) &
      result(are_equal)
    type(physical_dof_t), intent(in) :: first
    type(physical_dof_t), intent(in) :: second
    logical :: are_equal

    are_equal = first%node_id == second%node_id .and. &
      first%dof_type == second%dof_type
  end function are_physical_dofs_equal

end module tms_dof_types
