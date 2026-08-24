module tms_dof_map
  use tms_torsional_node, only : torsional_node_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    count_torsional_dofs, get_torsional_node_count, get_torsional_node, &
    validate_torsional_system
  implicit none
  private

  !> Bir fiziksel torsional düğüm kimliği ile matematiksel denklem kimliği
  !! arasındaki tek eşleme kaydını taşır.
  type, public :: dof_map_entry_t
    !> Topolojideki pozitif fiziksel düğüm kimliği [-].
    integer :: node_id = 0
    !> Aktif denklem kimliği [-]. Kısıtlı düğüm için sıfır, serbest düğüm için
    !! kesintisiz 1..n_active aralığındadır.
    integer :: equation_id = -1
  end type dof_map_entry_t

  !> Fiziksel düğüm kimliklerini solver denklem numaralarından ayıran harita.
  !! Kayıtlar private tutulur; böylece kısıt ve numaralandırma invariantları
  !! yalnız doğrulanmış başlatma yordamıyla kurulabilir.
  type, public :: dof_map_t
    private
    type(dof_map_entry_t), allocatable :: entries(:)
    integer :: active_dof_count = 0
  end type dof_map_t

  public :: initialize_dof_map
  public :: validate_dof_map
  public :: validate_dof_map_for_system
  public :: lookup_equation_id
  public :: get_active_dof_count
  public :: get_dof_map_entry_count
  public :: get_dof_map_entry

contains

  !> Genel torsional sistem için deterministik aktif DOF haritası oluşturur.
  !!
  !! Fiziksel açıklama: Kısıtlı dönel düğümler indirgenmiş denklem takımından
  !! çıkarılır; tüm fiziksel düğümler izlenebilirlik için haritada tutulur.
  !! Matematiksel açıklama: Ekleme sırasındaki serbest düğümlere kesintisiz
  !! equation_id=1..n verilir, kısıtlı düğümlere equation_id=0 atanır.
  !! Girdi: Doğrulanmış torsional_system_t. Çıktı: node_id ve boyutsuz
  !! equation_id kayıtlarını taşıyan dof_map_t.
  !! Varsayımlar ve geçerlilik: Her düğüm tek torsional DOF taşır. Sıfır kimliği
  !! yalnız homojen kısıt işaretidir; eksik düğüm anlamına gelmez.
  pure subroutine initialize_dof_map(mapping, system)
    type(dof_map_t), intent(out) :: mapping
    type(torsional_system_t), intent(in) :: system

    type(torsional_node_t) :: node
    integer :: node_count
    integer :: node_index

    call validate_torsional_system(system)
    node_count = get_torsional_node_count(system)
    allocate(mapping%entries(node_count))
    mapping%active_dof_count = 0

    do node_index = 1, node_count
      node = get_torsional_node(system, node_index)
      mapping%entries(node_index)%node_id = node%id

      if (node%constrained) then
        mapping%entries(node_index)%equation_id = 0
      else
        mapping%active_dof_count = mapping%active_dof_count + 1
        mapping%entries(node_index)%equation_id = mapping%active_dof_count
      end if
    end do

    call validate_dof_map(mapping)
  end subroutine initialize_dof_map

  !> DOF haritasının kimlik, benzersizlik ve süreklilik invariantlarını sınar.
  !! Matematiksel açıklama: node_id değerleri pozitif ve benzersiz; pozitif
  !! equation_id değerleri benzersiz ve tam olarak 1..n_active kümesidir.
  !! equation_id=0 kısıtlı düğüm için ayrılmıştır.
  !! Girdi dof_map_t, çıktı yoktur; geçersiz harita error stop ile reddedilir.
  pure subroutine validate_dof_map(mapping)
    type(dof_map_t), intent(in) :: mapping

    integer :: first_index
    integer :: second_index
    integer :: equation_id
    integer :: positive_equation_count

    if (.not. allocated(mapping%entries)) then
      error stop "DOF haritası kullanılmadan önce başlatılmalıdır."
    end if

    if (size(mapping%entries) == 0) then
      error stop "DOF haritası en az bir fiziksel düğüm içermelidir."
    end if

    if (mapping%active_dof_count < 0) then
      error stop "Aktif DOF sayısı negatif olamaz."
    end if

    positive_equation_count = 0
    do first_index = 1, size(mapping%entries)
      if (mapping%entries(first_index)%node_id <= 0) then
        error stop "DOF haritasındaki fiziksel düğüm kimliği pozitif olmalıdır."
      end if

      if (mapping%entries(first_index)%equation_id < 0 .or. &
          mapping%entries(first_index)%equation_id > &
          mapping%active_dof_count) then
        error stop "DOF haritasındaki denklem kimliği geçersiz."
      end if

      if (mapping%entries(first_index)%equation_id > 0) then
        positive_equation_count = positive_equation_count + 1
      end if

      do second_index = first_index + 1, size(mapping%entries)
        if (mapping%entries(first_index)%node_id == &
            mapping%entries(second_index)%node_id) then
          error stop "DOF haritasındaki düğüm kimlikleri benzersiz olmalıdır."
        end if

        if (mapping%entries(first_index)%equation_id > 0 .and. &
            mapping%entries(first_index)%equation_id == &
            mapping%entries(second_index)%equation_id) then
          error stop "Aktif denklem kimlikleri benzersiz olmalıdır."
        end if
      end do
    end do

    if (positive_equation_count /= mapping%active_dof_count) then
      error stop "DOF haritasındaki aktif denklem sayısı tutarsız."
    end if

    do equation_id = 1, mapping%active_dof_count
      if (.not. any(mapping%entries%equation_id == equation_id)) then
        error stop "Aktif denklem kimlikleri kesintisiz olmalıdır."
      end if
    end do
  end subroutine validate_dof_map

  !> DOF haritasının belirli fiziksel sistemle bire bir uyumunu doğrular.
  !! Matematiksel açıklama: Düğüm kümeleri aynıdır; constrained=.true. yalnız
  !! equation_id=0, serbest düğümler yalnız pozitif denklem kimliği taşır.
  !! Girdiler dof_map_t ve torsional_system_t, çıktı yoktur. Eski veya başka
  !! sisteme ait harita error stop ile reddedilir.
  pure subroutine validate_dof_map_for_system(mapping, system)
    type(dof_map_t), intent(in) :: mapping
    type(torsional_system_t), intent(in) :: system

    type(torsional_node_t) :: node
    integer :: entry_index
    integer :: node_index

    call validate_dof_map(mapping)
    call validate_torsional_system(system)

    if (size(mapping%entries) /= get_torsional_node_count(system)) then
      error stop "DOF haritası ile sistemin fiziksel düğüm sayısı farklı."
    end if

    if (mapping%active_dof_count /= count_torsional_dofs(system)) then
      error stop "DOF haritası ile sistemin aktif DOF sayısı farklı."
    end if

    do node_index = 1, get_torsional_node_count(system)
      node = get_torsional_node(system, node_index)
      entry_index = find_dof_map_entry_index(mapping, node%id)
      if (entry_index == 0) then
        error stop "Sistem düğümü DOF haritasında bulunamadı."
      end if

      if (node%constrained .neqv. &
          (mapping%entries(entry_index)%equation_id == 0)) then
        error stop "Düğüm kısıtı ile DOF haritası denklem kimliği uyumsuz."
      end if
    end do
  end subroutine validate_dof_map_for_system

  !> Fiziksel düğüm kimliğine karşılık gelen matematiksel denklem kimliğini
  !! döndürür. Girdi pozitif node_id, çıktı boyutsuz equation_id değeridir.
  !! Kısıtlı düğüm sıfır döndürür; bulunamayan düğüm error stop üretir ve
  !! böylece kısıtlı düğüm ile eksik düğüm birbirine karışmaz.
  pure function lookup_equation_id(mapping, node_id) result(equation_id)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer :: equation_id

    integer :: entry_index

    if (.not. allocated(mapping%entries)) then
      error stop "DOF haritası kullanılmadan önce başlatılmalıdır."
    end if

    entry_index = find_dof_map_entry_index(mapping, node_id)
    if (entry_index == 0) then
      error stop "Fiziksel düğüm kimliği DOF haritasında bulunamadı."
    end if

    equation_id = mapping%entries(entry_index)%equation_id
  end function lookup_equation_id

  !> Haritadaki boyutsuz aktif denklem sayısını döndürür.
  pure function get_active_dof_count(mapping) result(active_dof_count)
    type(dof_map_t), intent(in) :: mapping
    integer :: active_dof_count

    call validate_dof_map(mapping)
    active_dof_count = mapping%active_dof_count
  end function get_active_dof_count

  !> Haritada tutulan fiziksel düğüm kaydı sayısını döndürür.
  pure function get_dof_map_entry_count(mapping) result(entry_count)
    type(dof_map_t), intent(in) :: mapping
    integer :: entry_count

    call validate_dof_map(mapping)
    entry_count = size(mapping%entries)
  end function get_dof_map_entry_count

  !> Harita kaydının ekleme sırasındaki bağımsız bir kopyasını döndürür.
  !! Girdiler dof_map_t ve bir tabanlı boyutsuz indeks, çıktı node/equation
  !! kimliklerini taşıyan dof_map_entry_t değeridir.
  pure function get_dof_map_entry(mapping, index) result(entry)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: index
    type(dof_map_entry_t) :: entry

    call validate_dof_map(mapping)
    if (index < 1 .or. index > size(mapping%entries)) then
      error stop "DOF haritası kayıt indeksi geçersiz."
    end if

    entry = mapping%entries(index)
  end function get_dof_map_entry

  !> Fiziksel düğüm kimliğinin haritadaki bir tabanlı indeksini bulur.
  !! Bulunamayan kimlik için sıfır döndürür. Bu özel topoloji yardımcı yordamı
  !! fiziksel veya matematiksel katsayı hesabı yapmaz.
  pure function find_dof_map_entry_index(mapping, node_id) result(entry_index)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer :: entry_index

    integer :: current_index

    entry_index = 0
    if (.not. allocated(mapping%entries)) return

    do current_index = 1, size(mapping%entries)
      if (mapping%entries(current_index)%node_id == node_id) then
        entry_index = current_index
        return
      end if
    end do
  end function find_dof_map_entry_index

end module tms_dof_map
