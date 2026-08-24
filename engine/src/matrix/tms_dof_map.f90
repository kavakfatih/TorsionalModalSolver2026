module tms_dof_map
  use tms_dof_types, only : physical_dof_t, TORSIONAL_ROTATION, &
    validate_physical_dof
  use tms_torsional_node, only : torsional_node_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    count_torsional_dofs, get_torsional_node_count, get_torsional_node, &
    validate_torsional_system
  implicit none
  private

  !> Tek bir fiziksel torsional DOF'un full ve active denklem kimliklerini
  !! geriye uyumlu kayıt biçiminde birlikte taşır.
  type, public :: dof_map_entry_t
    !> Topolojideki pozitif fiziksel düğüm kimliği [-].
    integer :: node_id = 0

    !> Geriye uyumlu denklem kimliği [-]. initialize_dof_map ile active
    !! equation ID, initialize_full_dof_map ile full equation ID değeridir.
    integer :: equation_id = -1

    !> Fiziksel hareket bileşeninin tür kimliği [-]. Mevcut torsional modelde
    !! yalnız TORSIONAL_ROTATION desteklenir.
    integer :: dof_type = TORSIONAL_ROTATION

    !> Constraint durumundan bağımsız fiziksel DOF sıra kimliği [-].
    integer :: physical_dof_id = 0

    !> Constraint durumundan bağımsız full equation kimliği [-]; 1..n_full.
    integer :: full_equation_id = 0

    !> İndirgenmiş active equation kimliği [-]. Legacy haritada constrained
    !! DOF için sıfırdır; full haritada tüm değerler 1..n_full aralığındadır.
    integer :: active_equation_id = -1
  end type dof_map_entry_t

  !> Fiziksel DOF, full equation ve active equation uzayları arasındaki
  !! doğrulanmış eşlemeyi taşır.
  !!
  !! Kayıtlar private tutulur. Full mod constraint-blind assembly için tüm
  !! fiziksel DOF'ları korur; legacy mod mevcut constrained=.true. davranışını
  !! active_equation_id=0 ile sürdürür.
  type, public :: dof_map_t
    private
    type(dof_map_entry_t), allocatable :: entries(:)
    integer :: physical_dof_count = 0
    integer :: full_equation_count = 0
    integer :: active_dof_count = 0
    logical :: full_map_mode = .false.
  end type dof_map_t

  public :: initialize_dof_map
  public :: initialize_full_dof_map
  public :: validate_dof_map
  public :: validate_dof_map_for_system
  public :: lookup_equation_id
  public :: lookup_physical_dof
  public :: lookup_physical_dof_id
  public :: lookup_full_equation_id
  public :: lookup_active_equation_id
  public :: get_physical_dof_count
  public :: get_full_equation_count
  public :: get_active_equation_count
  public :: get_active_dof_count
  public :: get_retained_full_equation_indices
  public :: get_dof_map_entry_count
  public :: get_dof_map_entry
  public :: get_dof_map_entries
  public :: is_full_dof_map
  public :: is_full_equation_map

contains

  !> Genel torsional sistem için geriye uyumlu active DOF haritası oluşturur.
  !!
  !! Fiziksel açıklama: Kısıtlı dönel düğümler indirgenmiş denklem takımından
  !! çıkarılır; tüm fiziksel DOF'lar full kimlikleriyle izlenebilir kalır.
  !! Matematiksel açıklama: Ekleme sırasındaki tüm DOF'lara constraint'ten
  !! bağımsız full_equation_id=1..n_full verilir. Serbest DOF'lar kesintisiz
  !! active_equation_id=1..n_active, kısıtlı DOF'lar sıfır alır. Legacy
  !! equation_id, active_equation_id ile aynı tutulur.
  !! Girdi: Doğrulanmış torsional_system_t. Çıktı: Boyutsuz kimlikler taşıyan
  !! dof_map_t. Varsayım: Her düğüm tek TORSIONAL_ROTATION DOF'u taşır.
  pure subroutine initialize_dof_map(mapping, system)
    type(dof_map_t), intent(out) :: mapping
    type(torsional_system_t), intent(in) :: system

    call initialize_dof_map_internal(mapping, system, .false.)
  end subroutine initialize_dof_map

  !> Constraint-blind full equation haritasını oluşturur.
  !!
  !! Fiziksel açıklama: Full sistem, kinematik constraint uygulanmadan önce
  !! her fiziksel torsional DOF'u içerir.
  !! Matematiksel açıklama: physical_dof_id, full_equation_id,
  !! active_equation_id ve geriye uyumlu equation_id alanlarının tümü ekleme
  !! sırasında 1..n_full olarak atanır. Böylece full matrix logical indexi
  !! full_equation_id ile aynıdır.
  !! Girdi: Doğrulanmış torsional_system_t. Çıktı: Constraint durumundan
  !! bağımsız, boyutsuz kimlikler taşıyan dof_map_t.
  pure subroutine initialize_full_dof_map(mapping, system)
    type(dof_map_t), intent(out) :: mapping
    type(torsional_system_t), intent(in) :: system

    call initialize_dof_map_internal(mapping, system, .true.)
  end subroutine initialize_full_dof_map

  !> İstenen moda göre ortak fiziksel/full/active DOF kayıtlarını başlatır.
  !! Girdiler sistem ve boyutsuz full_mode seçimidir; çıktı doğrulanmış map'tir.
  pure subroutine initialize_dof_map_internal(mapping, system, full_mode)
    type(dof_map_t), intent(out) :: mapping
    type(torsional_system_t), intent(in) :: system
    logical, intent(in) :: full_mode

    type(torsional_node_t) :: node
    integer :: node_count
    integer :: node_index

    call validate_torsional_system(system)
    node_count = get_torsional_node_count(system)
    allocate(mapping%entries(node_count))
    mapping%physical_dof_count = node_count
    mapping%full_equation_count = node_count
    mapping%active_dof_count = 0
    mapping%full_map_mode = full_mode

    do node_index = 1, node_count
      node = get_torsional_node(system, node_index)
      mapping%entries(node_index)%node_id = node%id
      mapping%entries(node_index)%dof_type = TORSIONAL_ROTATION
      mapping%entries(node_index)%physical_dof_id = node_index
      mapping%entries(node_index)%full_equation_id = node_index

      if (full_mode .or. .not. node%constrained) then
        mapping%active_dof_count = mapping%active_dof_count + 1
        mapping%entries(node_index)%active_equation_id = &
          mapping%active_dof_count
      else
        mapping%entries(node_index)%active_equation_id = 0
      end if

      mapping%entries(node_index)%equation_id = &
        mapping%entries(node_index)%active_equation_id
    end do

    call validate_dof_map(mapping)
  end subroutine initialize_dof_map_internal

  !> DOF haritasının fiziksel, full ve active numaralandırma invariantlarını
  !! doğrular.
  !!
  !! Matematiksel açıklama: Fiziksel ve full kimlikler benzersiz ve kesintisiz
  !! 1..n_full; pozitif active kimlikler 1..n_active kümesidir. Legacy
  !! equation_id her zaman active_equation_id ile aynıdır. Full modda bütün
  !! active kimlikler full kimliklerle bire birdir.
  !! Girdi dof_map_t, çıktı yoktur; geçersiz harita error stop ile reddedilir.
  pure subroutine validate_dof_map(mapping)
    type(dof_map_t), intent(in) :: mapping

    type(physical_dof_t) :: physical_dof
    integer :: active_equation_id
    integer :: first_index
    integer :: full_equation_id
    integer :: physical_dof_id
    integer :: positive_active_count
    integer :: second_index

    if (.not. allocated(mapping%entries)) then
      error stop "DOF haritası kullanılmadan önce başlatılmalıdır."
    end if

    if (size(mapping%entries) == 0) then
      error stop "DOF haritası en az bir fiziksel DOF içermelidir."
    end if

    if (mapping%physical_dof_count /= size(mapping%entries) .or. &
        mapping%full_equation_count /= size(mapping%entries)) then
      error stop "Fiziksel DOF veya full equation sayısı kayıtlarla tutarsız."
    end if

    if (mapping%active_dof_count < 0 .or. &
        mapping%active_dof_count > mapping%full_equation_count) then
      error stop "Active DOF sayısı geçersiz."
    end if

    positive_active_count = 0
    do first_index = 1, size(mapping%entries)
      physical_dof = physical_dof_t( &
        node_id=mapping%entries(first_index)%node_id, &
        dof_type=mapping%entries(first_index)%dof_type)
      call validate_physical_dof(physical_dof)

      if (mapping%entries(first_index)%physical_dof_id < 1 .or. &
          mapping%entries(first_index)%physical_dof_id > &
          mapping%physical_dof_count) then
        error stop "Fiziksel DOF kimliği geçersiz."
      end if

      if (mapping%entries(first_index)%full_equation_id < 1 .or. &
          mapping%entries(first_index)%full_equation_id > &
          mapping%full_equation_count) then
        error stop "Full equation kimliği geçersiz."
      end if

      if (mapping%entries(first_index)%active_equation_id < 0 .or. &
          mapping%entries(first_index)%active_equation_id > &
          mapping%active_dof_count) then
        error stop "Active equation kimliği geçersiz."
      end if

      if (mapping%entries(first_index)%equation_id /= &
          mapping%entries(first_index)%active_equation_id) then
        error stop "Legacy ve active equation kimlikleri tutarsız."
      end if

      if (mapping%entries(first_index)%active_equation_id > 0) then
        positive_active_count = positive_active_count + 1
      end if

      do second_index = first_index + 1, size(mapping%entries)
        if (mapping%entries(first_index)%node_id == &
            mapping%entries(second_index)%node_id .and. &
            mapping%entries(first_index)%dof_type == &
            mapping%entries(second_index)%dof_type) then
          error stop "Fiziksel DOF anahtarları benzersiz olmalıdır."
        end if

        if (mapping%entries(first_index)%physical_dof_id == &
            mapping%entries(second_index)%physical_dof_id) then
          error stop "Fiziksel DOF kimlikleri benzersiz olmalıdır."
        end if

        if (mapping%entries(first_index)%full_equation_id == &
            mapping%entries(second_index)%full_equation_id) then
          error stop "Full equation kimlikleri benzersiz olmalıdır."
        end if

        if (mapping%entries(first_index)%active_equation_id > 0 .and. &
            mapping%entries(first_index)%active_equation_id == &
            mapping%entries(second_index)%active_equation_id) then
          error stop "Active equation kimlikleri benzersiz olmalıdır."
        end if
      end do
    end do

    if (positive_active_count /= mapping%active_dof_count) then
      error stop "DOF haritasındaki active equation sayısı tutarsız."
    end if

    do physical_dof_id = 1, mapping%physical_dof_count
      if (.not. any(mapping%entries%physical_dof_id == physical_dof_id)) then
        error stop "Fiziksel DOF kimlikleri kesintisiz olmalıdır."
      end if
    end do

    do full_equation_id = 1, mapping%full_equation_count
      if (.not. any( &
          mapping%entries%full_equation_id == full_equation_id)) then
        error stop "Full equation kimlikleri kesintisiz olmalıdır."
      end if
    end do

    do active_equation_id = 1, mapping%active_dof_count
      if (.not. any( &
          mapping%entries%active_equation_id == active_equation_id)) then
        error stop "Active equation kimlikleri kesintisiz olmalıdır."
      end if
    end do

    if (mapping%full_map_mode) then
      if (mapping%active_dof_count /= mapping%full_equation_count .or. &
          any(mapping%entries%active_equation_id /= &
          mapping%entries%full_equation_id)) then
        error stop "Full DOF haritası bütün full equation değerlerini korumalıdır."
      end if
    end if
  end subroutine validate_dof_map

  !> DOF haritasının belirli fiziksel sistemle bire bir uyumunu doğrular.
  !!
  !! Matematiksel açıklama: Düğüm ve fiziksel DOF kümeleri aynıdır; full
  !! kimlikler sistem ekleme sırasından bağımsız constraint değişiminde sabit
  !! kalır. Legacy modda constrained=.true. yalnız active_equation_id=0 ile
  !! eşleşir; full mod constraint durumunu numaralandırmaya katmaz.
  !! Girdiler dof_map_t ve torsional_system_t, çıktı yoktur.
  pure subroutine validate_dof_map_for_system(mapping, system)
    type(dof_map_t), intent(in) :: mapping
    type(torsional_system_t), intent(in) :: system

    type(torsional_node_t) :: node
    integer :: entry_index
    integer :: node_index

    call validate_dof_map(mapping)
    call validate_torsional_system(system)

    if (mapping%physical_dof_count /= get_torsional_node_count(system) .or. &
        mapping%full_equation_count /= get_torsional_node_count(system)) then
      error stop "DOF haritası ile sistemin fiziksel DOF sayısı farklı."
    end if

    if (mapping%full_map_mode) then
      if (mapping%active_dof_count /= get_torsional_node_count(system)) then
        error stop "Full DOF haritası bütün fiziksel DOF'ları içermelidir."
      end if
    else
      if (mapping%active_dof_count /= count_torsional_dofs(system)) then
        error stop "DOF haritası ile sistemin active DOF sayısı farklı."
      end if
    end if

    do node_index = 1, get_torsional_node_count(system)
      node = get_torsional_node(system, node_index)
      entry_index = find_dof_map_entry_index(mapping, node%id)
      if (entry_index == 0) then
        error stop "Sistem düğümü DOF haritasında bulunamadı."
      end if

      if (mapping%entries(entry_index)%physical_dof_id /= node_index .or. &
          mapping%entries(entry_index)%full_equation_id /= node_index) then
        error stop "Sistem ekleme sırası ile full DOF numaralandırması uyumsuz."
      end if

      if (.not. mapping%full_map_mode) then
        if (node%constrained .neqv. &
            (mapping%entries(entry_index)%active_equation_id == 0)) then
          error stop "Düğüm kısıtı ile active equation kimliği uyumsuz."
        end if
      end if
    end do
  end subroutine validate_dof_map_for_system

  !> Fiziksel düğüm kimliğinin geriye uyumlu equation kimliğini döndürür.
  !! initialize_dof_map ile kısıtlı düğüm için sıfır ve serbest düğüm için
  !! active kimlik; initialize_full_dof_map ile pozitif full kimlik döner.
  pure function lookup_equation_id(mapping, node_id, dof_type) &
      result(equation_id)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in), optional :: dof_type
    integer :: equation_id

    equation_id = lookup_active_equation_id(mapping, node_id, dof_type)
  end function lookup_equation_id

  !> Düğüm kimliğine ait fiziksel DOF anahtarını döndürür.
  !! Girdi ve çıktı kimlikleri boyutsuzdur; bulunamayan düğüm reddedilir.
  pure function lookup_physical_dof(mapping, node_id, dof_type) &
      result(physical_dof)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in), optional :: dof_type
    type(physical_dof_t) :: physical_dof

    integer :: entry_index

    entry_index = require_dof_map_entry_index(mapping, node_id, dof_type)
    physical_dof = physical_dof_t( &
      node_id=mapping%entries(entry_index)%node_id, &
      dof_type=mapping%entries(entry_index)%dof_type)
  end function lookup_physical_dof

  !> Düğüm kimliğine ait constraint-independent fiziksel DOF kimliğini döndürür.
  pure function lookup_physical_dof_id(mapping, node_id, dof_type) &
      result(physical_dof_id)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in), optional :: dof_type
    integer :: physical_dof_id

    integer :: entry_index

    entry_index = require_dof_map_entry_index(mapping, node_id, dof_type)
    physical_dof_id = mapping%entries(entry_index)%physical_dof_id
  end function lookup_physical_dof_id

  !> Düğüm kimliğine ait pozitif constraint-independent full equation
  !! kimliğini döndürür. Girdi ve çıktı kimlikleri boyutsuzdur.
  pure function lookup_full_equation_id(mapping, node_id, dof_type) &
      result(full_equation_id)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in), optional :: dof_type
    integer :: full_equation_id

    integer :: entry_index

    entry_index = require_dof_map_entry_index(mapping, node_id, dof_type)
    full_equation_id = mapping%entries(entry_index)%full_equation_id
  end function lookup_full_equation_id

  !> Düğüm kimliğine ait active equation kimliğini döndürür.
  !! Legacy haritada constrained DOF sıfır, full haritada bütün DOF'lar pozitif
  !! değer taşır. Girdi ve çıktı kimlikleri boyutsuzdur.
  pure function lookup_active_equation_id(mapping, node_id, dof_type) &
      result(active_equation_id)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in), optional :: dof_type
    integer :: active_equation_id

    integer :: entry_index

    entry_index = require_dof_map_entry_index(mapping, node_id, dof_type)
    active_equation_id = mapping%entries(entry_index)%active_equation_id
  end function lookup_active_equation_id

  !> Haritadaki constraint-independent fiziksel DOF sayısını döndürür.
  pure function get_physical_dof_count(mapping) result(physical_dof_count)
    type(dof_map_t), intent(in) :: mapping
    integer :: physical_dof_count

    call validate_dof_map(mapping)
    physical_dof_count = mapping%physical_dof_count
  end function get_physical_dof_count

  !> Haritadaki constraint-independent full equation sayısını döndürür.
  pure function get_full_equation_count(mapping) result(full_equation_count)
    type(dof_map_t), intent(in) :: mapping
    integer :: full_equation_count

    call validate_dof_map(mapping)
    full_equation_count = mapping%full_equation_count
  end function get_full_equation_count

  !> Haritada tutulan pozitif active equation sayısını döndürür.
  pure function get_active_equation_count(mapping) &
      result(active_equation_count)
    type(dof_map_t), intent(in) :: mapping
    integer :: active_equation_count

    call validate_dof_map(mapping)
    active_equation_count = mapping%active_dof_count
  end function get_active_equation_count

  !> Haritadaki active DOF sayısını geriye uyumlu adla döndürür.
  pure function get_active_dof_count(mapping) result(active_dof_count)
    type(dof_map_t), intent(in) :: mapping
    integer :: active_dof_count

    active_dof_count = get_active_equation_count(mapping)
  end function get_active_dof_count

  !> Active equation sırasına karşılık gelen full matrix logical indekslerini
  !! döndürür.
  !!
  !! Matematiksel açıklama: Sonuç retained(active_id)=full_id bağıntısını
  !! sağlar ve P^T*A*P principal-submatrix indirgemesini tanımlar.
  !! Girdi dof_map_t, çıktı boyutsuz full equation indeks dizisidir. Tamamen
  !! kısıtlı legacy harita için sıfır uzunluklu dizi döner.
  pure function get_retained_full_equation_indices(mapping) result(indices)
    type(dof_map_t), intent(in) :: mapping
    integer, allocatable :: indices(:)

    integer :: entry_index

    call validate_dof_map(mapping)
    allocate(indices(mapping%active_dof_count))

    do entry_index = 1, size(mapping%entries)
      if (mapping%entries(entry_index)%active_equation_id > 0) then
        indices(mapping%entries(entry_index)%active_equation_id) = &
          mapping%entries(entry_index)%full_equation_id
      end if
    end do
  end function get_retained_full_equation_indices

  !> Haritada tutulan fiziksel DOF kaydı sayısını döndürür.
  pure function get_dof_map_entry_count(mapping) result(entry_count)
    type(dof_map_t), intent(in) :: mapping
    integer :: entry_count

    call validate_dof_map(mapping)
    entry_count = size(mapping%entries)
  end function get_dof_map_entry_count

  !> Harita kaydının ekleme sırasındaki bağımsız bir kopyasını döndürür.
  !! Girdi bir tabanlı boyutsuz indeks, çıktı bütün DOF/equation kimliklerini
  !! taşıyan dof_map_entry_t değeridir.
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

  !> Haritadaki bütün DOF kayıtlarının bağımsız kopyasını tek doğrulama ile
  !! döndürür. Çıktı fiziksel ekleme sıralı dof_map_entry_t dizisidir.
  !! Bu toplu erişim, büyük modellerde her kayıt için O(n^2) bütünlük
  !! doğrulamasının tekrarlanmasını önler; private depolamayı değiştirmez.
  pure function get_dof_map_entries(mapping) result(entries)
    type(dof_map_t), intent(in) :: mapping
    type(dof_map_entry_t), allocatable :: entries(:)

    call validate_dof_map(mapping)
    entries = mapping%entries
  end function get_dof_map_entries

  !> Haritanın constraint-blind full modda oluşturulup oluşturulmadığını
  !! boyutsuz mantıksal değer olarak döndürür.
  pure function is_full_dof_map(mapping) result(is_full)
    type(dof_map_t), intent(in) :: mapping
    logical :: is_full

    call validate_dof_map(mapping)
    is_full = mapping%full_map_mode
  end function is_full_dof_map

  !> Haritanın constraint-blind full equation modunda olduğunu geriye uyumlu
  !! olmayan açık equation terimiyle bildirir. is_full_dof_map ile eşdeğerdir.
  pure function is_full_equation_map(mapping) result(is_full)
    type(dof_map_t), intent(in) :: mapping
    logical :: is_full

    is_full = is_full_dof_map(mapping)
  end function is_full_equation_map

  !> Fiziksel düğüm kimliğini zorunlu olarak harita indeksine dönüştürür.
  !! Bulunamayan kimlik error stop üretir; kısıt bilgisiyle karışmaz.
  pure function require_dof_map_entry_index(mapping, node_id, dof_type) &
      result(entry_index)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in), optional :: dof_type
    integer :: entry_index

    if (.not. allocated(mapping%entries)) then
      error stop "DOF haritası kullanılmadan önce başlatılmalıdır."
    end if

    entry_index = find_dof_map_entry_index(mapping, node_id, dof_type)
    if (entry_index == 0) then
      error stop "Fiziksel düğüm kimliği DOF haritasında bulunamadı."
    end if
  end function require_dof_map_entry_index

  !> Fiziksel düğüm kimliğinin haritadaki bir tabanlı indeksini bulur.
  !! Bulunamayan kimlik için sıfır döndürür. Bu özel topoloji yardımcı yordamı
  !! fiziksel veya matematiksel katsayı hesabı yapmaz.
  pure function find_dof_map_entry_index(mapping, node_id, dof_type) &
      result(entry_index)
    type(dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in), optional :: dof_type
    integer :: entry_index

    integer :: current_index
    integer :: requested_dof_type

    requested_dof_type = TORSIONAL_ROTATION
    if (present(dof_type)) requested_dof_type = dof_type

    entry_index = 0
    if (.not. allocated(mapping%entries)) return

    do current_index = 1, size(mapping%entries)
      if (mapping%entries(current_index)%node_id == node_id .and. &
          mapping%entries(current_index)%dof_type == requested_dof_type) then
        entry_index = current_index
        return
      end if
    end do
  end function find_dof_map_entry_index

end module tms_dof_map
