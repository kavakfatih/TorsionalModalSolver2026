module tms_constraint_manager
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_dof_types, only : physical_dof_t, TORSIONAL_ROTATION, &
    validate_physical_dof, are_physical_dofs_equal
  use tms_constraint_types, only : constraint_t, FIXED_CONSTRAINT, &
    validate_constraint
  use tms_torsional_node, only : torsional_node_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    get_torsional_node_count, get_torsional_node, &
    validate_torsional_system
  use tms_dof_map, only : dof_map_t, dof_map_entry_t, &
    get_dof_map_entry_count, get_dof_map_entries, &
    get_full_map_physical_dof_count => get_physical_dof_count, &
    get_full_equation_count, is_full_equation_map, &
    validate_dof_map_for_system
  implicit none
  private

  !> Doğrulanmış kinematik constraint kayıtlarının sahibi olan yönetici.
  !! Depolama private tutularak constraint kimliği ve fiziksel DOF
  !! benzersizliği yalnız public ekleme yordamı üzerinden korunur.
  type, public :: constraint_manager_t
    private
    type(constraint_t), allocatable :: constraints(:)
  end type constraint_manager_t

  !> Bir fiziksel DOF'un constraint sonrası aktif denklem eşlemesini taşır.
  !! Bu kayıt, full sistem koordinatı ile reduced sistem koordinatını birlikte
  !! korur; böylece matris indirgeme ve sonuç recovery aynı eşlemeyi kullanır.
  type, public :: active_dof_map_entry_t
    !> Fiziksel koordinatın (node_id,dof_type) anahtarı [-].
    type(physical_dof_t) :: physical_dof

    !> Fiziksel DOF koleksiyonundaki benzersiz pozitif kimlik [-].
    integer :: physical_dof_id = 0

    !> Constraint'ten bağımsız full sistem denklem kimliği [-].
    integer :: full_equation_id = 0

    !> Reduced sistem denklem kimliği [-]. Kısıtlı DOF için sıfır, serbest
    !! DOF için kesintisiz 1..n_active aralığındadır.
    integer :: active_equation_id = -1

    !> Fiziksel koordinatın kinematik olarak prescribed olup olmadığı [-].
    logical :: constrained = .false.

    !> Toplam fiziksel durum recovery'sinde yeniden eklenecek açı [rad].
    !! Modal pertürbasyon recovery'sinde bu değer kullanılmaz.
    real(dp) :: prescribed_value = 0.0_dp

    !> Uygulanan constraint kaydının kimliği [-]; serbest DOF için sıfırdır.
    integer :: constraint_id = 0

    !> Uygulanan constraint türü [-]; serbest DOF için sıfırdır.
    integer :: constraint_type = 0
  end type active_dof_map_entry_t

  !> Full equation uzayından constraint sonrası aktif equation uzayına geçişi
  !! temsil eder. Kayıtlar private depolanır ve fiziksel numaralandırma,
  !! constraint metadata'sı ile recovery bilgisi birlikte doğrulanır.
  type, public :: active_dof_map_t
    private
    type(active_dof_map_entry_t), allocatable :: entries(:)
    integer :: active_equation_count = 0
  end type active_dof_map_t

  public :: initialize_constraint_manager
  public :: add_constraint
  public :: initialize_constraint_manager_from_system
  public :: validate_constraint_manager
  public :: get_constraint_count
  public :: get_constraint
  public :: generate_active_dof_map
  public :: validate_active_dof_map
  public :: get_active_equation_count
  public :: get_physical_dof_count
  public :: get_active_dof_map_entry
  public :: get_active_dof_map_entries
  public :: lookup_active_equation_id
  public :: lookup_full_equation_id_for_active_map
  public :: get_retained_full_equation_ids

contains

  !> Constraint yöneticisini geçerli ve boş bir koleksiyonla başlatır.
  !! Girdi yoktur; çıktı sıfır constraint içeren constraint_manager_t'dir.
  !! İşlem fiziksel denklem veya matris üretmez.
  pure subroutine initialize_constraint_manager(manager)
    type(constraint_manager_t), intent(out) :: manager

    allocate(manager%constraints(0))
  end subroutine initialize_constraint_manager

  !> Sistemde bulunan tek bir fiziksel DOF için doğrulanmış constraint ekler.
  !!
  !! Fiziksel açıklama: Constraint, bir torsional düğümün dönme koordinatına
  !! fixed veya prescribed açı sınır koşulu uygular. Düğümün eski
  !! constrained=.true. işareti varsa yalnız matching fixed-zero kayıt kabul
  !! edilir; açık manager diğer düğümler için yetkili constraint kaynağıdır.
  !! Matematiksel açıklama: constraint_id ve (node_id,dof_type) anahtarı ayrı
  !! ayrı benzersiz olmalıdır. Aynı fiziksel DOF'a iki koşul uygulanmaz.
  !! Girdiler: manager, constraint ve doğrulanmış torsional sistem. Constraint
  !! value alanı torsional dönme için [rad] birimindedir. Çıktı güncellenmiş
  !! manager'dır; eksik düğüm veya çakışma error stop ile reddedilir.
  pure subroutine add_constraint(manager, constraint, system)
    type(constraint_manager_t), intent(inout) :: manager
    type(constraint_t), intent(in) :: constraint
    type(torsional_system_t), intent(in) :: system

    type(torsional_node_t) :: node
    integer :: constraint_index
    integer :: node_index

    call require_initialized_manager(manager)
    call validate_torsional_system(system)
    call validate_constraint(constraint)

    node_index = find_system_node_index(system, constraint%node_id)
    if (node_index == 0) then
      error stop "Constraint sistemde bulunmayan bir düğüme uygulanamaz."
    end if

    node = get_torsional_node(system, node_index)
    if (node%constrained .and. &
        constraint%constraint_type /= FIXED_CONSTRAINT) then
      error stop &
        "Legacy kısıtlı düğüm matching fixed-zero constraint gerektirir."
    end if

    do constraint_index = 1, size(manager%constraints)
      if (manager%constraints(constraint_index)%constraint_id == &
          constraint%constraint_id) then
        error stop "Constraint kimlikleri benzersiz olmalıdır."
      end if

      if (manager%constraints(constraint_index)%node_id == &
          constraint%node_id .and. &
          manager%constraints(constraint_index)%dof_type == &
          constraint%dof_type) then
        error stop "Aynı fiziksel DOF için birden çok constraint tanımlanamaz."
      end if
    end do

    manager%constraints = [manager%constraints, constraint]
  end subroutine add_constraint

  !> Eski node%constrained işaretlerini explicit fixed-zero constraint
  !! kayıtlarına dönüştüren geriye uyumluluk adaptörünü oluşturur.
  !!
  !! Fiziksel açıklama: constrained=.true. olan her torsional dönme DOF'u
  !! theta=0 rad ile sabitlenir; serbest düğüme constraint eklenmez.
  !! Matematiksel açıklama: Constraint kimlikleri düğüm ekleme sırasına göre
  !! kesintisiz 1..n_constraint olarak atanır. Girdi torsional sistem, çıktı
  !! başlatılmış constraint_manager_t'dir. Prescribed nonzero değer türetilmez.
  pure subroutine initialize_constraint_manager_from_system(manager, system)
    type(constraint_manager_t), intent(out) :: manager
    type(torsional_system_t), intent(in) :: system

    type(constraint_t) :: legacy_constraint
    type(torsional_node_t) :: node
    integer :: constraint_id
    integer :: node_index

    call validate_torsional_system(system)
    call initialize_constraint_manager(manager)
    constraint_id = 0

    do node_index = 1, get_torsional_node_count(system)
      node = get_torsional_node(system, node_index)
      if (.not. node%constrained) cycle

      constraint_id = constraint_id + 1
      legacy_constraint = constraint_t( &
        constraint_id=constraint_id, node_id=node%id, &
        dof_type=TORSIONAL_ROTATION, value=0.0_dp, &
        constraint_type=FIXED_CONSTRAINT)
      call add_constraint(manager, legacy_constraint, system)
    end do
  end subroutine initialize_constraint_manager_from_system

  !> Manager koleksiyonunu verilen sistem ve legacy kısıt sözleşmesiyle
  !! birlikte doğrular.
  !! Matematiksel açıklama: Tüm kayıtlar geçerli, kimlikler ve fiziksel DOF
  !! anahtarları benzersiz, düğümler sistemde mevcut olmalıdır. Eski
  !! constrained=.true. düğümlerin her biri matching fixed-zero kayıt taşır.
  !! Girdiler manager ve torsional sistem; çıktı yoktur.
  pure subroutine validate_constraint_manager(manager, system)
    type(constraint_manager_t), intent(in) :: manager
    type(torsional_system_t), intent(in) :: system

    type(torsional_node_t) :: node
    integer :: constraint_index
    integer :: first_index
    integer :: node_index
    integer :: second_index

    call require_initialized_manager(manager)
    call validate_torsional_system(system)

    do first_index = 1, size(manager%constraints)
      call validate_constraint(manager%constraints(first_index))
      if (find_system_node_index( &
          system, manager%constraints(first_index)%node_id) == 0) then
        error stop "Constraint yöneticisi sisteme ait olmayan düğüm içeriyor."
      end if

      do second_index = first_index + 1, size(manager%constraints)
        if (manager%constraints(first_index)%constraint_id == &
            manager%constraints(second_index)%constraint_id) then
          error stop "Constraint kimlikleri benzersiz olmalıdır."
        end if

        if (manager%constraints(first_index)%node_id == &
            manager%constraints(second_index)%node_id .and. &
            manager%constraints(first_index)%dof_type == &
            manager%constraints(second_index)%dof_type) then
          error stop &
            "Aynı fiziksel DOF için birden çok constraint tanımlanamaz."
        end if
      end do
    end do

    do node_index = 1, get_torsional_node_count(system)
      node = get_torsional_node(system, node_index)
      if (.not. node%constrained) cycle

      constraint_index = find_constraint_index( &
        manager, node%id, TORSIONAL_ROTATION)
      if (constraint_index == 0) then
        error stop &
          "Legacy kısıtlı düğüm matching fixed-zero constraint gerektirir."
      end if
      if (manager%constraints(constraint_index)%constraint_type /= &
          FIXED_CONSTRAINT .or. &
          abs(manager%constraints(constraint_index)%value) > 0.0_dp) then
        error stop &
          "Legacy kısıtlı düğüm matching fixed-zero constraint gerektirir."
      end if
    end do
  end subroutine validate_constraint_manager

  !> Yöneticideki doğrulanmış constraint kaydı sayısını döndürür [-].
  pure function get_constraint_count(manager) result(constraint_count)
    type(constraint_manager_t), intent(in) :: manager
    integer :: constraint_count

    call require_initialized_manager(manager)
    constraint_count = size(manager%constraints)
  end function get_constraint_count

  !> Constraint kaydının ekleme sırasındaki bağımsız kopyasını döndürür.
  !! Girdiler manager ve bir tabanlı boyutsuz indeks, çıktı constraint_t'dir.
  pure function get_constraint(manager, index) result(constraint)
    type(constraint_manager_t), intent(in) :: manager
    integer, intent(in) :: index
    type(constraint_t) :: constraint

    call require_initialized_manager(manager)
    if (index < 1 .or. index > size(manager%constraints)) then
      error stop "Constraint koleksiyonu indeksi geçersiz."
    end if

    constraint = manager%constraints(index)
  end function get_constraint

  !> Constraint'ten bağımsız full DOF haritasından aktif denklem haritasını
  !! üretir.
  !!
  !! Fiziksel açıklama: Prescribed torsional koordinatlar reduced sistemden
  !! çıkarılır; tüm fiziksel koordinatlar full denklem ve recovery bilgileriyle
  !! korunur.
  !! Matematiksel açıklama: Full equation ID'ler değiştirilmez. Serbest DOF'lar
  !! full harita sırasıyla kesintisiz active=1..n_active, constrained DOF'lar
  !! active=0 alır. İleride matris indirgeme A_r=P^T*A*P eşdeğerini kullanır.
  !! Girdiler: constraint_manager_t, constraint-blind full dof_map_t ve aynı
  !! torsional_system_t. Çıktı active_dof_map_t'dir. Prescribed value [rad]
  !! yalnız state recovery metadata'sıdır; K/M veya RHS burada değiştirilmez.
  pure subroutine generate_active_dof_map( &
      active_mapping, manager, full_mapping, system)
    type(active_dof_map_t), intent(out) :: active_mapping
    type(constraint_manager_t), intent(in) :: manager
    type(dof_map_t), intent(in) :: full_mapping
    type(torsional_system_t), intent(in) :: system

    type(constraint_t) :: constraint
    type(dof_map_entry_t) :: full_entry
    type(dof_map_entry_t), allocatable :: full_entries(:)
    type(physical_dof_t) :: physical_dof
    integer :: constraint_index
    integer :: entry_count
    integer :: entry_index

    call validate_constraint_manager(manager, system)

    if (.not. is_full_equation_map(full_mapping)) then
      error stop "Aktif DOF haritası constraint-blind full harita gerektirir."
    end if

    ! Full equation satır sırası, assembly yapılan sistemin fiziksel düğüm
    ! sırasıyla bire bir aynı olmalıdır. Yalnız node kümesi ve boyut kontrolü,
    ! başka sıradaki bir sistem haritasının yanlış K/M satırını elemesine izin
    ! verebilir.
    call validate_dof_map_for_system(full_mapping, system)

    entry_count = get_dof_map_entry_count(full_mapping)
    if (entry_count /= get_full_map_physical_dof_count(full_mapping) .or. &
        entry_count /= get_full_equation_count(full_mapping) .or. &
        entry_count /= get_torsional_node_count(system)) then
      error stop "Full DOF haritası ile torsional sistem boyutları uyumsuz."
    end if

    allocate(active_mapping%entries(entry_count))
    active_mapping%active_equation_count = 0
    full_entries = get_dof_map_entries(full_mapping)

    do entry_index = 1, entry_count
      full_entry = full_entries(entry_index)
      physical_dof = physical_dof_t( &
        node_id=full_entry%node_id, dof_type=full_entry%dof_type)
      call validate_physical_dof(physical_dof)

      if (find_system_node_index(system, physical_dof%node_id) == 0) then
        error stop "Full DOF haritası sisteme ait olmayan düğüm içeriyor."
      end if

      active_mapping%entries(entry_index)%physical_dof = physical_dof
      active_mapping%entries(entry_index)%physical_dof_id = &
        full_entry%physical_dof_id
      active_mapping%entries(entry_index)%full_equation_id = &
        full_entry%full_equation_id

      constraint_index = find_constraint_index( &
        manager, physical_dof%node_id, physical_dof%dof_type)
      if (constraint_index > 0) then
        constraint = manager%constraints(constraint_index)
        active_mapping%entries(entry_index)%active_equation_id = 0
        active_mapping%entries(entry_index)%constrained = .true.
        active_mapping%entries(entry_index)%prescribed_value = &
          constraint%value
        active_mapping%entries(entry_index)%constraint_id = &
          constraint%constraint_id
        active_mapping%entries(entry_index)%constraint_type = &
          constraint%constraint_type
      else
        active_mapping%active_equation_count = &
          active_mapping%active_equation_count + 1
        active_mapping%entries(entry_index)%active_equation_id = &
          active_mapping%active_equation_count
        active_mapping%entries(entry_index)%constrained = .false.
        active_mapping%entries(entry_index)%prescribed_value = 0.0_dp
        active_mapping%entries(entry_index)%constraint_id = 0
        active_mapping%entries(entry_index)%constraint_type = 0
      end if
    end do

    call validate_active_dof_map(active_mapping)
  end subroutine generate_active_dof_map

  !> Aktif DOF haritasının numaralandırma, constraint ve recovery
  !! invariantlarını doğrular.
  !!
  !! Matematiksel açıklama: Fiziksel DOF ve full equation kimlikleri pozitif,
  !! benzersiz ve 1..n aralığında kesintisizdir. Active kimlikler kısıtlı DOF
  !! için 0, serbest DOF için benzersiz 1..n_active değeridir. Prescribed açı
  !! [rad] sonludur. Girdi active_dof_map_t, çıktı yoktur.
  pure subroutine validate_active_dof_map(mapping)
    type(active_dof_map_t), intent(in) :: mapping

    type(constraint_t) :: reconstructed_constraint
    integer :: active_count
    integer :: equation_id
    integer :: first_index
    integer :: second_index

    if (.not. allocated(mapping%entries)) then
      error stop "Aktif DOF haritası kullanılmadan önce üretilmelidir."
    end if
    if (size(mapping%entries) == 0) then
      error stop "Aktif DOF haritası en az bir fiziksel DOF içermelidir."
    end if
    if (mapping%active_equation_count < 0 .or. &
        mapping%active_equation_count > size(mapping%entries)) then
      error stop "Aktif denklem sayısı geçersiz."
    end if

    active_count = 0
    do first_index = 1, size(mapping%entries)
      call validate_physical_dof(mapping%entries(first_index)%physical_dof)

      if (mapping%entries(first_index)%physical_dof_id < 1 .or. &
          mapping%entries(first_index)%physical_dof_id > &
          size(mapping%entries)) then
        error stop "Aktif haritadaki fiziksel DOF kimliği geçersiz."
      end if
      if (mapping%entries(first_index)%full_equation_id < 1 .or. &
          mapping%entries(first_index)%full_equation_id > &
          size(mapping%entries)) then
        error stop "Aktif haritadaki full denklem kimliği geçersiz."
      end if

      if (mapping%entries(first_index)%constrained) then
        if (mapping%entries(first_index)%active_equation_id /= 0) then
          error stop "Kısıtlı DOF aktif denklem kimliği sıfır olmalıdır."
        end if

        reconstructed_constraint = constraint_t( &
          constraint_id=mapping%entries(first_index)%constraint_id, &
          node_id=mapping%entries(first_index)%physical_dof%node_id, &
          dof_type=mapping%entries(first_index)%physical_dof%dof_type, &
          value=mapping%entries(first_index)%prescribed_value, &
          constraint_type=mapping%entries(first_index)%constraint_type)
        call validate_constraint(reconstructed_constraint)
      else
        active_count = active_count + 1
        if (mapping%entries(first_index)%active_equation_id < 1 .or. &
            mapping%entries(first_index)%active_equation_id > &
            mapping%active_equation_count) then
          error stop "Serbest DOF aktif denklem kimliği geçersiz."
        end if
        if (mapping%entries(first_index)%constraint_id /= 0 .or. &
            mapping%entries(first_index)%constraint_type /= 0 .or. &
            abs(mapping%entries(first_index)%prescribed_value) > 0.0_dp) then
          error stop "Serbest DOF constraint metadata'sı taşımamalıdır."
        end if
      end if

      if (.not. ieee_is_finite( &
          mapping%entries(first_index)%prescribed_value)) then
        error stop "Aktif harita prescribed değeri sonlu olmalıdır."
      end if

      do second_index = first_index + 1, size(mapping%entries)
        if (are_physical_dofs_equal( &
            mapping%entries(first_index)%physical_dof, &
            mapping%entries(second_index)%physical_dof)) then
          error stop "Aktif haritadaki fiziksel DOF'lar benzersiz olmalıdır."
        end if
        if (mapping%entries(first_index)%physical_dof_id == &
            mapping%entries(second_index)%physical_dof_id) then
          error stop "Aktif haritadaki fiziksel DOF kimlikleri benzersiz olmalıdır."
        end if
        if (mapping%entries(first_index)%full_equation_id == &
            mapping%entries(second_index)%full_equation_id) then
          error stop "Aktif haritadaki full denklem kimlikleri benzersiz olmalıdır."
        end if
        if (.not. mapping%entries(first_index)%constrained .and. &
            .not. mapping%entries(second_index)%constrained .and. &
            mapping%entries(first_index)%active_equation_id == &
            mapping%entries(second_index)%active_equation_id) then
          error stop "Aktif denklem kimlikleri benzersiz olmalıdır."
        end if
        if (mapping%entries(first_index)%constrained .and. &
            mapping%entries(second_index)%constrained .and. &
            mapping%entries(first_index)%constraint_id == &
            mapping%entries(second_index)%constraint_id) then
          error stop "Aktif haritadaki constraint kimlikleri benzersiz olmalıdır."
        end if
      end do
    end do

    if (active_count /= mapping%active_equation_count) then
      error stop "Aktif haritanın serbest DOF sayısı tutarsız."
    end if

    do equation_id = 1, size(mapping%entries)
      if (.not. any( &
          mapping%entries%physical_dof_id == equation_id)) then
        error stop "Fiziksel DOF kimlikleri kesintisiz olmalıdır."
      end if
      if (.not. any( &
          mapping%entries%full_equation_id == equation_id)) then
        error stop "Full denklem kimlikleri kesintisiz olmalıdır."
      end if
    end do

    do equation_id = 1, mapping%active_equation_count
      if (.not. any( &
          mapping%entries%active_equation_id == equation_id)) then
        error stop "Aktif denklem kimlikleri kesintisiz olmalıdır."
      end if
    end do
  end subroutine validate_active_dof_map

  !> Reduced sistemdeki aktif denklem sayısını döndürür [-]. Tamamen
  !! kısıtlanmış geçerli sistem için sıfır döndürür.
  pure function get_active_equation_count(mapping) result(equation_count)
    type(active_dof_map_t), intent(in) :: mapping
    integer :: equation_count

    call validate_active_dof_map(mapping)
    equation_count = mapping%active_equation_count
  end function get_active_equation_count

  !> Aktif haritada korunan toplam fiziksel DOF sayısını döndürür [-].
  pure function get_physical_dof_count(mapping) result(dof_count)
    type(active_dof_map_t), intent(in) :: mapping
    integer :: dof_count

    call validate_active_dof_map(mapping)
    dof_count = size(mapping%entries)
  end function get_physical_dof_count

  !> Aktif harita kaydının fiziksel sıra indeksindeki bağımsız kopyasını
  !! döndürür. Girdi bir tabanlı boyutsuz indeks, çıktı mapping kaydıdır.
  pure function get_active_dof_map_entry(mapping, index) result(entry)
    type(active_dof_map_t), intent(in) :: mapping
    integer, intent(in) :: index
    type(active_dof_map_entry_t) :: entry

    call validate_active_dof_map(mapping)
    if (index < 1 .or. index > size(mapping%entries)) then
      error stop "Aktif DOF haritası kayıt indeksi geçersiz."
    end if

    entry = mapping%entries(index)
  end function get_active_dof_map_entry

  !> Bütün active-map kayıtlarının bağımsız kopyasını tek doğrulama ile
  !! döndürür. Çıktı physical_dof_id sıralı active_dof_map_entry_t dizisidir.
  !! Toplu erişim, result recovery sırasında her kayıt için bütün haritanın
  !! yeniden doğrulanmasını önler; constraint metadata'sını değiştirmez.
  pure function get_active_dof_map_entries(mapping) result(entries)
    type(active_dof_map_t), intent(in) :: mapping
    type(active_dof_map_entry_t), allocatable :: entries(:)

    call validate_active_dof_map(mapping)
    entries = mapping%entries
  end function get_active_dof_map_entries

  !> Fiziksel (node_id,dof_type) anahtarının active equation kimliğini
  !! döndürür [-]. Kısıtlı DOF sıfır döndürür; haritada olmayan DOF error stop
  !! üretir ve böylece eksik ile kısıtlı durumları ayırır.
  pure function lookup_active_equation_id(mapping, node_id, dof_type) &
      result(active_equation_id)
    type(active_dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in) :: dof_type
    integer :: active_equation_id

    integer :: entry_index

    call validate_active_dof_map(mapping)
    entry_index = find_active_map_entry_index(mapping, node_id, dof_type)
    if (entry_index == 0) then
      error stop "Fiziksel DOF aktif haritada bulunamadı."
    end if

    active_equation_id = &
      mapping%entries(entry_index)%active_equation_id
  end function lookup_active_equation_id

  !> Fiziksel (node_id,dof_type) anahtarının constraint'ten bağımsız full
  !! equation kimliğini döndürür [-]. Kısıtlı ve serbest DOF'lar için pozitiftir.
  pure function lookup_full_equation_id_for_active_map( &
      mapping, node_id, dof_type) result(full_equation_id)
    type(active_dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in) :: dof_type
    integer :: full_equation_id

    integer :: entry_index

    call validate_active_dof_map(mapping)
    entry_index = find_active_map_entry_index(mapping, node_id, dof_type)
    if (entry_index == 0) then
      error stop "Fiziksel DOF aktif haritada bulunamadı."
    end if

    full_equation_id = mapping%entries(entry_index)%full_equation_id
  end function lookup_full_equation_id_for_active_map

  !> Reduced matris sırasına karşılık gelen full equation kimliklerini döndürür.
  !! Matematiksel anlam: Sonuç active equation a için seçim/prolongation
  !! matrisindeki P(full_equation_ids(a),a)=1 konumudur. Girdi aktif harita,
  !! çıktı n_active uzunluklu boyutsuz tamsayı dizisidir; tüm DOF'lar kısıtlıysa
  !! sıfır uzunluklu dizi döner.
  pure function get_retained_full_equation_ids(mapping) &
      result(full_equation_ids)
    type(active_dof_map_t), intent(in) :: mapping
    integer, allocatable :: full_equation_ids(:)

    integer :: entry_index

    call validate_active_dof_map(mapping)
    allocate(full_equation_ids(mapping%active_equation_count))

    do entry_index = 1, size(mapping%entries)
      if (mapping%entries(entry_index)%constrained) cycle
      full_equation_ids( &
        mapping%entries(entry_index)%active_equation_id) = &
        mapping%entries(entry_index)%full_equation_id
    end do
  end function get_retained_full_equation_ids

  !> Manager depolamasının initialize edildiğini doğrular.
  pure subroutine require_initialized_manager(manager)
    type(constraint_manager_t), intent(in) :: manager

    if (.not. allocated(manager%constraints)) then
      error stop "Constraint yöneticisi kullanılmadan önce başlatılmalıdır."
    end if
  end subroutine require_initialized_manager

  !> Verilen constraint kimliğinin manager içindeki indeksini bulur.
  !! Bulunamayan kayıt için sıfır döndürür; fizik hesabı yapmaz.
  pure function find_constraint_index(manager, node_id, dof_type) &
      result(constraint_index)
    type(constraint_manager_t), intent(in) :: manager
    integer, intent(in) :: node_id
    integer, intent(in) :: dof_type
    integer :: constraint_index

    integer :: current_index

    constraint_index = 0
    if (.not. allocated(manager%constraints)) return

    do current_index = 1, size(manager%constraints)
      if (manager%constraints(current_index)%node_id == node_id .and. &
          manager%constraints(current_index)%dof_type == dof_type) then
        constraint_index = current_index
        return
      end if
    end do
  end function find_constraint_index

  !> Düğüm kimliğinin torsional sistem koleksiyonundaki indeksini bulur.
  !! Bulunamayan kimlik için sıfır döndürür; fizik hesabı yapmaz.
  pure function find_system_node_index(system, node_id) result(node_index)
    type(torsional_system_t), intent(in) :: system
    integer, intent(in) :: node_id
    integer :: node_index

    type(torsional_node_t) :: node
    integer :: current_index

    node_index = 0
    do current_index = 1, get_torsional_node_count(system)
      node = get_torsional_node(system, current_index)
      if (node%id == node_id) then
        node_index = current_index
        return
      end if
    end do
  end function find_system_node_index

  !> Fiziksel DOF anahtarının aktif haritadaki indeksini bulur.
  !! Bulunamayan anahtar için sıfır döndürür; fizik hesabı yapmaz.
  pure function find_active_map_entry_index(mapping, node_id, dof_type) &
      result(entry_index)
    type(active_dof_map_t), intent(in) :: mapping
    integer, intent(in) :: node_id
    integer, intent(in) :: dof_type
    integer :: entry_index

    integer :: current_index

    entry_index = 0
    if (.not. allocated(mapping%entries)) return

    do current_index = 1, size(mapping%entries)
      if (mapping%entries(current_index)%physical_dof%node_id == node_id .and. &
          mapping%entries(current_index)%physical_dof%dof_type == dof_type) then
        entry_index = current_index
        return
      end if
    end do
  end function find_active_map_entry_index

end module tms_constraint_manager
