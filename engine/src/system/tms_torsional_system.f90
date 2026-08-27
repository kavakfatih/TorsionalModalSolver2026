module tms_generalized_torsional_system
  use tms_torsional_node, only : torsional_node_t, validate_torsional_node
  use tms_torsional_element, only : torsional_element_t, &
    validate_torsional_element
  implicit none
  private

  !> Genel ayrık-parametreli torsional sistemin düğüm ve eleman koleksiyonunu
  !! taşır. Koleksiyonlar yalnız public yönetim yordamlarıyla değiştirilerek
  !! kimlik ve bağlantı bütünlüğü korunur.
  type, public :: torsional_system_t
    private
    type(torsional_node_t), allocatable :: nodes(:)
    type(torsional_element_t), allocatable :: elements(:)
  end type torsional_system_t

  public :: add_torsional_node
  public :: add_torsional_element
  public :: count_torsional_dofs
  public :: get_torsional_node_count
  public :: get_torsional_element_count
  public :: get_torsional_node
  public :: get_torsional_element
  public :: validate_torsional_system

contains

  !> Genel torsional sisteme doğrulanmış bir fiziksel düğüm ekler.
  !!
  !! Fiziksel açıklama: Düğüm, yığılmış polar atalet ve dönel sınır koşulunu
  !! sisteme dahil eder.
  !! Matematiksel açıklama: Yeni düğümün pozitif kimlik ve atalet koşulları
  !! doğrulanır; kimlik mevcut düğüm kümesinde benzersiz olmalıdır.
  !! Girdiler: Değiştirilecek torsional_system_t ve J [kg*m^2], theta_0 [rad]
  !! alanlarını taşıyan torsional_node_t değeridir. Çıktı, güncellenmiş sistem
  !! koleksiyonudur.
  !! Varsayımlar ve geçerlilik: Ekleme sırası korunur. Geçersiz veya yinelenen
  !! kimlik error stop ile reddedilir; matris veya çözüm üretilmez.
  pure subroutine add_torsional_node(system, node)
    type(torsional_system_t), intent(inout) :: system
    type(torsional_node_t), intent(in) :: node

    integer :: node_index

    call validate_torsional_node(node)

    node_index = find_torsional_node_index(system, node%id)
    if (node_index /= 0) then
      error stop "Torsional sistemde düğüm kimlikleri benzersiz olmalıdır."
    end if

    if (allocated(system%nodes)) then
      system%nodes = [system%nodes, node]
    else
      allocate(system%nodes(1))
      system%nodes(1) = node
    end if
  end subroutine add_torsional_node

  !> Genel torsional sisteme doğrulanmış bir bağlantı elemanı ekler.
  !!
  !! Fiziksel açıklama: Eleman, mevcut iki dönel düğüm arasındaki lineer
  !! depolama rijitliği, kayıp rijitliği ve viskoz sönüm bağlantısını temsil
  !! eder.
  !! Matematiksel açıklama: Eleman kimliği benzersiz, uç kimlikleri farklı ve
  !! sistemde mevcut, K' > 0, K'' >= 0 ve c >= 0 olmalıdır.
  !! Girdiler: Değiştirilecek torsional_system_t ile K' ve K'' [N*m/rad],
  !! c [N*m*s/rad] alanlarını taşıyan torsional_element_t değeridir. Çıktı,
  !! güncellenmiş sistem koleksiyonudur.
  !! Varsayımlar ve geçerlilik: Aynı düğüm çifti arasında farklı kimlikli
  !! paralel elemanlara izin verilir. Bu koleksiyon yordamı assembly yapmaz;
  !! global M/K ayrı tms_matrix_assembly modülünde oluşturulur.
  pure subroutine add_torsional_element(system, element)
    type(torsional_system_t), intent(inout) :: system
    type(torsional_element_t), intent(in) :: element

    integer :: element_index

    call validate_torsional_element(element)

    if (find_torsional_node_index(system, element%node_i_id) == 0 .or. &
        find_torsional_node_index(system, element%node_j_id) == 0) then
      error stop "Torsional eleman yalnız sistemde bulunan düğümleri bağlayabilir."
    end if

    element_index = find_torsional_element_index(system, element%id)
    if (element_index /= 0) then
      error stop "Torsional sistemde eleman kimlikleri benzersiz olmalıdır."
    end if

    if (allocated(system%elements)) then
      system%elements = [system%elements, element]
    else
      allocate(system%elements(1))
      system%elements(1) = element
    end if
  end subroutine add_torsional_element

  !> Sistemdeki aktif dönel serbestlik derecesi sayısını döndürür.
  !!
  !! Fiziksel açıklama: Sabitlenmiş düğümler dönemez ve aktif koordinat
  !! oluşturmaz; serbest düğümlerin her biri bir torsional DOF taşır.
  !! Matematiksel açıklama: n_dof = sum_i not(constrained_i).
  !! Girdi: torsional_system_t. Çıktı: Boyutsuz, sıfır veya pozitif DOF sayısı.
  !! Varsayımlar ve geçerlilik: Düğüm kimlikleri DOF sıra numarası değildir.
  !! Boş sistem için sonuç sıfırdır; sistem geçerliliği ayrı yordamla sınanır.
  pure function count_torsional_dofs(system) result(dof_count)
    type(torsional_system_t), intent(in) :: system
    integer :: dof_count

    integer :: node_index

    dof_count = 0
    if (.not. allocated(system%nodes)) return

    do node_index = 1, size(system%nodes)
      if (.not. system%nodes(node_index)%constrained) then
        dof_count = dof_count + 1
      end if
    end do
  end function count_torsional_dofs

  !> Sistemde kayıtlı düğüm sayısını döndürür.
  !! Girdi torsional_system_t, çıktı boyutsuz ve sıfır veya pozitif tam sayıdır.
  !! Koleksiyon tahsis edilmemişse sıfır döndürülür; fizik hesabı yapılmaz.
  pure function get_torsional_node_count(system) result(node_count)
    type(torsional_system_t), intent(in) :: system
    integer :: node_count

    if (allocated(system%nodes)) then
      node_count = size(system%nodes)
    else
      node_count = 0
    end if
  end function get_torsional_node_count

  !> Sistemde kayıtlı eleman sayısını döndürür.
  !! Girdi torsional_system_t, çıktı boyutsuz ve sıfır veya pozitif tam sayıdır.
  !! Koleksiyon tahsis edilmemişse sıfır döndürülür; fizik hesabı yapılmaz.
  pure function get_torsional_element_count(system) result(element_count)
    type(torsional_system_t), intent(in) :: system
    integer :: element_count

    if (allocated(system%elements)) then
      element_count = size(system%elements)
    else
      element_count = 0
    end if
  end function get_torsional_element_count

  !> Ekleme sırasındaki indeksle torsional düğümün bir kopyasını döndürür.
  !! Girdiler torsional_system_t ve bir tabanlı boyutsuz indeks, çıktı
  !! torsional_node_t değeridir. Geçersiz indeks error stop ile reddedilir.
  pure function get_torsional_node(system, index) result(node)
    type(torsional_system_t), intent(in) :: system
    integer, intent(in) :: index
    type(torsional_node_t) :: node

    if (.not. allocated(system%nodes)) then
      error stop "Torsional sistemde okunabilecek düğüm bulunmuyor."
    end if

    if (index < 1 .or. index > size(system%nodes)) then
      error stop "Torsional düğüm koleksiyonu indeksi geçersiz."
    end if

    node = system%nodes(index)
  end function get_torsional_node

  !> Ekleme sırasındaki indeksle torsional elemanın bir kopyasını döndürür.
  !! Girdiler torsional_system_t ve bir tabanlı boyutsuz indeks, çıktı
  !! torsional_element_t değeridir. Geçersiz indeks error stop ile reddedilir.
  pure function get_torsional_element(system, index) result(element)
    type(torsional_system_t), intent(in) :: system
    integer, intent(in) :: index
    type(torsional_element_t) :: element

    if (.not. allocated(system%elements)) then
      error stop "Torsional sistemde okunabilecek eleman bulunmuyor."
    end if

    if (index < 1 .or. index > size(system%elements)) then
      error stop "Torsional eleman koleksiyonu indeksi geçersiz."
    end if

    element = system%elements(index)
  end function get_torsional_element

  !> Genel torsional sistemin koleksiyon ve bağlantı bütünlüğünü doğrular.
  !!
  !! Fiziksel açıklama: Sistem en az bir geçerli yığılmış atalet düğümü taşır;
  !! her bağlantı mevcut iki farklı fiziksel düğüm arasında tanımlıdır.
  !! Matematiksel açıklama: Düğüm ve eleman kimlik kümeleri kendi içinde
  !! benzersizdir; her elemanın uçları düğüm kümesinin elemanlarıdır.
  !! Girdi: torsional_system_t. Çıktı üretmez; geçersiz topoloji error stop ile
  !! reddedilir.
  !! Varsayımlar ve geçerlilik: Tek düğümlü, elemasız; tamamen sabitlenmiş;
  !! paralel elemanlı veya ayrık alt sistemli topolojiler geçerli olabilir.
  !! Bu doğrulayıcı assembly, bağlantılılık veya rijit-cisim modu analizi
  !! yapmaz; global M/K ayrı tms_matrix_assembly modülünde oluşturulur.
  pure subroutine validate_torsional_system(system)
    type(torsional_system_t), intent(in) :: system

    integer :: first_index
    integer :: second_index

    if (.not. allocated(system%nodes)) then
      error stop "Torsional sistem en az bir düğüm içermelidir."
    end if

    if (size(system%nodes) == 0) then
      error stop "Torsional sistem en az bir düğüm içermelidir."
    end if

    do first_index = 1, size(system%nodes)
      call validate_torsional_node(system%nodes(first_index))
      do second_index = first_index + 1, size(system%nodes)
        if (system%nodes(first_index)%id == system%nodes(second_index)%id) then
          error stop "Torsional sistemde yinelenen düğüm kimliği bulundu."
        end if
      end do
    end do

    if (.not. allocated(system%elements)) return

    do first_index = 1, size(system%elements)
      call validate_torsional_element(system%elements(first_index))

      if (find_torsional_node_index( &
          system, system%elements(first_index)%node_i_id) == 0 .or. &
          find_torsional_node_index( &
          system, system%elements(first_index)%node_j_id) == 0) then
        error stop "Torsional sistem elemanında tanımsız düğüm başvurusu var."
      end if

      do second_index = first_index + 1, size(system%elements)
        if (system%elements(first_index)%id == &
            system%elements(second_index)%id) then
          error stop "Torsional sistemde yinelenen eleman kimliği bulundu."
        end if
      end do
    end do
  end subroutine validate_torsional_system

  !> Verilen düğüm kimliğinin koleksiyondaki bir tabanlı indeksini bulur.
  !! Girdiler sistem ve boyutsuz düğüm kimliğidir; bulunamazsa sıfır döndürür.
  !! Fizik veya matris hesabı yapmayan özel topoloji yardımcı yordamıdır.
  pure function find_torsional_node_index(system, node_id) result(node_index)
    type(torsional_system_t), intent(in) :: system
    integer, intent(in) :: node_id
    integer :: node_index

    integer :: current_index

    node_index = 0
    if (.not. allocated(system%nodes)) return

    do current_index = 1, size(system%nodes)
      if (system%nodes(current_index)%id == node_id) then
        node_index = current_index
        return
      end if
    end do
  end function find_torsional_node_index

  !> Verilen eleman kimliğinin koleksiyondaki bir tabanlı indeksini bulur.
  !! Girdiler sistem ve boyutsuz eleman kimliğidir; bulunamazsa sıfır döndürür.
  !! Fizik veya matris hesabı yapmayan özel topoloji yardımcı yordamıdır.
  pure function find_torsional_element_index(system, element_id) &
      result(element_index)
    type(torsional_system_t), intent(in) :: system
    integer, intent(in) :: element_id
    integer :: element_index

    integer :: current_index

    element_index = 0
    if (.not. allocated(system%elements)) return

    do current_index = 1, size(system%elements)
      if (system%elements(current_index)%id == element_id) then
        element_index = current_index
        return
      end if
    end do
  end function find_torsional_element_index

end module tms_generalized_torsional_system
