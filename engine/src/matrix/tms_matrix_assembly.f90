module tms_matrix_assembly
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t, &
    get_local_stiffness, get_local_loss_stiffness, get_local_damping
  use tms_generalized_torsional_system, only : torsional_system_t, &
    get_torsional_node_count, get_torsional_element_count, &
    get_torsional_node, get_torsional_element, validate_torsional_system
  use tms_dof_map, only : dof_map_t, validate_dof_map_for_system, &
    lookup_full_equation_id, get_full_equation_count, &
    get_retained_full_equation_indices
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    initialize_stiffness_matrix, add_local_stiffness
  use tms_mass_matrix, only : mass_matrix_t, initialize_mass_matrix, &
    add_nodal_inertia
  use tms_loss_stiffness_matrix, only : loss_stiffness_matrix_t, &
    initialize_loss_stiffness_matrix, add_local_loss_stiffness
  use tms_damping_matrix, only : damping_matrix_t, &
    initialize_damping_matrix, add_local_damping
  use tms_matrix_reduction, only : reduce_matrix
  implicit none
  private

  public :: assemble_stiffness
  public :: assemble_inertia
  public :: assemble_loss_stiffness
  public :: assemble_damping
  public :: assemble_full_stiffness
  public :: assemble_full_inertia
  public :: assemble_full_loss_stiffness
  public :: assemble_full_damping

contains

  !> Torsional elemanların rijitliğini geriye uyumlu active K matrisinde toplar.
  !!
  !! Fiziksel açıklama: Her lineer bağlantının geri çağırıcı moment katkısı,
  !! bağlı olduğu aktif dönel denklemlere eklenir.
  !! Matematiksel açıklama: Önce constraint-blind K_full oluşturulur, ardından
  !! retained full equation indeksleriyle K_active=P^T*K_full*P hesaplanır.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! n_active x n_active, [N*m/rad] birimli stiffness_matrix_t.
  !! Varsayımlar ve geçerlilik: Ortak global dönme ekseni ve aynı pozitif açı
  !! yönü kullanılır; ek koordinat dönüşümü yoktur. equation_id=0 satır/sütunu
  !! homojen theta=0 kısıtı olarak elenir. Bu yordam yalnız K' üretir; K'' ve C
  !! ayrı semantic assembly yordamlarıyla oluşturulur. Yük vektörü ve çözüm bu
  !! yordamın kapsamı dışındadır.
  pure function assemble_stiffness(system, mapping) result(global_stiffness)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(stiffness_matrix_t) :: global_stiffness

    type(stiffness_matrix_t) :: full_stiffness
    integer, allocatable :: retained_full_equation_indices(:)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    full_stiffness = assemble_full_stiffness(system, mapping)
    retained_full_equation_indices = &
      get_retained_full_equation_indices(mapping)
    global_stiffness = reduce_matrix( &
      full_stiffness, retained_full_equation_indices)
  end function assemble_stiffness

  !> Torsional elemanların lokal rijitlik katkılarını constraint-blind full K
  !! matrisinde toplar.
  !!
  !! Fiziksel açıklama: Her fiziksel torsional DOF, kinematik constraint
  !! durumundan bağımsız olarak full geri çağırıcı moment denkleminde tutulur.
  !! Matematiksel açıklama: Her eleman için node_id -> full_equation_id
  !! eşlemesi yapılır ve K_e=k[[1,-1],[-1,1]] katsayıları n_full x n_full
  !! K matrisine scatter-add ile eklenir.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! Katsayıları [N*m/rad] olan full stiffness_matrix_t.
  !! Varsayımlar ve geçerlilik: Constraint eliminasyonu ve prescribed değer
  !! etkisi uygulanmaz; sonuç yalnız fiziksel topoloji ve eleman rijitliğine
  !! bağlıdır.
  pure function assemble_full_stiffness(system, mapping) &
      result(full_stiffness)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(stiffness_matrix_t) :: full_stiffness

    type(torsional_element_t) :: element
    integer :: element_index
    integer :: full_equation_ids(2)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    call initialize_stiffness_matrix( &
      full_stiffness, get_full_equation_count(mapping))

    do element_index = 1, get_torsional_element_count(system)
      element = get_torsional_element(system, element_index)
      full_equation_ids = [ &
        lookup_full_equation_id(mapping, element%node_i_id), &
        lookup_full_equation_id(mapping, element%node_j_id)]
      call add_local_stiffness( &
        full_stiffness, full_equation_ids, get_local_stiffness(element))
    end do
  end function assemble_full_stiffness

  !> Torsional elemanların K'' katkılarını geriye uyumlu active uzayda
  !! toplar.
  !!
  !! Fiziksel açıklama: Her pasif elemanın harmonik enerji kaybı
  !! katsayıları, kısıt sonrası aktif torsional denklemlerde korunur.
  !! K'' viskoz C ile birleştirilmez.
  !! Matematiksel açıklama: Önce K''_full assemble edilir, sonra
  !! K''_active=P^T*K''_full*P principal submatrix seçimi uygulanır.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! n_active x n_active K'' [N*m/rad] loss_stiffness_matrix_t.
  !! Varsayımlar: Legacy active map equation_id=0 ile homojen constraint
  !! uygular; canonical dynamic yol full assembly ve explicit reduction'dır.
  pure function assemble_loss_stiffness(system, mapping) &
      result(global_loss_stiffness)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(loss_stiffness_matrix_t) :: global_loss_stiffness

    type(loss_stiffness_matrix_t) :: full_loss_stiffness
    integer, allocatable :: retained_full_equation_indices(:)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    full_loss_stiffness = assemble_full_loss_stiffness(system, mapping)
    retained_full_equation_indices = &
      get_retained_full_equation_indices(mapping)
    global_loss_stiffness = reduce_matrix( &
      full_loss_stiffness, retained_full_equation_indices)
  end function assemble_loss_stiffness

  !> Torsional elemanların K'' katkılarını constraint-blind full matriste
  !! toplar.
  !!
  !! Fiziksel açıklama: Her fiziksel torsional DOF, kinematik constraint
  !! durumundan bağımsız kayıp rijitliği denkleminde tutulur.
  !! Matematiksel açıklama: node_id -> full_equation_id eşlemesiyle
  !! K''_e=k''[[1,-1],[-1,1]] katkıları scatter-add edilir.
  !! Girdiler: Geçerli torsional_system_t ve full dof_map_t. Çıktı:
  !! n_full x n_full K''_full [N*m/rad].
  !! Varsayımlar: Constraint eliminasyonu, viskoz C ve frekans çarpanı bu
  !! yordamda uygulanmaz; eleman K'' değerleri frozen-property kabul edilir.
  pure function assemble_full_loss_stiffness(system, mapping) &
      result(full_loss_stiffness)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(loss_stiffness_matrix_t) :: full_loss_stiffness

    type(torsional_element_t) :: element
    integer :: element_index
    integer :: full_equation_ids(2)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    call initialize_loss_stiffness_matrix( &
      full_loss_stiffness, get_full_equation_count(mapping))

    do element_index = 1, get_torsional_element_count(system)
      element = get_torsional_element(system, element_index)
      full_equation_ids = [ &
        lookup_full_equation_id(mapping, element%node_i_id), &
        lookup_full_equation_id(mapping, element%node_j_id)]
      call add_local_loss_stiffness( &
        full_loss_stiffness, full_equation_ids, &
        get_local_loss_stiffness(element))
    end do
  end function assemble_full_loss_stiffness

  !> Torsional elemanların viskoz C katkılarını geriye uyumlu active
  !! uzayda toplar.
  !!
  !! Fiziksel açıklama: Bağıl açısal hıza karşı viskoz moment
  !! katsayıları kısıt sonrası aktif torsional denklemlerde korunur ve
  !! K'' kayıp rijitliğinden ayrı tutulur.
  !! Matematiksel açıklama: Önce C_full assemble edilir, sonra
  !! C_active=P^T*C_full*P principal submatrix seçimi uygulanır.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! n_active x n_active C [N*m*s/rad] damping_matrix_t.
  !! Varsayımlar: Legacy active map equation_id=0 ile homojen constraint
  !! uygular; canonical dynamic yol full assembly ve explicit reduction'dır.
  pure function assemble_damping(system, mapping) result(global_damping)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(damping_matrix_t) :: global_damping

    type(damping_matrix_t) :: full_damping
    integer, allocatable :: retained_full_equation_indices(:)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    full_damping = assemble_full_damping(system, mapping)
    retained_full_equation_indices = &
      get_retained_full_equation_indices(mapping)
    global_damping = reduce_matrix( &
      full_damping, retained_full_equation_indices)
  end function assemble_damping

  !> Torsional elemanların viskoz C katkılarını constraint-blind full
  !! matriste toplar.
  !!
  !! Fiziksel açıklama: Her fiziksel torsional DOF, kinematik constraint
  !! durumundan bağımsız viskoz sönüm denkleminde tutulur.
  !! Matematiksel açıklama: node_id -> full_equation_id eşlemesiyle
  !! C_e=c[[1,-1],[-1,1]] katkıları scatter-add edilir.
  !! Girdiler: Geçerli torsional_system_t ve full dof_map_t. Çıktı:
  !! n_full x n_full C_full [N*m*s/rad].
  !! Varsayımlar: Constraint eliminasyonu, K'' ve omega frekans çarpanı bu
  !! yordamda uygulanmaz; eleman c değerleri lineer ve frozen kabul edilir.
  pure function assemble_full_damping(system, mapping) result(full_damping)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(damping_matrix_t) :: full_damping

    type(torsional_element_t) :: element
    integer :: element_index
    integer :: full_equation_ids(2)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    call initialize_damping_matrix( &
      full_damping, get_full_equation_count(mapping))

    do element_index = 1, get_torsional_element_count(system)
      element = get_torsional_element(system, element_index)
      full_equation_ids = [ &
        lookup_full_equation_id(mapping, element%node_i_id), &
        lookup_full_equation_id(mapping, element%node_j_id)]
      call add_local_damping( &
        full_damping, full_equation_ids, get_local_damping(element))
    end do
  end function assemble_full_damping

  !> Torsional düğümlerde yığılmış polar ataletleri global M matrisinde toplar.
  !!
  !! Fiziksel açıklama: Her aktif dönel denklem, bağlı düğümün açısal ivmeye
  !! direncini temsil eden polar kütle ataleti J değerini taşır.
  !! Matematiksel açıklama: Önce constraint-blind M_full oluşturulur, ardından
  !! retained full equation indeksleriyle M_active=P^T*M_full*P hesaplanır.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! n_active x n_active, [kg*m^2] birimli mass_matrix_t.
  !! Varsayımlar ve geçerlilik: Ataletler düğümde yığılmıştır. Kısıtlı düğüm
  !! equation_id=0 ile indirgenmiş matristen çıkarılır. Eleman ataleti,
  !! tutarlı mass matrix, modal reduction ve eigen çözümü uygulanmaz.
  pure function assemble_inertia(system, mapping) result(global_mass)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(mass_matrix_t) :: global_mass

    type(mass_matrix_t) :: full_mass
    integer, allocatable :: retained_full_equation_indices(:)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    full_mass = assemble_full_inertia(system, mapping)
    retained_full_equation_indices = &
      get_retained_full_equation_indices(mapping)
    global_mass = reduce_matrix(full_mass, retained_full_equation_indices)
  end function assemble_inertia

  !> Düğümlerde yığılmış polar ataletleri constraint-blind full M matrisinde
  !! toplar.
  !!
  !! Fiziksel açıklama: Kısıtlı olsa bile her fiziksel torsional düğümün açısal
  !! ivmeye karşı polar atalet katkısı full hareket denkleminde korunur.
  !! Matematiksel açıklama: Her düğüm için
  !! M_full(full_equation_id,full_equation_id) += J uygulanır; bu modelde
  !! köşegen dışı katsayılar sıfırdır.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! Katsayıları [kg*m^2] olan full mass_matrix_t.
  !! Varsayımlar ve geçerlilik: Ataletler düğümde yığılmıştır; constraint
  !! eliminasyonu bu yordamda uygulanmaz.
  pure function assemble_full_inertia(system, mapping) result(full_mass)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(mass_matrix_t) :: full_mass

    type(torsional_node_t) :: node
    integer :: full_equation_id
    integer :: node_index

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    call initialize_mass_matrix( &
      full_mass, get_full_equation_count(mapping))

    do node_index = 1, get_torsional_node_count(system)
      node = get_torsional_node(system, node_index)
      full_equation_id = lookup_full_equation_id(mapping, node%id)
      call add_nodal_inertia( &
        full_mass, full_equation_id, node%polar_inertia_kg_m2)
    end do
  end function assemble_full_inertia

end module tms_matrix_assembly
