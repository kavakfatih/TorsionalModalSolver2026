module tms_matrix_assembly
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t, get_local_stiffness
  use tms_generalized_torsional_system, only : torsional_system_t, &
    get_torsional_node_count, get_torsional_element_count, &
    get_torsional_node, get_torsional_element, validate_torsional_system
  use tms_dof_map, only : dof_map_t, validate_dof_map_for_system, &
    lookup_equation_id, get_active_dof_count
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    initialize_stiffness_matrix, add_local_stiffness
  use tms_mass_matrix, only : mass_matrix_t, initialize_mass_matrix, &
    add_nodal_inertia
  implicit none
  private

  public :: assemble_stiffness
  public :: assemble_inertia

contains

  !> Torsional elemanların lokal rijitlik katkılarını global K matrisinde toplar.
  !!
  !! Fiziksel açıklama: Her lineer bağlantının geri çağırıcı moment katkısı,
  !! bağlı olduğu aktif dönel denklemlere eklenir.
  !! Matematiksel açıklama: Her eleman için node_id -> equation_id eşlemesi
  !! yapılır ve K_e=k[[1,-1],[-1,1]] katsayıları global K konumlarına saçılır.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! n_active x n_active, [N*m/rad] birimli stiffness_matrix_t.
  !! Varsayımlar ve geçerlilik: Ortak global dönme ekseni ve aynı pozitif açı
  !! yönü kullanılır; ek koordinat dönüşümü yoktur. equation_id=0 satır/sütunu
  !! homojen theta=0 kısıtı olarak elenir. Global C, K'', yük vektörü ve eigen
  !! çözümü bu yordamın kapsamı dışındadır.
  pure function assemble_stiffness(system, mapping) result(global_stiffness)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(stiffness_matrix_t) :: global_stiffness

    type(torsional_element_t) :: element
    integer :: element_index
    integer :: equation_ids(2)

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    call initialize_stiffness_matrix( &
      global_stiffness, get_active_dof_count(mapping))

    do element_index = 1, get_torsional_element_count(system)
      element = get_torsional_element(system, element_index)
      equation_ids = [ &
        lookup_equation_id(mapping, element%node_i_id), &
        lookup_equation_id(mapping, element%node_j_id)]
      call add_local_stiffness( &
        global_stiffness, equation_ids, get_local_stiffness(element))
    end do
  end function assemble_stiffness

  !> Torsional düğümlerde yığılmış polar ataletleri global M matrisinde toplar.
  !!
  !! Fiziksel açıklama: Her aktif dönel denklem, bağlı düğümün açısal ivmeye
  !! direncini temsil eden polar kütle ataleti J değerini taşır.
  !! Matematiksel açıklama: M(equation_id,equation_id) += J; bu ilk modelde
  !! köşegen dışı katsayılar sıfırdır.
  !! Girdiler: Geçerli torsional_system_t ve ona ait dof_map_t. Çıktı:
  !! n_active x n_active, [kg*m^2] birimli mass_matrix_t.
  !! Varsayımlar ve geçerlilik: Ataletler düğümde yığılmıştır. Kısıtlı düğüm
  !! equation_id=0 ile indirgenmiş matristen çıkarılır. Eleman ataleti,
  !! tutarlı mass matrix, modal reduction ve eigen çözümü uygulanmaz.
  pure function assemble_inertia(system, mapping) result(global_mass)
    type(torsional_system_t), intent(in) :: system
    type(dof_map_t), intent(in) :: mapping
    type(mass_matrix_t) :: global_mass

    type(torsional_node_t) :: node
    integer :: equation_id
    integer :: node_index

    call validate_torsional_system(system)
    call validate_dof_map_for_system(mapping, system)
    call initialize_mass_matrix(global_mass, get_active_dof_count(mapping))

    do node_index = 1, get_torsional_node_count(system)
      node = get_torsional_node(system, node_index)
      equation_id = lookup_equation_id(mapping, node%id)
      call add_nodal_inertia( &
        global_mass, equation_id, node%polar_inertia_kg_m2)
    end do
  end function assemble_inertia

end module tms_matrix_assembly
