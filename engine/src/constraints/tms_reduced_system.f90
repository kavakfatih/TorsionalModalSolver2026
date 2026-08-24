module tms_reduced_system
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_generalized_torsional_system, only : torsional_system_t
  use tms_dof_map, only : dof_map_t, initialize_full_dof_map
  use tms_constraint_manager, only : constraint_manager_t, &
    active_dof_map_t, active_dof_map_entry_t, generate_active_dof_map, &
    validate_active_dof_map, get_active_equation_count, &
    get_active_dof_map_entries, &
    get_retained_full_equation_ids
  use tms_stiffness_matrix, only : stiffness_matrix_t, &
    validate_stiffness_matrix, get_stiffness_matrix_size
  use tms_mass_matrix, only : mass_matrix_t, validate_mass_matrix, &
    get_mass_matrix_size
  use tms_matrix_assembly, only : assemble_full_stiffness, &
    assemble_full_inertia
  use tms_matrix_reduction, only : reduce_matrix
  implicit none
  private

  !> Constraint uygulanmış torsional sistemin solver girdilerini ve sonuç
  !! recovery bilgisini birlikte taşır. Reduced K [N*m/rad] ve reduced M
  !! [kg*m^2] matrislerinin satır sırası active_mapping ile tanımlanır.
  !! Depolama private tutulur; böylece ilerideki sparse matrix backend veya
  !! eigen solver değişikliği public veri bütünlüğünü bozamaz.
  type, public :: reduced_torsional_system_t
    private
    type(stiffness_matrix_t) :: stiffness
    type(mass_matrix_t) :: mass
    type(active_dof_map_t) :: active_mapping
  end type reduced_torsional_system_t

  public :: build_reduced_torsional_system
  public :: validate_reduced_torsional_system
  public :: get_reduced_stiffness
  public :: get_reduced_mass
  public :: get_reduced_active_dof_map
  public :: get_reduced_active_dof_count
  public :: recover_physical_state
  public :: recover_mode_shape

contains

  !> Full torsional sistemi constraint sonrası principal alt sisteme indirger.
  !!
  !! Fiziksel açıklama: Önce bütün fiziksel dönme DOF'ları için full K ve M
  !! oluşturulur; fixed veya prescribed DOF'lar daha sonra solver denklem
  !! takımından çıkarılır.
  !! Matematiksel açıklama: Boyutsuz seçim matrisi P açıkça depolanmadan
  !! K_r=P^T*K*P ve M_r=P^T*M*P principal alt matrisleri oluşturulur.
  !! Girdiler: Geçerli torsional_system_t ve constraint_manager_t. Çıktı:
  !! K_r [N*m/rad], M_r [kg*m^2] ve physical/full/active eşlemeyi taşıyan
  !! reduced_torsional_system_t.
  !! Varsayımlar ve geçerlilik: Direct elimination uygulanır. Prescribed
  !! sıfır olmayan değerlerin yük vektörü düzeltmesi bu sürümün kapsamında
  !! değildir; değer yalnız fiziksel durum recovery bilgisi olarak korunur.
  pure function build_reduced_torsional_system(system, manager) &
      result(reduced_system)
    type(torsional_system_t), intent(in) :: system
    type(constraint_manager_t), intent(in) :: manager
    type(reduced_torsional_system_t) :: reduced_system

    type(dof_map_t) :: full_mapping
    type(stiffness_matrix_t) :: full_stiffness
    type(mass_matrix_t) :: full_mass
    integer, allocatable :: retained_equation_ids(:)

    call initialize_full_dof_map(full_mapping, system)
    full_stiffness = assemble_full_stiffness(system, full_mapping)
    full_mass = assemble_full_inertia(system, full_mapping)

    call generate_active_dof_map( &
      reduced_system%active_mapping, manager, full_mapping, system)
    retained_equation_ids = get_retained_full_equation_ids( &
      reduced_system%active_mapping)

    reduced_system%stiffness = reduce_matrix( &
      full_stiffness, retained_equation_ids)
    reduced_system%mass = reduce_matrix(full_mass, retained_equation_ids)
    call validate_reduced_torsional_system(reduced_system)
  end function build_reduced_torsional_system

  !> Reduced sistemin matris ve active equation boyut bütünlüğünü doğrular.
  !! Matematiksel açıklama: K_r ve M_r kare olmalı ve boyutları active DOF
  !! sayısına eşit olmalıdır. Katsayı birimleri sırasıyla [N*m/rad] ve
  !! [kg*m^2]'dir. Geçersiz veri error stop ile reddedilir.
  pure subroutine validate_reduced_torsional_system(reduced_system)
    type(reduced_torsional_system_t), intent(in) :: reduced_system

    integer :: active_count

    call validate_active_dof_map(reduced_system%active_mapping)
    call validate_stiffness_matrix(reduced_system%stiffness)
    call validate_mass_matrix(reduced_system%mass)

    active_count = get_active_equation_count( &
      reduced_system%active_mapping)
    if (get_stiffness_matrix_size(reduced_system%stiffness) /= &
        active_count .or. &
        get_mass_matrix_size(reduced_system%mass) /= active_count) then
      error stop "Reduced K/M boyutu active equation sayısıyla uyumsuz."
    end if
  end subroutine validate_reduced_torsional_system

  !> Reduced torsional rijitlik matrisinin bağımsız kopyasını döndürür.
  !! Çıktı active equation sıralı K_r [N*m/rad] değeridir.
  pure function get_reduced_stiffness(reduced_system) result(stiffness)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    type(stiffness_matrix_t) :: stiffness

    call validate_reduced_torsional_system(reduced_system)
    stiffness = reduced_system%stiffness
  end function get_reduced_stiffness

  !> Reduced torsional atalet matrisinin bağımsız kopyasını döndürür.
  !! Çıktı active equation sıralı M_r [kg*m^2] değeridir.
  pure function get_reduced_mass(reduced_system) result(mass)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    type(mass_matrix_t) :: mass

    call validate_reduced_torsional_system(reduced_system)
    mass = reduced_system%mass
  end function get_reduced_mass

  !> Physical DOF, full equation ve active equation ilişkisini taşıyan
  !! haritanın bağımsız kopyasını döndürür. Kimlikler boyutsuzdur.
  pure function get_reduced_active_dof_map(reduced_system) result(mapping)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    type(active_dof_map_t) :: mapping

    call validate_reduced_torsional_system(reduced_system)
    mapping = reduced_system%active_mapping
  end function get_reduced_active_dof_map

  !> Reduced sistemdeki kesintisiz active equation sayısını [-] döndürür.
  pure function get_reduced_active_dof_count(reduced_system) &
      result(active_count)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    integer :: active_count

    call validate_reduced_torsional_system(reduced_system)
    active_count = get_active_equation_count( &
      reduced_system%active_mapping)
  end function get_reduced_active_dof_count

  !> Reduced fiziksel durum vektörünü prescribed değerleri ekleyerek full
  !! physical DOF sırasına geri açar.
  !!
  !! Fiziksel açıklama: Active dönme değerleri [rad] korunur; constraint
  !! uygulanmış bileşenlere saklanan prescribed dönme [rad] yerleştirilir.
  !! Matematiksel açıklama: q=P*q_r+q_p. Girdi q_r uzunluğu n_active ve sonlu
  !! değerlerden oluşur. Çıktı physical_dof_id sıralı q vektörüdür.
  !! Varsayımlar ve geçerlilik: Bu yordam statik/kinematik durum recovery'si
  !! yapar; modal özvektör recovery'si için recover_mode_shape kullanılmalıdır.
  pure function recover_physical_state(reduced_system, reduced_state) &
      result(full_state)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    real(dp), intent(in) :: reduced_state(:)
    real(dp), allocatable :: full_state(:)

    type(active_dof_map_entry_t) :: entry
    type(active_dof_map_entry_t), allocatable :: entries(:)
    integer :: entry_index

    call validate_recovery_input(reduced_system, reduced_state)
    entries = get_active_dof_map_entries(reduced_system%active_mapping)
    allocate(full_state(size(entries)))
    full_state = 0.0_dp

    do entry_index = 1, size(full_state)
      entry = entries(entry_index)
      if (entry%constrained) then
        full_state(entry%physical_dof_id) = entry%prescribed_value
      else
        full_state(entry%physical_dof_id) = &
          reduced_state(entry%active_equation_id)
      end if
    end do
  end function recover_physical_state

  !> Reduced modal şekli constraint bileşenlerini sıfır yaparak full physical
  !! DOF sırasına geri açar.
  !!
  !! Fiziksel açıklama: Mode shape bir denge konumu etrafındaki değişimi
  !! temsil ettiğinden prescribed offset özvektöre eklenmez.
  !! Matematiksel açıklama: phi=P*phi_r; constrained phi bileşenleri sıfırdır.
  !! Girdi sonlu reduced modal genlikleri, çıktı physical_dof_id sıralı modal
  !! genliklerdir. Normalizasyon ve eigen çözümü bu yordamın kapsamı dışındadır.
  pure function recover_mode_shape(reduced_system, reduced_mode) &
      result(full_mode)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    real(dp), intent(in) :: reduced_mode(:)
    real(dp), allocatable :: full_mode(:)

    type(active_dof_map_entry_t) :: entry
    type(active_dof_map_entry_t), allocatable :: entries(:)
    integer :: entry_index

    call validate_recovery_input(reduced_system, reduced_mode)
    entries = get_active_dof_map_entries(reduced_system%active_mapping)
    allocate(full_mode(size(entries)))
    full_mode = 0.0_dp

    do entry_index = 1, size(full_mode)
      entry = entries(entry_index)
      if (.not. entry%constrained) then
        full_mode(entry%physical_dof_id) = &
          reduced_mode(entry%active_equation_id)
      end if
    end do
  end function recover_mode_shape

  !> Recovery girdisinin active equation boyutu ve sonluluğunu doğrular.
  !! Girdi reduced koordinat vektörüdür; birimi çağıran duruma göre [rad] veya
  !! modal genliktir. Geçersiz boyut, NaN veya sonsuz değer reddedilir.
  pure subroutine validate_recovery_input(reduced_system, reduced_values)
    type(reduced_torsional_system_t), intent(in) :: reduced_system
    real(dp), intent(in) :: reduced_values(:)

    call validate_reduced_torsional_system(reduced_system)
    if (size(reduced_values) /= get_active_equation_count( &
        reduced_system%active_mapping)) then
      error stop "Recovery vektörü active equation sayısıyla uyumsuz."
    end if

    if (.not. all(ieee_is_finite(reduced_values))) then
      error stop "Recovery vektörü yalnız sonlu değerler içermelidir."
    end if
  end subroutine validate_recovery_input

end module tms_reduced_system
