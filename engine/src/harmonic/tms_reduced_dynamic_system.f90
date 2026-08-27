module tms_reduced_dynamic_system
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_generalized_torsional_system, only : torsional_system_t
  use tms_dof_map, only : dof_map_t, initialize_full_dof_map
  use tms_constraint_manager, only : constraint_manager_t, &
    active_dof_map_t, active_dof_map_entry_t, &
    get_active_dof_map_entries, get_retained_full_equation_ids
  use tms_stiffness_matrix, only : stiffness_matrix_t
  use tms_mass_matrix, only : mass_matrix_t
  use tms_loss_stiffness_matrix, only : loss_stiffness_matrix_t, &
    validate_loss_stiffness_matrix, get_loss_stiffness_matrix_size
  use tms_damping_matrix, only : damping_matrix_t, &
    validate_damping_matrix, get_damping_matrix_size
  use tms_matrix_assembly, only : assemble_full_loss_stiffness, &
    assemble_full_damping
  use tms_matrix_reduction, only : reduce_matrix
  use tms_reduced_system, only : reduced_torsional_system_t, &
    build_reduced_torsional_system, validate_reduced_torsional_system, &
    get_reduced_stiffness, get_reduced_mass, &
    get_reduced_active_dof_map, get_reduced_active_dof_count
  implicit none
  private

  !> Frequency-domain torsional analiz için constraint uygulanmış reel
  !! sistem matrislerini ve recovery bağlamını birlikte taşır.
  !!
  !! Fiziksel anlam: Depolama rijitliği K' [N*m/rad], kayıp rijitliği
  !! K'' [N*m/rad], viskoz sönüm C [N*m*s/rad] ve polar atalet M
  !! [kg*m^2] ayrı constitutive kanallardır.
  !! Matematiksel anlam: Aynı active-equation sırasındaki matrisler ileride
  !! Z=K'-omega^2*M+i*(K''+omega*C) oluşturmak için kullanılır.
  !! Depolama private tutulur. Mevcut V0.5 K'/M reduced-system nesnesi
  !! bileşim yoluyla korunur; modal public API değiştirilmez.
  type, public :: reduced_dynamic_torsional_system_t
    private
    type(reduced_torsional_system_t) :: base_system
    type(loss_stiffness_matrix_t) :: loss_stiffness
    type(damping_matrix_t) :: damping
  end type reduced_dynamic_torsional_system_t

  public :: build_reduced_dynamic_torsional_system
  public :: validate_reduced_dynamic_torsional_system
  public :: get_reduced_dynamic_stiffness
  public :: get_reduced_dynamic_loss_stiffness
  public :: get_reduced_dynamic_damping
  public :: get_reduced_dynamic_mass
  public :: get_reduced_dynamic_active_dof_map
  public :: get_reduced_dynamic_active_dof_count
  public :: recover_harmonic_response

contains

  !> Full K', K'', C ve M matrislerini aynı constraint haritasıyla indirger.
  !!
  !! Fiziksel açıklama: Tüm fiziksel torsional DOF katkıları constraint
  !! uygulanmadan assemble edilir. K', K'', C ve M birbirine dönüştürülmeden
  !! aynı active torsional koordinatlara seçilir.
  !! Matematiksel açıklama: retained indeksleriyle A_r=P^T*A_full*P,
  !! A elemanı {K',K'',C,M} kümesindedir. Mevcut reduced_torsional_system_t
  !! K'/M ve active map için authoritative temel olarak yeniden kullanılır.
  !! Girdiler: SI birimli torsional_system_t ve constraint_manager_t.
  !! Çıktı: Ortak active sıralı reduced_dynamic_torsional_system_t.
  !! Varsayımlar ve geçerlilik: Direct elimination, lineer küçük genlik
  !! ve frozen-property kabul edilir. Prescribed dinamik hareket ve RHS
  !! correction uygulanmaz; K'' viskoz c'ye dönüştürülmez.
  pure function build_reduced_dynamic_torsional_system(system, manager) &
      result(reduced_system)
    type(torsional_system_t), intent(in) :: system
    type(constraint_manager_t), intent(in) :: manager
    type(reduced_dynamic_torsional_system_t) :: reduced_system

    type(active_dof_map_t) :: active_mapping
    type(damping_matrix_t) :: full_damping
    type(dof_map_t) :: full_mapping
    type(loss_stiffness_matrix_t) :: full_loss_stiffness
    integer, allocatable :: retained_equation_ids(:)

    reduced_system%base_system = &
      build_reduced_torsional_system(system, manager)
    active_mapping = &
      get_reduced_active_dof_map(reduced_system%base_system)
    retained_equation_ids = &
      get_retained_full_equation_ids(active_mapping)

    call initialize_full_dof_map(full_mapping, system)
    full_loss_stiffness = &
      assemble_full_loss_stiffness(system, full_mapping)
    full_damping = assemble_full_damping(system, full_mapping)

    reduced_system%loss_stiffness = reduce_matrix( &
      full_loss_stiffness, retained_equation_ids)
    reduced_system%damping = reduce_matrix( &
      full_damping, retained_equation_ids)
    call validate_reduced_dynamic_torsional_system(reduced_system)
  end function build_reduced_dynamic_torsional_system

  !> Reduced dynamic sistemin ortak active-equation boyutunu doğrular.
  !! Fiziksel açıklama: K', K'', C ve M aynı torsional koordinat takımına
  !! ait olmalıdır; aksi durumda dinamik denge denklemi kurulamaz.
  !! Matematiksel açıklama: Dört kare matrisin mertebesi active DOF sayısı
  !! n ile aynı olmalıdır. Birimler sırasıyla [N*m/rad], [N*m/rad],
  !! [N*m*s/rad] ve [kg*m^2]'dir.
  !! Girdi reduced_dynamic_torsional_system_t; çıktı yoktur. Geçersiz
  !! depolama veya boyut error stop ile reddedilir. n=0 reduction geçerlidir.
  pure subroutine validate_reduced_dynamic_torsional_system(reduced_system)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system

    integer :: active_count

    call validate_reduced_torsional_system(reduced_system%base_system)
    call validate_loss_stiffness_matrix(reduced_system%loss_stiffness)
    call validate_damping_matrix(reduced_system%damping)

    active_count = &
      get_reduced_active_dof_count(reduced_system%base_system)
    if (get_loss_stiffness_matrix_size( &
        reduced_system%loss_stiffness) /= active_count .or. &
        get_damping_matrix_size(reduced_system%damping) /= active_count) then
      error stop "Reduced K'', C ve K'/M boyutları birbiriyle uyumsuz."
    end if
  end subroutine validate_reduced_dynamic_torsional_system

  !> Reduced depolama rijitliği K'_r [N*m/rad] bağımsız kopyasını verir.
  pure function get_reduced_dynamic_stiffness(reduced_system) &
      result(stiffness)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system
    type(stiffness_matrix_t) :: stiffness

    call validate_reduced_dynamic_torsional_system(reduced_system)
    stiffness = get_reduced_stiffness(reduced_system%base_system)
  end function get_reduced_dynamic_stiffness

  !> Reduced kayıp rijitliği K''_r [N*m/rad] bağımsız kopyasını verir.
  pure function get_reduced_dynamic_loss_stiffness(reduced_system) &
      result(loss_stiffness)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system
    type(loss_stiffness_matrix_t) :: loss_stiffness

    call validate_reduced_dynamic_torsional_system(reduced_system)
    loss_stiffness = reduced_system%loss_stiffness
  end function get_reduced_dynamic_loss_stiffness

  !> Reduced viskoz sönüm C_r [N*m*s/rad] bağımsız kopyasını verir.
  pure function get_reduced_dynamic_damping(reduced_system) result(damping)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system
    type(damping_matrix_t) :: damping

    call validate_reduced_dynamic_torsional_system(reduced_system)
    damping = reduced_system%damping
  end function get_reduced_dynamic_damping

  !> Reduced polar atalet M_r [kg*m^2] bağımsız kopyasını verir.
  pure function get_reduced_dynamic_mass(reduced_system) result(mass)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system
    type(mass_matrix_t) :: mass

    call validate_reduced_dynamic_torsional_system(reduced_system)
    mass = get_reduced_mass(reduced_system%base_system)
  end function get_reduced_dynamic_mass

  !> Physical/full/active torsional DOF eşlemesinin bağımsız kopyasını
  !! verir. Kimlikler boyutsuzdur ve dört reduced matrisin ortak sırasını
  !! tanımlar.
  pure function get_reduced_dynamic_active_dof_map(reduced_system) &
      result(mapping)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system
    type(active_dof_map_t) :: mapping

    call validate_reduced_dynamic_torsional_system(reduced_system)
    mapping = get_reduced_active_dof_map(reduced_system%base_system)
  end function get_reduced_dynamic_active_dof_map

  !> Reduced dynamic sistemdeki active torsional denklem sayısını [-] verir.
  pure function get_reduced_dynamic_active_dof_count(reduced_system) &
      result(active_count)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system
    integer :: active_count

    call validate_reduced_dynamic_torsional_system(reduced_system)
    active_count = &
      get_reduced_active_dof_count(reduced_system%base_system)
  end function get_reduced_dynamic_active_dof_count

  !> Reduced kompleks harmonik dönme cevabını fiziksel DOF uzayına açar.
  !!
  !! Fiziksel açıklama: Complex phasor theta_hat, peak angular displacement
  !! [rad] genliğidir. Kısıtlı DOF'lar V0.6 homojen harmonik perturbation
  !! modelinde sıfır genlik taşır.
  !! Matematiksel açıklama: theta_hat=P*theta_hat_r. Bu recovery
  !! q=P*q_r+q_p statik bağıntısındaki q_p prescribed offset'ini eklemez.
  !! Girdi: n_active uzunluklu, reel ve sanal bileşenleri sonlu kompleks
  !! reduced_response [rad]. Çıktı: physical_dof_id sıralı kompleks
  !! full_response [rad].
  !! Varsayımlar ve geçerlilik: Harmonic convention exp(+i*omega*t) ve peak
  !! amplitude'dur. Zamanla değişen prescribed angle bu kapsamda değildir.
  pure function recover_harmonic_response(reduced_system, reduced_response) &
      result(full_response)
    type(reduced_dynamic_torsional_system_t), intent(in) :: reduced_system
    complex(dp), intent(in) :: reduced_response(:)
    complex(dp), allocatable :: full_response(:)

    type(active_dof_map_entry_t) :: entry
    type(active_dof_map_entry_t), allocatable :: entries(:)
    integer :: entry_index

    call validate_reduced_dynamic_torsional_system(reduced_system)
    if (size(reduced_response) /= &
        get_reduced_dynamic_active_dof_count(reduced_system)) then
      error stop "Harmonik recovery vektörü active DOF sayısıyla uyumsuz."
    end if

    if (.not. all(ieee_is_finite(real(reduced_response, kind=dp))) .or. &
        .not. all(ieee_is_finite(aimag(reduced_response)))) then
      error stop "Harmonik recovery vektörü yalnız sonlu değerler içermelidir."
    end if

    entries = get_active_dof_map_entries( &
      get_reduced_dynamic_active_dof_map(reduced_system))
    allocate(full_response(size(entries)))
    full_response = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do entry_index = 1, size(entries)
      entry = entries(entry_index)
      if (.not. entry%constrained) then
        full_response(entry%physical_dof_id) = &
          reduced_response(entry%active_equation_id)
      end if
    end do
  end function recover_harmonic_response

end module tms_reduced_dynamic_system
