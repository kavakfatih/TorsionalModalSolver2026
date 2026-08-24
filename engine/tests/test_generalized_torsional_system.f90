program test_generalized_torsional_system
  use tms_kinds, only : dp
  use tms_torsional_node, only : torsional_node_t
  use tms_torsional_element, only : torsional_element_t
  use tms_generalized_torsional_system, only : torsional_system_t, &
    add_torsional_node, add_torsional_element, count_torsional_dofs, &
    get_torsional_node_count, get_torsional_element_count, &
    get_torsional_node, get_torsional_element, validate_torsional_system
  use tms_torsional_system, only : two_inertia_tvd_system_t, &
    build_generalized_two_inertia_system
  implicit none

  real(dp), parameter :: tolerance = 1.0e-12_dp

  type(torsional_system_t) :: system
  type(torsional_system_t) :: free_free_system
  type(torsional_system_t) :: fixed_hub_system
  type(torsional_node_t) :: node
  type(torsional_element_t) :: element
  type(two_inertia_tvd_system_t) :: two_inertia_system
  character(len=32) :: validation_case

  ! CTest koleksiyon ve bağlantı bütünlüğü hatalarını ayrı WILL_FAIL vakaları
  ! olarak doğrudan üretim yönetim yordamları üzerinden sınar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  call add_torsional_node(system, torsional_node_t( &
    id=10, polar_inertia_kg_m2=0.15_dp, initial_angle_rad=0.01_dp, &
    constrained=.false.))
  call add_torsional_node(system, torsional_node_t( &
    id=20, polar_inertia_kg_m2=0.25_dp, initial_angle_rad=0.0_dp, &
    constrained=.true.))
  call add_torsional_element(system, torsional_element_t( &
    id=30, node_i_id=10, node_j_id=20, &
    stiffness_nm_per_rad=2500.0_dp, damping_nms_per_rad=5.0_dp))
  call validate_torsional_system(system)

  if (get_torsional_node_count(system) /= 2) then
    error stop "Genel sistem düğüm sayısı doğru değil."
  end if

  if (get_torsional_element_count(system) /= 1) then
    error stop "Genel sistem eleman sayısı doğru değil."
  end if

  if (count_torsional_dofs(system) /= 1) then
    error stop "Sabitlenmemiş düğümlerden aktif DOF sayısı doğru bulunamadı."
  end if

  node = get_torsional_node(system, 1)
  element = get_torsional_element(system, 1)
  if (node%id /= 10 .or. abs(node%polar_inertia_kg_m2 - 0.15_dp) > &
      tolerance) then
    error stop "Genel sistem düğüm koleksiyonu doğru okunamadı."
  end if
  if (element%node_i_id /= 10 .or. element%node_j_id /= 20 .or. &
      abs(element%stiffness_nm_per_rad - 2500.0_dp) > tolerance) then
    error stop "Genel sistem eleman koleksiyonu doğru okunamadı."
  end if

  ! Benchmark 004 iki-ataletli referansı yeni genel topolojiye yalnız mevcut
  ! üretim dönüşümüyle aktarılır; test herhangi bir M/K veya eigen hesabını
  ! yeniden kurmaz.
  two_inertia_system = two_inertia_tvd_system_t( &
    hub_polar_inertia_kg_m2=0.10_dp, &
    ring_polar_inertia_kg_m2=0.20_dp, &
    storage_stiffness_nm_per_rad=1000.0_dp, &
    loss_stiffness_nm_per_rad=100.0_dp, &
    loss_factor=0.1_dp, &
    material_reference_frequency_hz=100.0_dp, &
    material_temperature_k=293.15_dp)

  free_free_system = build_generalized_two_inertia_system(two_inertia_system)
  call validate_torsional_system(free_free_system)

  if (get_torsional_node_count(free_free_system) /= 2 .or. &
      get_torsional_element_count(free_free_system) /= 1) then
    error stop "İki ataletli genel topoloji iki düğüm ve bir eleman içermiyor."
  end if
  if (count_torsional_dofs(free_free_system) /= 2) then
    error stop "Serbest-serbest genel topolojinin DOF sayısı iki değil."
  end if

  node = get_torsional_node(free_free_system, 1)
  if (node%id /= 1 .or. node%constrained .or. &
      abs(node%polar_inertia_kg_m2 - 0.10_dp) > tolerance) then
    error stop "Göbek ataleti genel serbest-serbest düğüme aktarılamadı."
  end if
  node = get_torsional_node(free_free_system, 2)
  if (node%id /= 2 .or. node%constrained .or. &
      abs(node%polar_inertia_kg_m2 - 0.20_dp) > tolerance) then
    error stop "Halka ataleti genel serbest-serbest düğüme aktarılamadı."
  end if

  element = get_torsional_element(free_free_system, 1)
  if (element%node_i_id /= 1 .or. element%node_j_id /= 2) then
    error stop "Elastomer genel elemanının düğüm bağlantısı doğru değil."
  end if
  if (abs(element%stiffness_nm_per_rad - 1000.0_dp) > tolerance) then
    error stop "K' genel torsional eleman rijitliğine aktarılamadı."
  end if
  if (abs(element%damping_nms_per_rad) > tolerance) then
    error stop "K'' yanlış biçimde viskoz sönüm alanına aktarıldı."
  end if

  fixed_hub_system = build_generalized_two_inertia_system( &
    two_inertia_system, hub_constrained=.true.)
  call validate_torsional_system(fixed_hub_system)
  if (count_torsional_dofs(fixed_hub_system) /= 1) then
    error stop "Fixed-hub genel topolojinin DOF sayısı bir değil."
  end if
  node = get_torsional_node(fixed_hub_system, 1)
  if (.not. node%constrained) then
    error stop "Fixed-hub dönüşümünde göbek düğümü sabitlenmedi."
  end if

  print *, "Genel torsional sistem ve iki-atalet uyumluluğu doğrulandı."

contains

  !> Sistem yönetim yordamlarının boş sistem, yinelenen kimlik ve tanımsız
  !! bağlantı uçlarını reddettiğini sınar.
  !! Girdiler: J [kg*m^2], K [N*m/rad], c [N*m*s/rad]. Çıktı üretmez;
  !! seçilen geçersiz topoloji error stop oluşturmalıdır.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_system_t) :: invalid_system
    type(torsional_node_t) :: valid_node_i
    type(torsional_node_t) :: valid_node_j
    type(torsional_element_t) :: valid_element

    valid_node_i = torsional_node_t( &
      id=1, polar_inertia_kg_m2=0.1_dp, constrained=.false.)
    valid_node_j = torsional_node_t( &
      id=2, polar_inertia_kg_m2=0.2_dp, constrained=.false.)
    valid_element = torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=1000.0_dp, damping_nms_per_rad=0.0_dp)

    select case (case_name)
      case ("empty_system")
        call validate_torsional_system(invalid_system)
      case ("duplicate_node_id")
        call add_torsional_node(invalid_system, valid_node_i)
        valid_node_j%id = valid_node_i%id
        call add_torsional_node(invalid_system, valid_node_j)
      case ("missing_endpoint")
        call add_torsional_node(invalid_system, valid_node_i)
        call add_torsional_element(invalid_system, valid_element)
      case ("duplicate_element_id")
        call add_torsional_node(invalid_system, valid_node_i)
        call add_torsional_node(invalid_system, valid_node_j)
        call add_torsional_element(invalid_system, valid_element)
        valid_element%node_i_id = 2
        valid_element%node_j_id = 1
        call add_torsional_element(invalid_system, valid_element)
      case default
        error stop "Bilinmeyen geçersiz genel torsional sistem testi istendi."
    end select

    print *, "Geçersiz genel sistem beklenmedik biçimde kabul edildi: ", &
      case_name
  end subroutine exercise_invalid_case

end program test_generalized_torsional_system
