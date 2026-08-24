program test_torsional_node
  use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
  use tms_kinds, only : dp
  use tms_torsional_node, only : torsional_node_t, validate_torsional_node
  implicit none

  real(dp), parameter :: tolerance = 1.0e-12_dp

  type(torsional_node_t) :: node
  character(len=32) :: validation_case

  ! CTest geçersiz düğüm alanlarını ayrı WILL_FAIL vakaları olarak sınar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  node = torsional_node_t( &
    id=7, &
    polar_inertia_kg_m2=0.125_dp, &
    initial_angle_rad=0.02_dp, &
    constrained=.true.)
  call validate_torsional_node(node)

  if (node%id /= 7) then
    error stop "Torsional düğüm kimliği doğru saklanmadı."
  end if

  if (abs(node%polar_inertia_kg_m2 - 0.125_dp) > tolerance) then
    error stop "Torsional düğüm polar ataleti doğru saklanmadı."
  end if

  if (abs(node%initial_angle_rad - 0.02_dp) > tolerance) then
    error stop "Torsional düğüm başlangıç açısı doğru saklanmadı."
  end if

  if (.not. node%constrained) then
    error stop "Torsional düğüm sabitlenmişlik bilgisi doğru saklanmadı."
  end if

  print *, "Genel torsional düğüm veri yapısı doğrulandı."

contains

  !> Geçersiz düğüm kimliği, ataleti veya başlangıç açısının üretim
  !! doğrulayıcısı tarafından reddedildiğini sınar.
  !! Girdiler: Kimlik boyutsuz, J [kg*m^2], theta_0 [rad]. Çıktı üretmez;
  !! seçilen geçersiz durum error stop oluşturmalıdır.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_node_t) :: invalid_node

    invalid_node = torsional_node_t( &
      id=1, polar_inertia_kg_m2=0.1_dp, initial_angle_rad=0.0_dp)

    select case (case_name)
      case ("nonpositive_id")
        invalid_node%id = 0
      case ("nonpositive_inertia")
        invalid_node%polar_inertia_kg_m2 = 0.0_dp
      case ("nonfinite_angle")
        invalid_node%initial_angle_rad = &
          ieee_value(0.0_dp, ieee_quiet_nan)
      case default
        error stop "Bilinmeyen geçersiz torsional düğüm testi istendi."
    end select

    call validate_torsional_node(invalid_node)
    print *, "Geçersiz düğüm beklenmedik biçimde kabul edildi: ", case_name
  end subroutine exercise_invalid_case

end program test_torsional_node
