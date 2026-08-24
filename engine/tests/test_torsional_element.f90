program test_torsional_element
  use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
  use tms_kinds, only : dp
  use tms_torsional_element, only : torsional_element_t, &
    validate_torsional_element
  implicit none

  real(dp), parameter :: tolerance = 1.0e-12_dp

  type(torsional_element_t) :: element
  character(len=32) :: validation_case

  ! CTest geçersiz eleman alanlarını ayrı WILL_FAIL vakaları olarak sınar.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case))
    stop 0
  end if

  element = torsional_element_t( &
    id=3, &
    node_i_id=10, &
    node_j_id=20, &
    stiffness_nm_per_rad=5000.0_dp, &
    damping_nms_per_rad=12.5_dp)
  call validate_torsional_element(element)

  if (element%id /= 3) then
    error stop "Torsional eleman kimliği doğru saklanmadı."
  end if

  if (element%node_i_id /= 10 .or. element%node_j_id /= 20) then
    error stop "Torsional eleman düğüm bağlantısı doğru saklanmadı."
  end if

  if (abs(element%stiffness_nm_per_rad - 5000.0_dp) > tolerance) then
    error stop "Torsional eleman rijitliği doğru saklanmadı."
  end if

  if (abs(element%damping_nms_per_rad - 12.5_dp) > tolerance) then
    error stop "Torsional eleman viskoz sönümü doğru saklanmadı."
  end if

  print *, "Genel torsional eleman veri yapısı doğrulandı."

contains

  !> Geçersiz eleman kimliği, bağlantısı, rijitliği veya sönümünün üretim
  !! doğrulayıcısı tarafından reddedildiğini sınar.
  !! Girdiler: Kimlikler boyutsuz, K [N*m/rad], c [N*m*s/rad]. Çıktı üretmez;
  !! seçilen geçersiz durum error stop oluşturmalıdır.
  subroutine exercise_invalid_case(case_name)
    character(len=*), intent(in) :: case_name

    type(torsional_element_t) :: invalid_element

    invalid_element = torsional_element_t( &
      id=1, node_i_id=1, node_j_id=2, &
      stiffness_nm_per_rad=1000.0_dp, damping_nms_per_rad=0.0_dp)

    select case (case_name)
      case ("nonpositive_id")
        invalid_element%id = 0
      case ("nonpositive_node_id")
        invalid_element%node_i_id = 0
      case ("self_connection")
        invalid_element%node_j_id = invalid_element%node_i_id
      case ("nonpositive_stiffness")
        invalid_element%stiffness_nm_per_rad = 0.0_dp
      case ("negative_damping")
        invalid_element%damping_nms_per_rad = -1.0_dp
      case ("nonfinite_stiffness")
        invalid_element%stiffness_nm_per_rad = &
          ieee_value(0.0_dp, ieee_positive_inf)
      case default
        error stop "Bilinmeyen geçersiz torsional eleman testi istendi."
    end select

    call validate_torsional_element(invalid_element)
    print *, "Geçersiz eleman beklenmedik biçimde kabul edildi: ", case_name
  end subroutine exercise_invalid_case

end program test_torsional_element
