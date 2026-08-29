module tms_tts_diagnostics
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_vgp_point_t, &
    tts_cole_cole_point_t, is_runtime_export_usable
  implicit none
  private

  public :: build_tts_vgp_cloud
  public :: build_tts_cole_cole_cloud

contains

  !> Her VALID complex measured point için Van Gurp-Palmen cloud üretir.
  !! |G*|=hypot(G',G'') [Pa] ve delta=atan2(G'',G') [rad] kullanılır.
  !! Frequency shifting gerekmez; temperature ve source point provenance
  !! korunur. Universal TRS score veya PASS threshold üretilmez.
  pure function build_tts_vgp_cloud(family) result(points)
    type(tts_material_family_t), intent(in) :: family
    type(tts_vgp_point_t), allocatable :: points(:)

    integer :: count_value
    integer :: i
    integer :: j

    count_value = count_runtime_valid_points(family)
    allocate(points(count_value))
    count_value = 0
    do i = 1, size(family%isotherms)
      do j = 1, size(family%isotherms(i)%points)
        if (.not. is_runtime_export_usable( &
            family%isotherms(i)%points(j))) cycle
        count_value = count_value + 1
        points(count_value)%source_isotherm_index = i
        points(count_value)%source_point_index = j
        points(count_value)%source_temperature_k = &
          family%isotherms(i)%temperature_k
        points(count_value)%complex_modulus_magnitude_pa = hypot( &
          family%isotherms(i)%points(j)%storage_modulus_pa, &
          family%isotherms(i)%points(j)%loss_modulus_pa)
        points(count_value)%phase_angle_rad = atan2( &
          family%isotherms(i)%points(j)%loss_modulus_pa, &
          family%isotherms(i)%points(j)%storage_modulus_pa)
      end do
    end do
  end function build_tts_vgp_cloud

  !> Her VALID complex measured point için linear-axis Cole-Cole cloud üretir:
  !! ordinate G'' [Pa], abscissa G' [Pa]. Temperature/source provenance
  !! korunur; smoothing veya tek scalar acceptance metriği uygulanmaz.
  pure function build_tts_cole_cole_cloud(family) result(points)
    type(tts_material_family_t), intent(in) :: family
    type(tts_cole_cole_point_t), allocatable :: points(:)

    integer :: count_value
    integer :: i
    integer :: j

    count_value = count_runtime_valid_points(family)
    allocate(points(count_value))
    count_value = 0
    do i = 1, size(family%isotherms)
      do j = 1, size(family%isotherms(i)%points)
        if (.not. is_runtime_export_usable( &
            family%isotherms(i)%points(j))) cycle
        count_value = count_value + 1
        points(count_value)%source_isotherm_index = i
        points(count_value)%source_point_index = j
        points(count_value)%source_temperature_k = &
          family%isotherms(i)%temperature_k
        points(count_value)%storage_modulus_pa = &
          family%isotherms(i)%points(j)%storage_modulus_pa
        points(count_value)%loss_modulus_pa = &
          family%isotherms(i)%points(j)%loss_modulus_pa
      end do
    end do
  end function build_tts_cole_cole_cloud

  pure function count_runtime_valid_points(family) result(count_value)
    type(tts_material_family_t), intent(in) :: family
    integer :: count_value
    integer :: i
    integer :: j

    count_value = 0
    do i = 1, size(family%isotherms)
      do j = 1, size(family%isotherms(i)%points)
        if (is_runtime_export_usable( &
            family%isotherms(i)%points(j))) count_value = count_value + 1
      end do
    end do
  end function count_runtime_valid_points

end module tms_tts_diagnostics
