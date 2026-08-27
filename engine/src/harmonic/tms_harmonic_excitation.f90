module tms_harmonic_excitation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_dof_types, only : TORSIONAL_ROTATION, is_supported_dof_type
  use tms_constraint_manager, only : active_dof_map_t, &
    validate_active_dof_map, get_active_equation_count, &
    lookup_active_equation_id
  implicit none
  private

  !> Tek bir fiziksel torsional DOF'a uygulanan complex peak torque
  !! excitation kaydını taşır.
  type, public :: harmonic_excitation_t
    !> Uyarımın uygulandığı pozitif fiziksel düğüm kimliği [-].
    integer :: node_id = 0

    !> Uyarılan fiziksel DOF türü [-]. V0.6 yalnız torsional rotation destekler.
    integer :: dof_type = TORSIONAL_ROTATION

    !> exp(+i*omega*t) convention'ındaki complex PEAK torque genliği [N*m].
    !! Gerçek ve sanal bileşenler sırasıyla in-phase ve quadrature kanallardır;
    !! bu değer RMS değildir.
    complex(dp) :: torque_amplitude_nm = cmplx(0.0_dp, 0.0_dp, kind=dp)
  end type harmonic_excitation_t

  public :: validate_harmonic_excitation
  public :: assemble_harmonic_load_vector

contains

  !> Harmonic nodal torque kaydının fiziksel ve sayısal geçerliliğini sınar.
  !! Girdide node_id pozitif, DOF türü desteklenen torsional rotation ve torque
  !! genliğinin [N*m] gerçek/sanal bileşenleri sonlu olmalıdır. Çıktı üretmez;
  !! geçersiz veri error stop ile reddedilir. Faz derece veya RMS dönüşümü bu
  !! core veri sözleşmesinin kapsamında değildir.
  pure subroutine validate_harmonic_excitation(excitation)
    type(harmonic_excitation_t), intent(in) :: excitation

    if (excitation%node_id <= 0) then
      error stop "Harmonic excitation düğüm kimliği pozitif olmalıdır."
    end if
    if (.not. is_supported_dof_type(excitation%dof_type)) then
      error stop "Harmonic excitation DOF türü desteklenmiyor."
    end if
    if (.not. ieee_is_finite(real(excitation%torque_amplitude_nm, dp)) .or. &
        .not. ieee_is_finite(aimag(excitation%torque_amplitude_nm))) then
      error stop "Harmonic torque genliğinin complex bileşenleri sonlu olmalıdır."
    end if
  end subroutine validate_harmonic_excitation

  !> Fiziksel nodal torque katkılarını reduced active load vektöründe toplar.
  !!
  !! Fiziksel açıklama: Aynı DOF'a uygulanan birden çok complex peak torque
  !! lineer süperpozisyonla toplanır.
  !! Matematiksel açıklama: T_r(eq) <- T_r(eq)+T_hat_j scatter-add işlemi,
  !! active equation eşlemesi üzerinden uygulanır.
  !! Girdiler: Active DOF haritası ve [N*m] birimli, boş olmayan excitation
  !! dizisi. Çıktı: n_active uzunluklu complex(dp) T_hat_r [N*m] vektörü.
  !! Varsayımlar ve geçerlilik: Tamamen constrained hedef reaction torque
  !! niteliğindedir ve sessizce atılmaz; clean diagnostic ile reddedilir.
  !! Bilinmeyen fiziksel DOF, desteklenmeyen tür veya toplama taşması reddedilir.
  pure function assemble_harmonic_load_vector(mapping, excitations) &
      result(load_vector)
    type(active_dof_map_t), intent(in) :: mapping
    type(harmonic_excitation_t), intent(in) :: excitations(:)
    complex(dp), allocatable :: load_vector(:)

    complex(dp) :: updated_torque
    integer :: active_equation_id
    integer :: excitation_index

    call validate_active_dof_map(mapping)
    if (size(excitations) == 0) then
      error stop "Harmonic analysis en az bir excitation gerektirir."
    end if

    allocate(load_vector(get_active_equation_count(mapping)))
    load_vector = cmplx(0.0_dp, 0.0_dp, kind=dp)

    do excitation_index = 1, size(excitations)
      call validate_harmonic_excitation(excitations(excitation_index))
      active_equation_id = lookup_active_equation_id( &
        mapping, excitations(excitation_index)%node_id, &
        excitations(excitation_index)%dof_type)

      if (active_equation_id == 0) then
        error stop &
          "Harmonic torque tamamen constrained bir DOF'a uygulanamaz."
      end if

      updated_torque = load_vector(active_equation_id) + &
        excitations(excitation_index)%torque_amplitude_nm
      if (.not. ieee_is_finite(real(updated_torque, dp)) .or. &
          .not. ieee_is_finite(aimag(updated_torque))) then
        error stop "Harmonic torque scatter-add sonucu sonlu olmalıdır."
      end if
      load_vector(active_equation_id) = updated_torque
    end do
  end function assemble_harmonic_load_vector

end module tms_harmonic_excitation
