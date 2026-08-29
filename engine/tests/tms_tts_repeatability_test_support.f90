module tms_tts_repeatability_test_support
  use tms_kinds, only : dp
  use tms_constants, only : universal_gas_constant_j_per_mol_k
  use tms_tts_types, only : tts_material_family_t, &
    tts_pair_shift_result_t, TTS_IDENTIFICATION_SUCCESS
  use tms_tts_identification, only : identify_tts_master_curve
  use tms_tts_repeatability_types, only : tts_repeatability_campaign_t
  use tms_tts_test_support, only : make_exact_trs_family
  implicit none
  private

  public :: make_repeatability_campaign
  public :: make_arrhenius_shifts
  public :: make_wlf_shifts
  public :: make_linear_temperature_shifts
  public :: reverse_pair_order_and_orientation
  public :: replace_campaign_shift_evidence

contains

  !> Exact horizontal-TRS synthetic family'yi gerçek V0.8.1 production
  !! identification API'sinden geçirerek complete repeatability campaign
  !! fixture üretir. Sıcaklık [K], shift boyutsuzdur. Synthetic değerler
  !! yalnız numerical test truth'udur; material kabul limiti değildir.
  function make_repeatability_campaign( &
      temperatures_k, physical_shifts, reference_temperature_k, &
      campaign_identifier, replicate_basis) result(campaign)
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: physical_shifts(:)
    real(dp), intent(in) :: reference_temperature_k
    character(len=*), intent(in) :: campaign_identifier
    integer, intent(in) :: replicate_basis
    type(tts_repeatability_campaign_t) :: campaign

    type(tts_material_family_t) :: family
    character(len=256) :: reference_identifier
    character(len=32) :: index_text
    integer :: i

    if (size(temperatures_k) /= size(physical_shifts)) then
      error stop "Repeatability fixture temperature/shift boyutları eşit olmalı."
    end if
    family = make_exact_trs_family(temperatures_k, physical_shifts)
    family%family_identifier = "REPEATABILITY-"//trim(campaign_identifier)
    family%common_state%source_metadata = &
      "Synthetic campaign source: "//trim(campaign_identifier)
    reference_identifier = ""
    do i = 1, size(family%isotherms)
      write(index_text, '(I0)') i
      family%isotherms(i)%specimen_identifier = &
        trim(campaign_identifier)//"-SPECIMEN-"//trim(index_text)
      family%isotherms(i)%source_identifier = &
        trim(campaign_identifier)//"-SOURCE-"//trim(index_text)
      if (abs(temperatures_k(i) - reference_temperature_k) <= &
          64.0_dp*epsilon(1.0_dp)*max(abs(temperatures_k(i)), &
            abs(reference_temperature_k), 1.0_dp)) then
        reference_identifier = family%isotherms(i)%isotherm_identifier
      end if
    end do
    if (len_trim(reference_identifier) == 0) then
      error stop "Synthetic campaign measured reference içermiyor."
    end if

    campaign%campaign_identifier = campaign_identifier
    campaign%replicate_basis = replicate_basis
    campaign%laboratory_identifier = "TMS26-SYNTHETIC-LAB"
    campaign%operator_identifier = "TMS26-SYNTHETIC-OPERATOR"
    campaign%instrument_identifier = "TMS26-SYNTHETIC-DMA"
    campaign%test_protocol_identifier = "TMS26-SYNTHETIC-PROTOCOL"
    campaign%calibration_reference = "SYNTHETIC-CALIBRATION-TRACE"
    campaign%run_identifier = "RUN-"//trim(campaign_identifier)
    campaign%test_date_metadata = "2026-08-29"
    campaign%identification = identify_tts_master_curve( &
      family, trim(reference_identifier))
    if (campaign%identification%status /= TTS_IDENTIFICATION_SUCCESS) then
      error stop "Synthetic complete V0.8.1 campaign üretilemedi."
    end if
  end function make_repeatability_campaign

  !> Arrhenius synthetic truth üretir:
  !! s=Ea_app/(R ln(10))*(1/T-1/T_ref). T/T_ref [K], Ea_app [J/mol],
  !! çıktı s=log10(a_T) boyutsuzdur. Chemical-aging modeli değildir.
  pure function make_arrhenius_shifts( &
      temperatures_k, reference_temperature_k, apparent_energy_j_per_mol) &
      result(shifts)
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: apparent_energy_j_per_mol
    real(dp), allocatable :: shifts(:)
    real(dp) :: beta_k

    beta_k = apparent_energy_j_per_mol / &
      (universal_gas_constant_j_per_mol_k*log(10.0_dp))
    shifts = beta_k*(1.0_dp/temperatures_k - &
      1.0_dp/reference_temperature_k)
  end function make_arrhenius_shifts

  !> WLF synthetic truth üretir: s=-C1*dT/(C2+dT). T/C2 [K], C1
  !! boyutsuz ve çıktı s=log10(a_T) boyutsuzdur. Pole-safe fixture seçimi
  !! caller sorumluluğundadır.
  pure function make_wlf_shifts( &
      temperatures_k, reference_temperature_k, c1, c2_k) result(shifts)
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: c1
    real(dp), intent(in) :: c2_k
    real(dp), allocatable :: shifts(:)
    real(dp), allocatable :: delta_temperature_k(:)

    delta_temperature_k = temperatures_k - reference_temperature_k
    shifts = -c1*delta_temperature_k/(c2_k + delta_temperature_k)
  end function make_wlf_shifts

  !> Large-C2 WLF limitini temsil eden linear-temperature synthetic truth
  !! üretir: s=-p*(T-T_ref). p [1/K], T [K], s boyutsuzdur.
  pure function make_linear_temperature_shifts( &
      temperatures_k, reference_temperature_k, slope_per_k) result(shifts)
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: reference_temperature_k
    real(dp), intent(in) :: slope_per_k
    real(dp), allocatable :: shifts(:)

    shifts = -slope_per_k*(temperatures_k - reference_temperature_k)
  end function make_linear_temperature_shifts

  !> Pair storage sırasını ters çevirir ve her pair orientation'ını swap eder.
  !! delta_s_new=-delta_s_old ile aynı physical key/value korunur. Bu helper,
  !! production matching'in array index ve orientation'a bağımlı olmadığını
  !! sınamak için kullanılır.
  subroutine reverse_pair_order_and_orientation(campaign)
    type(tts_repeatability_campaign_t), intent(inout) :: campaign
    type(tts_pair_shift_result_t), allocatable :: reversed(:)
    integer :: i
    integer :: old_reference_index
    integer :: source_index

    allocate(reversed(size(campaign%identification%pair_shift_results)))
    do i = 1, size(reversed)
      source_index = size(reversed) - i + 1
      reversed(i) = campaign%identification%pair_shift_results(source_index)
      old_reference_index = reversed(i)%reference_isotherm_index
      reversed(i)%reference_isotherm_index = reversed(i)%moving_isotherm_index
      reversed(i)%moving_isotherm_index = old_reference_index
      reversed(i)%delta_s = -reversed(i)%delta_s
    end do
    campaign%identification%pair_shift_results = reversed
  end subroutine reverse_pair_order_and_orientation

  !> Targeted failure fixture için complete campaign'in empirical s(T) ve
  !! adjacent pair delta_s evidence'ını aynı finite physical key tablosuyla
  !! değiştirir. Production analysis çağrılmaz; input sıcaklık [K], absolute
  !! shift boyutsuzdur. Master data yalnız test scaffold'udur.
  subroutine replace_campaign_shift_evidence( &
      campaign, temperatures_k, absolute_shifts)
    type(tts_repeatability_campaign_t), intent(inout) :: campaign
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: absolute_shifts(:)

    real(dp) :: moving_shift
    real(dp) :: reference_anchor_shift
    real(dp) :: reference_shift
    integer :: i
    integer :: moving_index
    integer :: reference_index

    if (size(temperatures_k) /= size(absolute_shifts)) then
      error stop "Replacement temperature/shift boyutları eşit olmalıdır."
    end if
    reference_anchor_shift = lookup_shift(temperatures_k, absolute_shifts, &
      campaign%identification%reference_temperature_k)
    do i = 1, size(campaign%identification%empirical_shifts)
      reference_shift = lookup_shift(temperatures_k, absolute_shifts, &
        campaign%identification%empirical_shifts(i)%temperature_k)
      campaign%identification%empirical_shifts(i)%log10_a_t = &
        reference_shift - reference_anchor_shift
      campaign%identification%empirical_shifts(i)%a_t = 10.0_dp**( &
        reference_shift - reference_anchor_shift)
    end do
    do i = 1, size(campaign%identification%pair_shift_results)
      reference_index = campaign%identification%pair_shift_results(i) &
        %reference_isotherm_index
      moving_index = campaign%identification%pair_shift_results(i) &
        %moving_isotherm_index
      reference_shift = lookup_shift(temperatures_k, absolute_shifts, &
        campaign%identification%source_family%isotherms( &
          reference_index)%temperature_k)
      moving_shift = lookup_shift(temperatures_k, absolute_shifts, &
        campaign%identification%source_family%isotherms( &
          moving_index)%temperature_k)
      campaign%identification%pair_shift_results(i)%delta_s = &
        moving_shift - reference_shift
    end do
  end subroutine replace_campaign_shift_evidence

  pure function lookup_shift(temperatures_k, shifts, temperature_k) &
      result(value)
    real(dp), intent(in) :: temperatures_k(:)
    real(dp), intent(in) :: shifts(:)
    real(dp), intent(in) :: temperature_k
    real(dp) :: value
    integer :: i

    value = huge(1.0_dp)
    do i = 1, size(temperatures_k)
      if (abs(temperatures_k(i) - temperature_k) <= &
          64.0_dp*epsilon(1.0_dp)*max(abs(temperatures_k(i)), &
            abs(temperature_k), 1.0_dp)) then
        value = shifts(i)
        return
      end if
    end do
    error stop "Synthetic shift lookup measured temperature bulamadı."
  end function lookup_shift

end module tms_tts_repeatability_test_support
