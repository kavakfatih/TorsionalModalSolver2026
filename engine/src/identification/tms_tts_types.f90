module tms_tts_types
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: tts_text_length = 256

  integer, parameter, public :: MEASUREMENT_VALID = 1
  integer, parameter, public :: BELOW_RELIABLE_FLOOR = 2
  integer, parameter, public :: MEASUREMENT_UNAVAILABLE = 3
  integer, parameter, public :: MEASUREMENT_REJECTED = 4

  integer, parameter, public :: TTS_DEFORMATION_MODE_SHEAR = 1

  integer, parameter, public :: TTS_CHANNEL_STORAGE = 1
  integer, parameter, public :: TTS_CHANNEL_LOSS = 2
  integer, parameter, public :: TTS_CHANNEL_JOINT = 3
  integer, parameter, public :: SHIFT_FROM_JOINT = 1
  integer, parameter, public :: SHIFT_FROM_STORAGE_ONLY = 2

  integer, parameter, public :: PAIR_SHIFT_SUCCESS = 0
  integer, parameter, public :: PAIR_SHIFT_STORAGE_ONLY = 1
  integer, parameter, public :: PAIR_SHIFT_INVALID_INPUT = 2
  integer, parameter, public :: PAIR_SHIFT_NO_OVERLAP = 3
  integer, parameter, public :: PAIR_SHIFT_INSUFFICIENT_SUPPORT = 4
  integer, parameter, public :: PAIR_SHIFT_NO_INTERIOR_MINIMUM = 5
  integer, parameter, public :: PAIR_SHIFT_OPTIMIZATION_FAILED = 6

  integer, parameter, public :: TTS_IDENTIFICATION_SUCCESS = 0
  integer, parameter, public :: TTS_IDENTIFICATION_INVALID_INPUT = 1
  integer, parameter, public :: TTS_IDENTIFICATION_INCONSISTENT_STATE = 2
  integer, parameter, public :: TTS_IDENTIFICATION_REFERENCE_NOT_FOUND = 3
  integer, parameter, public :: TTS_IDENTIFICATION_INSUFFICIENT_ISOTHERMS = 4
  integer, parameter, public :: TTS_IDENTIFICATION_INSUFFICIENT_SUPPORT = 5
  integer, parameter, public :: TTS_IDENTIFICATION_CHAIN_BROKEN = 6
  integer, parameter, public :: TTS_IDENTIFICATION_OPTIMIZATION_FAILED = 7
  integer, parameter, public :: TTS_IDENTIFICATION_MASTER_CURVE_FAILED = 8
  integer, parameter, public :: TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED = 9

  integer, parameter, public :: DEFAULT_TTS_COARSE_SCAN_POINT_COUNT = 65

  !> Tek DMA/dynamic-shear ölçüm noktasını ve channel-specific quality
  !! semantiğini taşır. Frekans [Hz], G' ve G'' [Pa], quality değerleri
  !! explicit enum sabitleridir. G''=0, VALID olduğunda fiziksel olarak geçerli
  !! fakat log-loss objective açısından kullanılamayan passive bir noktadır.
  type, public :: tts_measurement_point_t
    real(dp) :: frequency_hz = 0.0_dp
    real(dp) :: storage_modulus_pa = 0.0_dp
    real(dp) :: loss_modulus_pa = 0.0_dp
    integer :: storage_quality = MEASUREMENT_UNAVAILABLE
    integer :: loss_quality = MEASUREMENT_UNAVAILABLE
  end type tts_measurement_point_t

  !> Aynı sıcaklıktaki ölçülmüş isotherm'i taşır. Nokta dizisi caller'dan
  !! deep-copy edilir; sıralama authoritative input sırasıdır ve sessiz sort
  !! uygulanmaz. Sıcaklık mutlak Kelvin'dir.
  type, public :: tts_isotherm_t
    character(len=tts_text_length) :: isotherm_identifier = ""
    real(dp) :: temperature_k = 0.0_dp
    character(len=tts_text_length) :: specimen_identifier = ""
    character(len=tts_text_length) :: source_identifier = ""
    type(tts_measurement_point_t), allocatable :: points(:)
  end type tts_isotherm_t

  !> Bir TTS family içindeki bütün isotherm'lerin ortak ve authoritative test
  !! durumudur. Boyutsuz dynamic strain amplitude ve static prestrain her
  !! isotherm'de tekrar edilmez; floating-point equality ile family kararı
  !! verilmez. Batch/material-state tek kimlikle sabitlenir.
  type, public :: tts_common_test_state_t
    character(len=tts_text_length) :: material_identifier = ""
    character(len=tts_text_length) :: batch_state_identifier = ""
    real(dp) :: dynamic_strain_amplitude_ratio = 0.0_dp
    real(dp) :: static_prestrain_ratio = 0.0_dp
    integer :: deformation_mode = 0
    character(len=tts_text_length) :: conditioning_description = ""
    character(len=tts_text_length) :: test_method = ""
    character(len=tts_text_length) :: source_metadata = ""
  end type tts_common_test_state_t

  !> Ortak deney durumu altında measured isotherm koleksiyonudur. V0.8.1
  !! reference implementation tek batch/material-state family kullanır.
  type, public :: tts_material_family_t
    character(len=tts_text_length) :: family_identifier = ""
    type(tts_common_test_state_t) :: common_state
    type(tts_isotherm_t), allocatable :: isotherms(:)
  end type tts_material_family_t

  !> Programı durdurmadan dönen experimental validation sonucudur.
  type, public :: tts_validation_result_t
    integer :: status = TTS_IDENTIFICATION_INVALID_INPUT
    logical :: valid = .false.
    character(len=tts_text_length) :: message = ""
  end type tts_validation_result_t

  !> Tek channel'ın contiguous VALID ve positive-log-usable parçasıdır.
  !! x=log10(f/Hz), y=log10(G/Pa); source indeksleri provenance içindir.
  type, public :: tts_log_segment_t
    integer :: channel = 0
    real(dp), allocatable :: x(:)
    real(dp), allocatable :: y(:)
    integer, allocatable :: source_point_indices(:)
  end type tts_log_segment_t

  !> Belirli relative shift'teki exact piecewise-linear objective sonucudur.
  !! Objective log10-modulus residual karesidir; overlap boyutları decade'dir.
  type, public :: tts_pair_objective_evaluation_t
    logical :: valid = .false.
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: overlap_width_decades = 0.0_dp
    real(dp) :: overlap_fraction = 0.0_dp
    integer :: interpolation_interval_count = 0
  end type tts_pair_objective_evaluation_t

  !> Adjacent-pair shift aramasının numerical konfigürasyonudur. Tolerance
  !! s=log10(a_T) boyutsuz koordinatında numerical stopping ölçütüdür.
  type, public :: tts_pair_shift_configuration_t
    integer :: coarse_scan_point_count = DEFAULT_TTS_COARSE_SCAN_POINT_COUNT
    real(dp) :: absolute_tolerance = 8.0_dp*sqrt(epsilon(1.0_dp))
    real(dp) :: relative_tolerance = 8.0_dp*sqrt(epsilon(1.0_dp))
    integer :: maximum_iterations = 200
  end type tts_pair_shift_configuration_t

  !> İki measured isotherm arasındaki horizontal relative shift ve TRS
  !! kanıtlarını taşır. delta_s, moving eğrisine x_shifted=x+delta_s olarak
  !! uygulanır. Başarı yalnız matematiksel minimum bulunduğunu ifade eder.
  type, public :: tts_pair_shift_result_t
    integer :: status = PAIR_SHIFT_INVALID_INPUT
    character(len=tts_text_length) :: reference_isotherm_identifier = ""
    character(len=tts_text_length) :: moving_isotherm_identifier = ""
    integer :: reference_isotherm_index = 0
    integer :: moving_isotherm_index = 0
    integer :: production_channel = 0
    logical :: shift_available = .false.
    logical :: joint_shift_available = .false.
    logical :: storage_shift_available = .false.
    logical :: loss_shift_available = .false.
    real(dp) :: delta_s = 0.0_dp
    real(dp) :: delta_s_joint = 0.0_dp
    real(dp) :: delta_s_storage = 0.0_dp
    real(dp) :: delta_s_loss = 0.0_dp
    real(dp) :: storage_loss_shift_discrepancy = 0.0_dp
    real(dp) :: objective_minimum = huge(1.0_dp)
    real(dp) :: overlap_width_decades = 0.0_dp
    real(dp) :: overlap_fraction = 0.0_dp
    real(dp) :: storage_overlap_width_decades = 0.0_dp
    real(dp) :: loss_overlap_width_decades = 0.0_dp
    real(dp) :: objective_curvature = 0.0_dp
    integer :: iteration_count = 0
    integer :: evaluation_count = 0
  end type tts_pair_shift_result_t

  !> Reference-anchored empirical temperature shift tablosunun tek satırıdır.
  !! a_T boyutsuz, s=log10(a_T), sıcaklık K'dir.
  type, public :: tts_empirical_shift_t
    integer :: source_isotherm_index = 0
    character(len=tts_text_length) :: source_isotherm_identifier = ""
    real(dp) :: temperature_k = 0.0_dp
    real(dp) :: log10_a_t = 0.0_dp
    real(dp) :: a_t = 1.0_dp
  end type tts_empirical_shift_t

  !> Temperature-sorted adjacent linklerden reference-anchored absolute shift
  !! tablosu üretiminin ara sonucudur. Bir zorunlu link çözülemezse complete
  !! table unavailable olur; non-adjacent bridge denenmez.
  type, public :: tts_shift_chain_result_t
    integer :: status = TTS_IDENTIFICATION_INVALID_INPUT
    logical :: available = .false.
    integer :: reference_isotherm_index = 0
    type(tts_pair_shift_result_t), allocatable :: pair_shift_results(:)
    type(tts_empirical_shift_t), allocatable :: empirical_shifts(:)
  end type tts_shift_chain_result_t

  !> Original ölçümün değerlerini değiştirmeden shifted master-coordinate
  !! provenance'ını taşır. reduced frequency [Hz], shift boyutsuzdur.
  type, public :: tts_master_cloud_point_t
    integer :: source_isotherm_index = 0
    integer :: source_point_index = 0
    character(len=tts_text_length) :: source_isotherm_identifier = ""
    character(len=tts_text_length) :: specimen_identifier = ""
    character(len=tts_text_length) :: source_identifier = ""
    real(dp) :: source_temperature_k = 0.0_dp
    real(dp) :: source_frequency_hz = 0.0_dp
    real(dp) :: log10_a_t = 0.0_dp
    real(dp) :: reduced_frequency_hz = 0.0_dp
    real(dp) :: storage_modulus_pa = 0.0_dp
    real(dp) :: loss_modulus_pa = 0.0_dp
    integer :: storage_quality = MEASUREMENT_UNAVAILABLE
    integer :: loss_quality = MEASUREMENT_UNAVAILABLE
    logical :: contributes_to_validation = .false.
    logical :: contributes_to_runtime_extension = .false.
  end type tts_master_cloud_point_t

  !> Strictly increasing, single-valued solver master tablosunun bir satırıdır.
  !! Yalnız VALID G'>0 ve VALID G''>=0 ölçümleri girebilir; averaging yoktur.
  type, public :: tts_runtime_master_point_t
    real(dp) :: reduced_frequency_hz = 0.0_dp
    real(dp) :: storage_modulus_pa = 0.0_dp
    real(dp) :: loss_modulus_pa = 0.0_dp
    integer :: source_isotherm_index = 0
    integer :: source_point_index = 0
    character(len=tts_text_length) :: source_isotherm_identifier = ""
  end type tts_runtime_master_point_t

  !> İki farklı source bölgesi arasındaki master stitching tanısıdır. Gap
  !! log10(f_r) decade cinsindedir; overlap varsa log-modulus mismatch saklanır.
  type, public :: tts_master_boundary_diagnostic_t
    integer :: left_runtime_point_index = 0
    integer :: right_runtime_point_index = 0
    real(dp) :: boundary_gap_decades = 0.0_dp
    logical :: has_overlap = .false.
    real(dp) :: storage_log10_mismatch = 0.0_dp
    real(dp) :: loss_log10_mismatch = 0.0_dp
  end type tts_master_boundary_diagnostic_t

  !> Van Gurp-Palmen cloud noktasıdır. |G*| [Pa], delta [rad]; frekans shift'i
  !! gerekmez ve measured source provenance korunur.
  type, public :: tts_vgp_point_t
    integer :: source_isotherm_index = 0
    integer :: source_point_index = 0
    real(dp) :: source_temperature_k = 0.0_dp
    real(dp) :: complex_modulus_magnitude_pa = 0.0_dp
    real(dp) :: phase_angle_rad = 0.0_dp
  end type tts_vgp_point_t

  !> Cole-Cole cloud noktasıdır. Primary gösterim linear-axis G'' [Pa] vs
  !! G' [Pa] değerleridir; temperature/source provenance korunur.
  type, public :: tts_cole_cole_point_t
    integer :: source_isotherm_index = 0
    integer :: source_point_index = 0
    real(dp) :: source_temperature_k = 0.0_dp
    real(dp) :: storage_modulus_pa = 0.0_dp
    real(dp) :: loss_modulus_pa = 0.0_dp
  end type tts_cole_cole_point_t

  type, public :: tts_trs_diagnostics_t
    type(tts_vgp_point_t), allocatable :: vgp_points(:)
    type(tts_cole_cole_point_t), allocatable :: cole_cole_points(:)
    type(tts_master_boundary_diagnostic_t), allocatable :: boundaries(:)
    logical :: full_complex_pair_support = .false.
  end type tts_trs_diagnostics_t

  !> Top-level identification sonucudur. Input family deep-copy edildiğinden
  !! caller'ın sonraki mutation'ı sonucu değiştirmez. SUCCESS yalnız
  !! matematiksel construction başarısıdır; universal TRS PASS anlamı taşımaz.
  type, public :: tts_identification_result_t
    integer :: status = TTS_IDENTIFICATION_INVALID_INPUT
    character(len=tts_text_length) :: message = ""
    integer :: reference_isotherm_index = 0
    character(len=tts_text_length) :: reference_isotherm_identifier = ""
    real(dp) :: reference_temperature_k = 0.0_dp
    type(tts_material_family_t) :: source_family
    type(tts_pair_shift_result_t), allocatable :: pair_shift_results(:)
    type(tts_empirical_shift_t), allocatable :: empirical_shifts(:)
    type(tts_master_cloud_point_t), allocatable :: master_cloud(:)
    type(tts_runtime_master_point_t), allocatable :: runtime_master_table(:)
    type(tts_trs_diagnostics_t) :: diagnostics
    logical :: shift_chain_available = .false.
    logical :: master_cloud_available = .false.
    logical :: runtime_export_ready = .false.
  end type tts_identification_result_t

  public :: validate_tts_material_family
  public :: is_measurement_quality_known
  public :: is_storage_log_usable
  public :: is_loss_log_usable
  public :: is_runtime_export_usable
  public :: are_tts_values_machine_equivalent

contains

  !> Quality integer'ının tanımlı explicit ölçüm semantiğine ait olduğunu
  !! döndürür. Magic numeric quality değerlerini engeller.
  pure elemental function is_measurement_quality_known(quality) result(known)
    integer, intent(in) :: quality
    logical :: known

    known = quality == MEASUREMENT_VALID .or. &
      quality == BELOW_RELIABLE_FLOOR .or. &
      quality == MEASUREMENT_UNAVAILABLE .or. &
      quality == MEASUREMENT_REJECTED
  end function is_measurement_quality_known

  !> Storage log-objective kullanımı için G'>0, finite ve VALID koşulunu
  !! uygular. Girdi Pa, çıktı boyutsuz logical'dır.
  pure elemental function is_storage_log_usable(point) result(usable)
    type(tts_measurement_point_t), intent(in) :: point
    logical :: usable

    usable = point%storage_quality == MEASUREMENT_VALID .and. &
      ieee_is_finite(point%storage_modulus_pa) .and. &
      point%storage_modulus_pa > 0.0_dp
  end function is_storage_log_usable

  !> Loss log-objective kullanımı için G''>0, finite ve VALID koşulunu
  !! uygular. VALID G''=0 bu objective'te false döner; epsilon eklenmez.
  pure elemental function is_loss_log_usable(point) result(usable)
    type(tts_measurement_point_t), intent(in) :: point
    logical :: usable

    usable = point%loss_quality == MEASUREMENT_VALID .and. &
      ieee_is_finite(point%loss_modulus_pa) .and. &
      point%loss_modulus_pa > 0.0_dp
  end function is_loss_log_usable

  !> Solver-ready runtime export için iki channel'ın authoritative koşulunu
  !! uygular: VALID G'>0 ve VALID G''>=0. G''=0 geçerlidir.
  pure elemental function is_runtime_export_usable(point) result(usable)
    type(tts_measurement_point_t), intent(in) :: point
    logical :: usable

    usable = is_storage_log_usable(point) .and. &
      point%loss_quality == MEASUREMENT_VALID .and. &
      ieee_is_finite(point%loss_modulus_pa) .and. &
      point%loss_modulus_pa >= 0.0_dp
  end function is_runtime_export_usable

  !> İki sonlu değerin yalnız representation düzeyindeki eşdeğerliğini sınar.
  !! Bu tolerans deney uncertainty'si veya physical acceptance değildir.
  pure elemental function are_tts_values_machine_equivalent(a, b) &
      result(equivalent)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: equivalent
    real(dp) :: scale

    scale = max(abs(a), abs(b), tiny(1.0_dp))
    equivalent = ieee_is_finite(a) .and. ieee_is_finite(b) .and. &
      abs(a - b) <= 64.0_dp*epsilon(1.0_dp)*scale
  end function are_tts_values_machine_equivalent

  !> Family-level common state, measured isotherm ve point quality
  !! sözleşmelerini doğrular. Frekanslar [Hz] sonlu, pozitif ve strictly
  !! increasing; sıcaklıklar [K] sonlu/pozitif olmalıdır. Authoritative VALID
  !! G'>0 ve G''>=0 koşulları uygulanır. Input sırası değiştirilmez ve validation
  !! failure error-stop yerine açık result/status döndürür.
  pure function validate_tts_material_family(family) result(validation)
    type(tts_material_family_t), intent(in) :: family
    type(tts_validation_result_t) :: validation

    integer :: i
    integer :: j
    integer :: point_index

    validation%status = TTS_IDENTIFICATION_INVALID_INPUT
    if (len_trim(family%family_identifier) == 0 .or. &
        len_trim(family%common_state%material_identifier) == 0 .or. &
        len_trim(family%common_state%batch_state_identifier) == 0) then
      validation%message = "Family/material/batch kimlikleri boş olamaz."
      return
    end if
    if (.not. ieee_is_finite( &
        family%common_state%dynamic_strain_amplitude_ratio) .or. &
        family%common_state%dynamic_strain_amplitude_ratio <= 0.0_dp .or. &
        .not. ieee_is_finite(family%common_state%static_prestrain_ratio) .or. &
        family%common_state%static_prestrain_ratio < 0.0_dp) then
      validation%status = TTS_IDENTIFICATION_INCONSISTENT_STATE
      validation%message = "Common strain amplitude/prestrain geçersiz."
      return
    end if
    if (family%common_state%deformation_mode /= &
        TTS_DEFORMATION_MODE_SHEAR) then
      validation%status = TTS_IDENTIFICATION_INCONSISTENT_STATE
      validation%message = "V0.8.1 yalnız common dynamic-shear state destekler."
      return
    end if
    if (len_trim(family%common_state%conditioning_description) == 0 .or. &
        len_trim(family%common_state%test_method) == 0 .or. &
        len_trim(family%common_state%source_metadata) == 0) then
      validation%status = TTS_IDENTIFICATION_INCONSISTENT_STATE
      validation%message = "Common conditioning/method/source metadata eksik."
      return
    end if
    if (.not. allocated(family%isotherms)) then
      validation%status = TTS_IDENTIFICATION_INSUFFICIENT_ISOTHERMS
      validation%message = "Identification en az iki measured isotherm ister."
      return
    end if
    if (size(family%isotherms) < 2) then
      validation%status = TTS_IDENTIFICATION_INSUFFICIENT_ISOTHERMS
      validation%message = "Identification en az iki measured isotherm ister."
      return
    end if

    do i = 1, size(family%isotherms)
      if (len_trim(family%isotherms(i)%isotherm_identifier) == 0 .or. &
          len_trim(family%isotherms(i)%specimen_identifier) == 0 .or. &
          len_trim(family%isotherms(i)%source_identifier) == 0) then
        validation%message = "Isotherm/specimen/source kimliği eksik."
        return
      end if
      if (.not. ieee_is_finite(family%isotherms(i)%temperature_k) .or. &
          family%isotherms(i)%temperature_k <= 0.0_dp) then
        validation%message = "Isotherm sıcaklığı sonlu ve pozitif K olmalıdır."
        return
      end if
      if (.not. allocated(family%isotherms(i)%points)) then
        validation%message = "Her isotherm en az iki measured point ister."
        return
      end if
      if (size(family%isotherms(i)%points) < 2) then
        validation%message = "Her isotherm en az iki measured point ister."
        return
      end if

      do j = 1, i - 1
        if (trim(family%isotherms(i)%isotherm_identifier) == &
            trim(family%isotherms(j)%isotherm_identifier)) then
          validation%message = "Isotherm kimlikleri unique olmalıdır."
          return
        end if
        if (are_tts_values_machine_equivalent( &
            family%isotherms(i)%temperature_k, &
            family%isotherms(j)%temperature_k)) then
          validation%message = "Reference implementation duplicate T içeremez."
          return
        end if
      end do

      do point_index = 1, size(family%isotherms(i)%points)
        if (.not. ieee_is_finite( &
            family%isotherms(i)%points(point_index)%frequency_hz) .or. &
            family%isotherms(i)%points(point_index)%frequency_hz <= 0.0_dp) then
          validation%message = "Measured frequency sonlu ve pozitif Hz olmalıdır."
          return
        end if
        if (point_index > 1) then
          if (family%isotherms(i)%points(point_index)%frequency_hz <= &
              family%isotherms(i)%points(point_index - 1)%frequency_hz .or. &
              are_tts_values_machine_equivalent( &
                family%isotherms(i)%points(point_index)%frequency_hz, &
                family%isotherms(i)%points(point_index - 1)%frequency_hz)) then
            validation%message = &
              "Measured frequencies strictly increasing olmalıdır."
            return
          end if
        end if
        if (.not. is_measurement_quality_known( &
            family%isotherms(i)%points(point_index)%storage_quality) .or. &
            .not. is_measurement_quality_known( &
              family%isotherms(i)%points(point_index)%loss_quality)) then
          validation%message = "Tanımsız measurement quality değeri."
          return
        end if
        if (family%isotherms(i)%points(point_index)%storage_quality == &
            MEASUREMENT_VALID .and. .not. is_storage_log_usable( &
              family%isotherms(i)%points(point_index))) then
          validation%message = "VALID storage modulus sonlu ve pozitif olmalıdır."
          return
        end if
        if (family%isotherms(i)%points(point_index)%loss_quality == &
            MEASUREMENT_VALID) then
          if (.not. ieee_is_finite( &
              family%isotherms(i)%points(point_index)%loss_modulus_pa) .or. &
              family%isotherms(i)%points(point_index)%loss_modulus_pa < &
                0.0_dp) then
            validation%message = &
              "VALID loss modulus sonlu ve negatif olmayan Pa olmalıdır."
            return
          end if
        end if
      end do
    end do

    validation%status = TTS_IDENTIFICATION_SUCCESS
    validation%valid = .true.
    validation%message = "Experimental TTS family doğrulandı."
  end function validate_tts_material_family

end module tms_tts_types
