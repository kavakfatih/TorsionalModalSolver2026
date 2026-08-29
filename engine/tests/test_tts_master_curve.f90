program test_tts_master_curve
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tms_kinds, only : dp
  use tms_tts_types, only : tts_material_family_t, tts_shift_chain_result_t, &
    tts_empirical_shift_t, tts_master_cloud_point_t, &
    tts_runtime_master_point_t, tts_master_boundary_diagnostic_t, &
    BELOW_RELIABLE_FLOOR, MEASUREMENT_VALID, MEASUREMENT_UNAVAILABLE, &
    TTS_IDENTIFICATION_RUNTIME_DOMAIN_GAP, &
    TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED
  use tms_tts_shift_chain, only : build_tts_shift_chain
  use tms_tts_master_curve, only : build_tts_master_experimental_cloud, &
    stitch_tts_runtime_master_table
  use tms_tts_test_support, only : make_exact_trs_family, &
    populate_isotherm_from_log_grid, assert_true, assert_close
  implicit none

  real(dp), parameter :: temperatures(3) = &
    [273.15_dp, 293.15_dp, 313.15_dp]
  real(dp), parameter :: shifts(3) = [1.0_dp, 0.0_dp, -1.0_dp]
  real(dp), parameter :: narrow_grid(3) = [-1.0_dp, 0.0_dp, 1.0_dp]
  real(dp), parameter :: wide_grid(5) = &
    [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
  real(dp), parameter :: low_extension_grid(5) = &
    [-3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp]
  type(tts_material_family_t) :: family
  type(tts_shift_chain_result_t) :: chain
  type(tts_empirical_shift_t), allocatable :: empirical(:)
  type(tts_master_cloud_point_t), allocatable :: cloud(:)
  type(tts_runtime_master_point_t), allocatable :: table(:)
  type(tts_master_boundary_diagnostic_t), allocatable :: boundaries(:)
  logical :: success
  logical :: found
  integer :: i
  integer :: stitch_status

  family = make_exact_trs_family(temperatures, shifts)
  chain = build_tts_shift_chain(family, "ISO-2")
  call assert_true(chain%available, "Master test shift chain kurulamadı.")
  call build_tts_master_experimental_cloud(family, chain%empirical_shifts, &
    cloud, success)
  call assert_true(success .and. size(cloud) == 21, &
    "Experimental cloud bütün original points'i korumadı.")
  call assert_true(cloud(1)%source_isotherm_index == 1 .and. &
    cloud(1)%source_point_index == 1 .and. &
    trim(cloud(1)%specimen_identifier) == "SPECIMEN-1", &
    "Master cloud source/specimen provenance kayboldu.")
  call assert_close(cloud(1)%reduced_frequency_hz, &
    10.0_dp*cloud(1)%source_frequency_hz, 1.0e-6_dp, &
    "Master cloud f_r=10^(x+s) koordinatı hatalı.")

  call stitch_tts_runtime_master_table(family, 2, chain%empirical_shifts, &
    cloud, table, boundaries, success)
  call assert_true(success .and. size(table) > &
    size(family%isotherms(2)%points), &
    "Low/high reference-centered extension oluşmadı.")
  do i = 2, size(table)
    call assert_true(table(i)%reduced_frequency_hz > &
      table(i - 1)%reduced_frequency_hz, &
      "Runtime master table strictly increasing değil.")
  end do
  found = .false.
  do i = 1, size(table)
    if (abs(log10(table(i)%reduced_frequency_hz)) < 1.0e-12_dp) then
      found = table(i)%source_isotherm_index == 2
    end if
  end do
  call assert_true(found, &
    "Duplicate reduced frequency'de reference priority korunmadı.")
  call assert_true(size(boundaries) >= 2, &
    "Master source transition boundary diagnostics üretilmedi.")
  do i = 1, size(boundaries)
    call assert_true(boundaries(i)%boundary_gap_decades > 0.0_dp, &
      "Boundary log-frequency gap tanısı pozitif değil.")
    if (boundaries(i)%has_overlap) then
      call assert_close(boundaries(i)%storage_log10_mismatch, 0.0_dp, &
        2.0e-12_dp, "Exact TRS boundary G' mismatch sıfır değil.")
      call assert_close(boundaries(i)%loss_log10_mismatch, 0.0_dp, &
        2.0e-12_dp, "Exact TRS boundary G'' mismatch sıfır değil.")
    end if
  end do

  ! Tek adjacent isotherm reference aralığını iki tarafta genişletebilmelidir.
  family = make_exact_trs_family([293.15_dp, 303.15_dp], [0.0_dp, 0.0_dp])
  call populate_isotherm_from_log_grid(family%isotherms(1), narrow_grid, &
    0.0_dp, 0.25_dp, 0.15_dp)
  call populate_isotherm_from_log_grid(family%isotherms(2), wide_grid, &
    0.0_dp, 0.25_dp, 0.15_dp)
  empirical = make_zero_shift_table(family)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success)
  call assert_true(success .and. size(table) == 5 .and. &
    table(1)%source_isotherm_index == 2 .and. &
    table(5)%source_isotherm_index == 2, &
    "Both-side extension aynı isotherm'den eklenemedi.")

  ! Yeni isotherm range genişletmiyorsa error olmamalı ve overlap points
  ! runtime table'da duplicate edilmemelidir.
  call populate_isotherm_from_log_grid(family%isotherms(2), narrow_grid, &
    0.0_dp, 0.25_dp, 0.15_dp)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success)
  call assert_true(success .and. size(table) == 3 .and. &
    all(table%source_isotherm_index == 1), &
    "No-extension isotherm duplicate runtime points ekledi.")

  ! Aynı low-side extension'i sunan iki curve'de reference'a temperature
  ! olarak daha yakın source kazanmalıdır.
  family = make_exact_trs_family(temperatures, [0.0_dp, 0.0_dp, 0.0_dp])
  call populate_isotherm_from_log_grid(family%isotherms(2), narrow_grid, &
    0.0_dp, 0.25_dp, 0.15_dp)
  call populate_isotherm_from_log_grid(family%isotherms(1), &
    low_extension_grid, 0.0_dp, 0.25_dp, 0.15_dp)
  call populate_isotherm_from_log_grid(family%isotherms(3), &
    low_extension_grid, 0.0_dp, 0.25_dp, 0.15_dp)
  empirical = make_zero_shift_table(family)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call stitch_tts_runtime_master_table(family, 2, empirical, cloud, table, &
    boundaries, success)
  call assert_true(success .and. size(table) == 5 .and. &
    table(1)%source_isotherm_index == 1 .and. &
    table(2)%source_isotherm_index == 1, &
    "Nearest-temperature duplicate priority uygulanmadı.")

  ! BELOW_RELIABLE_FLOOR numeric value taşımış olsa da authoritative runtime
  ! extension'a giremez. Aynı VALID G''=0 ise girebilir.
  family%isotherms(1)%points(1)%loss_quality = BELOW_RELIABLE_FLOOR
  family%isotherms(1)%points(2)%loss_modulus_pa = 0.0_dp
  family%isotherms(1)%points(2)%loss_quality = MEASUREMENT_VALID
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call stitch_tts_runtime_master_table(family, 2, empirical, cloud, table, &
    boundaries, success)
  found = .false.
  do i = 1, size(table)
    if (table(i)%source_isotherm_index == 1 .and. &
        table(i)%source_point_index == 1) found = .true.
  end do
  call assert_true(.not. found, &
    "Unreliable loss-quality point runtime table'a sızdı.")

  ! İki source'ta da aynı internal measurement-quality hole varsa final
  ! runtime table bu aralığı V0.8.0 interpolation'ına açmamalıdır. Experimental
  ! cloud korunur, fakat unsupported continuous solver domain reddedilir.
  family = make_exact_trs_family( &
    [293.15_dp, 313.15_dp], [0.0_dp, 0.0_dp])
  family%isotherms(1)%points(4)%loss_quality = BELOW_RELIABLE_FLOOR
  family%isotherms(2)%points(4)%loss_quality = BELOW_RELIABLE_FLOOR
  empirical = make_zero_shift_table(family)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call assert_true(success .and. size(cloud) == 14, &
    "Internal gap experimental cloud provenance'ını kaybetti.")
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success, stitch_status)
  call assert_true(.not. success, &
    "Unsupported internal quality gap runtime'a açıldı.")
  call assert_true(stitch_status == TTS_IDENTIFICATION_RUNTIME_DOMAIN_GAP, &
    "Internal quality gap explicit status üretmedi.")
  call assert_true(.not. allocated(table), &
    "Unsupported gap için unsafe runtime table bırakıldı.")

  ! Edge'deki unavailable ölçümler solver domain'ini ilk contiguous usable
  ! noktaya daraltabilir; internal hole olmadığı için export başarısız olmaz.
  family = make_exact_trs_family( &
    [293.15_dp, 313.15_dp], [0.0_dp, 0.0_dp])
  do i = 1, 2
    family%isotherms(1)%points(i)%loss_quality = MEASUREMENT_UNAVAILABLE
    family%isotherms(2)%points(i)%loss_quality = MEASUREMENT_UNAVAILABLE
  end do
  empirical = make_zero_shift_table(family)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success, stitch_status)
  call assert_true(success, "Edge quality gap usable domain'e daraltılamadı.")
  call assert_close(table(1)%reduced_frequency_hz, &
    family%isotherms(1)%points(3)%frequency_hz, 1.0e-14_dp, &
    "Runtime minimum frequency ilk usable edge point'ten başlamadı.")

  ! Reference'taki internal hole, ikinci isotherm'in shifted adjacent valid
  ! interval'ları tarafından gerçekten örtülüyorsa global union continuous'tur.
  family = make_exact_trs_family( &
    [293.15_dp, 313.15_dp], [0.0_dp, 0.0_dp])
  family%isotherms(1)%points(4)%loss_quality = MEASUREMENT_UNAVAILABLE
  empirical = make_zero_shift_table(family)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success, stitch_status)
  call assert_true(success, &
    "Başka isotherm'in valid coverage bridge'i kabul edilmedi.")

  ! VALID G''=0 loss-log objective'e giremez fakat adjacent runtime coverage
  ! interval'ını kesmez ve zero-loss measured point solver table'da korunur.
  family = make_exact_trs_family( &
    [293.15_dp, 313.15_dp], [0.0_dp, 0.0_dp])
  family%isotherms(1)%points(4)%loss_modulus_pa = 0.0_dp
  family%isotherms(2)%points(4)%loss_modulus_pa = 0.0_dp
  empirical = make_zero_shift_table(family)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success, stitch_status)
  call assert_true(success .and. any(table%loss_modulus_pa == 0.0_dp), &
    "VALID zero-loss adjacent runtime coverage'dan çıkarıldı.")

  ! Extreme shift'ler direct stitching API'sinde fake huge/tiny frequency
  ! üretmemeli; overflow, underflow ve nonfinite shift clean failure olmalıdır.
  family = make_exact_trs_family( &
    [293.15_dp, 313.15_dp], [0.0_dp, 0.0_dp])
  empirical = make_zero_shift_table(family)
  call build_tts_master_experimental_cloud(family, empirical, cloud, success)
  empirical(1)%log10_a_t = log10(huge(1.0_dp)) + 1.0_dp
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success, stitch_status)
  call assert_true(.not. success .and. &
    stitch_status == TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED, &
    "Reduced-frequency overflow clean stitching failure üretmedi.")
  call assert_true(.not. allocated(table), &
    "Overflow fake huge runtime table bıraktı.")

  empirical = make_zero_shift_table(family)
  empirical(1)%log10_a_t = log10(tiny(1.0_dp)) - 1.0_dp
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success, stitch_status)
  call assert_true(.not. success .and. &
    stitch_status == TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED, &
    "Reduced-frequency underflow clean stitching failure üretmedi.")
  call assert_true(.not. allocated(table), &
    "Underflow fake tiny runtime table bıraktı.")

  empirical = make_zero_shift_table(family)
  empirical(1)%log10_a_t = ieee_value(1.0_dp, ieee_quiet_nan)
  call stitch_tts_runtime_master_table(family, 1, empirical, cloud, table, &
    boundaries, success, stitch_status)
  call assert_true(.not. success .and. &
    stitch_status == TTS_IDENTIFICATION_RUNTIME_EXPORT_FAILED, &
    "Nonfinite empirical shift clean stitching failure üretmedi.")

  print *, "V0.8.1 master cloud ve deterministic stitching doğrulandı."

contains

  function make_zero_shift_table(source_family) result(result_shifts)
    type(tts_material_family_t), intent(in) :: source_family
    type(tts_empirical_shift_t), allocatable :: result_shifts(:)
    integer :: index_value

    allocate(result_shifts(size(source_family%isotherms)))
    do index_value = 1, size(source_family%isotherms)
      result_shifts(index_value)%source_isotherm_index = index_value
      result_shifts(index_value)%source_isotherm_identifier = &
        source_family%isotherms(index_value)%isotherm_identifier
      result_shifts(index_value)%temperature_k = &
        source_family%isotherms(index_value)%temperature_k
      result_shifts(index_value)%log10_a_t = 0.0_dp
      result_shifts(index_value)%a_t = 1.0_dp
    end do
  end function make_zero_shift_table

end program test_tts_master_curve
