module tms_material_aware_harmonic_response
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_dynamic_material_metadata, only : &
    validate_dynamic_material_metadata
  use tms_dynamic_modulus_provider, only : LINEAR_FREQUENCY, &
    LINEAR_LOG_FREQUENCY, are_machine_equivalent
  use tms_harmonic_response, only : harmonic_response_t, &
    validate_harmonic_response, get_harmonic_frequency_count, &
    get_harmonic_frequencies_hz
  use tms_material_state_trace, only : material_binding_trace_t, &
    material_state_trace_t
  implicit none
  private

  !> V0.6 harmonic response ile V0.7 material/dataset izlerini bileşim yoluyla
  !! taşır. Mevcut harmonic_response_t API'sini büyütmez veya değiştirmez.
  type, public :: material_aware_harmonic_response_t
    private
    type(harmonic_response_t) :: harmonic_response
    type(material_binding_trace_t), allocatable :: binding_traces(:)
    type(material_state_trace_t), allocatable :: state_traces(:, :)
  end type material_aware_harmonic_response_t

  public :: create_material_aware_harmonic_response
  public :: validate_material_aware_harmonic_response
  public :: get_base_harmonic_response
  public :: get_material_binding_count
  public :: get_material_binding_traces
  public :: get_material_state_traces
  public :: get_material_state_trace

contains

  !> Harmonic çözümü, binding-level metadata ve binding x frequency malzeme
  !! izleriyle tek immutable sonuç nesnesinde birleştirir.
  !! Girdi/çıktı: Harmonic response [rad ve solver diagnostics], C_theta [m^3],
  !! G'/G'' [Pa], K'/K'' [N*m/rad], f [Hz], T [K] kayıtlarıdır.
  !! Varsayım: Material evaluation bütün sweep için solver çağrısından önce
  !! tamamlanmıştır; bu nedenle singular çözüm sütunlarının da trace'i vardır.
  pure function create_material_aware_harmonic_response( &
      harmonic_response, binding_traces, state_traces) result(response)
    type(harmonic_response_t), intent(in) :: harmonic_response
    type(material_binding_trace_t), intent(in) :: binding_traces(:)
    type(material_state_trace_t), intent(in) :: state_traces(:, :)
    type(material_aware_harmonic_response_t) :: response

    response%harmonic_response = harmonic_response
    response%binding_traces = binding_traces
    response%state_traces = state_traces
    call validate_material_aware_harmonic_response(response)
  end function create_material_aware_harmonic_response

  !> Wrapper sonuç boyutlarını, dataset trace'ini ve constitutive passivity
  !! invariantlarını doğrular. Her binding/frequency kaydı harmonic grid ile
  !! eşleşmeli; G'>0, G''>=0, K'>0, K''>=0 ve tan(delta)>=0 olmalıdır.
  pure subroutine validate_material_aware_harmonic_response(response)
    type(material_aware_harmonic_response_t), intent(in) :: response

    real(dp), allocatable :: frequencies_hz(:)
    integer :: binding_index
    integer :: frequency_index

    call validate_harmonic_response(response%harmonic_response)
    if (.not. allocated(response%binding_traces) .or. &
        .not. allocated(response%state_traces)) then
      error stop "Material-aware harmonic response trace içermelidir."
    end if
    if (size(response%binding_traces) == 0) then
      error stop "Material-aware response en az bir dynamic binding gerektirir."
    end if
    if (size(response%state_traces, 1) /= &
        size(response%binding_traces) .or. &
        size(response%state_traces, 2) /= &
        get_harmonic_frequency_count(response%harmonic_response)) then
      error stop "Material-aware response trace boyutları harmonic grid ile uyumsuz."
    end if

    frequencies_hz = get_harmonic_frequencies_hz( &
      response%harmonic_response)
    do binding_index = 1, size(response%binding_traces)
      if (response%binding_traces(binding_index)%element_id <= 0 .or. &
          .not. ieee_is_finite( &
            response%binding_traces(binding_index)%geometry_factor_m3) .or. &
          response%binding_traces(binding_index)%geometry_factor_m3 <= 0.0_dp) then
        error stop "Material binding trace kimliği veya C_theta değeri geçersiz."
      end if
      call validate_dynamic_material_metadata( &
        response%binding_traces(binding_index)%metadata)

      do frequency_index = 1, size(frequencies_hz)
        call validate_state_trace( &
          response%state_traces(binding_index, frequency_index), &
          response%binding_traces(binding_index)%element_id, &
          frequencies_hz(frequency_index), &
          response%binding_traces(binding_index)%geometry_factor_m3)
      end do
    end do
  end subroutine validate_material_aware_harmonic_response

  !> İç V0.6 harmonic_response_t nesnesinin bağımsız kopyasını döndürür.
  pure function get_base_harmonic_response(response) result(harmonic_response)
    type(material_aware_harmonic_response_t), intent(in) :: response
    type(harmonic_response_t) :: harmonic_response

    call validate_material_aware_harmonic_response(response)
    harmonic_response = response%harmonic_response
  end function get_base_harmonic_response

  !> Dynamic material binding sayısını [-] döndürür.
  pure function get_material_binding_count(response) result(binding_count)
    type(material_aware_harmonic_response_t), intent(in) :: response
    integer :: binding_count

    call validate_material_aware_harmonic_response(response)
    binding_count = size(response%binding_traces)
  end function get_material_binding_count

  !> Dataset-level binding trace dizisinin bağımsız kopyasını döndürür.
  pure function get_material_binding_traces(response) result(traces)
    type(material_aware_harmonic_response_t), intent(in) :: response
    type(material_binding_trace_t), allocatable :: traces(:)

    call validate_material_aware_harmonic_response(response)
    traces = response%binding_traces
  end function get_material_binding_traces

  !> Binding x frequency sıralı material state trace matrisinin bağımsız
  !! kopyasını döndürür.
  pure function get_material_state_traces(response) result(traces)
    type(material_aware_harmonic_response_t), intent(in) :: response
    type(material_state_trace_t), allocatable :: traces(:, :)

    call validate_material_aware_harmonic_response(response)
    traces = response%state_traces
  end function get_material_state_traces

  !> Bir tabanlı binding/frequency indeksindeki material trace'i döndürür.
  pure function get_material_state_trace( &
      response, binding_index, frequency_index) result(trace)
    type(material_aware_harmonic_response_t), intent(in) :: response
    integer, intent(in) :: binding_index
    integer, intent(in) :: frequency_index
    type(material_state_trace_t) :: trace

    call validate_material_aware_harmonic_response(response)
    if (binding_index < 1 .or. &
        binding_index > size(response%state_traces, 1) .or. &
        frequency_index < 1 .or. &
        frequency_index > size(response%state_traces, 2)) then
      error stop "Material state trace indeksi sınır dışında."
    end if
    trace = response%state_traces(binding_index, frequency_index)
  end function get_material_state_trace

  !> Tek trace kaydının harmonic frequency, binding kimliği, SI sonluluğu,
  !! interpolation bracket'i ve passive G*/K* koşullarını doğrular.
  pure subroutine validate_state_trace( &
      trace, element_id, frequency_hz, geometry_factor_m3)
    type(material_state_trace_t), intent(in) :: trace
    integer, intent(in) :: element_id
    real(dp), intent(in) :: frequency_hz
    real(dp), intent(in) :: geometry_factor_m3

    real(dp) :: expected_loss_factor
    real(dp) :: expected_loss_stiffness
    real(dp) :: expected_storage_stiffness

    if (trace%element_id /= element_id .or. &
        .not. are_machine_equivalent(trace%frequency_hz, frequency_hz)) then
      error stop "Material state trace eleman veya frequency eşlemesi hatalı."
    end if
    if (.not. ieee_is_finite(trace%temperature_k) .or. &
        trace%temperature_k <= 0.0_dp .or. &
        .not. ieee_is_finite(trace%storage_modulus_pa) .or. &
        trace%storage_modulus_pa <= 0.0_dp .or. &
        .not. ieee_is_finite(trace%loss_modulus_pa) .or. &
        trace%loss_modulus_pa < 0.0_dp .or. &
        .not. ieee_is_finite(trace%storage_stiffness_nm_per_rad) .or. &
        trace%storage_stiffness_nm_per_rad <= 0.0_dp .or. &
        .not. ieee_is_finite(trace%loss_stiffness_nm_per_rad) .or. &
        trace%loss_stiffness_nm_per_rad < 0.0_dp .or. &
        .not. ieee_is_finite(trace%loss_factor) .or. &
        trace%loss_factor < 0.0_dp) then
      error stop "Material state trace passive ve sonlu SI değerleri sağlamıyor."
    end if
    if (trace%interpolation_policy /= LINEAR_FREQUENCY .and. &
        trace%interpolation_policy /= LINEAR_LOG_FREQUENCY) then
      error stop "Material state trace interpolation policy geçersiz."
    end if
    if (.not. ieee_is_finite(trace%lower_frequency_hz) .or. &
        .not. ieee_is_finite(trace%upper_frequency_hz) .or. &
        .not. ieee_is_finite(trace%interpolation_alpha) .or. &
        trace%lower_frequency_hz <= 0.0_dp .or. &
        trace%upper_frequency_hz < trace%lower_frequency_hz .or. &
        trace%interpolation_alpha < 0.0_dp .or. &
        trace%interpolation_alpha > 1.0_dp) then
      error stop "Material state trace interpolation bracket'i geçersiz."
    end if

    expected_storage_stiffness = geometry_factor_m3*trace%storage_modulus_pa
    expected_loss_stiffness = geometry_factor_m3*trace%loss_modulus_pa
    expected_loss_factor = trace%loss_modulus_pa/trace%storage_modulus_pa
    if (.not. are_trace_values_close( &
        trace%storage_stiffness_nm_per_rad, expected_storage_stiffness) .or. &
        .not. are_trace_values_close( &
          trace%loss_stiffness_nm_per_rad, expected_loss_stiffness) .or. &
        .not. are_trace_values_close(trace%loss_factor, expected_loss_factor)) then
      error stop "Material state trace G* -> K* mapping invariantı hatalı."
    end if

    if (trace%exact_table_point) then
      if (.not. are_machine_equivalent( &
          trace%lower_frequency_hz, trace%upper_frequency_hz) .or. &
          abs(trace%interpolation_alpha) > 0.0_dp) then
        error stop "Exact material trace bracket/alpha semantiği hatalı."
      end if
    else
      if (trace%lower_frequency_hz >= trace%upper_frequency_hz .or. &
          trace%frequency_hz <= trace%lower_frequency_hz .or. &
          trace%frequency_hz >= trace%upper_frequency_hz) then
        error stop "Interpolated material trace bracket semantiği hatalı."
      end if
    end if
  end subroutine validate_state_trace

  !> Aynı SI birimindeki trace ve bağımsız türetilmiş değeri yalnız
  !! floating-point roundoff düzeyinde karşılaştırır.
  pure elemental function are_trace_values_close(a, b) result(is_close)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    logical :: is_close

    is_close = abs(a-b) <= &
      256.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(a), abs(b))
  end function are_trace_values_close

end module tms_material_aware_harmonic_response
