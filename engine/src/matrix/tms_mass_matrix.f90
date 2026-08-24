module tms_mass_matrix
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_matrix_types, only : dense_matrix_t, initialize_dense_matrix, &
    validate_dense_matrix, add_dense_matrix_entry, get_dense_matrix_entry, &
    get_dense_matrix_values, get_dense_matrix_row_count, &
    get_dense_matrix_column_count
  implicit none
  private

  !> Full veya active torsional denklem kümesinin global dönel atalet
  !! matrisini taşır.
  !! Adı geleneksel hareket denklemi gösterimine uygun mass matrix olsa da
  !! katsayıları polar kütle ataleti J [kg*m^2] değerleridir.
  type, public :: mass_matrix_t
    private
    type(dense_matrix_t) :: storage
  end type mass_matrix_t

  public :: initialize_mass_matrix
  public :: validate_mass_matrix
  public :: add_nodal_inertia
  public :: get_mass_matrix_value
  public :: get_mass_matrix_values
  public :: get_mass_matrix_size
  public :: extract_mass_principal_submatrix

contains

  !> Global dönel atalet matrisini sıfır olarak başlatır.
  !! Girdi negatif olmayan logical equation sayısıdır [-]. Çıktı n x n,
  !! [kg*m^2] birimli mass_matrix_t değeridir; sıfır DOF 0x0 üretir.
  pure subroutine initialize_mass_matrix(matrix, dof_count)
    type(mass_matrix_t), intent(out) :: matrix
    integer, intent(in) :: dof_count

    if (dof_count < 0) then
      error stop "Global atalet matrisi DOF sayısı negatif olamaz."
    end if

    call initialize_dense_matrix(matrix%storage, dof_count, dof_count)
  end subroutine initialize_mass_matrix

  !> Global dönel atalet matrisinin kare ve başlatılmış olduğunu sınar.
  pure subroutine validate_mass_matrix(matrix)
    type(mass_matrix_t), intent(in) :: matrix

    call validate_dense_matrix(matrix%storage)
    if (get_dense_matrix_row_count(matrix%storage) /= &
        get_dense_matrix_column_count(matrix%storage)) then
      error stop "Global atalet matrisi kare olmalıdır."
    end if
  end subroutine validate_mass_matrix

  !> Düğümde yığılmış polar kütle ataletini global matris köşegenine ekler.
  !!
  !! Fiziksel model: M(equation_id,equation_id) += J. J [kg*m^2], açısal
  !! ivmeye karşı dönel eylemsizliği temsil eder.
  !! Girdiler: Boyutsuz denklem kimliği ve sonlu pozitif J [kg*m^2]. Çıktı:
  !! Güncellenmiş mass_matrix_t. equation_id=0 homojen kısıtlı DOF olduğu için
  !! indirgenmiş matrise eklenmez. Eleman ataleti ve tutarlı M uygulanmaz.
  pure subroutine add_nodal_inertia(matrix, equation_id, polar_inertia_kg_m2)
    type(mass_matrix_t), intent(inout) :: matrix
    integer, intent(in) :: equation_id
    real(dp), intent(in) :: polar_inertia_kg_m2

    integer :: matrix_size

    call validate_mass_matrix(matrix)
    matrix_size = get_mass_matrix_size(matrix)

    if (equation_id < 0 .or. equation_id > matrix_size) then
      error stop "Düğüm ataletinin denklem kimliği geçersiz."
    end if

    if (.not. ieee_is_finite(polar_inertia_kg_m2) .or. &
        polar_inertia_kg_m2 <= 0.0_dp) then
      error stop "Düğüm polar ataleti sonlu ve pozitif olmalıdır."
    end if

    if (equation_id == 0) return
    call add_dense_matrix_entry( &
      matrix%storage, equation_id, equation_id, polar_inertia_kg_m2)
  end subroutine add_nodal_inertia

  !> Global dönel atalet matrisinin tek katsayısını [kg*m^2] döndürür.
  pure function get_mass_matrix_value(matrix, row, column) result(value)
    type(mass_matrix_t), intent(in) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: column
    real(dp) :: value

    call validate_mass_matrix(matrix)
    value = get_dense_matrix_entry(matrix%storage, row, column)
  end function get_mass_matrix_value

  !> Global dönel atalet katsayılarının [kg*m^2] bağımsız kopyasını döndürür.
  pure function get_mass_matrix_values(matrix) result(values)
    type(mass_matrix_t), intent(in) :: matrix
    real(dp), allocatable :: values(:, :)

    call validate_mass_matrix(matrix)
    values = get_dense_matrix_values(matrix%storage)
  end function get_mass_matrix_values

  !> Global dönel atalet matrisinin boyutsuz satır/sütun sayısını döndürür.
  pure function get_mass_matrix_size(matrix) result(matrix_size)
    type(mass_matrix_t), intent(in) :: matrix
    integer :: matrix_size

    call validate_mass_matrix(matrix)
    matrix_size = get_dense_matrix_row_count(matrix%storage)
  end function get_mass_matrix_size

  !> Seçilen logical equation indekslerinin principal atalet alt matrisini
  !! storage ayrıntısını dışarı açmadan üretir.
  !!
  !! Fiziksel açıklama: Korunacak torsional serbestlik derecelerinin polar
  !! atalet katkıları full denklem takımından seçilir.
  !! Matematiksel açıklama: B(a,b)=A(indices(a),indices(b)); bu seçim
  !! P^T*A*P indirgemesinin katsayı biçimidir.
  !! Girdiler: Katsayıları [kg*m^2] olan matrix ve benzersiz, bir tabanlı,
  !! boyutsuz logical equation indeksleri. Çıktı: Aynı birimde n x n
  !! mass_matrix_t. Sıfır uzunluklu indeks dizisi 0x0 matris üretir.
  !! Varsayımlar ve geçerlilik: İndeks sırası çıktı denklem sırasıdır; tüm
  !! indeksler kaynak matris sınırlarında ve benzersiz olmalıdır.
  pure function extract_mass_principal_submatrix(matrix, indices) &
      result(submatrix)
    type(mass_matrix_t), intent(in) :: matrix
    integer, intent(in) :: indices(:)
    type(mass_matrix_t) :: submatrix

    integer :: column
    integer :: first_index
    integer :: matrix_size
    integer :: row
    integer :: second_index

    call validate_mass_matrix(matrix)
    matrix_size = get_mass_matrix_size(matrix)

    do first_index = 1, size(indices)
      if (indices(first_index) < 1 .or. &
          indices(first_index) > matrix_size) then
        error stop "Atalet alt matrisi indeksi kaynak matris sınırları dışında."
      end if

      do second_index = first_index + 1, size(indices)
        if (indices(first_index) == indices(second_index)) then
          error stop "Atalet alt matrisi indeksleri benzersiz olmalıdır."
        end if
      end do
    end do

    call initialize_mass_matrix(submatrix, size(indices))
    do row = 1, size(indices)
      do column = 1, size(indices)
        call add_dense_matrix_entry( &
          submatrix%storage, row, column, &
          get_dense_matrix_entry( &
          matrix%storage, indices(row), indices(column)))
      end do
    end do
  end function extract_mass_principal_submatrix

end module tms_mass_matrix
