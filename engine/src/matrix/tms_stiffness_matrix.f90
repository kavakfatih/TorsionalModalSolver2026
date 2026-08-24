module tms_stiffness_matrix
  use tms_kinds, only : dp
  use tms_local_matrix, only : local_matrix_2x2
  use tms_matrix_types, only : dense_matrix_t, initialize_dense_matrix, &
    validate_dense_matrix, add_dense_matrix_entry, get_dense_matrix_entry, &
    get_dense_matrix_values, get_dense_matrix_row_count, &
    get_dense_matrix_column_count
  implicit none
  private

  !> Aktif torsional denklemler için global rijitlik matrisini taşır.
  !! Katsayılar [N*m/rad] SI birimindedir. Dense depolama private olduğundan
  !! gelecekte sparse gerçekleştirim üst seviye assembly API'sini bozmadan
  !! eklenebilir.
  type, public :: stiffness_matrix_t
    private
    type(dense_matrix_t) :: storage
  end type stiffness_matrix_t

  public :: initialize_stiffness_matrix
  public :: validate_stiffness_matrix
  public :: add_local_stiffness
  public :: get_stiffness_matrix_value
  public :: get_stiffness_matrix_values
  public :: get_stiffness_matrix_size

contains

  !> Global torsional rijitlik matrisini sıfır olarak başlatır.
  !! Girdi negatif olmayan aktif DOF sayısıdır [-]. Çıktı n_dof x n_dof,
  !! [N*m/rad] birimli stiffness_matrix_t değeridir; sıfır DOF 0x0 üretir.
  pure subroutine initialize_stiffness_matrix(matrix, dof_count)
    type(stiffness_matrix_t), intent(out) :: matrix
    integer, intent(in) :: dof_count

    if (dof_count < 0) then
      error stop "Global rijitlik matrisi DOF sayısı negatif olamaz."
    end if

    call initialize_dense_matrix(matrix%storage, dof_count, dof_count)
  end subroutine initialize_stiffness_matrix

  !> Global torsional rijitlik matrisinin kare ve başlatılmış olduğunu sınar.
  pure subroutine validate_stiffness_matrix(matrix)
    type(stiffness_matrix_t), intent(in) :: matrix

    call validate_dense_matrix(matrix%storage)
    if (get_dense_matrix_row_count(matrix%storage) /= &
        get_dense_matrix_column_count(matrix%storage)) then
      error stop "Global rijitlik matrisi kare olmalıdır."
    end if
  end subroutine validate_stiffness_matrix

  !> İki uçlu elemanın lokal rijitlik katkısını global matrise toplar.
  !!
  !! Matematiksel model: Yerel a,b katsayısı, equation_ids(a),
  !! equation_ids(b) global konumuna eklenir. Sıfır denklem kimliği homojen
  !! kısıtlı DOF'tur ve ilgili satır/sütun indirgenmiş matris dışında bırakılır.
  !! Girdiler: [N*m/rad] local_matrix_2x2 ve boyutsuz iki denklem kimliği.
  !! Çıktı: Güncellenmiş [N*m/rad] global stiffness_matrix_t.
  !! Varsayımlar: Yerel sıra [node_i,node_j], global koordinat dönüşümü yoktur.
  pure subroutine add_local_stiffness(matrix, equation_ids, local_stiffness)
    type(stiffness_matrix_t), intent(inout) :: matrix
    integer, intent(in) :: equation_ids(2)
    type(local_matrix_2x2), intent(in) :: local_stiffness

    integer :: global_column
    integer :: global_row
    integer :: local_column
    integer :: local_row
    integer :: matrix_size

    call validate_stiffness_matrix(matrix)
    matrix_size = get_stiffness_matrix_size(matrix)

    if (any(equation_ids < 0) .or. any(equation_ids > matrix_size)) then
      error stop "Lokal rijitlik katkısının denklem kimliği geçersiz."
    end if

    do local_row = 1, 2
      global_row = equation_ids(local_row)
      if (global_row == 0) cycle

      do local_column = 1, 2
        global_column = equation_ids(local_column)
        if (global_column == 0) cycle

        call add_dense_matrix_entry( &
          matrix%storage, global_row, global_column, &
          local_stiffness%value(local_row, local_column))
      end do
    end do
  end subroutine add_local_stiffness

  !> Global torsional rijitlik matrisinin tek katsayısını [N*m/rad] döndürür.
  pure function get_stiffness_matrix_value(matrix, row, column) result(value)
    type(stiffness_matrix_t), intent(in) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: column
    real(dp) :: value

    call validate_stiffness_matrix(matrix)
    value = get_dense_matrix_entry(matrix%storage, row, column)
  end function get_stiffness_matrix_value

  !> Global torsional rijitlik katsayılarının [N*m/rad] bağımsız kopyasını
  !! döndürür.
  pure function get_stiffness_matrix_values(matrix) result(values)
    type(stiffness_matrix_t), intent(in) :: matrix
    real(dp), allocatable :: values(:, :)

    call validate_stiffness_matrix(matrix)
    values = get_dense_matrix_values(matrix%storage)
  end function get_stiffness_matrix_values

  !> Global rijitlik matrisinin boyutsuz satır/sütun sayısını döndürür.
  pure function get_stiffness_matrix_size(matrix) result(matrix_size)
    type(stiffness_matrix_t), intent(in) :: matrix
    integer :: matrix_size

    call validate_stiffness_matrix(matrix)
    matrix_size = get_dense_matrix_row_count(matrix%storage)
  end function get_stiffness_matrix_size

end module tms_stiffness_matrix
