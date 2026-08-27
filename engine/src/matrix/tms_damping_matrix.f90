module tms_damping_matrix
  use tms_kinds, only : dp
  use tms_local_matrix, only : local_matrix_2x2
  use tms_matrix_types, only : dense_matrix_t, initialize_dense_matrix, &
    validate_dense_matrix, add_dense_matrix_entry, get_dense_matrix_entry, &
    get_dense_matrix_values, get_dense_matrix_row_count, &
    get_dense_matrix_column_count
  implicit none
  private

  !> Full veya reduced torsional denklem kümesinin viskoz sönüm matrisini
  !! taşır.
  !! Fiziksel anlam: C, bağıl açısal hızla orantılı viskoz moment
  !! katsayılarını temsil eder ve K'' kayıp rijitliğinden ayrıdır.
  !! Matematiksel anlam: Kompleks dinamik rijitliğe i*omega*C katkısı
  !! verir. Katsayılar [N*m*s/rad] birimindedir. Depolama private tutulur.
  type, public :: damping_matrix_t
    private
    type(dense_matrix_t) :: storage
  end type damping_matrix_t

  public :: initialize_damping_matrix
  public :: validate_damping_matrix
  public :: add_local_damping
  public :: get_damping_matrix_value
  public :: get_damping_matrix_values
  public :: get_damping_matrix_size
  public :: extract_damping_principal_submatrix

contains

  !> Global viskoz sönüm matrisini sıfır olarak başlatır.
  !! Girdi negatif olmayan denklem sayısıdır [-]. Çıktı n x n,
  !! katsayıları [N*m*s/rad] olan damping_matrix_t'dir. Sıfır DOF,
  !! tamamen kısıtlanmış sistem için geçerli 0x0 matris üretir.
  pure subroutine initialize_damping_matrix(matrix, dof_count)
    type(damping_matrix_t), intent(out) :: matrix
    integer, intent(in) :: dof_count

    if (dof_count < 0) then
      error stop "Global viskoz sönüm matrisi DOF sayısı negatif olamaz."
    end if

    call initialize_dense_matrix(matrix%storage, dof_count, dof_count)
  end subroutine initialize_damping_matrix

  !> Viskoz sönüm matrisinin başlatılmış ve kare olduğunu doğrular.
  !! Girdi katsayıları [N*m*s/rad] olan matristir; geçersiz depolama error
  !! stop ile reddedilir. Eleman pasifliği assembly kaynağında doğrulanır.
  pure subroutine validate_damping_matrix(matrix)
    type(damping_matrix_t), intent(in) :: matrix

    call validate_dense_matrix(matrix%storage)
    if (get_dense_matrix_row_count(matrix%storage) /= &
        get_dense_matrix_column_count(matrix%storage)) then
      error stop "Global viskoz sönüm matrisi kare olmalıdır."
    end if
  end subroutine validate_damping_matrix

  !> İki uçlu elemanın lokal viskoz sönüm katkısını toplar.
  !! Fiziksel açıklama: Elemanın bağıl açısal hıza karşı moment
  !! katkısı, bağlı torsional denklemlerin C katsayılarına eklenir.
  !! Matematiksel açıklama: Yerel (a,b) katsayısı
  !! C(equation_ids(a),equation_ids(b)) konumuna scatter-add edilir.
  !! Girdiler [N*m*s/rad] local_matrix_2x2 ve iki boyutsuz denklem kimliğidir.
  !! Çıktı güncellenmiş damping_matrix_t'dir. Sıfır denklem kimliği
  !! homojen kısıtlı DOF olarak atlanır.
  pure subroutine add_local_damping(matrix, equation_ids, local_damping)
    type(damping_matrix_t), intent(inout) :: matrix
    integer, intent(in) :: equation_ids(2)
    type(local_matrix_2x2), intent(in) :: local_damping

    integer :: global_column
    integer :: global_row
    integer :: local_column
    integer :: local_row
    integer :: matrix_size

    call validate_damping_matrix(matrix)
    matrix_size = get_damping_matrix_size(matrix)

    if (any(equation_ids < 0) .or. any(equation_ids > matrix_size)) then
      error stop "Lokal viskoz sönüm katkısının denklem kimliği geçersiz."
    end if

    do local_row = 1, 2
      global_row = equation_ids(local_row)
      if (global_row == 0) cycle

      do local_column = 1, 2
        global_column = equation_ids(local_column)
        if (global_column == 0) cycle

        call add_dense_matrix_entry( &
          matrix%storage, global_row, global_column, &
          local_damping%value(local_row, local_column))
      end do
    end do
  end subroutine add_local_damping

  !> Global C matrisinin tek katsayısını [N*m*s/rad] döndürür.
  pure function get_damping_matrix_value(matrix, row, column) result(value)
    type(damping_matrix_t), intent(in) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: column
    real(dp) :: value

    call validate_damping_matrix(matrix)
    value = get_dense_matrix_entry(matrix%storage, row, column)
  end function get_damping_matrix_value

  !> Global C katsayılarının [N*m*s/rad] bağımsız kopyasını döndürür.
  !! Dönen dizi iç private depolamayı değiştiremez.
  pure function get_damping_matrix_values(matrix) result(values)
    type(damping_matrix_t), intent(in) :: matrix
    real(dp), allocatable :: values(:, :)

    call validate_damping_matrix(matrix)
    values = get_dense_matrix_values(matrix%storage)
  end function get_damping_matrix_values

  !> Global C matrisinin boyutsuz satır/sütun sayısını döndürür.
  pure function get_damping_matrix_size(matrix) result(matrix_size)
    type(damping_matrix_t), intent(in) :: matrix
    integer :: matrix_size

    call validate_damping_matrix(matrix)
    matrix_size = get_dense_matrix_row_count(matrix%storage)
  end function get_damping_matrix_size

  !> Seçilen full-equation indekslerinin principal C alt matrisini üretir.
  !! Fiziksel açıklama: Kısıtlı torsional DOF'ların viskoz sönüm
  !! katsayıları active denklem takımından elenir.
  !! Matematiksel açıklama: B(a,b)=A(indices(a),indices(b)), yani
  !! C_r=P^T*C_full*P. Girdi C_full [N*m*s/rad] ve benzersiz bir tabanlı
  !! indekslerdir [-]; çıktı aynı birimde reduced matristir.
  !! Boş indeks dizisi 0x0 matris üretir; sınır dışı veya yinelenen
  !! indeksler error stop ile reddedilir.
  pure function extract_damping_principal_submatrix(matrix, indices) &
      result(submatrix)
    type(damping_matrix_t), intent(in) :: matrix
    integer, intent(in) :: indices(:)
    type(damping_matrix_t) :: submatrix

    integer :: column
    integer :: first_index
    integer :: matrix_size
    integer :: row
    integer :: second_index

    call validate_damping_matrix(matrix)
    matrix_size = get_damping_matrix_size(matrix)

    do first_index = 1, size(indices)
      if (indices(first_index) < 1 .or. &
          indices(first_index) > matrix_size) then
        error stop "Viskoz sönüm alt matrisi indeksi sınır dışında."
      end if

      do second_index = first_index + 1, size(indices)
        if (indices(first_index) == indices(second_index)) then
          error stop "Viskoz sönüm alt matrisi indeksleri benzersiz olmalıdır."
        end if
      end do
    end do

    call initialize_damping_matrix(submatrix, size(indices))
    do row = 1, size(indices)
      do column = 1, size(indices)
        call add_dense_matrix_entry( &
          submatrix%storage, row, column, &
          get_dense_matrix_entry( &
          matrix%storage, indices(row), indices(column)))
      end do
    end do
  end function extract_damping_principal_submatrix

end module tms_damping_matrix
