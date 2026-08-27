module tms_loss_stiffness_matrix
  use tms_kinds, only : dp
  use tms_local_matrix, only : local_matrix_2x2
  use tms_matrix_types, only : dense_matrix_t, initialize_dense_matrix, &
    validate_dense_matrix, add_dense_matrix_entry, get_dense_matrix_entry, &
    get_dense_matrix_values, get_dense_matrix_row_count, &
    get_dense_matrix_column_count
  implicit none
  private

  !> Full veya reduced torsional denklem kümesinin kayıp rijitliği
  !! matrisini taşır.
  !! Fiziksel anlam: K'' harmonik çevrimde yapısal/elastomerik enerji
  !! kaybını temsil eder ve viskoz C matrisinden ayrıdır.
  !! Matematiksel anlam: Kompleks dinamik rijitliğin sanal kısmına
  !! doğrudan i*K'' katkısı verir. Katsayılar [N*m/rad] birimindedir.
  !! Depolama private tutularak fiziksel semantik dense ayrıntısından ayrılır.
  type, public :: loss_stiffness_matrix_t
    private
    type(dense_matrix_t) :: storage
  end type loss_stiffness_matrix_t

  public :: initialize_loss_stiffness_matrix
  public :: validate_loss_stiffness_matrix
  public :: add_local_loss_stiffness
  public :: get_loss_stiffness_matrix_value
  public :: get_loss_stiffness_matrix_values
  public :: get_loss_stiffness_matrix_size
  public :: extract_loss_stiffness_principal_submatrix

contains

  !> Global kayıp rijitliği matrisini sıfır olarak başlatır.
  !! Girdi negatif olmayan denklem sayısıdır [-]. Çıktı n x n,
  !! katsayıları [N*m/rad] olan loss_stiffness_matrix_t'dir. Sıfır DOF,
  !! tamamen kısıtlanmış sistem için geçerli 0x0 matris üretir.
  pure subroutine initialize_loss_stiffness_matrix(matrix, dof_count)
    type(loss_stiffness_matrix_t), intent(out) :: matrix
    integer, intent(in) :: dof_count

    if (dof_count < 0) then
      error stop "Global kayıp rijitliği matrisi DOF sayısı negatif olamaz."
    end if

    call initialize_dense_matrix(matrix%storage, dof_count, dof_count)
  end subroutine initialize_loss_stiffness_matrix

  !> Kayıp rijitliği matrisinin başlatılmış ve kare olduğunu doğrular.
  !! Girdi katsayıları [N*m/rad] olan matristir; geçersiz depolama error
  !! stop ile reddedilir. Eleman pasifliği assembly kaynağında doğrulanır.
  pure subroutine validate_loss_stiffness_matrix(matrix)
    type(loss_stiffness_matrix_t), intent(in) :: matrix

    call validate_dense_matrix(matrix%storage)
    if (get_dense_matrix_row_count(matrix%storage) /= &
        get_dense_matrix_column_count(matrix%storage)) then
      error stop "Global kayıp rijitliği matrisi kare olmalıdır."
    end if
  end subroutine validate_loss_stiffness_matrix

  !> İki uçlu elemanın lokal kayıp rijitliği katkısını toplar.
  !! Fiziksel açıklama: Elemanın harmonik enerji kaybı, bağlı olduğu
  !! torsional denklemlerin K'' katsayılarına eklenir.
  !! Matematiksel açıklama: Yerel (a,b) katsayısı
  !! K''(equation_ids(a),equation_ids(b)) konumuna scatter-add edilir.
  !! Girdiler [N*m/rad] local_matrix_2x2 ve iki boyutsuz denklem kimliğidir.
  !! Çıktı güncellenmiş loss_stiffness_matrix_t'dir. Sıfır denklem
  !! kimliği homojen kısıtlı DOF olarak atlanır.
  pure subroutine add_local_loss_stiffness( &
      matrix, equation_ids, local_loss_stiffness)
    type(loss_stiffness_matrix_t), intent(inout) :: matrix
    integer, intent(in) :: equation_ids(2)
    type(local_matrix_2x2), intent(in) :: local_loss_stiffness

    integer :: global_column
    integer :: global_row
    integer :: local_column
    integer :: local_row
    integer :: matrix_size

    call validate_loss_stiffness_matrix(matrix)
    matrix_size = get_loss_stiffness_matrix_size(matrix)

    if (any(equation_ids < 0) .or. any(equation_ids > matrix_size)) then
      error stop "Lokal kayıp rijitliği katkısının denklem kimliği geçersiz."
    end if

    do local_row = 1, 2
      global_row = equation_ids(local_row)
      if (global_row == 0) cycle

      do local_column = 1, 2
        global_column = equation_ids(local_column)
        if (global_column == 0) cycle

        call add_dense_matrix_entry( &
          matrix%storage, global_row, global_column, &
          local_loss_stiffness%value(local_row, local_column))
      end do
    end do
  end subroutine add_local_loss_stiffness

  !> Global K'' matrisinin tek katsayısını [N*m/rad] döndürür.
  pure function get_loss_stiffness_matrix_value( &
      matrix, row, column) result(value)
    type(loss_stiffness_matrix_t), intent(in) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: column
    real(dp) :: value

    call validate_loss_stiffness_matrix(matrix)
    value = get_dense_matrix_entry(matrix%storage, row, column)
  end function get_loss_stiffness_matrix_value

  !> Global K'' katsayılarının [N*m/rad] bağımsız kopyasını döndürür.
  !! Dönen dizi iç private depolamayı değiştiremez.
  pure function get_loss_stiffness_matrix_values(matrix) result(values)
    type(loss_stiffness_matrix_t), intent(in) :: matrix
    real(dp), allocatable :: values(:, :)

    call validate_loss_stiffness_matrix(matrix)
    values = get_dense_matrix_values(matrix%storage)
  end function get_loss_stiffness_matrix_values

  !> Global K'' matrisinin boyutsuz satır/sütun sayısını döndürür.
  pure function get_loss_stiffness_matrix_size(matrix) result(matrix_size)
    type(loss_stiffness_matrix_t), intent(in) :: matrix
    integer :: matrix_size

    call validate_loss_stiffness_matrix(matrix)
    matrix_size = get_dense_matrix_row_count(matrix%storage)
  end function get_loss_stiffness_matrix_size

  !> Seçilen full-equation indekslerinin principal K'' alt matrisini üretir.
  !! Fiziksel açıklama: Kısıtlı torsional DOF'ların kayıp rijitliği
  !! katsayıları active denklem takımından elenir.
  !! Matematiksel açıklama: B(a,b)=A(indices(a),indices(b)), yani
  !! K''_r=P^T*K''_full*P. Girdi K''_full [N*m/rad] ve benzersiz bir tabanlı
  !! indekslerdir [-]; çıktı aynı birimde reduced matristir.
  !! Boş indeks dizisi 0x0 matris üretir; sınır dışı veya yinelenen
  !! indeksler error stop ile reddedilir.
  pure function extract_loss_stiffness_principal_submatrix( &
      matrix, indices) result(submatrix)
    type(loss_stiffness_matrix_t), intent(in) :: matrix
    integer, intent(in) :: indices(:)
    type(loss_stiffness_matrix_t) :: submatrix

    integer :: column
    integer :: first_index
    integer :: matrix_size
    integer :: row
    integer :: second_index

    call validate_loss_stiffness_matrix(matrix)
    matrix_size = get_loss_stiffness_matrix_size(matrix)

    do first_index = 1, size(indices)
      if (indices(first_index) < 1 .or. &
          indices(first_index) > matrix_size) then
        error stop "Kayıp rijitliği alt matrisi indeksi sınır dışında."
      end if

      do second_index = first_index + 1, size(indices)
        if (indices(first_index) == indices(second_index)) then
          error stop "Kayıp rijitliği alt matrisi indeksleri benzersiz olmalıdır."
        end if
      end do
    end do

    call initialize_loss_stiffness_matrix(submatrix, size(indices))
    do row = 1, size(indices)
      do column = 1, size(indices)
        call add_dense_matrix_entry( &
          submatrix%storage, row, column, &
          get_dense_matrix_entry( &
          matrix%storage, indices(row), indices(column)))
      end do
    end do
  end function extract_loss_stiffness_principal_submatrix

end module tms_loss_stiffness_matrix
