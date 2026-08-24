module tms_matrix_types
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  implicit none
  private

  !> Genel amaçlı yoğun matris depolamasını temsil eder.
  !!
  !! Matematiksel anlam: value(i,j), i satırı ile j sütununun katsayısıdır.
  !! Fiziksel birim bu genel türde kodlanmaz; matrisi kullanan üst katmanın
  !! sözleşmesi tarafından belirlenir.
  !! Depolama private tutularak gelecekte aynı üst seviye API arkasında sparse
  !! bir gerçekleştirim kullanılabilmesine olanak sağlanır.
  type, public :: dense_matrix_t
    private
    real(dp), allocatable :: values(:, :)
  end type dense_matrix_t

  public :: initialize_dense_matrix
  public :: validate_dense_matrix
  public :: add_dense_matrix_entry
  public :: get_dense_matrix_entry
  public :: get_dense_matrix_values
  public :: get_dense_matrix_row_count
  public :: get_dense_matrix_column_count

contains

  !> Yoğun matrisi verilen satır ve sütun sayısıyla sıfır olarak başlatır.
  !! Girdiler boyutsuz ve negatif olmayan satır/sütun sayılarıdır. Çıktı,
  !! yeniden tahsis edilip sıfırlanmış dense_matrix_t değeridir. Sıfır boyut,
  !! tamamen kısıtlanmış bir sistemin 0x0 indirgenmiş matrisi için geçerlidir.
  pure subroutine initialize_dense_matrix(matrix, row_count, column_count)
    type(dense_matrix_t), intent(out) :: matrix
    integer, intent(in) :: row_count
    integer, intent(in) :: column_count

    if (row_count < 0 .or. column_count < 0) then
      error stop "Matris boyutları negatif olamaz."
    end if

    allocate(matrix%values(row_count, column_count))
    matrix%values = 0.0_dp
  end subroutine initialize_dense_matrix

  !> Yoğun matris depolamasının başlatıldığını doğrular.
  !! Girdi dense_matrix_t, çıktı yoktur; tahsis edilmemiş depolama error stop
  !! ile reddedilir. Fiziksel katsayıların birimi üst katmana aittir.
  pure subroutine validate_dense_matrix(matrix)
    type(dense_matrix_t), intent(in) :: matrix

    if (.not. allocated(matrix%values)) then
      error stop "Yoğun matris kullanılmadan önce başlatılmalıdır."
    end if
  end subroutine validate_dense_matrix

  !> Yoğun matrisin tek katsayısına sonlu bir katkı ekler.
  !! Matematiksel model: A(row,column) <- A(row,column) + contribution.
  !! Girdiler bir tabanlı boyutsuz indeksler ve üst katmanın SI birimindeki
  !! sonlu katkısıdır. Çıktı güncellenmiş matristir; boyut aşımı ve sonlu
  !! olmayan değerler error stop ile reddedilir.
  pure subroutine add_dense_matrix_entry(matrix, row, column, contribution)
    type(dense_matrix_t), intent(inout) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: column
    real(dp), intent(in) :: contribution

    real(dp) :: updated_value

    call validate_dense_matrix(matrix)
    call validate_dense_matrix_index(matrix, row, column)

    if (.not. ieee_is_finite(contribution)) then
      error stop "Matris katkısı sonlu olmalıdır."
    end if

    updated_value = matrix%values(row, column) + contribution
    if (.not. ieee_is_finite(updated_value)) then
      error stop "Matris katsayısı sonlu sayı aralığını aşamaz."
    end if

    matrix%values(row, column) = updated_value
  end subroutine add_dense_matrix_entry

  !> Yoğun matrisin bir katsayısını döndürür.
  !! Girdiler matris ile bir tabanlı boyutsuz satır/sütun indeksleridir.
  !! Çıktının fiziksel birimi üst katmanın matris sözleşmesinden gelir.
  !! Geçersiz indeks error stop ile reddedilir.
  pure function get_dense_matrix_entry(matrix, row, column) result(value)
    type(dense_matrix_t), intent(in) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: column
    real(dp) :: value

    call validate_dense_matrix(matrix)
    call validate_dense_matrix_index(matrix, row, column)
    value = matrix%values(row, column)
  end function get_dense_matrix_entry

  !> Yoğun matris katsayılarının bağımsız bir kopyasını döndürür.
  !! Girdi dense_matrix_t, çıktı aynı boyut ve üst katman SI birimine sahip
  !! allocatable gerçek dizidir. Dönen kopya iç depolamayı değiştiremez.
  pure function get_dense_matrix_values(matrix) result(values)
    type(dense_matrix_t), intent(in) :: matrix
    real(dp), allocatable :: values(:, :)

    call validate_dense_matrix(matrix)
    values = matrix%values
  end function get_dense_matrix_values

  !> Yoğun matrisin boyutsuz satır sayısını döndürür.
  pure function get_dense_matrix_row_count(matrix) result(row_count)
    type(dense_matrix_t), intent(in) :: matrix
    integer :: row_count

    call validate_dense_matrix(matrix)
    row_count = size(matrix%values, 1)
  end function get_dense_matrix_row_count

  !> Yoğun matrisin boyutsuz sütun sayısını döndürür.
  pure function get_dense_matrix_column_count(matrix) result(column_count)
    type(dense_matrix_t), intent(in) :: matrix
    integer :: column_count

    call validate_dense_matrix(matrix)
    column_count = size(matrix%values, 2)
  end function get_dense_matrix_column_count

  !> Bir tabanlı yoğun matris indeksinin depolama sınırlarında olduğunu sınar.
  !! Girdiler boyutsuz satır/sütun indeksleridir; çıktı yoktur. Bu özel yardımcı
  !! fiziksel hesap yapmaz ve sınır dışı erişimi error stop ile reddeder.
  pure subroutine validate_dense_matrix_index(matrix, row, column)
    type(dense_matrix_t), intent(in) :: matrix
    integer, intent(in) :: row
    integer, intent(in) :: column

    if (row < 1 .or. row > size(matrix%values, 1) .or. &
        column < 1 .or. column > size(matrix%values, 2)) then
      error stop "Matris satır veya sütun indeksi geçersiz."
    end if
  end subroutine validate_dense_matrix_index

end module tms_matrix_types
