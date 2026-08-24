module tms_constraint_types
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use tms_kinds, only : dp
  use tms_dof_types, only : is_supported_dof_type
  implicit none
  private

  !> Dönel serbestlik derecesini sıfır değerde sabitleyen constraint türü [-].
  integer, parameter, public :: FIXED_CONSTRAINT = 1

  !> Dönel serbestlik derecesine sonlu bir açı değeri atayan constraint türü [-].
  integer, parameter, public :: PRESCRIBED_VALUE_CONSTRAINT = 2

  !> Tek bir fiziksel serbestlik derecesine uygulanan kinematik constraint'i
  !! temsil eder. V0.4.0 yalnız torsional dönme için fixed ve prescribed değer
  !! koşullarını destekler; kayıt yapısı gelecekte yeni DOF türlerine açıktır.
  type, public :: constraint_t
    !> Constraint kaydının benzersiz pozitif kimliği [-].
    integer :: constraint_id = 0

    !> Constraint'in uygulandığı fiziksel düğümün pozitif kimliği [-].
    integer :: node_id = 0

    !> Kısıtlanan fiziksel serbestlik derecesinin anlamlı tür kimliği [-].
    integer :: dof_type = 0

    !> Prescribed torsional dönme değeri theta_bar [rad]. Matematiksel olarak
    !! q = P*q_r + q_bar ilişkisindeki kısıtlı koordinat bileşenidir. Modal
    !! pertürbasyon vektörüne eklenmez.
    real(dp) :: value = 0.0_dp

    !> Fixed veya prescribed değer davranışını seçen constraint türü [-].
    integer :: constraint_type = 0
  end type constraint_t

  public :: validate_constraint

contains

  !> Constraint veri kaydının topolojik, fiziksel ve sayısal önkoşullarını
  !! doğrular.
  !!
  !! Fiziksel açıklama: V0.4.0'da yalnız ortak eksen çevresindeki torsional
  !! dönme kısıtlanabilir. Fixed koşul theta=0 rad, prescribed koşul ise sonlu
  !! theta=theta_bar rad değerini temsil eder.
  !! Matematiksel açıklama: Kimlikler pozitiftir, DOF türü
  !! TORSIONAL_ROTATION'dır, value sonludur ve constraint türü desteklenen
  !! kümededir. FIXED_CONSTRAINT için value tam olarak sıfır olmalıdır.
  !! Girdi: constraint_t. Çıktı yoktur; geçersiz kayıt error stop ile
  !! reddedilir. Bu yordam matris indirgeme veya yük vektörü oluşturmaz.
  pure subroutine validate_constraint(constraint)
    type(constraint_t), intent(in) :: constraint

    if (constraint%constraint_id <= 0) then
      error stop "Constraint kimliği pozitif olmalıdır."
    end if

    if (constraint%node_id <= 0) then
      error stop "Constraint düğüm kimliği pozitif olmalıdır."
    end if

    if (.not. is_supported_dof_type(constraint%dof_type)) then
      error stop "Constraint serbestlik derecesi türü desteklenmiyor."
    end if

    if (.not. ieee_is_finite(constraint%value)) then
      error stop "Constraint prescribed değeri sonlu olmalıdır."
    end if

    select case (constraint%constraint_type)
      case (FIXED_CONSTRAINT)
        if (abs(constraint%value) > 0.0_dp) then
          error stop "Fixed constraint değeri sıfır radyan olmalıdır."
        end if
      case (PRESCRIBED_VALUE_CONSTRAINT)
        continue
      case default
        error stop "Constraint türü desteklenmiyor."
    end select
  end subroutine validate_constraint

end module tms_constraint_types
