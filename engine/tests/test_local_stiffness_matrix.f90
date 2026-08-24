program test_local_stiffness_matrix
  use tms_kinds, only : dp
  use tms_local_matrix, only : local_matrix_2x2
  use tms_torsional_element, only : torsional_element_t, &
    calculate_local_stiffness
  implicit none

  real(dp), parameter :: stiffness = 100.0_dp
  real(dp), parameter :: tolerance = 1.0e-12_dp * stiffness

  type(torsional_element_t) :: element
  type(local_matrix_2x2) :: local_stiffness
  real(dp) :: expected(2, 2)
  real(dp) :: theta(2)
  real(dp) :: quadratic_form
  character(len=32) :: validation_case

  element = torsional_element_t( &
    id=1, &
    node_i_id=10, &
    node_j_id=20, &
    stiffness_nm_per_rad=stiffness, &
    damping_nms_per_rad=0.0_dp)

  ! Negatif rijitlik, üretim matris yordamı üzerinden ayrı WILL_FAIL testiyle
  ! reddedilir; diğer eleman alanları geçerli bırakılır.
  if (command_argument_count() > 0) then
    call get_command_argument(1, validation_case)
    call exercise_invalid_case(trim(validation_case), element)
    stop 0
  end if

  local_stiffness = calculate_local_stiffness(element)
  expected = reshape([ &
    100.0_dp, -100.0_dp, &
    -100.0_dp, 100.0_dp], [2, 2])

  ! Bilinen k = 100 N*m/rad değeri için üretim yordamının dört lokal matris
  ! katsayısı doğrudan analitik referansla karşılaştırılır.
  if (any(abs(local_stiffness%value - expected) > tolerance)) then
    error stop "Lokal torsional rijitlik matrisi analitik referansla uyuşmuyor."
  end if

  ! Ke(1,2) = Ke(2,1) simetrisi, lineer elastik bağlantının karşılıklılık ve
  ! konservatif enerji özelliğini doğrular.
  if (abs(local_stiffness%value(1, 2) - &
      local_stiffness%value(2, 1)) > tolerance) then
    error stop "Lokal torsional rijitlik matrisi simetrik değil."
  end if

  ! Her satırın sıfır toplamı, iki uca aynı açı eklendiğinde iç momentin
  ! değişmediğini ve ortak rijit-cisim dönmesinin serbest olduğunu doğrular.
  if (abs(sum(local_stiffness%value(1, :))) > tolerance .or. &
      abs(sum(local_stiffness%value(2, :))) > tolerance) then
    error stop "Lokal torsional rijitlik matrisinin satır toplamı sıfır değil."
  end if

  ! theta^T*Ke*theta >= 0 kontrolü, üretim formülünü testte yeniden kurmadan
  ! elemanın negatif elastik enerji üretemediğini sınar.
  theta = [0.30_dp, -0.20_dp]
  quadratic_form = dot_product( &
    theta, matmul(local_stiffness%value, theta))
  if (quadratic_form <= tolerance) then
    error stop "Bağıl dönme için lokal matris pozitif enerji üretmedi."
  end if

  ! Ortak dönme vektörü lokal matrisin rijit-cisim null modudur; bu durumda
  ! theta^T*Ke*theta sıfır olmalıdır.
  theta = [0.25_dp, 0.25_dp]
  quadratic_form = dot_product( &
    theta, matmul(local_stiffness%value, theta))
  if (abs(quadratic_form) > tolerance) then
    error stop "Ortak dönme, lokal torsional elemanda yapay enerji üretti."
  end if

  print *, "Lokal torsional eleman rijitlik matrisi doğrulandı."

contains

  !> Lokal rijitlik üretim yordamının geçersiz fiziksel girdiyi reddettiğini
  !! sınar.
  !! Girdi: Geçerli topoloji alanlarına sahip torsional eleman ve test vakası.
  !! Negatif k [N*m/rad], pasif elemanda negatif enerjiye yol açacağı için
  !! error stop üretmelidir. Çıktı yoktur.
  subroutine exercise_invalid_case(case_name, valid_element)
    character(len=*), intent(in) :: case_name
    type(torsional_element_t), intent(in) :: valid_element

    type(torsional_element_t) :: invalid_element
    type(local_matrix_2x2) :: invalid_matrix

    invalid_element = valid_element

    select case (case_name)
      case ("negative_stiffness")
        invalid_element%stiffness_nm_per_rad = -100.0_dp
      case default
        error stop "Bilinmeyen geçersiz lokal rijitlik testi istendi."
    end select

    invalid_matrix = calculate_local_stiffness(invalid_element)
    print *, "Negatif rijitlik beklenmedik biçimde kabul edildi: ", &
      invalid_matrix%value(1, 1)
  end subroutine exercise_invalid_case

end program test_local_stiffness_matrix
