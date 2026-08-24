# Karar 0008: Global Matris Assembly Tasarımı

- Durum: Kabul edildi
- Tarih: 2026-08-24

## Bağlam

Karar 0006 fiziksel node-element topolojisini, Karar 0007 ise elemanın bağımsız
2x2 lokal rijitlik katkısını tanımladı. Global modal denkleme geçmeden önce
fiziksel düğüm kimliklerinin matematiksel denklem numaralarından ayrılması ve
lokal katkıların doğrulanabilir global M/K matrislerine toplanması gerekir.

Fiziksel kimlikler kullanıcı ve model topolojisi için kararlı etiketlerdir;
matris indeksleri ise kısıtlara ve solver politikasına göre değişebilir. İki
kavramın aynı tamsayı olarak kullanılması sınır koşulu ve yeniden numaralandırma
hatalarına yol açar.

## Karar

- `dof_map_t`, tüm fiziksel node ID'lerini tutacak; serbest düğümlere eklenme
  sırasıyla kesintisiz `1..n_active`, kısıtlı düğümlere `0` equation ID verecek.
- Eksik node araması hata olacak; `0` yalnız kısıtlı DOF anlamına gelecektir.
- Harita ile sistem node kümesi ve kısıt durumu assembly öncesi çapraz
  doğrulanacaktır.
- Tamamen kısıtlı sistem 0 aktif DOF ve 0x0 M/K matrisleriyle desteklenir.
- `dense_matrix_t` ilk depolama gerçekleştirimidir. Allocatable katsayı dizisi
  private tutulacak; boyut ve değer erişimi yordamlar üzerinden yapılacaktır.
- Global `stiffness_matrix_t` ve `mass_matrix_t`, dense ayrıntısını private
  tutan fiziksel üst seviye türler olacaktır. Bu sınır gelecekte sparse
  depolama eklenmesini mümkün kılar.
- Global K assembly, `get_local_stiffness` sonucunu açık
  `node_id -> equation_id -> scatter-add` zinciriyle toplayacaktır.
- Global M assembly, aktif düğümlerin `J [kg·m²]` değerlerini diagonal olarak
  toplayacaktır.
- Denklem kimliği sıfır olan satır ve sütunlar atlanarak homojen sıfır dönme
  kısıtının indirgenmiş matrisi oluşturulacaktır.
- Mevcut `calculate_local_stiffness` API'si korunacak; yeni
  `get_local_stiffness` sarmalayıcısı standart eleman katkı arayüzünü kuracaktır.
- Global C, K'', viskoz damping, eleman ataleti, sparse matris, LAPACK, eigen
  solver, nonlinear solver ve FEM mesh V0.3.0 kapsamına alınmayacaktır.

## Sonuçlar

- Fiziksel node ID artık hiçbir global matris yordamında doğrudan indeks değildir.
- Fiziksel kısıt bilgisi düğümde korunur; bunun `0/1..n` denklem gösterimi ve
  numaralandırma politikası ayrı, test edilebilir DOF haritasında türetilir.
- Lokal eleman fiziği sistem boyutundan bağımsız kalırken assembly yalnız
  denklem eşleme ve katsayı toplama sorumluluğunu taşır.
- Serbest global K için simetri, sıfır satır toplamı ve rijit-cisim null modu;
  global M için diagonal polar atalet yapısı analitik testlerle korunur.
- Dense depolama performans ve sadelik için yeterlidir; üst seviye API private
  depolama sayesinde gelecekteki sparse gerçekleştirimle değiştirilebilir.
- Sıfırdan farklı prescribed dönmeler ileride RHS düzeltmesi gerektirir.
- Modal reduction bu tam-uzay DOF haritasını bozmak yerine ayrı bir baz dönüşüm
  katmanı olarak eklenmelidir.
