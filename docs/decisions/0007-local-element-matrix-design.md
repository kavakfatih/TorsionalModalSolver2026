# Karar 0007: Lokal Eleman Matrisi Tasarımı

- Durum: Kabul edildi
- Tarih: 2026-08-24

## Bağlam

Karar 0006 genel torsional node-element topolojisini ve gelecekteki global M/K
formülasyonunun veri temelini oluşturdu. Global assembly katmanından önce her
lineer elemanın kendi fiziksel rijitlik katkısını bağımsız, saf ve doğrudan test
edilebilir biçimde üretmesi gerekir. Bu katkının topoloji koleksiyonuna veya
global DOF numaralandırmasına bağımlı olmaması istenir.

## Karar

- `tms_local_matrix` modülü, yalnız `tms_kinds` modülüne bağlı sabit boyutlu
  `local_matrix_2x2` veri taşıyıcısını tanımlayacaktır.
- Genel taşıyıcı kendi birimini kodlamaz. Fiziksel birim üretici yordamın
  sözleşmesinden gelir; torsional rijitlik sonucu için `N·m/rad` kullanılır.
- `tms_torsional_element`, matrisi kullanacak ve saf
  `calculate_local_stiffness` yordamını sağlayacaktır.
- Yerel koordinat sırası `[theta_i, theta_j]` olacak ve yordam
  `Ke = k[[1,-1],[-1,1]]` katkısını üretecektir.
- Mevcut eleman doğrulayıcısı çağrılarak sonlu ve pozitif `k` zorunlu tutulur.
- Bağımlılık yönü `tms_kinds -> tms_local_matrix -> tms_torsional_element ->
  tms_generalized_torsional_system` olacaktır.
- V0.2.4 global DOF eşlemesi, matrix assembly, M/C matrisi veya eigen çözümü
  sağlamayacaktır.

## Sonuçlar

- Lokal fizik, sistem boyutundan ve koleksiyon sırasından bağımsız test edilir.
- Matris simetrisi, sıfır satır toplamı, pozitif yarı-tanımlılık ve rijit-cisim
  null modu assembly katmanı eklenmeden doğrulanır.
- Sabit 2x2 tür, küçük eleman katkısının boyutunu açık tutar ve ham dizilerin
  API boyunca farklı anlamlarla kullanılmasını sınırlar.
- Gelecekteki global assembly, eleman uç kimliklerini DOF indekslerine eşlemek
  ve lokal katsayıları global K matrisine toplamakla sorumlu olacaktır.
- Damping alanı bu reel rijitlik matrisine katılmaz; K'' veya viskoz c için
  ayrı fiziksel model ve karar gerekir.
