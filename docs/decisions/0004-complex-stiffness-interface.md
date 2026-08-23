# Karar 0004: Kompleks Rijitlik Arayüzü

- Durum: Kabul edildi
- Tarih: 2026-08-23

## Bağlam

V0.2.0 dinamik elastomer modeli G' ve G'' değerlerini saklar ancak bu değerleri
torsional solver'a bağlamaz. Yeni çözümün mevcut skaler rijitlik API'sini
bozmadan kompleks rijitlik bileşenlerini, kayıp faktörünü ve çalışma noktası
bilgisini açıkça taşıması gerekir.

## Karar

- Mevcut `calculate_torsional_stiffness` yordamı ve skaler sonucu korunacaktır.
- Dinamik hesap ayrı `tms_dynamic_torsional_stiffness` modülünde tutulacaktır.
- Sonuç, Fortran `complex` skaleri yerine alanları açıkça adlandırılmış
  `complex_torsional_stiffness_t` türüyle taşınacaktır.
- Annüler `Jp` hesabı `tms_geometry` içinde ortak bir saf yordam olacak ve iki
  rijitlik solver'ı tarafından kullanılacaktır.
- V0.2.1 dinamik hesabı `dynamic_rubber_material_t` içindeki geriye uyumlu tek
  çalışma noktası alanlarını kullanacak; `frequency_points` seçimi yapmayacaktır.

## Sonuçlar

- K' ve K'' bileşenlerinin birimleri ve fiziksel anlamları API'de görünürdür.
- Statik solver kullanan mevcut kodun çağrı biçimi ve sayısal davranışı değişmez.
- Geometrik moment formülü iki solver arasında yinelenmez.
- Çoklu çalışma noktası seçimi ve interpolasyon için ileride ayrı bir API kararı
  gerekir.
- Kompleks doğal frekans, eigen çözümü veya frekans cevabı bu arayüzün kapsamında
  değildir.
