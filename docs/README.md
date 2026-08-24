# TMS26 Dokümantasyonu

Bu dizin, TMS26 hesap motorunun teknik tasarımını, matematiksel temelini,
fiziksel model kabullerini, analitik doğrulama kanıtlarını, geliştirme
süreçlerini ve kalıcı kararlarını içerir.

## Dizinler

- [`architecture/`](architecture/): modül sınırları ve veri akışı
- [`mathematics/`](mathematics/): denklemler, birimler ve sayısal varsayımlar
- [`physics/`](physics/): fiziksel modellerin anlamı ve geçerlilik kapsamı
- [`validation/`](validation/): analitik referanslar ve doğrulama sonuçları
- [`development/`](development/): build, test ve katkı süreçleri
- [`decisions/`](decisions/): mimari karar kayıtları

Kod ile belgenin tutarlılığı aynı değişiklik içinde korunmalıdır. Fiziksel veya
matematiksel bir hesap değiştirildiğinde ilgili matematik ve fizik belgesi,
testler ve gerekiyorsa benchmark sonuçları birlikte güncellenir.
