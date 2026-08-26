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

## V0.5.0 modal solver belgeleri

- [`architecture/V0.5_modal_eigen_solver.md`](architecture/V0.5_modal_eigen_solver.md):
  reduced system, solver facade, DSYGV backend ve physical recovery katmanları
- [`mathematics/generalized_modal_eigenproblem.md`](mathematics/generalized_modal_eigenproblem.md):
  generalized eigenproblem, rigid-mode toleransı, normalizasyon ve residual
- [`validation/modal_eigen_validation.md`](validation/modal_eigen_validation.md):
  analitik modeller, repeated/sign-invariant testler ve kalite toleransları
- [`decisions/0010-generalized-eigen-solver-backend.md`](decisions/0010-generalized-eigen-solver-backend.md):
  dense DSYGV reference backend ve gelecek sparse/Lanczos-family facade kararı

Kod ile belgenin tutarlılığı aynı değişiklik içinde korunmalıdır. Fiziksel veya
matematiksel bir hesap değiştirildiğinde ilgili matematik ve fizik belgesi,
testler ve gerekiyorsa benchmark sonuçları birlikte güncellenir.
