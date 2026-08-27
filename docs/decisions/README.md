# Mimari Karar Kayıtları

Bu dizin, TMS26 için önemli ve kalıcı teknik kararları numaralı kayıtlar halinde
tutar. Her kayıt bağlam, karar, sonuçlar, durum ve tarih bilgilerini içermelidir.

## Mevcut kararlar

- [`0001-internal-si-unit-system.md`](0001-internal-si-unit-system.md): çekirdek
  hesaplarda SI birim sistemi kullanımı
- [`0002-initial-torsional-physics-model.md`](0002-initial-torsional-physics-model.md):
  ilk lineer ve tek serbestlik dereceli torsional fizik kapsamı
- [`0003-dynamic-elastomer-data-model.md`](0003-dynamic-elastomer-data-model.md):
  dinamik modül çalışma noktaları ve geriye uyumlu malzeme genişletmesi
- [`0004-complex-stiffness-interface.md`](0004-complex-stiffness-interface.md):
  kompleks rijitlik sonuç türü, ortak geometri hesabı ve çalışma noktası seçimi
- [`0005-bonded-annular-rubber-torsion-model.md`](0005-bonded-annular-rubber-torsion-model.md):
  TVD elastomer rijitliği için tam bağlı eş merkezli silindirik burç modeli
- [`0006-generalized-torsional-topology.md`](0006-generalized-torsional-topology.md):
  genel düğüm-eleman topolojisi, mevcut iki-atalet API uyumluluğu ve K''/viskoz
  sönüm ayrımı
- [`0007-local-element-matrix-design.md`](0007-local-element-matrix-design.md):
  sabit boyutlu lokal matris taşıyıcısı, eleman bağımlılığı ve global assembly
  sınırı
- [`0008-global-matrix-assembly-design.md`](0008-global-matrix-assembly-design.md):
  fiziksel node/denklem ayrımı, dense global M/K assembly ve depolama soyutlaması
- [`0009-constraint-reduction-architecture.md`](0009-constraint-reduction-architecture.md):
  Physical DOF, constraint'ten bağımsız tam denklem numaralandırması, ayrı aktif
  harita, direct elimination ve result recovery mimarisi
- [`0010-generalized-eigen-solver-backend.md`](0010-generalized-eigen-solver-backend.md):
  LP64 LAPACK DSYGV dense reference backend'i, backend-neutral modal facade ve
  gelecekteki sparse/Lanczos-family backend genişleme sınırı
- [`0011-frequency-domain-complex-solver.md`](0011-frequency-domain-complex-solver.md):
  `exp(+i*omega*t)` harmonic convention, K'/K''/C ayrımı, complex-symmetric
  ZSYSVX reference backend'i ve status-aware frequency-sweep sözleşmesi
- [`0012-tabulated-dynamic-material-provider.md`](0012-tabulated-dynamic-material-provider.md):
  primary G'/G'', measured-isotherm interpolation, provider/binding sınırı,
  harmonic-only dynamic override, prevalidation ve traceability kararı

Mevcut karar değiştirilecekse eski kayıt silinmez; yeni bir karar kaydıyla
önceki kararın yerini aldığı belirtilir.
