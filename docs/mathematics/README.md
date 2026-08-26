# Matematik Dokümantasyonu

Bu dizin, TMS26 hesaplarında kullanılan denklemleri, sembolleri, SI birimlerini,
varsayımları ve geçerlilik sınırlarını içerir.

## Mevcut belgeler

- [`unit-conversions.md`](unit-conversions.md): temel mühendislik birimi
  dönüşümleri
- [`torsional-physics-core.md`](torsional-physics-core.md): annüler halka
  ataleti, elastomer rijitliği ve ilk doğal frekans zinciri
- [`torsional_frequency_model.md`](torsional_frequency_model.md): tek serbestlik
  dereceli torsional hareket denklemi ve sönümsüz doğal frekans
- [`dynamic_elastomer_model.md`](dynamic_elastomer_model.md): kompleks dinamik
  kayma modülü, kompleks burulma rijitliği bağlantısı ve kayıp faktörü
- [`two_inertia_modal_model.md`](two_inertia_modal_model.md): fixed-hub ve
  serbest-serbest iki ataletli TVD sisteminin analitik frekansları ile mod
  şekilleri
- [`generalized_torsional_system_model.md`](generalized_torsional_system_model.md):
  düğüm ataleti, eleman rijitlik/sönüm bağlantıları, aktif DOF eşlemesi ve
  global M/K formülasyonunun veri temeli
- [`local_torsional_element_matrix.md`](local_torsional_element_matrix.md):
  iki uçlu lineer torsional elemanın lokal rijitlik matrisi, enerji türetimi ve
  global assembly öncesi invariantları
- [`global_matrix_assembly.md`](global_matrix_assembly.md): fiziksel düğüm ile
  denklem kimliği ayrımı, lokal-global rijitlik toplama ve dönel atalet matrisi
- [`constraint_reduction.md`](constraint_reduction.md): tam M/K sisteminden
  aktif Kr/Mr sistemine direct elimination, seçim matrisi, prescribed değer ve
  result recovery bağıntıları
- [`generalized_modal_eigenproblem.md`](generalized_modal_eigenproblem.md):
  `K_r phi=lambda M_r phi` problemi, rigid-mode toleransı, mass normalization,
  relative residual, ortogonallik ve repeated eigenspace sözleşmesi

Fiziksel veya matematiksel bir yordam değiştirildiğinde denklem, girdi/çıktı
birimleri, varsayımlar ve test referansları aynı değişiklikte güncellenmelidir.
