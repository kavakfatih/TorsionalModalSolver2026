# TMS26 Dokümantasyonu

Bu dizin, TMS26 hesap motorunun teknik tasarımını, matematiksel temelini,
fiziksel model kabullerini, analitik doğrulama kanıtlarını, geliştirme
süreçlerini ve kalıcı kararlarını içerir.

## Dizinler

- [`architecture/`](architecture/): modül sınırları ve veri akışı
- [`mathematics/`](mathematics/): denklemler, birimler ve sayısal varsayımlar
- [`physics/`](physics/): fiziksel modellerin anlamı ve geçerlilik kapsamı
- [`materials/`](materials/): constitutive malzeme verisi ve deney izlenebilirliği
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

## V0.6.0 frequency-domain response belgeleri

- [`architecture/V0.6_frequency_domain_response.md`](architecture/V0.6_frequency_domain_response.md):
  full/reduced dynamic matrisler, complex solver facade, sweep ve recovery
- [`mathematics/harmonic_torsional_response.md`](mathematics/harmonic_torsional_response.md):
  harmonic equation, residual, FRF ve dissipated energy matematiği
- [`physics/torsional_damping_models.md`](physics/torsional_damping_models.md):
  depolama rijitliği, kayıp rijitliği ve viskoz sönümün fiziksel ayrımı
- [`validation/harmonic_response_validation.md`](validation/harmonic_response_validation.md):
  analitik, matris, solver-status, excitation ve passivity doğrulama kapsamı
- [`decisions/0011-frequency-domain-complex-solver.md`](decisions/0011-frequency-domain-complex-solver.md):
  complex-symmetric LP64 ZSYSVX reference backend kararı

## V0.7.0 dynamic material provider belgeleri

- [`architecture/V0.7_dynamic_material_provider.md`](architecture/V0.7_dynamic_material_provider.md):
  provider, binding, preparation, material-aware assembly ve result trace
- [`materials/tabulated_dynamic_elastomer_material.md`](materials/tabulated_dynamic_elastomer_material.md):
  canonical SI dataset, metadata, operating state ve DMA→TVD transfer sınırı
- [`mathematics/dynamic_modulus_interpolation.md`](mathematics/dynamic_modulus_interpolation.md):
  iki interpolation policy, machine tolerance, isotherm ve no-extrapolation
- [`validation/dynamic_material_provider_validation.md`](validation/dynamic_material_provider_validation.md):
  analytical provider/mapping/harmonic/trace doğrulama kapıları
- [`decisions/0012-tabulated-dynamic-material-provider.md`](decisions/0012-tabulated-dynamic-material-provider.md):
  harmonic-only tabulated provider ve V0.8 genişleme sınırı

## V0.8.0 thermorheological runtime belgeleri

- [`architecture/V0.8_thermorheological_runtime.md`](architecture/V0.8_thermorheological_runtime.md):
  shift provider bileşimi, physical/reduced coordinate ayrımı, mevcut binding
  ve harmonic solver reuse mimarisi
- [`materials/thermorheological_dynamic_elastomer.md`](materials/thermorheological_dynamic_elastomer.md):
  reference master curve, operating-state metadata, horizontal TTS ve DMA→TVD
  aktarım sınırları
- [`mathematics/temperature_shift_functions.md`](mathematics/temperature_shift_functions.md):
  canonical `a_T`, WLF, Arrhenius, tabulated log-shift ve log-space reduced
  frequency denklemleri
- [`validation/thermorheological_runtime_validation.md`](validation/thermorheological_runtime_validation.md):
  shift physics, dual-domain, material, harmonic, trace ve regression kapıları
- [`decisions/0013-thermorheological-runtime-convention.md`](decisions/0013-thermorheological-runtime-convention.md):
  `a_T=tau(T)/tau(T_ref)`, `f_r=a_Tf` ve horizontal-only runtime kararı
- [`development/V0.8.0_development_report.md`](development/V0.8.0_development_report.md):
  kalıcı architecture/technical/validation/CI geliştirme raporu

Kod ile belgenin tutarlılığı aynı değişiklik içinde korunmalıdır. Fiziksel veya
matematiksel bir hesap değiştirildiğinde ilgili matematik ve fizik belgesi,
testler ve gerekiyorsa benchmark sonuçları birlikte güncellenir.
