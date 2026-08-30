# TMS26 Doğrulama Belgeleri

Bu dizin, hesap motorunun üretim yordamlarını bağımsız analitik referanslar ve
matematiksel invariantlarla karşılaştıran doğrulama kanıtlarını içerir.

- [`torsional_validation.md`](torsional_validation.md): V0.3.1 torsional
  foundation için lokal eleman, DOF eşleme, global M/K assembly, matris kalite
  ölçütleri ve iki ataletli analitik benchmark doğrulaması
- [`numerical_robustness_validation.md`](numerical_robustness_validation.md):
  V0.3.2 için sonlu girdi sözleşmesi, toleranslı geometri ara yüzleri, uç ölçek
  numerik kararlılığı ve fixed-DOF indirgenmiş matris regresyonları
- [`modal_eigen_validation.md`](modal_eigen_validation.md): V0.5.0 generalized
  modal solver için analitik eigenvalue/frekans, mass normalization, residual,
  orthogonality, sign/repeated eigenspace ve physical recovery doğrulaması
- [`harmonic_response_validation.md`](harmonic_response_validation.md): V0.6.0
  K'/K''/C/M assembly, complex-symmetric ZSYSVX, status/diagnostics, 1-DOF ve
  iki-atalet analitik cevapları, excitation, FRF ve passivity doğrulaması
- [`dynamic_material_provider_validation.md`](dynamic_material_provider_validation.md):
  V0.7 provider data-quality/interpolation, G*→K* mapping, mixed/multiple
  material-aware harmonic zincir, prevalidation ve singular-point trace kanıtı
- [`thermorheological_runtime_validation.md`](thermorheological_runtime_validation.md):
  V0.8 tabulated/WLF/Arrhenius shift physics, log-space reduced frequency,
  dual-domain material evaluation, mevcut harmonic API entegrasyonu, trace ve
  V0.1–V0.7 regression kapıları
- [`V0.8.1_tts_validation.md`](V0.8.1_tts_validation.md): experimental quality,
  Brent safety, exact/non-TRS/weak synthetic evidence, stitching matrix,
  provenance ve V0.8.0 provider round-trip kapıları
- [`V0.8.2_shift_law_validation.md`](V0.8.2_shift_law_validation.md): adjacent
  pair Arrhenius/WLF recovery, reference invariance, large-C2 identifiability,
  LOTO, immutable empirical result ve parametric runtime round-trip kapıları
- [`V0.8.3_repeatability_validation.md`](V0.8.3_repeatability_validation.md):
  sample statistics, canonical campaign mapping, pseudoreplication prevention,
  cohort availability ve deterministic coupled cluster-bootstrap kapıları
- [`V0.8.4_uncertainty_validation.md`](V0.8.4_uncertainty_validation.md):
  uncertainty propagation/gap semantics, exact weighted/Huber integralleri,
  heteroscedastic/robust fixtures ve input immutability kapıları
- [`V0.8.5_covariance_validation.md`](V0.8.5_covariance_validation.md): polar/log
  covariance propagation, SPD/gap/interpolation, independent integral oracle,
  correlation yönü ve support/covariance decomposition kapıları

Buradaki çalışmalar **analytical code verification** niteliğindedir. Deneysel
ölçüm veya saha verisiyle fiziksel model korelasyonu ayrıca yürütülmelidir.
