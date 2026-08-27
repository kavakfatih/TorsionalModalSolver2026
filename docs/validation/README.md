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

Buradaki çalışmalar **analytical code verification** niteliğindedir. Deneysel
ölçüm veya saha verisiyle fiziksel model korelasyonu ayrıca yürütülmelidir.
