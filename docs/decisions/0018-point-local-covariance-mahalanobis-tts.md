# ADR 0018: Point-Local Covariance ve Mahalanobis TTS Sensitivity

- Durum: Kabul edildi
- Tarih: 2026-08-31

## Bağlam

V0.8.4 storage ve loss standard uncertainty'lerini diagonal residual variance
olarak kullanır. Aynı DMA point'inde iki modulus aynı amplitude/phase ve
instrument chain'inden türeyebileceği için cross-channel covariance sıfır
olmayabilir. Covariance effect'i support daralmasıyla karıştırmadan incelemek ve
V0.8.1 authoritative result'ını korumak gerekir.

## Karar

1. Canonical external covariance linear physical G-space 2x2 matrix ve birimi
   Pa²'dir; log covariance dış input yapılmaz.
2. Log covariance `Sigma_y=D Sigma_G D^T` first-order Jacobian propagation ile
   üretilir ve provenance flag'i taşır.
3. Covariance yalnız aynı physical `(temperature, frequency)` noktasında
   `G'↔G''` bağımlılığını modeller.
4. Cross-isotherm ve cross-frequency covariance modellenmez; V0.8.5 full GLS
   olarak adlandırılmaz.
5. Covariance matrix pozitif diyagonal, pozitif determinant ve machine-safe
   reciprocal condition ile SPD olmalıdır. Singular/near-singular açık status
   üretir.
6. Automatic jitter, epsilon-I, eigenvalue clipping, shrinkage ve pseudoinverse
   uygulanmaz.
7. Missing covariance support gap'idir; gap bridge ve extrapolation yoktur.
8. Valid segmentte full matrix elemanları `log10(f)` üzerinde lineer
   interpolate edilir; rho veya SD ayrı interpolate edilmez.
9. Mahalanobis yalnız storage, positive loss, iki uncertainty ve covariance
   kesişimi olan bivariate common support `O_B` üzerinde çalışır.
10. Loss unavailable ise storage-only Mahalanobis fallback yoktur; V0.8.4
    storage-only sonucu ayrı evidence olarak kalabilir.
11. `J_M=(2|O_B|)^-1 integral r^T Sigma_r^-1 r dx` analitik/grid-free
    cubic-over-quadratic interval integraliyle değerlendirilir.
12. Aynı `O_B` üzerinde off-diagonal sıfırlanmış matched-support diagonal
    control zorunludur.
13. `delta_support`, `delta_covariance` ve `delta_total` ayrı raporlanır;
    threshold veya automatic preferred shift yoktur.
14. V0.8.1 empirical baseline authoritative ve immutable kalır. V0.8.4
    original weighted result independent evidence olarak yeniden kullanılır.
15. API hem V0.8.4 uncertainty hem V0.8.5 covariance alırsa matrix
    diyagonalleri `u_G²` ile machine-equivalent olmalıdır; uyuşmazlık sessizce
    çözülmez.
16. Covariance-aware Huber V0.8.5'e eklenmez; önce covariance physics bağımsız
    doğrulanır.
17. V0.8.6 common-mode/correlated-field uncertainty, cross-frequency,
    cross-isotherm, temperature-axis ve low-rank nuisance-source modellerine
    ayrılır.

## Sonuçlar

Support restriction ile covariance ellipse geometry'sinin shift etkileri ayrı
gözlenebilir. Diagonal limit mevcut scalar weighted integral ile uyumludur ve
runtime solver bağımlılığı oluşmaz. Bunun karşılığında input covariance'nın
explicit measurement evidence/provenance ile sağlanması gerekir ve point-local
model shared systematic effects'i açıklamaz.

## Reddedilen seçenekler

- Universal/hard-coded rubber correlation coefficient,
- V0.8.3 scatter veya V0.8.4 residual'dan automatic covariance estimation,
- singular matrix'i regularize ederek sessizce kabul,
- covariance support yokken storage-only Mahalanobis,
- arbitrary dense resampling,
- complete DMA spectrum için dense full GLS,
- Mahalanobis sonucu ile empirical master/runtime tablolarını overwrite etmek.
