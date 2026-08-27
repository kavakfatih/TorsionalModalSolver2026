# Dynamic Material Provider Doğrulaması

## Kapsam

V0.7 doğrulaması data/provider, annular property mapping, harmonic coupling,
traceability ve V0.1–V0.6 regression kapılarını kapsar. Test toleransları deney
belirsizliği değil IEEE double-precision ve bağımsız analitik karşılaştırma
amaçlıdır.

## Gate A — Data / Provider

- exact stored point: iki interpolation policy için stored G'/G'' aynen,
- one-ULP frequency ve temperature: machine-equivalent exact match,
- `LINEAR_FREQUENCY`: known midpoint `alpha=0.5`,
- `LINEAR_LOG_FREQUENCY`: 10–100 Hz log midpoint `sqrt(1000)` Hz,
- no extrapolation: alt/üst domain iki policy için reddedilir,
- en az iki point, strictly increasing f, finite/positive G', finite/nonnegative
  G'', geçerli isotherm ve metadata,
- negative G'' için clipping olmadan failure,
- different point temperature ve physically different query temperature için
  failure,
- provider input/getter kopyası mutation'ına karşı immutability.

## Gate B — Physical Property Mapping

Bağımsız annular geometri hesabıyla

\[
C_\theta=\frac{4\pi Lr_i^2r_o^2}{r_o^2-r_i^2},\quad
K'=C_\theta G',\quad K''=C_\theta G''
\]

doğrulanır. `K''/K'=G''/G'`, passivity ve nonnegative dissipated energy
kontrol edilir. TENSILE mode direct torsional binding olarak reddedilir.

## Gate C — Harmonic Coupling

Fixed-hub 1-DOF referansı:

\[
\hat\theta(f)=\frac{\hat T}
{K'(f)-\omega^2J+i[K''(f)+\omega c]}
\]

ile interpolate→map→assemble→ZSYSVX zinciri karşılaştırılır. Stored nominal
K'/K'' deliberately farklı seçilerek material-aware override ve no-double-
counting kanıtlanır. Aynı sistemin V0.6 frozen analizi nominal değerlerle
ayrıca doğrulanır; c'nin iki yolda da ayrı kaldığı sınanır.

İki free DOF'lu seri zincir, bir dynamic + bir constant eleman ve iki ayrı
dynamic provider için bağımsız complex 2×2 inverse ile doğrulanır.

## Gate D — Traceability

Solved noktada dataset ID, C_theta, G'/G'', K'/K'', tan(delta), policy,
bracket ve alpha solver state ile eşleşir. Exact resonance ile singular Z
üretilen noktada harmonic response unavailable olurken material trace mevcut
kalır.

## Hata ve prevalidation testleri

Duplicate binding, unknown element, invalid annular geometry, uninitialized
provider, incompatible mode, operating temperature mismatch, empty binding ve
provider kapsamı dışına çıkan full sweep ayrı expected-failure süreçlerinde
reddedilir. Full sweep evaluation solver loop'undan önce yapıldığı için domain
hatası partial solve başlatmaz.

## Regression

Mevcut DSYGV modal ve frozen ZSYSVX harmonic testleri değiştirilmeden CTest
içinde çalışır. V0.7 yeni external numerical dependency eklemez ve aynı LP64
LAPACK backend'lerini kullanır.
