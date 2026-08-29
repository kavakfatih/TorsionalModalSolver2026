# Dynamic Modulus Measurement Uncertainty

## Canonical quantity

V0.8.4 dış girdisi pointwise **standard uncertainty** `u_G [Pa]` değeridir.
Storage ve loss modulus için availability ve provenance ayrı taşınır. Expanded
uncertainty `U`, coverage factor bilinmeden standard uncertainty değildir;
relative accuracy/tolerance specification da otomatik `u_G` sayılmaz. Gerekli
dönüşüm gelecekte importer sorumluluğudur, core ambiguous conversion yapmaz.

Core içinde `y=log10(G)` için first-order propagation:

```text
u_log10_G = u_G/(G ln(10))
```

olarak uygulanır. JCGM 100/GUM ve NIST TN 1297'deki sensitivity-coefficient
yaklaşımı terminoloji çerçevesidir; TMS26 bu sürüm için bir standarda uygunluk
veya accreditation iddiasında bulunmaz.

## Provenance ve fiziksel anahtar

Uncertainty overlay, authoritative V0.8.1 measurement family'yi mutate etmez.
Kayıtlar array index yerine unique physical `(temperature K, frequency Hz)`
anahtarıyla eşlenir. Source kind ve serbest metadata; Type A, Type B, combined
standard veya repeat-measurement kökenini kaydedebilir.

V0.8.3 independent-specimen campaign scatter'ı otomatik pointwise measurement
uncertainty'ye çevrilmez. Campaign scatter aynı anda measurement, specimen,
material ve process variation içerebilir. Measurement repeatability de tek
başına material/process variability değildir.

## Quality ve support

Weighted storage point için `quality=VALID`, `G'>0`, sonlu `u_G'>0`; weighted
loss point için `quality=VALID`, `G''>0`, sonlu `u_G''>0` gerekir. `u_G=0`,
negative veya nonfinite değer reddedilir; epsilon yoktur. Missing uncertainty
support gap'tir. `G''=0` passive runtime değeri olabilir fakat log-loss
objective'e girmez.

## Covariance sınırı

V0.8.4 residual variance'ı diagonal toplamla kurar. Şunlar modellenmez:

- `G'`/`G''` cross-channel covariance,
- compared isotherm'ler arası covariance,
- shared amplitude/phase ve common-mode calibration covariance.

JCGM 102 multiple-output propagation ve covariance için, JCGM 101 ise
linearization yetersiz olduğunda Monte Carlo için gelecek bağlamdır. Bunlar
V0.8.5'e bırakılmıştır.

## Test-method ve literatür konumu

ASTM D5992 ve ISO 4664-1 forced-vibration dynamic property bağlamını; ASTM
E2254, E2425, E1867 ve E3301 dynamic mechanical analysis/verification
terminolojisini destekleyen bağlam kaynaklarıdır. Kullanılan fixture'lar bu
standartlara conformance iddiası değildir ve copyrighted prosedür metni burada
yeniden üretilmez.

BOTTS (Soft Matter, 2024, DOI `10.1039/D4SM00798K`) güncel TTS fitting bağlamı
sağlar; ancak farklı fitting/interpolation formulation kullanır ve TMS26
piecewise-linear weighted objective'inin birebir doğrulaması değildir.
Measurement-uncertainty weighting genel WLS/metrology ilkeleri ve güncel TTS
literatürüyle motive edilir; TMS26 kendi deterministic formulation'ını korur.

## DMA ile ürün validation ayrımı

DMA uncertainty, bonded TVD product uncertainty değildir. Angle/torque ölçüm
uncertainty'si, torsional stiffness correlation uncertainty'si ve torsional
fatigue/strength uncertainty'si ayrı test zincirleridir. İyi weighted TTS
agreement, crank-pulley ürün validation'ının yerine geçmez.

## Bağlam kaynakları

- [BIPM/JCGM Guides in Metrology](https://www.bipm.org/en/web/guest/publications/guides):
  JCGM 100, JCGM 101 ve JCGM 102 yayın ailesi
- [NIST Technical Note 1297](https://www.nist.gov/pml/nist-technical-note-1297):
  standard/combined/expanded uncertainty ve first-order propagation bağlamı
- [NIST Engineering Statistics Handbook — Weighted Least Squares](https://itl.nist.gov/div898/handbook/pmd/section1/pmd143.htm):
  değişken precision altında inverse-variance weighting bağlamı
- [ASTM D5992](https://store.astm.org/standards/d5992): vulkanize rubber için
  vibratory dynamic testing kapsamı
- [ISO 4664-1:2022](https://www.iso.org/standard/80465.html): vulkanize veya
  termoplastik rubber dynamic-property genel guidance kapsamı
- [ASTM thermal-analysis standards](https://store.astm.org/products-services/standards-and-publications/standards/thermal-analysis-standards.html):
  E2254, E2425, E1867 ve E3301 DMA calibration/conformance bağlamı
- [BOTTS, Soft Matter 2024](https://pubs.rsc.org/en/content/articlelanding/2024/sm/d4sm00798k):
  broadband TTS ve noise-weighted fitting literatür bağlamı
