# Experimental Repeatability ve Uncertainty

## Mühendislik bağlamı

V0.8.3, aynı laboratuvar/test sistemi/prosedürü ve kontrol edilen yakın
measurement conditions altında yürütülen complete DMA/TTS campaign'lerinin
**intralaboratory repeatability** evidence'ını üretir. Repeatability precision
ile ilgilidir; reproducibility, accuracy veya bias sonucu değildir. Düşük
scatter, measurement'ın doğru veya tarafsız olduğunu kanıtlamaz.

Tek-laboratuvar campaign'lerinden interlaboratory reproducibility iddiası
üretilmez. Aynı specimen rerun'ları thermal/preconditioning/fixture history
nedeniyle independent specimen campaign sayılmaz.

## Fiziksel nedensellik zinciri

```text
raw material / compound state
  -> rubber process / cure state
  -> repeated DMA dynamic measurements
  -> V0.8.1 empirical TTS
  -> V0.8.2 WLF / Arrhenius identification
  -> V0.8.3 repeatability / uncertainty
  -> G'(f,T), G''(f,T)
  -> K'(f,T), K''(f,T)
  -> TVD harmonic response
  -> product dynamic stiffness / angle-torque / fatigue correlation
  -> DOE / sensitivity
  -> material and process tolerances
  -> raw-material specification / process window
  -> production validation
```

V0.8.3 yalnız DMA/master-curve evidence katmanındadır. DMA repeatability,
bonded TVD product validation yerine geçmez; product-level correlation ayrıca
geliştirilecektir.

## Campaign ve provenance

Her campaign unique kimlik ve açık replicate basis taşır. Traceability için
laboratory, operator, instrument, test protocol, calibration reference, run,
test-date metadata'sı ile V0.8.1 specimen/source provenance'ı korunur.

Physics-critical common state şunları kapsar:

- material ve batch-state identifier,
- deformation mode,
- dynamic strain amplitude ve static prestrain,
- conditioning description ve test method.

Specimen/source/campaign kimliklerinin farklı olması tek başına incompatible
state değildir. Calibration metadata yazılım acceptance threshold'u değil,
traceability evidence'ıdır.

## Bağımsızlık ve descriptive policy

Complete independent specimen campaign statistical sample unit'tir. Aynı
campaign'deki çok sayıdaki frequency point veya isotherm `n` değerini artırmaz.
Default policy aynı-specimen rerun'ları descriptive statistics'e dahil eder,
fakat independent cluster bootstrap population'ından çıkarır. Tüm provenance
raporda kalır; hiçbir rerun gizlenmez veya silinmez.

Median, MAD, scaled MAD ve mean–median farkı robust descriptive diagnostics'tir.
Bunlar automatic outlier rejection, trimming, Huber weighting veya
authoritative empirical shift değişikliği yapmaz.

## Shift ve model yorum sınırları

- Measurement repeatability, thermorheological simplicity (TRS) validity
  değildir.
- Bootstrap confidence interval, engineering acceptance tolerance değildir.
- Common reference'taki yapısal sıfır, zero measurement uncertainty değildir.
- `Ea_app` scatter, chemical-aging activation-energy uncertainty değildir.
- WLF parameter scatter, WLF'nin doğru physical model olduğunun kanıtı
  değildir.
- DMA master-curve evidence, bonded TVD ürün doğrulaması değildir.

Arrhenius/WLF automatic winner seçilmez. Empirical V0.8.1 table reference model
olarak kalır ve V0.8.2 fit status/identifiability semantics'i korunur.

## Standard ve literatür konumu

Aşağıdaki yayınlar terminology, test provenance ve statistical-method context
için referanstır; telifli yöntem metni kopyalanmaz ve conformance iddiası
üretilmez:

- ASTM D5992-96(2024), vulcanized rubber dynamic vibratory testing,
- ISO 4664-1:2022, rubber dynamic-properties general guidance,
- ASTM E177, precision ve bias terminology,
- ASTM E691, interlaboratory study concepts için yalnız reference context,
- ASTM E2254, E2425, E1867 ve uygulanabildiğinde E3301, DMA calibration/
  conformance provenance,
- Sheridan, Zauscher ve Brinson, BOTTS, *Soft Matter* (2024),
- TTS package literature,
- NIST Engineering Statistics Handbook, bootstrap ve MAD methodology.

TMS26 V0.8.3, ASTM E691 compliant bir interlaboratory study implementation'ı
olduğunu iddia etmez. Listedeki standardlar software pass/fail limitleri
değildir.

## Kapsam dışı

Measurement uncertainty'den weighted V0.8.1 pair objective, Huber/Tukey
robust shifting, full covariance, BCa/Bayesian inference, vertical shift,
Payne/amplitude/prestrain normalization, self-heating, smoothing, Prony veya
product-level response uncertainty bu sürümde yoktur.
