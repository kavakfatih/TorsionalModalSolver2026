# TTS Experimental Data Model

## Family-Level Common State

Bir TTS family, aynı material/compound ve aynı physical test state altında
ölçülmüş isotherm'lerden oluşur. Authoritative common state yalnız
`tts_common_test_state_t` içinde bir kez tutulur:

- material/compound identifier,
- batch/material-state identifier,
- dynamic strain amplitude ratio `gamma_a` [-],
- static prestrain ratio `gamma_0` [-],
- deformation mode,
- conditioning/precycling description,
- test method ve source metadata.

Amplitude veya prestrain her isotherm'de tekrar edilmez. Family üyeliği
machine-epsilon karşılaştırmasıyla kurulmaz. V0.8.1 farklı batch'leri tek
family içinde sessizce birleştirmez. Specimen identifier ise isotherm'e özgü
olabilir.

## Isotherm

`tts_isotherm_t` şu alanları taşır:

- unique isotherm identifier,
- absolute temperature [K],
- specimen identifier,
- source identifier,
- authoritative measured points.

Frekanslar [Hz] sonlu, pozitif ve strictly increasing olmalıdır. Duplicate
veya unordered input reddedilir; silent sort yapılmaz. Reference isotherm
public API'de explicit measured kimlikle seçilir. Synthetic/interpolated
reference oluşturulmaz.

## Measurement Point ve Quality

Her `tts_measurement_point_t` frekans [Hz], `G'` [Pa], `G''` [Pa] ve iki ayrı
quality değeri taşır:

- `MEASUREMENT_VALID`,
- `BELOW_RELIABLE_FLOOR`,
- `MEASUREMENT_UNAVAILABLE`,
- `MEASUREMENT_REJECTED`.

Magic numeric quality kullanılmaz. `VALID` storage için `G'>0`; `VALID` loss
için `G''>=0` ve sonluluk zorunludur. Floor/unavailable/rejected point numeric
değer taşısa bile authoritative runtime table'a alınmaz.

## Zero Loss Semantiği

`G''=0, quality=VALID` passive ve fiziksel olarak geçerli olabilir:

- experimental input validation'dan geçer,
- runtime master table'a girebilir,
- VGP/Cole-Cole cloud'unda korunur,
- fakat `log10(G'')` tanımsız olduğundan loss shift objective'inde kullanılmaz.

Zero loss için epsilon, nearest value veya artificial floor eklenmez.

## Contiguous Segments

Storage objective yalnız contiguous `VALID, G'>0` koşularını; loss objective
yalnız contiguous `VALID, G''>0` koşularını kullanır. Aradaki tek bir invalid
nokta segmenti böler. Örneğin `VALID, VALID, FLOOR, VALID, VALID` dizisi iki
ayrı interpolation segmentidir; floor noktasının üzerinden interpolation
yapılmaz.

## Provenance ve Immutability

Experimental master cloud her original noktada isotherm, specimen, source,
temperature, physical frequency, applied shift, reduced frequency, iki modulus
ve iki quality değerini korur. Runtime table ayrı, solver-specific distilled
representation'dır. Caller inputunu identification sonrası değiştirmek sonuç
ve cloud'u değiştirmez.

## Kapsam Sınırı

V0.8.1 unit guessing, CSV/Excel importer, uncertainty weighting, repeated-test
statistics, amplitude/prestrain interpolation, vertical shift veya
self-heating içermez. Bu alanlar metadata'da açık state olarak taşınır fakat
fit edilmez.
