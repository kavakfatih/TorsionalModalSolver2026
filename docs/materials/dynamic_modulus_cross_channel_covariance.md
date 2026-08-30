# Dynamic Modulus Cross-Channel Covariance

## Measurement semantics

V0.8.5 external covariance, aynı physical `(temperature, frequency)` DMA
noktasında linear modulus space'de tutulur:

```text
Var(G') [Pa²]    Cov(G',G'') [Pa²]
Cov(G',G'')      Var(G'')    [Pa²]
```

Bu quantity specimen-to-specimen variability, process variability, temperature
uncertainty veya TVD product covariance ile aynı değildir. Covariance explicit
measurement evidence veya magnitude/phase gibi açık bir measurement modelinden
gelmelidir. V0.8.3 scatter'ı, adjacent-frequency farkı, residual veya assumed
rubber `rho` otomatik covariance yapılmaz.

## Provenance

Her point aşağıdaki metadata'yı taşıyabilir:

- source kind: direct, repeat measurement, magnitude/phase propagation veya
  calibration model,
- measurement method,
- instrument identifier,
- calibration reference,
- source identifier ve açıklama.

Metadata'nın varlığı modeli doğrulamaz; yalnız hangi evidence'ın kullanıldığını
izlenebilir yapar. `covariance_available=false` gerçek support gap'tir.

## Standards position

[ASTM D5992](https://store.astm.org/d5992-96r18.html) ve
[ISO 4664-1:2022](https://www.iso.org/standard/80465.html) rubber dynamic
property measurement bağlamını verir. ASTM thermal-analysis kataloğundaki
[E2254, E2425, E1867 ve E3301](https://store.astm.org/products-services/standards-and-publications/standards/thermal-analysis-standards.html)
sırasıyla storage-modulus calibration, loss-modulus conformance ve DMA
temperature calibration provenance'ı için ilişkilidir. Bu standartların hiçbiri
TMS26 için universal `Cov(G',G'')` veya correlation coefficient sağlamaz; böyle
bir katsayı hard-code edilmez ve telifli standart metni yeniden üretilmez.

BOTTS çalışması broadband DMA verisinde signal-to-noise, linear error
propagation ve noise-weighted shifting'in pratik önemini gösterir; çalışma
TMS26 point-local covariance katsayısı vermez. Kaynak:
[Sheridan, Zauscher ve Brinson, Soft Matter 2024](https://pubs.rsc.org/en/content/articlelanding/2024/sm/d4sm00798k).

## Aynı nokta sınırı

V0.8.5 şunları kapsamaz:

- farklı frequency noktaları arasındaki covariance,
- farklı temperature/isotherm'ler arasındaki covariance,
- common instrument gain veya phase offset,
- geometry-factor ve fixture-compliance shared uncertainty,
- temperature ve frequency axes uncertainty.

Bu etkiler point-local 2x2 matrix'e gömülmemelidir. Özellikle temperature
uncertainty bağımsız ekseni ve TTS shift'i etkiler; `Sigma_G` içine eklenmez.

## DMA → ürün sınırı

```text
raw material / process
  -> measured G'(f,T), G''(f,T) and same-point covariance
  -> material TTS sensitivity
  -> dynamic torsional stiffness
  -> TVD response
  -> product test / CAE correlation
```

İlk covariance katmanı son ürün angle, torque, fatigue veya strength
covariance'ı değildir. Ürün seviyesine propagation; geometri, bonding,
manufacturing ve boundary-condition uncertainty kaynaklarını ayrıca gerektirir.
