# ADR 0014: Experimental Master-Curve Identification

- Durum: Kabul edildi
- Tarih: 2026-08-29

## Bağlam

V0.8.0, önceden doğrulanmış reference master curve ve temperature-shift
provider'ını runtime'da birleştirir. Measured DMA/dynamic-shear isotherm'lerden
bu iki authoritative runtime girdisini üretmek için ayrı, izlenebilir ve
failure-aware bir identification katmanı gerekir.

## Karar

1. Experimental identification ile runtime material katmanları ayrılır.
   Runtime provider ve harmonic API'leri fitting logic bilmez.
2. Public API explicit measured reference isotherm ister; automatic veya
   synthetic reference oluşturulmaz.
3. Canonical koordinatlar `x=log10(f)`, `s=log10(a_T)` ve `x_r=x+s` olarak
   V0.8.0 convention'ını korur.
4. Production strategy reference'tan iki yönde ilerleyen adjacent-pair
   shifting'dir. Global multi-isotherm optimizer kullanılmaz.
5. Objective, contiguous quality-valid segmentlerde exact piecewise-linear
   normalized L2 integralidir. Fixed grid sampling authoritative değildir.
6. Storage ve loss objective'leri ayrı normalize edilir; ikisi destekliyse
   eşit ağırlıklı joint objective kullanılır. Positive-loss support yoksa
   storage-only sonuç açıkça işaretlenir.
7. Feasible shift domain measured curve domain'lerinden türetilir. Default 65
   noktalı deterministic coarse scan yapılır.
8. Brent yalnız üç valid noktada gerçek interior minimum bracket'ı varsa
   çağrılır. Boundary minimum başarı sayılmaz.
9. `G''=0,VALID` physical/runtime için geçerli, log-loss objective için
   kullanılamaz; epsilon clamping yoktur.
10. Pair result overlap, objective, curvature, evaluation/iteration sayısı ve
    storage-loss shift discrepancy taşır. Universal TRS threshold yoktur.
11. Bütün shifted measured points ve provenance experimental master cloud'da
    korunur. Bu cloud research/validation truth'tur.
12. Solver runtime table ayrı, reference-centered ve strictly increasing bir
    distilled representation'dır. Averaging, smoothing ve spline uygulanmaz.
13. Duplicate reduced-frequency priority'si reference, reference'a daha yakın
    temperature ve daha uzak temperature sırasıdır.
14. VGP ve Cole-Cole outputs point cloud olarak sağlanır; tek scalar PASS/FAIL'e
    indirgenmez.
15. Empirical shift table mevcut V0.8.0 tabulated shift provider'ına, stitched
    table mevcut tabulated modulus provider'ına doğrudan aktarılır.

## V0.8.1'de Bilinçli Olarak Alınmayan Kararlar

- WLF `C1/C2` fitting ve Arrhenius `E_a` fitting yoktur.
- Uncertainty weighting, covariance, confidence interval ve bootstrap yoktur.
- Huber/robust objective veya outlier optimizer yoktur.
- Automatic reference selection yoktur.
- Vertical shift `b_T`, smoothing, PCHIP/spline ve global optimizer yoktur.
- Universal overlap, curvature, channel discrepancy veya TRS acceptance eşiği
  yoktur.

## Sonuçlar

Experimental failure bütün process'i `error stop` ile sonlandırmaz; status ve
availability döner. Bir adjacent link çözülmezse chain kırılır ve complete
runtime export üretilmez. Runtime/harmonic kodu stabil kalırken identification
algoritmaları bağımsız geliştirilebilir ve test edilebilir.

V0.8.2; parametric shift fitting, uncertainty-aware identification, robust
sensitivities, advanced non-adjacent consistency ve repeated-test statistics
için mevcut result/configuration extension noktalarını kullanacaktır.
