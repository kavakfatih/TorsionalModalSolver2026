# ADR 0016: Repeatability ve Complete-Campaign Cluster Bootstrap

- Durum: Kabul edildi
- Tarih: 2026-08-29

## Bağlam

V0.8.1 complete DMA/TTS campaign başına authoritative empirical shift ve
master curve; V0.8.2 aynı evidence'dan additive Arrhenius/WLF fit üretir.
Repeated experiment sonuçlarının scatter'ı değerlendirilirken tek campaign
içindeki frequency, isotherm veya adjacent-pair noktalarını independent sample
saymak yapay confidence ve pseudoreplication üretir. Same-specimen rerun'ları
da specimen/process bağımsızlığını temsil etmeyebilir.

## Karar

1. Statistical sampling unit complete independent DMA/TTS campaign'dir.
2. Frequency point, isotherm ve adjacent pair replicate sayılmaz;
   pseudoreplication yasaktır.
3. `INDEPENDENT_SPECIMEN_CAMPAIGN`, `SAME_SPECIMEN_RERUN` ve unspecified basis
   explicit taşınır. Same-specimen rerun default descriptive evidence'a
   girebilir fakat independent bootstrap population'ına girmez.
4. Campaign temperature ve pair sonuçları array index ile değil canonical
   physical Kelvin key ve normalized pair orientation ile eşlenir.
5. Farklı original references,
   `s_common(T)=s_original(T)-s_original(T_common)` ile ortak measured
   reference'a taşınır; input mutate edilmez.
6. Common reference'taki sıfır structural anchor'dır, measurement uncertainty
   kanıtı değildir. Result bu semantic'i explicit taşır.
7. Mean, `n-1` sample SD, SE, median, MAD, scaled MAD, min/max ve
   `|mean-median|` descriptive statistics olarak raporlanır. Signed shift için
   coefficient of variation yoktur.
8. Uncertainty yöntemi nonparametric complete-campaign cluster bootstrap'tir.
   Sampling replacement ile yapılır.
9. Tek campaign-index draw planı bütün adjacent/absolute/Arrhenius/WLF
   niceliklerinde paylaşılır; campaign-level dependence korunur.
10. Derived result availability population draw'ından sonra uygulanır. Bir
    quantity/draw içinde iki usable value yoksa draw unavailable olur; zero,
    NaN veya `huge` placeholder kullanılmaz.
11. RNG explicit seed'li portable Park–Miller/Schrage yordamıdır. Global
    `random_number` state'i kullanılmaz.
12. Primary target cohort mean, interval convention Hyndman–Fan Type-7
    percentile CI'dır. BCa uygulanmaz.
13. V0.8.1 empirical result ve V0.8.2 parameter fit authoritative kalır;
    V0.8.3 yalnız additive evidence üretir.
14. Huber/uncertainty-weighted shifting, automatic outlier deletion ve
    automatic WLF/Arrhenius winner V0.8.3'e alınmaz.
15. Universal material/engineering pass-fail threshold tanımlanmaz.

## Sonuçlar

Study sonucu descriptive ve independent-bootstrap availability'sini ayırır;
total/independent/rerun campaign sayıları ile gerçek bootstrap population
kimlikleri izlenebilir kalır. Partial WLF/Arrhenius availability ortak cluster
planını bozmaz. Deterministic RNG ve quantile convention platformlar arası
regresyonu mümkün kılar.

Bu karar repeatability'yi reproducibility, accuracy, bias veya TRS validity ile
eşitlemez. Bootstrap confidence interval engineering acceptance tolerance
değildir. DMA evidence bonded TVD product validation yerine geçmez.

## Kapsam dışında

- full covariance, generalized least squares, BCa ve Bayesian uncertainty,
- pointwise uncertainty-weighted veya Huber V0.8.1 pair identification,
- outlier deletion, trimming ve automatic model selection,
- vertical shift, smoothing, Prony/Maxwell production fit,
- interlaboratory reproducibility ve product-level uncertainty propagation.

Bu başlıklar V0.8.4 veya daha sonraki ayrı karar kayıtları gerektirir.
