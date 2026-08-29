# ADR 0017: Measurement-Uncertainty Weighted ve Robust TTS Sensitivity

- Durum: Kabul edildi
- Tarih: 2026-08-30

## Bağlam

V0.8.1 measured isotherm pair'lerini unweighted exact piecewise-linear L2 ile
shift eder ve empirical sonuçları authoritative kabul eder. V0.8.3 complete
campaign repeatability/bootstrap scatter'ı üretir; bu scatter pointwise
measurement uncertainty ile aynı nicelik değildir. Frekansa bağlı ölçüm
kalitesinin shift sensitivity üzerindeki etkisini, baseline'ı değiştirmeden ve
grid sampling hatası eklemeden incelemek gerekir.

## Karar

1. V0.8.1 empirical unweighted-L2 pair shift authoritative baseline kalır.
2. V0.8.4 sonucu additive sensitivity evidence'dır; baseline, empirical shift,
   master cloud veya runtime table otomatik değiştirilmez.
3. Canonical external uncertainty quantity pointwise standard uncertainty
   `u_G [Pa]` değeridir. Expanded uncertainty veya relative specification
   ambiguous biçimde dönüştürülmez.
4. `y=log10(G)` için first-order propagation
   `u_log10_G=u_G/(G ln(10))` kullanılır.
5. `u_G<=0` ve nonfinite değer reddedilir; epsilon veya infinite weight yoktur.
6. Uncertainty overlay, source ölçüme unique physical `(T,f)` anahtarıyla
   eşlenir. Array index matching yapılmaz.
7. Missing uncertainty contiguous support'u böler. Gap bridge ve extrapolation
   yoktur.
8. Valid segment içinde `v_y=u_log10_G^2`, log-frequency ekseninde
   piecewise-linear interpolate edilir. Bu deterministic numerical policy'dir.
9. Pair residual variance diagonal toplamdır: `v_r=v_i+v_j`. Cross-channel,
   cross-isotherm ve common-mode covariance V0.8.4'te modellenmez.
10. Weighted objective, her merged interval'de `r^2/v` analitik integrali ve
    gerçek total overlap-width normalization kullanır. Arbitrary dense
    resampling yoktur.
11. Joint objective, channel'lar ayrı normalize edildikten sonra eşit
    `0.5/0.5` ağırlık kullanır; point count ağırlığı yoktur.
12. Optimization mevcut deterministic coarse scan, strict interior bracket ve
    bounded Brent policy'sini korur. Boundary-only success yoktur.
13. Robust objective raw residual'a değil standardized `z=r/sqrt(v)` değerine
    Huber loss uygular.
14. Default `c=1.345` configurable'dır ve rubber acceptance limiti değildir.
15. Huber regime crossing'leri quadratic denklemle bulunur; quadratic ve tail
    bölgeleri analitik, grid-free entegre edilir.
16. Huber hiçbir measurement point'i silmez. Tail point invalid measurement
    sayılmaz ve automatic pass/fail uygulanmaz.
17. Weighted ve Huber shift'ler authoritative baseline'ın yerine otomatik
    geçirilmez; yalnız signed sensitivity delta'ları raporlanır.
18. V0.8.3 independent-specimen scatter'ı otomatik `u_G` yapılmaz.
19. Covariance-aware measurement uncertainty, Mahalanobis/GLS ve multiple-output
    propagation V0.8.5 kapsamına ayrılır.

## Sonuçlar

Pointwise measurement uncertainty kaynak/provenance ile izlenir; heteroscedastic
ve robust influence deterministic fixture'larla görülebilir. Disjoint support
measure doğru normalize edilir ve varying variance altında standardized
residual'ın yanlışlıkla lineer varsayılması engellenir. Bunun karşılığında
covariance içermeyen sonuç conditional sensitivity evidence olarak yorumlanır.

## Kapsam dışında

- full covariance, generalized least squares ve Mahalanobis objective,
- JCGM 101 Monte Carlo ve JCGM 102 multiple-output implementation,
- automatic outlier deletion veya bad-point rejection,
- baseline/master/runtime replacement,
- vertical shift, smoothing, spline/Prony veya global multi-isotherm optimizer,
- bonded TVD angle/torque/fatigue uncertainty propagation ve product validation.
