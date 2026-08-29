# ADR 0015: Parametric Shift-Law Identification

- Durum: Kabul edildi
- Tarih: 2026-08-29

## Bağlam

V0.8.1 measured DMA isotherm'lerinden authoritative adjacent relative shifts,
reference-anchored empirical shift table ve solver-ready master curve üretir.
CAE runtime'da daha kompakt WLF veya Arrhenius law kullanmak istenebilir;
ancak parametric model experimental evidence'ı gizlememeli veya yeniden
hizalamamalıdır.

## Karar

1. Parametric identification V0.8.1 sonucunu tüketen additive katmandır.
2. Primary fit observations adjacent `delta_s=s_j-s_i` sonuçlarıdır;
   cumulative absolute shifts independent observations sayılmaz.
3. Default weights equal'dır. Pair curvature, overlap ve residual uncertainty
   veya inverse variance olarak yorumlanmaz.
4. Arrhenius `Ea_app` analytical one-parameter least-squares ile çözülür;
   nonlinear optimizer kullanılmaz.
5. WLF fixed C2 için C1 analytical profillenir; yalnız pole-safe C2 mevcut
   generic Brent minimizer ile çözülür.
6. WLF large-C2 linear limitinde `p=C1/C2` ve `q=1/C2` diagnostics saklanır.
   Interior identifiable C1/C2 yoksa result poorly identified olur.
7. Fit, residual validation, predictive validation ve WLF parameter
   identifiability ayrı availability alanlarıdır.
8. LOTO, shared-isotherm leakage'i azaltmak için primary predictive
   diagnostic'tir; minimum data yoksa mathematical fit reddedilmez.
9. Automatic WLF/Arrhenius best-model kararı üretilmez. Empirical table
   authoritative kalır.
10. Parametric runtime export mevcut V0.8.0 providers'ı reuse eder ve yalnız
    measured calibrated temperature domain'inde geçerlidir.
11. Poorly identified WLF runtime'a export edilmez.

## Sonuçlar

Empirical ve parametric evidence birlikte izlenebilir kalır. Runtime material
ve harmonic katmanları fitting logic import etmez. Reference-temperature
reparameterization physics'i değiştirmez. Parameter uncertainty, covariance,
robust fitting ve automatic model selection gelecekte ayrı karar gerektirir.

## Kapsam dışında

- global multi-isotherm optimizer ve graph refinement,
- covariance/uncertainty weighting, bootstrap ve robust objective,
- vertical shift, smoothing, Prony veya production Maxwell fit,
- temperature extrapolation ve automatic Tg/law selection.
