# TTS Repeatability ve Cluster Bootstrap Matematiği

## Semboller ve birimler

Complete independent campaign sayısı `n`, bootstrap draw sayısı `B` ile
gösterilir. `s=log10(a_T)` ve adjacent `delta_s` boyutsuzdur. Arrhenius
`beta` birimi K, `Ea_app` birimi J/mol; WLF `C1` boyutsuz, `C2` K ve
`p=C1/C2`, `q=1/C2` birimi 1/K'dir.

Sampling unit complete campaign'dir. Aynı campaign içindeki frequency,
isotherm ve adjacent-pair noktaları ayrı istatistiksel replicate değildir.

## Sample mean, SD ve SE

Sonlu scalar değerler `x_i` için arithmetic mean:

```text
x_bar = (1/n) sum_i x_i
```

Sample standard deviation, population SD değil, `n-1` paydalı tahmindir:

```text
s_x = sqrt(sum_i (x_i-x_bar)^2 / (n-1))
```

Mean standard error:

```text
SE(x_bar) = s_x / sqrt(n)
```

Hesap Welford güncellemesiyle yapılır. `n=0` için bütün sonuç unavailable,
`n=1` için mean/median/min/max available fakat sample SD ve SE unavailable'dır.
Tek değeri SD=0 kabul ederek yapay certainty üretilmez. `n>=2` sabit veri için
SD=0 geçerli descriptive sonuçtur.

## Median ve MAD

Median, artan sıralı örneklemin orta sıra istatistiği; çift `n` için iki orta
değerin ortalamasıdır. Median absolute deviation:

```text
MAD = median_i |x_i - median(x)|
```

Normal dağılımla tutarlı robust scale diagnostic'i:

```text
scaled_MAD = MAD / 0.6744897501960817
```

`scaled_MAD`, ordinary standard deviation diye etiketlenmez. Ayrıca
`|mean-median|`, minimum ve maximum saklanır. Median/MAD yalnız robust
descriptive evidence'dır; authoritative TTS sonucunu değiştirmez, outlier
silmez ve downweight etmez. Signed veya sıfıra yakın `delta_s` için
coefficient of variation hesaplanmaz.

## Common-reference normalization

Campaign `r` için:

```text
s_common,r(T) = s_original,r(T) - s_original,r(T_common)
```

Adjacent fark:

```text
delta_s,ij = s(T_j) - s(T_i)
```

referans seçiminden bağımsızdır. Ortak referansta normalization tanımı gereği
`s_common,r(T_common)=0` olur. Bu **structural reference anchor**'dır; ölçüm
uncertainty'sinin sıfır olduğu anlamına gelmez.

V0.8.3, per-temperature marginal statistics verir. Absolute shift
temperature'ları arasındaki full covariance matrisi hesaplanmaz.

## Nonparametric complete-campaign cluster bootstrap

Eligible bağımsız population `C={c_1,...,c_R}` olsun. Her draw'da `R` campaign
kimliği replacement ile seçilir:

```text
C_b* = {c_i1, c_i2, ..., c_iR}
```

Tek bir draw planı aşağıdaki bütün scalar niceliklerde paylaşılır:

- adjacent `delta_s` mean,
- common-reference `s(T)` mean,
- Arrhenius `beta` ve `Ea_app` mean,
- WLF `C1`, `C2`, `p` ve `q` mean.

Bu ortak plan, örneğin campaign düzeyinde `B=10A` ise her available draw'da
`mean(B)=10 mean(A)` coupling'ini korur. Pair veya isotherm bazında bağımsız
resampling pseudoreplication olur ve yasaktır.

## Partial result availability

Arrhenius/WLF availability, draw population'ını değiştirmez. Complete campaign
draw'ı üretildikten sonra ilgili quantity için usable değerler seçilir. Cohort
mean uncertainty'si için her draw'da en az iki usable sampled value gerekir.
Daha azı varsa o draw unavailable olur.

Her interval:

- requested bootstrap draw count,
- valid draw count,
- unavailable draw count,
- confidence level,
- seed,
- lower/upper ve availability

alanlarını taşır. Unavailable değerler zero/NaN placeholder olarak quantile'a
eklenmez. Bir same-specimen-only study descriptive sonuç verebilir fakat
independent cluster bootstrap interval vermez.

## Deterministic RNG

Yerel Park–Miller minimal-standard recurrence kullanılır:

```text
z_(k+1) = 16807 z_k mod 2147483647
```

Schrage ayrıştırması çarpımdaki tanımsız signed integer overflow'u engeller.
Seed explicit'tir; default `20260803`, draw count `1000`, confidence level
`0.95` değeridir. Bunlar engineering acceptance threshold değildir. RNG
cryptographic amaç taşımaz ve global Fortran `random_number` state'ini
kullanmaz.

## Hyndman–Fan Type-7 quantile

Artan sıralı `m` valid bootstrap değeri `x_(1)...x_(m)` ve `0<=p<=1` için:

```text
h = 1 + (m-1)p
j = floor(h)
gamma = h-j
Q_7(p) = (1-gamma)x_(j) + gamma x_(j+1)
```

Uçlarda ilk/son sıra istatistiği kullanılır. Confidence level `CL` için:

```text
alpha = 1-CL
CI = [Q_7(alpha/2), Q_7(1-alpha/2)]
```

Bu interval nonparametric percentile bootstrap cohort-mean interval'ıdır.
BCa değildir ve engineering acceptance tolerance olarak yorumlanmaz.

## Parametric cohort sınırı

Arrhenius cohort yalnız valid `fit_available` sonuçlardan `beta` ve `Ea_app`
toplar. `Ea_app`, chemical aging activation energy değildir. WLF parameter
statistics yalnız `fit_available` ve `parameter_identifiable` sonuçlardan
oluşur. Poorly identified ve invalid campaign sayıları korunur; C1/C2/p/q
population'ına placeholder girmez. V0.8.2 `fit_tts_shift_laws()` authoritative
fit yordamıdır; V0.8.3 denklemleri yeniden çözmez.
