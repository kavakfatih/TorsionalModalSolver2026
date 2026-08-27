# Dinamik Modül Interpolation Matematiği

## Primary değerler

V0.7 tablosunun primary büyüklükleri `G'(f)` ve `G''(f)` [Pa]'dır. Loss factor

\[
\tan\delta(f)=\frac{G''(f)}{G'(f)}
\]

interpolation sonrasında türetilir; bağımsız interpolate edilmez.

## LINEAR_FREQUENCY

Default policy `LINEAR_FREQUENCY`'dir. `f_1<f<f_2` için

\[
\alpha=\frac{f-f_1}{f_2-f_1},\qquad
G(f)=(1-\alpha)G_1+\alpha G_2
\]

bağıntısı ayrı ayrı `G'` ve `G''` üzerine uygulanır.

## LINEAR_LOG_FREQUENCY

Seçimlik `LINEAR_LOG_FREQUENCY` policy'sinde

\[
x=\log_{10}f,\qquad
\alpha=\frac{x-x_1}{x_2-x_1},\qquad
G(f)=(1-\alpha)G_1+\alpha G_2.
\]

Yalnız frequency axis logaritmiktir; modulus'un logaritması alınmaz. V0.7
log-log, spline, PCHIP, polynomial, smoothing veya regression fitting içermez.

## Exact point ve no-extrapolation

Stored frequency eşleşmesi bitwise equality ile değil,

\[
|a-b|\le 64\,\epsilon\,\max(|a|,|b|,\mathrm{tiny})
\]

ölçeğinde yalnız floating-point representation toleransıyla belirlenir. Bu
tolerans deney belirsizliği veya fiziksel frequency bandwidth değildir. Exact
eşleşmede stored `G'` ve `G''` aynen döner.

`f<f_min` veya `f>f_max` için endpoint hold, nearest value veya extrapolation
yoktur; provider domain hatası üretir.

## Measured-isotherm sıcaklık kuralı

Query operating temperature dataset temperature ile aynı machine-scale
tolerans içinde olmalıdır. Temperature interpolation, nearest isotherm,
ortalama, WLF, Arrhenius, reduced frequency ve TTS uygulanmaz. Bu tolerans
`±0.5 °C` gibi fiziksel test belirsizliği anlamına gelmez.

## Passivity

İki geçerli endpoint arasında `0<=alpha<=1` olduğundan convex interpolation
`G'>0` ve `G''>=0` koşullarını korur. Bonded-annular `C_theta>0` için

\[
K'=C_\theta G'>0,\qquad K''=C_\theta G''\geq0
\]

olur ve `K''/K'=G''/G'` eşitliği korunur.

## Causality sınırı

Parçalı interpolate edilmiş `G'(f),G''(f)` ampirik bir frequency-domain
gösterimdir. Bu yöntem kendi başına causal time-domain viscoelastic
constitutive law veya Kramers–Kronig consistency garanti etmez. V0.7 transient
malzeme modeli veya Prony/generalized-Maxwell fit'i değildir.
