# TTS Shift Identification Mathematics

## Coordinate and Sign Convention

V0.8.0 convention değişmeden korunur:

```text
a_T(T) = tau(T) / tau(T_ref)
s(T)   = log10(a_T)
x      = log10(f / 1 Hz)
x_r    = x + s
f_r    = a_T f
```

İki isotherm arasındaki relative shift `delta_s`, moving curve'e
`x_j,shifted=x_j+delta_s` olarak uygulanır. Absolute shift chain için
`s_j=s_i+delta_s` kullanılır. Reference ölçümde `s_ref=0`, `a_T=1` tamdır.

## Feasible Shift Domain

Measured log-frequency domain'leri `i:[L_i,U_i]` ve `j:[L_j,U_j]` ise gerçek
overlap için açık domain:

```text
L_i - U_j < delta_s < U_i - L_j
```

Arama sınırları doğrudan measured domain'den gelir. `[-20,20]` gibi arbitrary
physical bounds kullanılmaz. Açık uçlar yalnız machine-scale iç margin ile
scan edilir.

## Channel Coordinates

Storage ve positive-loss ordinatları:

```text
y'(x)  = log10(G'(x) / 1 Pa)
y''(x) = log10(G''(x) / 1 Pa)
```

Storage için `VALID,G'>0`; loss için `VALID,G''>0` gerekir. `VALID,G''=0`
log-objective'e alınmaz ve epsilon ile değiştirilmez. Her quality gap ayrı
contiguous piecewise-linear segment üretir.

## Exact Piecewise-Linear L2 Residual

İki valid segmentin shifted overlap knot'ları birleştirilir. Bir
`[x_0,x_1]` subinterval'ında residual lineer ise:

```text
h  = x_1 - x_0
r0 = r(x_0)
r1 = r(x_1)

Integral r(x)^2 dx = h/3 (r0^2 + r0 r1 + r1^2)
```

Bu exact integral authoritative implementation'dır; fixed grid sampling
kullanılmaz. Tek nokta overlap gerçek interpolation interval'i olmadığı için
mathematical support sayılmaz.

Channel objective kendi toplam valid overlap measure'ıyla normalize edilir:

```text
J_c = sum(Integral r_c^2 dx) / sum(overlap width)
```

Her iki channel support taşıyorsa production objective:

```text
J_joint = 0.5 (J_storage + J_loss)
```

G' support yeterli fakat positive-valid G'' interval support yoksa production
shift storage-only olarak açıkça işaretlenir. Bu sonuç full complex-modulus TRS
kanıtı değildir.

## Overlap Diagnostics

`overlap_width_decades`, valid interval genişliklerinin toplamıdır.
`overlap_fraction`, bu toplamın iki curve'ün toplam valid segment measure'ından
küçük olana oranıdır. Joint sonuçta iki channel değerinin küçüğü raporlanır.
Universal minimum-overlap threshold uygulanmaz.

## Coarse Scan ve Brent Bracket

Feasible interval default 65 deterministic noktada taranır. Her scan point'in
objective availability'si değerden ayrı tutulur. Brent yalnız ardışık üç valid
noktada aşağıdaki interior minimum varsa çağrılır:

```text
s[k-1] < s[k] < s[k+1]
J[k-1] > J[k] < J[k+1]
```

Floating-point plateau gürültüsünün sahte bracket üretmemesi için objective
drop, `64 epsilon max(1,|J|)` numerical ölçeğini de aşmalıdır. Bu eşik
experimental confidence veya TRS acceptance eşiği değildir. Interior bracket
yoksa boundary minimum başarı kabul edilmez.

Generic Brent implementation safeguarded parabolic interpolation ve
golden-section adımlarını kullanır. Default absolute/relative stopping ölçeği:

```text
8 sqrt(epsilon(dp))
```

Binary64 için yaklaşık `1.2e-7` boyutsuz shift mertebesidir. Machine precision'e
kadar gereksiz iterasyon yapılmaz.

## Identifiability Curvature

Optimum yakınında finite-difference evidence:

```text
H ≈ [J(s+h) - 2J(s) + J(s-h)] / h^2
```

Low curvature weak identifiability göstergesidir; V0.8.1 universal `H`
threshold veya PASS/FAIL üretmez. Ayrı storage/loss optimumları bulunduğunda:

```text
D_channel = |delta_s_storage - delta_s_loss|
```

değeri ayrıca tutulur.

## Shift Chain

Isotherm'ler yalnız index listesinde temperature-sorted edilir. Explicit
reference'tan colder ve hotter yönlere adjacent ilerlenir. Bir zorunlu link
çözülemezse `CHAIN_BROKEN` döner; non-adjacent bridge veya global optimizer
yoktur.

## Runtime Stitching

Shifted experimental cloud research/validation truth'tur. Solver table önce
reference'in valid complex noktalarını, sonra reference'a temperature olarak
en yakın curve'leri alır. Yalnız mevcut reduced-frequency range'ini gerçekten
genişleten noktalar eklenir. Overlap averaging, smoothing ve spline yoktur.
Machine-equivalent duplicate priority'si reference, daha yakın temperature ve
daha uzak temperature sırasıdır.
