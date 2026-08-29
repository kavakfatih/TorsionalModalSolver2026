# TTS Measurement-Uncertainty Weighted Objective

## Koordinatlar ve uncertainty propagation

TMS26 mevcut koordinatlarını korur:

```text
x = log10(f / Hz)
y = log10(G / Pa)
s = log10(a_T)
x_reduced = x + s
```

Pointwise modulus standard uncertainty `u_G [Pa]` için first-order sensitivity
coefficient sonucu:

```text
u_y = u_G / (G ln(10))
v_y = u_y^2
```

`u_y` boyutsuz log10-modulus standard uncertainty, `v_y` onun variance'ıdır.
Bu dönüşüm `G>0`, sonlu `u_G>0` gerektirir. Sıfır uncertainty infinite weight
olarak yorumlanmaz ve denominator'a epsilon eklenmez.

## Support segmentation ve interpolation

Measurement quality ile uncertainty availability kesişimi contiguous support
segmentlerini belirler. Missing uncertainty bir gap oluşturur. Her segment
içinde `x=log10(f)` ekseninde hem `y` hem de `v_y` endpoint değerlerini koruyan
piecewise-linear interpolation uygulanır. Variance interpolation fizik yasası
değil, deterministic TMS26 numerical policy'sidir. Gap bridge ve extrapolation
yoktur.

## Residual variance

Mevcut shift orientation ile pair residual:

```text
r(x;s) = y_i(x) - y_j(x-s)
```

V0.8.4 diagonal independence approximation'ı:

```text
v_r(x;s) = v_i(x) + v_j(x-s)
z(x;s)   = r(x;s) / sqrt(v_r(x;s))
```

`v_r` sonlu ve pozitif olmalıdır. Cross-channel, cross-isotherm ve common-mode
calibration covariance modellenmez.

## Channel objective

Gerçek ve gerekirse disjoint overlap support'u `O_k` için:

```text
J_W,k(s) = [sum integral_O r_k(x;s)^2 / v_r,k(x;s) dx] / |O_k|
|O_k|    = sum gerçek interval genişlikleri
```

Objective boyutsuzdur. Normalization point sayısı veya gap içeren convex hull
ile değil toplam valid overlap width [decade] ile yapılır. İki channel varsa:

```text
J_W = 0.5 J_W,storage + 0.5 J_W,loss
```

Bu nedenle daha çok sample taşıyan storage channel otomatik olarak daha büyük
ağırlık almaz. Positive loss support yoksa explicit storage-only sonuç mümkündür.

## Analitik interval integrali

Merged linear interval'de `t in [0,h]` için:

```text
r(t) = a + b t
v(t) = c + d t > 0
I_W  = integral_0^h (a+b t)^2/(c+d t) dt
```

Implementation endpoint'leri normalize edilmiş `q=t/h` koordinatına taşır ve
büyük variance endpoint'ini base seçer. Böylece `p=v_end/v_base-1` değeri
`[-1,0]` aralığında kalır. Gerekli analitik momentler:

```text
A_n(p) = integral_0^1 q^n/(1+p q) dq, n=0,1,2
I_W = h/v_base [r0^2 A0 + 2 r0 dr A1 + dr^2 A2]
```

Machine-small `p` için convergent power series, diğer durumda logarithmic
closed form/recurrence kullanılır. Bu seçim `d -> 0` limitindeki cancellation'ı
azaltır; arbitrary engineering threshold veya dense quadrature grid içermez.
Sabit variance limit'i:

```text
I_W = h/(3v) (r0^2 + r0 r1 + r1^2)
```

## Optimization ve invariants

Objective mevcut measured feasible domain, deterministic coarse scan, strict
interior bracket ve bounded Brent minimizer ile çözülür. Boundary-only ve flat
objective explicit failure'dır.

Bütün standard uncertainties ortak pozitif `lambda` ile çarpılırsa weighted
objective `1/lambda^2` ölçeklenir fakat minimizer değişmez. Bu invariant Huber
için zorlanmaz; çünkü standardized residual ve dolayısıyla Huber regime'i
`lambda` ile değişebilir.

Weighted shift, V0.8.1 authoritative unweighted shift'in replacement'ı değil
additive sensitivity evidence'ıdır. Düşük weighted residual tek başına
thermorheological simplicity kanıtı değildir.
