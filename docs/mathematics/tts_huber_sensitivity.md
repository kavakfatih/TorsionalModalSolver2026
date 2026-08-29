# Standardized Huber TTS Sensitivity

## Standardized residual ve loss

Huber loss raw log-modulus residual'a değil measurement uncertainty ile
standardize edilmiş

```text
z = r / sqrt(v_r)
```

değerine uygulanır:

```text
rho_c(z) = 0.5 z^2                 , |z| <= c
rho_c(z) = c |z| - 0.5 c^2        , |z| > c
```

Default `c=1.345` configurable'dır. Bu değer standardized Gaussian residual
bağlamındaki yaklaşık yüzde 95 asymptotic efficiency için classical default'tur;
rubber engineering acceptance limiti değildir. Valid standard uncertainty yoksa
Huber objective kurulmaz; raw residual için keyfi scale tahmin edilmez.

## Grid-free regime splitting

Her merged interval'de residual ve variance lineerdir:

```text
r(t) = a + b t
v(t) = c0 + d t > 0
```

Regime crossing noktaları:

```text
r(t)^2 - c^2 v(t) = 0
```

denkleminin interval içindeki gerçek kökleridir. Denklem en çok quadratic'tir.
Katsayılar ölçeklenir, cancellation-resistant quadratic root biçimi kullanılır,
machine-equivalent endpoint/duplicate kökler compact edilir ve near-tangent
discriminant yalnız floating-point tolerance içinde sıfıra yuvarlanır. Gerçekten
interval dışındaki kökler içeri clamp edilmez.

## Quadratic ve tail integralleri

`|z|<=c` bölgesinde:

```text
integral rho_c(z) dt = 0.5 integral r^2/v dt
```

olup weighted analitik integral yeniden kullanılır.

Tail bölgesinde residual işareti crossing olmadan sabittir:

```text
integral rho_c(z) dt =
  c integral |r|/sqrt(v) dt - 0.5 c^2 h
```

`integral r/sqrt(v)` için `u=sqrt(v)` substitution'ı analitik uygulanır.
Implementation büyük variance endpoint'ini base seçip stable momentleri

```text
B0 = 2/(sqrt(v_end/v_base)+1)
B1 = 2(sqrt(v_end/v_base)+2) /
     [3(sqrt(v_end/v_base)+1)^2]
```

kullanır. Bu biçim constant-variance limitinde de süreklidir. Sabit arbitrary
quadrature grid yoktur.

## Diagnostics ve yorum

Her channel ve production joint sonuç için total valid overlap, quadratic width,
tail width, tail fraction ve RMS standardized residual saklanır. Tail fraction
için hard-coded pass/fail eşiği yoktur.

Huber tail'e giren bir nokta invalid measurement veya otomatik silinecek
outlier anlamına gelmez. Hiçbir point silinmez; yalnız tail loss'un influence
büyümesi quadratic L2'ye göre sınırlandırılır. Robust shift yine V0.8.1
authoritative shift'in yerine geçmez ve düşük Huber objective TRS/product
validation sonucu değildir.
