# TTS Parametric Shift-Law Identification

## Ortak gözlem modeli

TMS26 `s(T)=log10(a_T)` ve `f_r=a_Tf` convention'ını kullanır. Parametric fit
girdisi adjacent relative shift'tir:

```text
delta_s_ij = s(T_j) - s(T_i)
```

Sıcaklıklar Kelvin, `s` ve `delta_s` boyutsuzdur. Equal-weight residual:

```text
r_ij = delta_s_measured,ij - delta_s_predicted,ij
```

Cumulative absolute `s(T)` değerleri independent observations değildir ve fit
denklemine girmez.

## Arrhenius

Apparent activation energy modeli:

```text
s(T) = Ea_app / (R ln(10)) * (1/T - 1/T_ref)
beta = Ea_app / (R ln(10))
x_ij = 1/T_j - 1/T_i
delta_s_ij = beta x_ij
```

`R=8.31446261815324 J/(mol K)`, `Ea_app [J/mol]`, `beta [K]` birimindedir.
Equal-weight analytical çözüm:

```text
beta = sum(x_ij delta_s_ij) / sum(x_ij^2)
Ea_app = R ln(10) beta
```

Kod, `x` değerlerini `max(abs(x))` ile normalize ederek `sum(x^2)` underflow
riskini azaltır. Conventional runtime modeli için `Ea_app>0` gerekir.
`Ea_app`, reference temperature'dan bağımsızdır; yeni reference yalnız absolute
shift sıfırını değiştirir:

```text
s_new(T) = s_old(T) - s_old(T_ref,new)
```

`Ea_app` chemical aging reaction activation energy'si olarak yorumlanmaz.

## WLF profiled fit

```text
s(T) = -C1 (T-T_ref) / (C2 + T-T_ref)
g(T,C2) = -(T-T_ref) / (C2 + T-T_ref)
delta_g_ij = g(T_j,C2) - g(T_i,C2)
```

Fixed C2 için C1 analytical profillenir:

```text
C1(C2) = sum(delta_g_ij delta_s_ij) / sum(delta_g_ij^2)
J(C2) = mean((delta_s_ij - C1(C2) delta_g_ij)^2)
```

Yalnız `J(C2)` mevcut bounded Brent minimizer ile minimize edilir. `C1`
boyutsuz, `C2 [K]` ve conventional admissibility `C1>0`, `C2>0` koşuludur.

Measured domain boyunca pole yasaktır:

```text
C2 + T - T_ref > 0
C2 > max(0, T_ref-Tmin) + machine_margin
```

İlk C2 scale measured temperature span'den türetilir. Upper sample iki kat
genişletilir ve en çok 64 expansion yapılır. Bunlar numerical algorithm
constants'tır; material acceptance threshold değildir.

## WLF identifiability

```text
p = C1/C2
q = 1/C2
s(T) = -p delta_T / (1 + q delta_T)
```

`C2` büyüdükçe `q->0` ve `s(T)≈-p delta_T` olur. Profile objective large-C2
limitine kaçıyorsa veya finite optimum linear-limit objective'e yalnız machine
scale kadar iyileşme getiriyorsa `p` saklanır; `C1/C2` parameter çifti
`parameter_identifiable=false` olur. Böyle bir sonuç düşük residual taşısa da
WLF runtime export için yeterli değildir.

Reference `T_ref'=T_ref+d` olduğunda exact dönüşüm:

```text
C2' = C2 + d
C1' = C1 C2 / (C2+d)
C1' C2' = C1 C2
```

## Residual ve predictive diagnostics

Her iki model pair RMSE, maximum absolute residual ve mean residual taşır.
Residual degree-of-freedom Arrhenius için `n-1`, WLF için `n-2`'dir.

En az beş unique temperature varsa LOTO uygulanır. Held-out `T_k` içeren bütün
pair'ler training set'ten çıkarılır; kalan pair'lerle model yeniden fit edilir
ve çıkarılan relative shifts tahmin edilir. Minimum dataset'te LOTO mümkün
değilse fit reddedilmez, yalnız `predictive_validation_available=false` kalır.

Bu ölçütler numerical/model evidence'dır. Universal rubber TRS acceptance
threshold'u veya automatic WLF/Arrhenius winner üretmez.
