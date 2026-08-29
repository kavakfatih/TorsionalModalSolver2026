# Parametric Temperature Shift Laws

## Empirical ve parametrik model ayrımı

V0.8.1 empirical `temperature, log10(a_T)` tablosu ölçümden türetilen
authoritative davranıştır. V0.8.2 Arrhenius ve WLF sonuçları bu davranışın
parametric approximations'ıdır. Fit başarılı olsa bile experimental shift
tablosu değiştirilmez ve model otomatik olarak daha doğru ilan edilmez.

Arrhenius sonucu `Ea_app` apparent activation energy olarak adlandırılır. Bu
değer horizontal dynamic-modulus shifting için etkili bir parametredir;
chemical aging, oxidation veya cure reaction activation energy'si değildir.

WLF sonucu C1/C2 yanında `p=C1/C2` ve `q=1/C2` tanılarını taşır. Large-C2
limitinde yalnız p belirlenebiliyorsa ayrı C1/C2 material constants gibi
raporlanmaz ve runtime export kapalı kalır.

## Material-state sözleşmesi

Fit yalnız aynı V0.8.1 common test state içindeki isotherm family için
geçerlidir:

- material ve compound/batch state,
- cure/process ve conditioning state,
- dynamic shear strain amplitude,
- static shear prestrain,
- deformation mode ve test method.

Aynı law farklı strain amplitude, prestrain, batch veya conditioning family'ye
otomatik transfer edilmez.

## Fiziksel kapsam

Model horizontal thermorheological shifting ile `G*(f,T)` davranışını temsil
eder. Gerçek filled rubber genel olarak aşağıdaki daha geniş bağımlılığa sahip
olabilir:

```text
G* = G*(f, T, dynamic strain amplitude, static prestrain,
        compound/batch state, conditioning, ...)
```

İyi WLF veya Arrhenius residual'ı Payne effect'i, strain-amplitude dependence,
static-prestrain dependence, self heating veya product-level TVD correlation
problemini çözmez.

## Runtime sınırı

Parametric runtime export yalnız measured calibrated `Tmin<=T<=Tmax`
domain'inde geçerlidir. Temperature extrapolation yoktur. Master dynamic
modulus table V0.8.1'den gelir; yalnız temperature-shift provider explicit
olarak empirical tabulated, Arrhenius veya identifiable WLF seçilir.

## Mühendislik zinciri

```text
raw material / compound state
  -> rubber process / cure state
  -> DMA dynamic measurements
  -> V0.8.1 empirical TTS
  -> V0.8.2 WLF / Arrhenius identification
  -> G'(f,T), G''(f,T)
  -> K'(f,T), K''(f,T)
  -> TVD harmonic response
  -> product dynamic stiffness / angle-torque / fatigue correlation
  -> material/process tolerances
```

V0.8.2 yalnız parametric temperature-shift identification katmanını kapsar;
product correlation veya tolerance calibration yapmaz.
