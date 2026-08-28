# Thermorheological Runtime Doğrulaması

## Kapsam ve kanıt yaklaşımı

V0.8.0 doğrulaması shift physics, reduced-coordinate hesabı, dynamic material
evaluation, mevcut V0.7 harmonic entegrasyonu, trace ve V0.1–V0.7 regression
kapılarını birlikte kapsar. Toleranslar deney belirsizliği değil IEEE
double-precision code verification içindir. Analitik referanslar production
formülünü çağırmadan bağımsız hesaplanır.

## Gate A — Shift Physics

### WLF

- `T=T_ref`: `log10(a_T)=0`, `a_T=1`, `f_r=f` identity,
- positive `C1,C2` ile `T>T_ref -> a_T<1 -> f_r<f`,
- positive `C1,C2` ile `T<T_ref -> a_T>1 -> f_r>f`,
- hand-calculated WLF `log10(a_T)` ve reduced frequency,
- `C1<=0`, `C2<=0`, invalid/nonfinite `T_ref`, `T_min` veya `T_max`,
  reference temperature'ın domain dışında kalması ve invalid query rejection,
- `D(T)=C2+(T-T_ref)` paydasının validated domain içinde sıfıra ulaşması veya
  pole'un diğer tarafına geçmesi halinde construction failure.

### Arrhenius

- `T=T_ref` identity,
- hand-calculated
  `ln(a_T)=Ea/R(1/T-1/T_ref)` ile provider sonucu,
- positive `Ea` için doğru sıcaklık yönü,
- `Ea<=0`, `T<=0`, invalid/nonfinite reference/domain ve domain dışı query
  rejection.

### Tabulated `log10(a_T)`

- exact temperature point ve reference point'te exact zero shift,
- iki nokta arasında `log10(a_T)` linear midpoint interpolation,
- machine-equivalent exact-temperature davranışı,
- en az iki nokta, finite/positive/strictly increasing temperature,
- duplicate ve unordered temperature ile NaN/Inf input rejection,
- `log10(a_T)` sonluluk kontrolü,
- her iki domain ucunda no extrapolation.

Measured shift değerlerinin monotonic olması bir acceptance testi değildir.
Yalnız temperature axis strictly increasing olmalıdır.

Concrete modelden bağımsız ortak provider sınırı; operating/query ve reference
sıcaklık eşleşmesini, model kind'i, `log10(a_T)` ile `a_T` tutarlılığını ve
tabulated/analytical bracket semantiğini lookup öncesinde doğrular. Kasıtlı
tutarsız test double ve NaN bracket girdileri clean rejection üretmelidir.

## Gate B — Reduced Frequency

Her model için:

\[
\log_{10}(f_r)=\log_{10}(f)+\log_{10}(a_T),
\qquad f_r=a_T f
\]

bağıntısı bağımsız referansla doğrulanır. Büyük positive/negative synthetic
shift'lerde implementation önce log-space master-curve domain kontrolü yapmalı;
geçersiz reduced point için ara overflow/underflow üretmeden temiz tanı
vermelidir.

Shift temperature domain'i içinde kalan fakat `f_r` master curve dışında olan
bir sorgu reddedilir. Örneğin `f=1000 Hz`, `a_T=1e-4` ile `f_r=0.1 Hz`, master
curve domain'i `1..1e6 Hz` ise material evaluation başarısız olmalıdır.

## Gate C — Material Evaluation

### Reference-temperature consistency

Master curve metadata `T_ref_master` ile shift provider
`T_ref_shift` machine-equivalent değilse thermorheological provider
construction başarısız olmalıdır.

### `T_ref` identity regression

Thermorheological provider `T=T_ref` için aynı V0.7 master-curve provider'ın
aynı `f` sorgusuyla aynı `G'` ve `G''` sonuçlarını üretmelidir.

### Physical ve reduced koordinat ayrımı

`f=100 Hz`, `a_T=0.1` örneğinde:

- physical frequency `100 Hz`,
- reduced/lookup frequency `10 Hz`,
- returned modulus frequency `100 Hz`,
- returned modulus temperature operating `T`,
- interpolation bracket 10 Hz çevresindeki master-curve bracket'ı

olmalıdır. Physical 100 Hz bracket olarak kullanılmamalıdır.

### Passivity

Validated master curve'de `G'>0`, `G''>=0` için bütün geçerli temperature ve
frequency sorgularında output `G'>0`, `G''>=0`; bonded-annular mapping
sonrasında `K'>0`, `K''>=0` kalmalıdır.

## Gate D — Harmonic Integration

Full fixed-hub 1-DOF test zinciri bağımsız olarak hesaplanır:

```text
T -> a_T -> f_r -> master-curve interpolation
  -> G',G'' -> K',K'' -> Z -> theta
```

Üretim testi özellikle şu mevcut zinciri kullanmalıdır:

```text
thermorheological provider
  -> dynamic_torsional_property_binding_t
  -> analyze_material_aware_harmonic_response()
```

Yeni harmonic orchestration veya solver facade kullanılmamalıdır. Aynı modelde
bir element WLF, başka bir element Arrhenius provider kullanabilmeli; her biri
kendi `T_ref`, domain ve shift factor'ını uygulamalıdır.

## Gate E — Trace

Trace aşağıdaki authoritative durumları birbirine karıştırmadan taşımalıdır:

- physical `f` ve operating `T`,
- reference temperature,
- shift model kind, `log10(a_T)` ve `a_T`,
- reduced/lookup `f_r`,
- master-curve interpolation policy, bracket ve alpha,
- `G'`, `G''`, `K'` ve `K''`.

Interpolated point invariantı lookup axis üzerinde sınanır:

\[
f_{lower,lookup}<f_{lookup}<f_{upper,lookup}.
\]

Valid material evaluation sonrasında harmonic `Z` singular olabilir. Bu
durumda harmonic status `SINGULAR` olur; material trace kaybolmaz. Unshifted
V0.7 trace için `temperature_shift_applied=false`, `log10(a_T)=0`, `a_T=1` ve
`lookup_frequency=physical_frequency` geriye uyumluluğu korunur.

## Gate F — Regression

Aşağıdaki grupların tamamı aynı CTest çalışmasında geçmelidir:

- V0.1–V0.4 core, assembly ve constraint testleri,
- V0.5 DSYGV modal testleri,
- V0.6 frozen harmonic ve ZSYSVX testleri,
- V0.7 tabulated provider, binding ve material-aware harmonic testleri,
- V0.8 shift/provider/harmonic/trace testleri.

V0.8 testleri V0.5 modal material semantics'ini frequency-dependent hale
getirmez ve V0.6 frozen harmonic sonucunu değiştirmez.

## Gate G — Platformlar

Definition of Done clean Debug Ninja build ve tüm CTest takımını gerektirir.
Push edilen commit için macOS GNU Fortran ve Windows MinGW64 GNU Fortran GitHub
Actions işleri ayrı ayrı başarılı olmadan V0.8.0 tamamlanmış sayılmaz.

## Validation sınırı

Bu kanıtlar analytical code verification'dır. Master-curve identification,
TRS acceptance ve gerçek DMA/TVD experimental correlation V0.8.0'ın test
kapsamı değildir.
