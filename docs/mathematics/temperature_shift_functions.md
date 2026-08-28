# Temperature Shift Functions

## Canonical TMS26 convention

TMS26 bütün temperature-shift provider'larında tek convention kullanır:

\[
a_T(T)=\frac{\tau(T)}{\tau(T_{ref})},\qquad a_T(T_{ref})=1,
\]

\[
t_r=\frac{t}{a_T},\qquad f_r=a_T f.
\]

Authoritative computational quantity:

\[
s(T)=\log_{10}(a_T)
\]

olup `a_T=10^s` türetilmiş büyüklüktür. `s` shift'i decade cinsinden doğal
temsil eder ve reduced-frequency domain kontrolünün log-space'te yapılmasını
sağlar.

Tipik thermally activated polymer/rubber için `T>T_ref` olduğunda relaxation
hızlanır, `a_T<1` ve `f_r<f`; `T<T_ref` olduğunda `a_T>1` ve `f_r>f` beklenir.
Bu yön bir evrensel data-quality monotonicity kuralı değildir; WLF ve
Arrhenius analytical sanity testidir.

## WLF modeli

TMS26 WLF denklemi:

\[
\log_{10}(a_T)=
-\frac{C_1(T-T_{ref})}{C_2+(T-T_{ref})}.
\]

Canonical birimler ve kısıtlar:

| Büyüklük | Anlam | Birim / kısıt |
|---|---|---|
| `T` | operating temperature | K, `T>0` |
| `T_ref` | reference temperature | K, `T_ref>0` |
| `C1` | WLF katsayısı | boyutsuz, `C1>0` |
| `C2` | WLF sıcaklık katsayısı | K, `C2>0` |

Provider açık `T_min,T_max` domain'i taşır ve
`T_min<=T_ref<=T_max` olmalıdır. Payda:

\[
D(T)=C_2+(T-T_{ref})
\]

validated domain boyunca strictly positive kalmalıdır. Domain WLF pole'una
ulaşamaz veya pole'un diğer tarafına geçemez. `T=T_ref` için pay sıfırdır,
dolayısıyla `s=0`, `a_T=1` ve `f_r=f` elde edilir.

ANSYS ile parameter alışverişinde yalnız formül açıkça eşleştirildikten sonra
dönüşüm yapılmalıdır. TMS26 core içinde ANSYS'e ait shift quantity adı veya
işaret convention'ı kullanılmaz.

## Arrhenius modeli

TMS26 Arrhenius denklemi:

\[
\ln(a_T)=\frac{E_a}{R}
\left(\frac{1}{T}-\frac{1}{T_{ref}}\right),
\]

ve eşdeğer decimal-log biçimi:

\[
\log_{10}(a_T)=\frac{E_a}{2.303R}
\left(\frac{1}{T}-\frac{1}{T_{ref}}\right).
\]

| Büyüklük | Anlam | Birim / kısıt |
|---|---|---|
| `E_a` | activation energy | J/mol, `E_a>0` |
| `R` | universal gas constant | J/(mol K), `R>0` |
| `T`, `T_ref` | mutlak sıcaklık | K, pozitif ve sonlu |

Arrhenius provider da explicit `T_min,T_max` validated domain taşır.
Denklemin domain dışında sonlu bir sayı vermesi extrapolation izni değildir.
`T=T_ref` için identity, positive `E_a` ile sıcaklık yönü analitik olarak
doğrulanır.

## Tabulated `log10(a_T)` modeli

Input dizileri `temperature_k(:)` ve `log10_a_t(:)` ile reference temperature
taşır. En az iki nokta; sonlu, positive ve strictly increasing sıcaklık;
duplicate olmayan T; sonlu `s` gerekir. Reference temperature dataset içinde
explicit bulunmalı ve bu noktada `s=0` machine tolerance içinde olmalıdır.

`T_1<=T<=T_2` için:

\[
\alpha=\frac{T-T_1}{T_2-T_1},
\]

\[
s(T)=(1-\alpha)s_1+\alpha s_2,
\qquad a_T=10^{s(T)}.
\]

Interpolation doğrudan `a_T` üzerinde yapılmaz. Measured shift values için
monotonicity zorlanmaz. `T<T_min` veya `T>T_max` reddedilir; extrapolation,
endpoint hold ve nearest-temperature seçimi yoktur.

## Log-space reduced-frequency hesabı

Doğrudan `10^s*f` hesaplamak büyük `|s|` durumunda gereksiz ara
overflow/underflow üretebilir. TMS26 önce:

\[
\log_{10}(f_r)=\log_{10}(f)+s
\]

hesabını yapar ve sonucu master-curve sınırlarıyla karşılaştırır:

\[
\log_{10}(f_{r,min})\leq\log_{10}(f_r)
\leq\log_{10}(f_{r,max}).
\]

Yalnız domain geçerliyse `f_r=10^{\log_{10}(f_r)}` oluşturulur. Böylece iki
zorunlu domain ayrı kalır:

\[
T_{min}\leq T\leq T_{max},
\]

\[
f_{r,min}\leq a_T(T)f\leq f_{r,max}.
\]

Bir domain'in geçmesi diğerinin yerine geçmez. Hiçbir adım temperature veya
master-curve frequency extrapolation uygulamaz.

## Model sınırları

Bu denklemler horizontal-only `b_T=1`, lineer viskoelastik, küçük dinamik
genlikli ve sabit amplitude/prestrain durum içindir. Vertical modulus shift,
self-heating, thermal expansion, temperature-dependent density/inertia,
nonlinear harmonic iteration, combined WLF–Arrhenius ve transient
viscoelasticity dahil değildir.
