# Torsional Sönüm Modelleri

## Fiziksel kanallar

TMS26 V0.6.0, elastik depolama, yapısal/viskoelastik kayıp ve viskoz sönümü
üç ayrı constitutive quantity olarak temsil eder.

### Depolama rijitliği K'

`K' [N·m/rad]`, dönme ile faz içi geri çağırıcı moment bileşenidir. Bir çevrim
içinde depolanan elastik enerjiyi geri verebilen kısmı temsil eder ve kararlı
lineer modelde `K'>0` olmalıdır.

### Kayıp rijitliği K''

`K'' [N·m/rad]`, kompleks rijitliğin dönmeye göre faz dışı bileşenidir:

\[
K^*=K'+iK''.
\]

Frozen-property harmonic noktada elastomerin çevrimsel enerji kaybını temsil
eder. Passive modelde `K''>=0` olmalıdır.

### Viskoz sönüm c

`c [N·m·s/rad]`, bağıl açısal hızla orantılı moment katsayısıdır:

\[
T_c=c(\dot\theta_i-\dot\theta_j).
\]

`e^{+i omega t}` altında kompleks contribution `i omega c` olur. Passive
modelde `c>=0` olmalıdır.

## K'' ile c neden ayrı tutulur?

`K''` ve `c` farklı birimlere, fiziksel kökene ve frekans davranışına sahiptir.
Tek bir frekansta:

\[
iK''=i\omega(K''/\omega)
\]

yazılabilmesi, `K''/omega` değerinin genel bir constitutive viscous coefficient
olduğu anlamına gelmez. TMS26 bu nedenle:

- `c=K''` ataması yapmaz,
- public fizik modelinde otomatik `c=K''/omega` dönüşümü yapmaz,
- iki büyüklüğü generic bir damping matrix içinde karıştırmaz.

Elemanın toplam faz dışı dynamic coefficient değeri yalnız çözüm frekansında:

\[
D_e(\omega)=K''_e+\omega c_e
\]

olarak cebirsel biçimde birleşir.

## Harmonic convention ve passive faz

TMS26:

\[
\theta(t)=\operatorname{Re}\{\hat\theta e^{+i\omega t}\}
\]

konvansiyonunu ve peak complex amplitude değerlerini kullanır. Bu seçimle
passive `K''` ve `c` dynamic stiffness'in pozitif sanal kısmına katkı verir:

\[
Z=K'-\omega^2M+i(K''+\omega C).
\]

Pozitif reel torque için cevap fazı resonance çevresinde negatife ilerler. Eksi
işaretli bir sanal damping terimi bu konvansiyonla passive modeli değil, enerji
üreten aktif davranışı göstereceğinden kabul edilmez.

## Enerji ve passivity

Eleman relative complex rotation değeri:

\[
\Delta\hat\theta_e=\hat\theta_i-\hat\theta_j
\]

ise peak amplitude konvansiyonunda ortalama kayıp gücü:

\[
P_{avg,e}=\frac{\omega}{2}
(K''_e+\omega c_e)|\Delta\hat\theta_e|^2
\]

ve çevrim başına kaybedilen enerji:

\[
E_{cycle,e}=\pi
(K''_e+\omega c_e)|\Delta\hat\theta_e|^2
\]

olur. Birimler sırasıyla watt `[W]` ve joule/çevrim `[J/cycle]` değeridir.
Passive elemanda sonuçlar negatif olamaz. Floating-point roundoff için
ölçeğe duyarlı tolerans kullanılabilir; anlamlı negatif enerji sessizce sıfıra
kırpılmaz.

`K''=0` ve `c=0` ise eleman enerji tüketmez. Ortak rigid rotation veya ortak
angular velocity de relative deformation üretmediği için iç kayıp sıfırdır.

## TVD yorumu

TVD elemanında `K'` hub ile inertia ring arasındaki geri çağırıcı momenti,
`K''` elastomer kaybını, `c` ise ayrıca tanımlanmış viskoz bağlantıyı temsil
eder. V0.6 derived sonuçları:

- hub/ring kompleks açıları,
- oriented relative angle,
- kompleks dynamic element torque,
- transmitted-torque magnitude,
- ortalama dissipated power,
- çevrim başına dissipated energy

üretebilir. Bu değerler ileride termal, fatigue, DMA ve test korelasyonu için
girdi sağlayabilir; V0.6 bu korelasyonları kendisi çözmez.

## Frozen-property sınırı

V0.6 sweep boyunca `K'`, `K''`, `c` ve `M` değerlerini sabit kabul eder.
Mevcut malzeme frequency-point depolamasından interpolation yapılmaz. Özellikle:

- `G'(f,T)` ve `G''(f,T)` provider yoktur,
- WLF/Arrhenius shift uygulanmaz,
- amplitude dependence ve Payne effect yoktur,
- `K''` time-domain transient modele taşınmaz.

Frequency-dependent malzeme provider V0.7 için ayrı geliştirme alanıdır.

## Geçerlilik varsayımları

- lineer ve küçük genlikli torsional davranış,
- sinusoidal steady state,
- sonlu ve pozitif excitation frequency,
- passive elemanlar,
- peak, RMS olmayan kompleks genlikler,
- frozen/current malzeme özellikleri.

Nonlinear harmonic, transient viskoelastisite, Prony state variables ve
temperature-amplitude coupled material davranışı bu modelin kapsamında
değildir.
