# V0.6 Harmonik Cevap Doğrulaması

## Kapsam ve kanıt düzeyi

Bu belge, direct/full-order, lineer, frequency-domain ve frozen-property
torsional harmonic solver'ın analytical code verification kapsamını tanımlar.
Deneysel TVD korelasyonu, material calibration veya solverlar arası ürün
eşdeğerliği iddiası değildir.

Tüm kompleks genlikler peak amplitude ve zaman konvansiyonu
`exp(+i*omega*t)` ile yorumlanır.

## Genel tolerans ilkeleri

- Moderate-scale analitik cevaplar: yaklaşık `1e-10` bağıl hata.
- Exact assembly/topology katsayıları: uygun olduğunda exact veya `1e-12`
  mutlak hata.
- Backend-independent relative residual: moderate problemlerde en çok
  `1e-10`.
- Simetri: problem ölçeği ve machine epsilon ile ölçeklenen tolerans.
- Input immutability: authoritative input kopyalarında exact eşitlik.
- `RCOND`, `FERR` ve `BERR`: exact platformlar arası değer zorlanmaz; çözülen
  noktada sonlu ve negatif olmayan değerler beklenir.

Assertion yordamları actual, expected ve tolerance değerlerinin sonluluğunu
kontrol etmelidir. Kompleks değerlerde reel ve sanal bileşenler ayrı IEEE
finite denetiminden geçmelidir.

## Eleman ve global matris regresyonları

### Lokal loss stiffness

İki düğümlü eleman için:

\[
K''_e=k''\begin{bmatrix}1&-1\\-1&1\end{bmatrix}
\]

exact katsayı, simetri, sıfır satır toplamı ve `[1,1]^T` rigid rotation
invariantlarıyla doğrulanır.

### Lokal viscous damping

\[
C_e=c\begin{bmatrix}1&-1\\-1&1\end{bmatrix}
\]

aynı invariantlarla sınanır. Ortak angular velocity iç moment üretmemelidir.

### Üç düğümlü zincir

İki farklı elemanlı üç düğümlü modelde `K'`, `K''`, `C` ve diagonal `M`
matrisleri bağımsız exact referanslarla karşılaştırılır. Dört matrisin Physical
DOF ve Equation ID sırası aynı olmalıdır. `K'`, `K''`, `C` ve `M` simetrisi
ayrı kontrol edilir.

### Constraint reduction

İlk düğümü fixed üçlü zincirde reduced matrisler, full matrislerin aynı
retained-index principal alt matrisleri olmalıdır:

\[
A_r=P^T A_{full}P,
\qquad A\in\{K',K'',C,M\}.
\]

## Dynamic stiffness doğrulaması

Seçilen moderate frequency için üretim builder sonucu:

\[
\operatorname{Re}Z=K'-\omega^2M,
\qquad
\operatorname{Im}Z=K''+\omega C
\]

bağıntılarıyla sınanır. Kontroller:

- `Z^T=Z`,
- kayıp bulunduğunda genel olarak `Z^H!=Z`,
- original `K'`, `K''`, `C` ve `M` değerlerinin değişmemesi,
- `K''` ve `C` katkılarının ayrı frequency scaling davranışı,
- sonlu olmayan dynamic-stiffness çıktısının reddi.

## Complex linear solver testleri

### Complex symmetric, Hermitian olmayan sistem

Sentetik non-Hermitian complex-symmetric matris, önceden seçilmiş kompleks
`x` ve `b=A*x` ile çözülür. Test:

- `A^T=A`,
- `A^H!=A`,
- hesaplanan `x` ile bilinen `x` eşleşmesi,
- küçük backend-independent residual

koşullarını doğrular. Böylece Hermitian solver'ın yanlışlıkla kullanılması
regresyonla engellenir.

### Exact singular

Örneğin `A=diag(1,0)` deterministic olarak singular'dır. Beklenen durum
`SINGULAR`, `RCOND=0` ve geçerli çözüm bulunmamasıdır. Bu analysis-state
sonucu crash veya `error stop` değildir; çözüm vektörü uydurulmamalıdır.

### Working-precision ill-conditioned

`A=diag(1,epsilon^2)` gibi nonsingular fakat güçlü ölçek farkı taşıyan
deterministic bir sistem `RCOND<epsilon` yolunu sınamak için kullanılır.
Beklenen durum `SOLVED_ILL_CONDITIONED` olup çözüm, `RCOND`, `FERR`, `BERR` ve
independent residual korunmalıdır. Exact `RCOND/FERR/BERR` değerleri farklı
LAPACK sağlayıcıları arasında zorlanmaz.

### Multiple RHS ve immutability

Low-level solver contract en az sentetik iki RHS ile boyut ve per-RHS
`FERR/BERR` eşleşmesini doğrular. `A` ve `B`, ZSYSVX'in çalışma kopyalarını
overwrite etmesine rağmen solve öncesi ve sonrası exact aynı kalmalıdır.

## Tek-DOF analitik doğrulamalar

Fixed-ground eleman ve tek aktif inertia için:

\[
\hat\theta=\frac{\hat T}
{k'-m\omega^2+i(k''+\omega c)}.
\]

Aşağıdaki bağımsız vakalar production assembly, reduction, dynamic builder,
complex solver, recovery ve result katmanını uçtan uca çalıştırır:

| Vaka | Fizik kanalları | Kontroller |
| --- | --- | --- |
| Viscous-only | `k''=0`, `c>0` | reel/sanal cevap, magnitude, phase, residual |
| Loss-only | `k''>0`, `c=0` | magnitude, phase, residual, dissipated energy |
| Combined | `k''>0`, `c>0` | iki katkının toplamı ve birbirine dönüştürülmemesi |
| Low-frequency | `k''=c=0`, küçük `f>0` | `theta_hat -> T_hat/k'` quasi-static limiti |

### Faz konvansiyonu

Pozitif reel torque ve passive damping ile response phase:

- resonance altında `0` ile `-pi/2` arasında,
- resonance çevresinde yaklaşık `-pi/2`,
- resonance üstünde `-pi/2` ile `-pi` arasında

olmalıdır. Test, asymptotik exact eşitlik yerine bu aralıkları ve sanal işaret
konvansiyonunu kilitler.

## İki ataletli TVD doğrulaması

### Fixed hub

Göbek constrained, halka aktif iken generalized harmonic cevap:

\[
\hat\theta_r=\frac{\hat T}
{K'-\omega^2J_r+i(K''+\omega c)}
\]

analitik referansla karşılaştırılır. Physical recovery sonrası hub kompleks
cevabı sıfır, ring cevabı reduced değerle aynı olmalıdır.

### K'' köprüsü

`two_inertia_tvd_system_t%loss_stiffness_nm_per_rad` değeri generalized
elemanın `loss_stiffness_nm_per_rad` alanına exact aktarılmalı;
`damping_nms_per_rad` alanına aktarılmamalıdır.

### Free-free finite frequency

Free-free sistem finite `f>0` değerinde çözülmelidir. Dengeli torque çifti için
relative response:

\[
\Delta\hat\theta=\frac{\hat T}
{K'-\omega^2J_{eq}+i(K''+\omega c)},
\qquad J_{eq}=\frac{J_hJ_r}{J_h+J_r}
\]

ile doğrulanır. V0.5'te rigid mode bulunması finite-frequency harmonic çözümü
otomatik olarak geçersiz kılmamalıdır.

### Modal-harmonic cross-validation

Çok küçük fakat pozitif damping ile harmonic response peak'i V0.5 natural
frequency çevresinde olmalıdır. Coarse ve refined explicit grid'lerde refined
peak location modal frekansa yaklaşmalıdır. Grid'in exact resonance noktasını
içermesi zorunlu tutulmaz.

## Excitation ve frequency-sweep doğrulaması

Geçerli yük testleri sonlu kompleks torque'u ve aynı active DOF'a birden çok
katkının scatter-add edilmesini kapsar. Ayrı expected-failure süreçleri:

- tanımsız node,
- desteklenmeyen DOF türü,
- constrained hedef,
- reel veya sanal bileşende NaN/+Inf/-Inf,
- scatter-add sırasında sonlu aralık taşması

vakalarını sınar.

Frequency array non-empty, sonlu, strictly positive ve strictly increasing
olmalıdır. Sıfır, negatif, NaN, +Inf, -Inf, unordered ve duplicate noktalar
ayrı test edilmelidir. Bir singular nokta bulunan üç noktalı sweep diğer iki
çözümü korumalıdır.

## Result, FRF ve kinematik türevler

Çözülen her noktada result frequency, status, `RCOND`, residual, `FERR/BERR`,
reduced response, physical response ve backend identity alanlarını aynı indeks
altında tutmalıdır. Magnitude/phase ile:

\[
\hat\Omega=i\omega\hat\theta,
\qquad \hat\alpha=-\omega^2\hat\theta
\]

yardımcıları doğrudan kompleks cebirle doğrulanır.

Tek tanımlı ve sıfır olmayan torque input için receptance, mobility ve
accelerance birimleri ve sayısal değerleri ayrı sınanır. Arbitrary multi-point
forced response'a generic FRF etiketi verilmediği API/test adlarında
korunmalıdır.

## Relative angle, torque ve passivity

`node_i -> node_j` orientation için relative angle ve:

\[
\hat T_e=[K'_e+i(K''_e+\omega c_e)]\Delta\hat\theta_e
\]

değeri analitik eleman referansıyla karşılaştırılır. Internal generalized-force
işaretleri `T_i=+T_e`, `T_j=-T_e` olmalıdır.

`K''>0` ve/veya `c>0` için `P_avg>=0`, `E_cycle>=0`; ikisi sıfırken exact sıfır
beklenir. Negative/nonfinite `K''` ve `c` üretim doğrulamasınca reddedilir.

## Analysis-state ve input-error ayrımı

- Exact singular veya working-precision ill-conditioned frequency point,
  result status'tur ve sweep'i öldürmez.
- Geçersiz frekans, excitation, matris boyutu ya da passive-property ihlali
  input/programming error'dır ve temiz `error stop` olabilir.
- Fully constrained active sistem LAPACK çağrısından önce anlamlı tanıyla
  reddedilir.

## Platform ve regression kapsamı

V0.1--V0.5 CTest takımının tamamı korunur. V0.6 testleri hem macOS Homebrew
LP64 LAPACK hem Windows MSYS2 LP64 OpenBLAS/LAPACK üzerinde çalıştırılır.
ZSYSVX symbol availability, explicit Fortran ABI ve runtime bağlantı böylece
iki CI platformunda gerçekten doğrulanır.

## Kapsam dışı doğrulamalar

Frequency-dependent material interpolation, nonlinear harmonic, transient
solver, mode-superposition, sparse/iterative solver, dynamic prescribed angle,
reaction torque ve phase unwrapping bu validation kapsamına dahil değildir.
