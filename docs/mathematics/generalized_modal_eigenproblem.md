# Genelleştirilmiş Torsional Modal Özdeğer Problemi

## Model

TMS26 V0.5.0, lineer ve sönümsüz torsional sistemin constraint sonrası
hareket denklemini çözer:

\[
M_r\,\ddot{q}_r + K_r\,q_r = 0
\]

Harmonik çözüm varsayımı `q_r(t)=phi exp(i omega t)` ile:

\[
K_r\,\phi_i = \lambda_i M_r\,\phi_i, \qquad
\lambda_i=\omega_i^2, \qquad
f_i=\frac{\omega_i}{2\pi}
\]

Burada `K_r` torsional rijitlik matrisi `[N·m/rad]`, `M_r` polar kütle
ataleti matrisi `[kg·m²]`, `lambda` kavramsal olarak `[1/s²]`, `omega`
`[rad/s]` ve `f` `[Hz]` birimindedir. Radyan boyutsuz kabul edilir.

## Matematiksel önkoşullar

- `K_r` ve `M_r` kare, eş boyutlu, sonlu ve simetrik olmalıdır.
- `M_r` symmetric positive definite (SPD) olmalıdır.
- `K_r` positive definite olmak zorunda değildir. Serbest sistemlerde
  positive semidefinite ve singular olabilir.
- Sıfır active DOF, geçerli bir reduction sonucu olsa da modal çözüm içermez;
  LAPACK çağrısından önce anlamlı bir tanı ile reddedilir.

V0.5.0, `M_r` için ikinci bir özel Cholesky algoritması yazmaz. SPD koşulunun
nihai doğrulaması LAPACK `DSYGV` faktörizasyonu ve `INFO` sözleşmesiyle yapılır.

## Reference dense backend

Dense reference backend, LAPACK `DSYGV` yordamını şu sabit sözleşmeyle çağırır:

- `ITYPE=1`: `K_r phi = lambda M_r phi`
- `JOBZ='V'`: eigenvalue ve eigenvector üretimi
- `UPLO='U'`: simetrik matrislerin üst üçgeni
- LP64 ABI: LAPACK `INTEGER` 32 bit

Optimal çalışma alanı `LWORK=-1` query ile öğrenilir. `DSYGV`, giriş
matrislerini değiştirdiği için query ve gerçek çözüm ayrı çalışma kopyalarında
yürütülür; sonuç doğrulaması her zaman değişmemiş original `K_r` ve `M_r`
ile yapılır.

## Rijit-cisim modu ve ölçeğe duyarlı tolerans

Sabit bir frekans eşiği kullanılmaz. Karakteristik özdeğer ölçeği:

\[
s_\lambda = \max\left(\max_i |\lambda_i|,
\frac{\lVert K_r\rVert_\infty}{\lVert M_r\rVert_\infty}\right)
\]

olarak alınır. Otomatik tolerans:

\[
\tau_\lambda = C_\mathrm{rigid}\,\epsilon\,s_\lambda
\]

biçimindedir. `C_rigid`, adı olan ve varsayılanı `O(10²)` mertebesindeki
sayısal çarpandır; fiziksel bir Hz sınırı değildir.

- `lambda < -tau_lambda`: kararsız/geçersiz fizik modeli
- `-tau_lambda <= lambda <= tau_lambda`: `RIGID_MODE`
- `lambda > tau_lambda`: `ELASTIC_MODE`

Tamamen sıfır rijitlikte ölçek sıfırdır; exact zero eigenvalue rijit moddur.
Birden fazla ayrık serbest bileşen birden fazla rijit mod üretebilir.

## Kütle normalizasyonu

Varsayılan mode-shape normalizasyonu:

\[
\phi_i^T M_r\phi_i=1
\]

şeklindedir. `DSYGV`, `ITYPE=1` için `M_r`-orthonormal eigenvector üretir.
TMS26 önce bu sonucu doğrular. Yalnız sayısal sapma tanımlı sınırı aşarsa,
pozitif kütle normu kullanılarak açık bir cleanup uygulanır.
Mass-normalized `phi` genliği genel olarak boyutsuz değildir; ölçeği
`[kg·m²]^{-1/2}` kütle metriğine bağlıdır ve fiziksel mode shape oranları
global işaret/ölçekten bağımsız yorumlanır.

Tüm mod matrisi `Phi` için global kalite ölçütü:

\[
E_M=\lVert\Phi^T M_r\Phi-I\rVert_F
\]

olarak saklanır. Test altyapısı ayrıca
`Phi^T K_r Phi ~= diag(lambda)` bağıntısını doğrular.

## Boyutsuz eigenpair residual

Her modun kalıntısı:

\[
r_i=K_r\phi_i-\lambda_i M_r\phi_i
\]

ve göreli residual değeri:

\[
\rho_i=\frac{\lVert r_i\rVert_2}
{(\lVert K_r\rVert_\infty+|\lambda_i|\lVert M_r\rVert_\infty)
\lVert\phi_i\rVert_2}
\]

biçimindedir. Payda ve pay birlikte sayısal olarak sıfırsa `rho_i=0`
tanımlanır; sıfır payda ile sıfır olmayan kalıntı başarılı sonuç sayılmaz.

## İşaret ve tekrarlı özdeğerler

`phi` ile `-phi` aynı fiziksel mode shape'tir. Testler doğrudan bileşen
eşitliği yerine sign-invariant modal correlation veya işaret hizalama kullanır.

Repeated veya clustered eigenvalue durumunda tekil eigenvector unique değildir.
Bu nedenle multiplicity, residual, `M`-orthogonality ve ilgili eigenspace
projector/overlap doğrulanır; belirli bir baz vektörü zorunlu tutulmaz.

## Fiziksel mode recovery

Constraint sonrası reduced mode:

\[
\phi=P\phi_r
\]

ile fiziksel DOF uzayına açılır. Constraint bileşenleri sıfırdır. Statik durum
recovery denklemindeki prescribed offset `q_p`, mode shape'e eklenmez.

## Fiziksel kapsam

Sonuç `linear`, `undamped` ve `frozen-property` modal çözümdür. Elastomer
`G'` veya `K'` hangi frekans-sıcaklık çalışma noktasından geldiyse çözüm o
değerleri sabit kabul eder. V0.5.0,
`G'(f) -> K(f) -> eigenfrequency` öz-tutarlı iterasyonunu çözmez.
