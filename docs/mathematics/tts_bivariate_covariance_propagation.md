# TTS Bivariate Covariance Propagation

## Canonical physical covariance

Tek physical DMA noktası için dış input:

\[
\boldsymbol\Sigma_G=
\begin{bmatrix}
\operatorname{Var}(G') & \operatorname{Cov}(G',G'')\\
\operatorname{Cov}(G',G'') & \operatorname{Var}(G'')
\end{bmatrix}\quad[\mathrm{Pa^2}].
\]

TMS26 universal rubber correlation coefficient varsaymaz. Matrix doğrudan
ölçüm evidence'ından veya açık bir measurement modelinden gelmelidir.

## Magnitude/phase measurement modeli

Complex modulus polar olarak ölçüldüğünde

\[
G^*=M e^{i\delta},\qquad G'=M\cos\delta,\qquad G''=M\sin\delta,
\]

burada `M` Pa, `delta` radyandır. Input covariance

\[
\boldsymbol\Sigma_{M\delta}=
\begin{bmatrix}
\operatorname{Var}(M)&\operatorname{Cov}(M,\delta)\\
\operatorname{Cov}(M,\delta)&\operatorname{Var}(\delta)
\end{bmatrix}
\]

ve Jacobian

\[
\mathbf J=
\begin{bmatrix}
\cos\delta&-M\sin\delta\\
\sin\delta&M\cos\delta
\end{bmatrix}
\]

ile first-order propagation

\[
\boldsymbol\Sigma_G=\mathbf J\boldsymbol\Sigma_{M\delta}\mathbf J^T
\]

uygulanır. `Cov(M,delta)=0` olsa bile ortak input'lar nedeniyle genel durumda
`Cov(G',G'')` sıfır değildir. Independent magnitude/phase, independent
storage/loss anlamına gelmez.

## Log-space propagation

TTS residual coordinates:

\[
\mathbf y=\begin{bmatrix}\log_{10}G'\\\log_{10}G''\end{bmatrix},\qquad
\mathbf D=\frac{1}{\ln 10}
\begin{bmatrix}1/G'&0\\0&1/G''\end{bmatrix}.
\]

`G'>0` ve `G''>0` için

\[
\boldsymbol\Sigma_y=\mathbf D\boldsymbol\Sigma_G\mathbf D^T
\]

boyutsuz covariance verir. `G''=0` passive runtime material açısından geçerli
olabilir fakat `log10(G'')` tanımsızdır; V0.8.5 support'una alınmaz ve epsilon
clamp uygulanmaz.

Pozitif diagonal Jacobian scaling altında correlation invariant'tır:

\[
\rho_{\log G',\log G''}=\rho_{G',G''}.
\]

Pa↔MPa dönüşümünde modulus, standard uncertainty ve covariance tutarlı
ölçeklenirse `Sigma_y` değişmez.

## First-order anlamı

Sonuçlar `first_order_covariance_propagation=.true.` taşır. Jacobian yöntemi
exact probability-distribution dönüşümü değildir. Büyük relative uncertainty
ve güçlü nonlinearity durumunda distribution propagation gerekebilir;
[JCGM 101:2008](https://www.bipm.org/en/doi/10.59161/jcgm101-2008) Monte Carlo
yaklaşımı gelecek bağlamıdır, V0.8.5 implementation'ı değildir.

Measurement model ve covariance propagation terminolojisi
[JCGM 100:2008](https://www.bipm.org/en/doi/10.59161/jcgm100-2008e),
[JCGM GUM-6:2020](https://www.bipm.org/fr/doi/10.59161/jcgmgum-6-2020) ve
multiple-output bağlamı için
[JCGM 102:2011](https://www.bipm.org/en/doi/10.59161/jcgm102-2011) ile uyumlu
konumlandırılmıştır. Bu referanslar TMS26'ya özel covariance katsayısı vermez.

## SPD ve conditioning

\[
v_s>0,\quad v_l>0,\quad \det\Sigma=v_sv_l-c^2>0,
\quad |\rho|<1
\]

zorunludur. Ölçeklenmiş 2x2 eigenvalue oranından elde edilen reciprocal
condition estimate yalnız machine-derived `sqrt(epsilon)` sınırıyla kontrol
edilir. Perfect ve materially near-perfect correlation açık status üretir.
Jitter, eigenvalue clipping, shrinkage veya pseudoinverse uygulanmaz.
