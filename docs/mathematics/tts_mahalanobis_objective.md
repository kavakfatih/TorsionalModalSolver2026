# TTS Mahalanobis Objective

## Residual ve covariance

V0.8.1 convention korunur: moving curve `x_shifted=x+s` ile taşınır,
`s=log10(a_T)` ve residual

\[
\mathbf r(x;s)=
\begin{bmatrix}
y'_i(x)-y'_j(x-s)\\
y''_i(x)-y''_j(x-s)
\end{bmatrix}.
\]

V0.8.5 cross-isotherm covariance modellemez. Bu nedenle

\[
\boldsymbol\Sigma_r(x;s)=\boldsymbol\Sigma_i(x)+
\boldsymbol\Sigma_j(x-s).
\]

Genel form `Sigma_i+Sigma_j-C_ij-C_ij^T` olurdu; bu sürümde `C_ij=0` kabul
edilir. Cross-frequency covariance da yoktur; dolayısıyla bu objective tüm DMA
dataset'i üzerinde GLS değildir.

## Squared Mahalanobis residual

\[
\boldsymbol\Sigma_r=\begin{bmatrix}v_s&c\\c&v_l\end{bmatrix}
\]

için

\[
d_M^2=\mathbf r^T\boldsymbol\Sigma_r^{-1}\mathbf r=
\frac{v_l r_s^2-2cr_sr_l+v_s r_l^2}{v_sv_l-c^2}.
\]

Pozitif correlation covariance ellipse'ının `[1,1]` yönünü `[1,-1]`
yönünden daha düşük mesafeli yapması, negatif correlation'da yönün tersine
dönmesi low-level fizik testidir. Point diagnostic generic inverse oluşturmaz;
2x2 Cholesky solve kullanır. Covariance-weighted quadratic formun temel
yorumu NIST'in [covariance matrix](https://www.itl.nist.gov/div898/handbook/pmc/section5/pmc541.htm)
ve [Mahalanobis/Hotelling distance](https://www.itl.nist.gov/div898/handbook/pmc/section3/pmc341.htm)
açıklamalarıyla uyumludur.

## Bivariate common support

`O_B`, storage ve positive-loss measurement quality, iki V0.8.4 uncertainty
channel'ı ve point-local covariance support'unun kesişimidir. Missing/invalid
covariance gap'tir; gap üzerinden bridge veya extrapolation yoktur. Mahalanobis
storage-only moda düşmez.

## Objective ve matched control

\[
J_M(s)=\frac{1}{2|O_B|}\int_{O_B}d_M^2\,dx.
\]

`1/2`, mevcut iki-channel equal-weight diagonal convention ile uyumluluk için
kasıtlıdır. Sample count değil, gerçek log-frequency overlap genişliği
normalize edilir.

Saf covariance etkisini support etkisinden ayırmak için aynı `O_B` üzerinde

\[
J_{D\cap}(s)=\frac{1}{2|O_B|}\int_{O_B}
\left(\frac{r_s^2}{v_s}+\frac{r_l^2}{v_l}\right)dx
\]

hesaplanır. `c=0` iken `J_M=J_Dcap` ve minimizer'lar eşittir. Original V0.8.4
support farklıysa `s_weighted_original` ile eşitlik zorunlu değildir.

## Grid-free interval integrali

Merged bir interval'de `q in [0,1]` için iki residual ve covariance matrix
elemanları lineerdir. Böylece

\[
Q_2(q)=a(q)b(q)-c(q)^2
\]

quadratic, adjugate numerator

\[
P_3(q)=b r_s^2-2cr_sr_l+a r_l^2
\]

en fazla cubic olur. Üretim integrali

\[
\frac{P_3}{Q_2}=P_1+\frac{R_1}{Q_2}
\]

polynomial division uygular; polynomial terimi doğrudan, remainder'ı
`ln(Q2)` ve discriminant'a göre `atan`, logarithmic veya repeated-root
primitive ile entegre eder. Nearly-linear determinant machine-scale limitinde
stable linear-denominator momentlerine geçilir. Bu numerical limit dalı fixed
quadrature değildir. Test oracle'ı üretim helper'ından bağımsız adaptive
Simpson integration kullanır.

SPD endpoint matrix'lerin convex entrywise interpolation'ı interval boyunca
SPD kalır. Yine de endpoint ve midpoint reciprocal conditioning savunmalı
kontrol edilir; singular denominator içinden integral alınmaz.

## Optimization ve invariants

Objective mevcut deterministic coarse scan, strict
`J[k-1]>J[k]<J[k+1]` interior bracket ve bounded Brent ile minimize edilir.
Flat veya boundary-only optimum explicit unavailable status'tur.

Ortak `Sigma -> lambda^2 Sigma` ölçeklemesi iki objective'i `1/lambda^2`
ölçekler fakat minimizer'ları değiştirmez. V0.8.5 covariance-aware Huber loss
uygulamaz; multivariate robust loss covariance geometrisi bağımsız
doğrulandıktan sonraki bir araştırma konusudur.
