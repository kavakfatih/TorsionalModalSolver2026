# Harmonik Torsional Cevap Matematiği

## Harmonik konvansiyon ve genlik

TMS26 V0.6.0 aşağıdaki zaman konvansiyonunu kullanır:

\[
\theta(t)=\operatorname{Re}\{\hat{\theta}\,e^{+i\omega t}\},\qquad
T(t)=\operatorname{Re}\{\hat{T}\,e^{+i\omega t}\}.
\]

`hat` ile gösterilen kompleks büyüklükler **peak amplitude** değerleridir;
RMS değildir. Buna göre:

\[
\dot{\theta}\rightarrow i\omega\hat{\theta},\qquad
\ddot{\theta}\rightarrow-\omega^2\hat{\theta},\qquad
\omega=2\pi f.
\]

Frekans `f` `[Hz]`, açısal frekans `omega` `[rad/s]` birimindedir. Harmonic
cycle tanımı ve kayıp yorumu nedeniyle solver noktaları sonlu ve `f>0` olmak
zorundadır.

## Reduced hareket denklemi

Storage stiffness `K'_r`, loss stiffness `K''_r`, viscous damping `C_r` ve
polar inertia `M_r` için reduced denklem:

\[
\left[K'_r-\omega^2M_r+i\left(K''_r+\omega C_r\right)\right]
\hat{\theta}_r=\hat{T}_r
\]

ve dynamic stiffness:

\[
Z_r(\omega)=K'_r-\omega^2M_r+i(K''_r+\omega C_r)
\]

olarak tanımlanır.

| Büyüklük | Fiziksel anlam | SI birimi |
| --- | --- | --- |
| `K'` | Depolama/geri çağırıcı torsional rijitlik | `N·m/rad` |
| `K''` | Çevrimsel kayıp torsional rijitliği | `N·m/rad` |
| `C` | Viskoz torsional sönüm matrisi | `N·m·s/rad` |
| `M` | Polar kütle ataleti matrisi | `kg·m²` |
| `theta_hat` | Kompleks açısal yer değiştirme | `rad` |
| `T_hat` | Kompleks peak torque | `N·m` |

Radyan SI cebrinde boyutsuz kabul edildiğinde bütün `Z theta_hat` terimleri
`N·m` birimindedir.

## Eleman katkıları

İki düğümlü lineer torsional elemanda yerel sıra
`[theta_i,theta_j]` ve bağlantı operatörü:

\[
L_e=\begin{bmatrix}1&-1\\-1&1\end{bmatrix}
\]

olmak üzere:

\[
K'_e=k'_eL_e,\qquad K''_e=k''_eL_e,\qquad C_e=c_eL_e.
\]

Üç matris simetriktir, satır toplamları sıfırdır ve ortak rigid rotation
`[1,1]^T` için iç katkı üretmez. Passive eleman varsayımı:

\[
k'_e>0,\qquad k''_e\ge0,\qquad c_e\ge0
\]

koşullarını gerektirir.

## Full assembly ve constraint reduction

Eleman scatter operatörü `A_e` ile:

\[
K'_{full}=\sum_e A_e^TK'_eA_e,
\quad K''_{full}=\sum_e A_e^TK''_eA_e,
\quad C_{full}=\sum_e A_e^TC_eA_e.
\]

Düğüm polar ataletleri `M_full` matrisine mevcut lumped-inertia yaklaşımıyla
eklenir. Direct-elimination seçim matrisi `P` için:

\[
K'_r=P^TK'_{full}P,\quad K''_r=P^TK''_{full}P,\quad
C_r=P^TC_{full}P,\quad M_r=P^TM_{full}P.
\]

Harmonic recovery:

\[
\hat{\theta}=P\hat{\theta}_r
\]

biçimindedir. Statik prescribed offset `q_p` bu perturbation cevabına
eklenmez.

## Kompleks simetri

Reel fizik matrisleri simetrik olduğunda:

\[
Z_r^T=Z_r.
\]

Fakat `K''_r+omega C_r` sıfır değilse çoğunlukla:

\[
Z_r^H\ne Z_r.
\]

Dolayısıyla problem complex symmetric sınıfındadır; Hermitian değildir.
Transpoz ile conjugate-transpose ayrımı solver seçiminin parçasıdır.

## Backend-independent relative residual

Çözülen her RHS için:

\[
r=Zx-b
\]

ve boyutsuz göreli kalıntı:

\[
\rho=\frac{\lVert Zx-b\rVert_2}
{\lVert Z\rVert_\infty\lVert x\rVert_2+\lVert b\rVert_2}
\]

kullanılır. Payda ile pay birlikte sayısal olarak sıfırsa `rho=0` tanımlanır.
Payda sıfır fakat pay sıfır değilse çözüm doğrulanmış kabul edilmez. Norm
hesapları overflow/underflow'a karşı ölçeklenmeli; `rho` sonlu ve negatif
olmamalıdır.

`RCOND`, `FERR` ve `BERR` LAPACK backend diagnostics değerleridir. `rho` ise
gelecekteki backend'lerle ortak TMS26 doğrulama ölçütüdür.

## Tek serbestlik dereceli analitik çözüm

Kayıp rijitliği time-domain'de anlık bir `iK''` terimi olarak yazılmaz. Frozen
complex-stiffness yaklaşımı doğrudan harmonic phasor denkleminde tanımlanır:

\[
\left[k'-m\omega^2+i(k''+\omega c)\right]\hat\theta=\hat T,
\qquad
\hat{\theta}=\frac{\hat{T}}
{k'-m\omega^2+i(k''+\omega c)}.
\]

Bu ifade viscous-only (`k''=0`), loss-only (`c=0`) ve combined damping
regresyonlarının bağımsız analitik referansıdır. Pozitif reel torque altında
`e^{+i omega t}` konvansiyonu passive sistem için cevabın fazını aşağıdaki
şekilde verir:

- rezonansın altında yaklaşık `0`,
- rezonans çevresinde yaklaşık `-pi/2`,
- rezonansın üstünde `-pi` yönünde.

## İki ataletli TVD referansları

Göbek fixed ve yalnız halka aktif ise:

\[
\hat{\theta}_r=\frac{\hat{T}}
{K'-\omega^2J_r+i(K''+\omega c)}.
\]

Free-free iki ataletli sistemde dengeli dış momentler
`[T_hat,-T_hat]^T` için eşdeğer atalet:

\[
J_{eq}=\frac{J_hJ_r}{J_h+J_r}
\]

ve relative angle:

\[
\Delta\hat{\theta}=\hat{\theta}_h-\hat{\theta}_r
=\frac{\hat{T}}
{K'-\omega^2J_{eq}+i(K''+\omega c)}.
\]

Finite `f>0` değerinde inertia terimi rigid-body koordinatını genel olarak
çözülebilir kılar. `f=0` bu harmonic solver'ın parçası değildir.

## Kinematik türevler ve phase

Kompleks açısal hız ve ivme:

\[
\hat{\Omega}=i\omega\hat{\theta},\qquad
\hat{\alpha}=-\omega^2\hat{\theta}
\]

olur. Kompleks bir sonuç için:

\[
|z|=\operatorname{abs}(z),\qquad
\varphi=\operatorname{atan2}(\operatorname{Im}z,\operatorname{Re}z).
\]

Core phase birimi radyandır. Phase unwrapping yapılmaz; tam sıfır genlikte
fazın fiziksel olarak tanımsız olduğu çağıran katmanda dikkate alınmalıdır.

## Harmonic response ve FRF

Genel `T_hat -> theta_hat` sonucu harmonic response'dur. Tek tanımlı ve sıfır
olmayan input torque channel için rotational FRF'ler:

\[
H_{\theta T}=\frac{\hat{\theta}}{\hat{T}}
\quad [rad/(N\,m)],
\]

\[
H_{\Omega T}=\frac{i\omega\hat{\theta}}{\hat{T}}
\quad [rad/(s\,N\,m)],
\]

\[
H_{\alpha T}=\frac{-\omega^2\hat{\theta}}{\hat{T}}
\quad [rad/(s^2\,N\,m)]
\]

olarak tanımlanır. Birden çok bağımsız eşzamanlı yükün toplam cevabını toplam
moment genliğine bölmek genel bir FRF tanımı değildir.

## Eleman relative angle, moment ve enerji kaybı

`node_i -> node_j` orientation için:

\[
\Delta\hat{\theta}_e=\hat{\theta}_i-\hat{\theta}_j,
\]

\[
\hat{T}_e=\left[k'_e+i(k''_e+\omega c_e)\right]
\Delta\hat{\theta}_e.
\]

Matrix-equilibrium iç kuvvet konvansiyonu `T_i=+T_e`, `T_j=-T_e` olarak
tanımlanır. Peak amplitude konvansiyonunda:

\[
P_{avg,e}=\frac{\omega}{2}(k''_e+\omega c_e)
|\Delta\hat{\theta}_e|^2\quad[W],
\]

\[
E_{cycle,e}=\pi(k''_e+\omega c_e)
|\Delta\hat{\theta}_e|^2\quad[J/cycle].
\]

Passive elemanda bu iki değer negatif olamaz. `k''=c=0` olduğunda her ikisi
de sıfırdır.

## Geçerlilik sınırları

Matrisler sweep boyunca sabit/frozen değerlerdir. V0.6; malzeme
interpolasyonu, nonlinear response, transient çözüm, mode superposition,
dynamic prescribed angle veya reaction torque çözmez. Constrained ve
sönümsüz bir sistemde `f -> 0+` limiti quasi-static compliance'a yaklaşabilir;
ancak `f=0` doğrudan solver girdisi değildir.
