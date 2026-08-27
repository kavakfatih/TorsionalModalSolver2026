# Karar 0011: Frekans Bölgesi Kompleks Çözücü

- Durum: Kabul edildi
- Tarih: 2026-08-27

## Bağlam

V0.5.0, lineer sönümsüz torsional modal problemi dense LP64 LAPACK DSYGV ile
çözmektedir. V0.6.0'ın amacı mode-superposition değil, sinusoidal torque altında
direct/full-order steady-state complex response üretmektir. Elastomer kayıp
rijitliği ile viskoz sönümün farklı fiziksel kanallar olarak korunması ve
singular frequency point'in bütün sweep'i sonlandırmaması gerekir.

## Karar

### Harmonic convention ve genlik

TMS26:

\[
\theta(t)=\operatorname{Re}\{\hat\theta e^{+i\omega t}\},\qquad
T(t)=\operatorname{Re}\{\hat T e^{+i\omega t}\}
\]

konvansiyonunu kullanacaktır. Kompleks genlikler peak amplitude olacak, RMS
olarak yorumlanmayacaktır.

### Dynamic stiffness

Reduced denklem:

\[
Z_r(\omega)\hat\theta_r=\hat T_r,
\qquad
Z_r=K'_r-\omega^2M_r+i(K''_r+\omega C_r)
\]

olacaktır. `K'' [N·m/rad]` ve `c [N·m·s/rad]` ayrı veri, matris ve assembly
kanalları olarak kalacaktır. Otomatik `c=K''` veya `c=K''/omega` dönüşümü
yapılmayacaktır.

### Matris sınıfı

`Z` complex symmetric (`Z^T=Z`) fakat genel olarak Hermitian olmayan
(`Z^H!=Z`) bir matristir. Hermitian `ZHESV` ve `ZHESVX` kullanılmayacaktır.

### Reference backend

Dense reference backend LAPACK `ZSYSVX` olacaktır:

- `FACT='N'`,
- `UPLO='U'`,
- explicit Fortran interface; `bind(C)` yok,
- `complex(dp)` ile double-complex ABI,
- default 32-bit integer ve `BLA_SIZEOF_INTEGER=4` LP64 sözleşmesi,
- `LWORK=-1` optimal workspace query,
- authoritative `A` ve `B` yerine internal çalışma kopyaları.

Tam mantıksal `Z`, residual ve diagnostics için korunacaktır; yalnız üst üçgen
seçimi backend faktörizasyon ayrıntısıdır.

### Katman sınırı

`harmonic_analysis`, ZSYSVX'i doğrudan çağırmayacaktır. Akış:

```text
harmonic_analysis
      |
complex_linear_solver facade
      |
LAPACK ZSYSVX backend
```

şeklinde olacaktır. Public harmonic API, LAPACK argumentlerini veya raw dense
storage'ı bilmeyecektir. Low-level contract birden çok RHS'ye açık tutulacak,
V0.6 public harmonic load case tek RHS ile başlayabilecektir.

### Status semantiği

Raw `INFO` public result'a sızdırılmayacaktır:

- `INFO=0` -> `SOLVED`,
- `1<=INFO<=N` -> `SINGULAR`; unique çözüm ve error bounds yoktur,
- `INFO=N+1` -> `SOLVED_ILL_CONDITIONED`; çözüm ve error bounds korunur,
- `INFO<0` -> backend/programming contract hatası.

Singular ve ill-conditioned noktalar analysis-state sonucudur. Singular nokta
sweep'i sonlandırmaz ve cevap vektörü uydurulmaz. `RCOND` yalnız numerical
conditioning diagnostic'idir; resonance detector olarak kullanılmaz.

### Diagnostics

ZSYSVX'in `RCOND`, per-RHS `FERR` ve `BERR` çıktıları korunacaktır. Ayrıca
backend-independent:

\[
\rho=\frac{\lVert Zx-b\rVert_2}
{\lVert Z\rVert_\infty\lVert x\rVert_2+\lVert b\rVert_2}
\]

relative residual değeri hesaplanacaktır.

### Constraint ve recovery

V0.4 direct-elimination mimarisi korunacaktır. Harmonic perturbation recovery:

\[
\hat\theta=P\hat\theta_r
\]

olacak; constrained bileşenler sıfır kalacaktır. Stored static prescribed
offset `q_p` eklenmeyecektir. Dynamic prescribed angle ve reaction recovery
V0.6 kapsamı dışındadır.

## Değerlendirilen alternatifler

### ZHESV / ZHESVX

Reddedildi; Hermitian problem içindir ve `Z` genel olarak Hermitian değildir.

### ZGESVX

Genel kompleks matrisleri çözebilir, fakat V0.6'nın complex-symmetric yapısını
doğrudan ifade etmez. Gelecekte doğrulama/fallback adayı olabilir; reference
backend olarak seçilmedi.

### Mode-superposition harmonic

V0.6 için reddedildi. Frequency-dependent gelecek malzeme gereksinimi ve
non-proportional `K''/C`, reel sönümsüz modal bazın damping'i her durumda
diagonalize etmesini engeller. Direct/full-order yol ilk reference çözüm olarak
daha açıktır.

### Sparse veya iterative complex solver

V0.6 foundation kapsamında değildir. Public facade gelecekte bu backend'lerin
eklenmesini engellemeyecek biçimde ayrılmıştır.

## Sonuçlar

- Mevcut `LAPACK::LAPACK` bağımlılığı yeterlidir; yeni external dependency
  eklenmez.
- V0.5 DSYGV modal yolu değişmeden kalır.
- Harmonic result sözleşmesi ve belgeleri çözümü linear, direct, full-order,
  frequency-domain, frozen-property ve peak-amplitude olarak açıkça tanımlar.
- Response peak, relative-angle peak, transmitted torque ve phase evolution
  mühendislik resonance yorumuna temel olabilir; düşük RCOND tek başına olamaz.
- Complex `dp` yalnız phasor/dynamic-stiffness cebrinde kullanılır. Semantic
  malzeme ve eleman verilerinde `K'`, `K''` ve `c` açık alanlar olarak kalır.

## Gelecek adayları

V0.6 sonrası ayrı benchmark ve karar gerektiren solver adayları:

- `ZSYSVXX`,
- rook-pivoting symmetric variant,
- `ZGESVX` general-complex fallback,
- sparse complex direct solver,
- iterative complex Krylov solver.

V0.7'de öncelikli fizik genişletmesi, frequency/temperature-dependent
`G'(f,T)` ve `G''(f,T)` provider'ının harmonic solver'a bağlanmasıdır.
