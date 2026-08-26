# V0.5 Modal Eigen Solver Doğrulaması

## Kapsam

Bu belge dense LAPACK DSYGV reference backend'i ile backend-neutral modal
katmanın analytical code verification kapsamını tanımlar. Deneysel modal test
ve malzeme model calibration çalışmasının yerine geçmez.

## Toleranslar

- Analitik pozitif eigenvalue ve frekans: `1e-10` bağıl hata
- Mass normalization ve sign-invariant mode correlation: `1e-10`
- Relative eigenpair residual: en çok `1e-10`
- Mass orthogonality Frobenius hatası: problem boyutuyla ölçeklenen `1e-10`
- Modal stiffness projection: `Phi^T K Phi ≈ diag(lambda)`, `1e-10` bağıl hata
- Input immutability: exact katsayı eşitliği
- Rigid-mode ayrımı: sabit Hz yerine machine epsilon ve spectral scale

## Analitik ve fizik regresyonları

| Model | Beklenen doğrulama |
|---|---|
| Fixed–k–J, 1 DOF | `lambda=k/J`, frekans, mass normalization, residual ve `[0,phi_r]` recovery |
| Free-free iki atalet | `lambda=[0,k(1/Jh+1/Jr)]`, analitik solver cross-check ve sign-invariant mode |
| Eş üç düğümlü zincir | `[0,k/J,3k/J]`, bir rigid ve iki elastic mode |
| İlk düğümü fixed üçlü zincir | `lambda=(k/J)(3±sqrt(5))/2`, frekanslar ve constrained physical bileşenin sıfır olması |
| İki identical oscillator | Repeated multiplicity, eigenspace, residual ve M-orthogonality |
| İki ayrık free-free alt sistem | Birden fazla `RIGID_MODE` |
| Fully constrained sistem | LAPACK çağrısından önce temiz `no active DOF` tanısı |

## Generalized eigenproblem contract testleri

Synthetic K/M matrisleri, production fizik tiplerine unsafe setter eklenmeden
`generalized_eigen_problem_t` factory üzerinden sınanır:

- kare olmayan ve eş boyutlu olmayan K/M,
- nonsymmetric K veya M,
- NaN, pozitif sonsuz ve negatif sonsuz katsayı,
- sıfır, negatif veya indefinite/non-SPD M,
- singular symmetric positive-semidefinite K'nin kabulü,
- anlamlı negatif eigenvalue'nun modal katmanda reddi,
- automatic tolerance içindeki küçük negatif değerin rigid kabulü,
- DSYGV öncesi ve sonrası original K/M immutability.

Beklenen-hata testleri yalnız nonzero çıkışa bağlı değildir. Non-SPD M,
fully constrained sistem ve anlamlı negatif özdeğer vakalarında CMake
sarmalayıcısı beklenen TMS26 diagnostic metnini de doğrular. Böylece raw LAPACK
`INFO` değerinin açıklamasız dışarı sızması başarı sayılmaz.

## Repeated eigenspace doğrulaması

Repeated eigenvalue için tekil eigenvector beklenmez. Hesaplanan ve beklenen
M-orthonormal bazlar arasındaki overlap matrisi kullanılarak alt uzayın aynı
olduğu doğrulanır. Bu ölçüt eigenvector işaretinden ve repeated alt uzay içindeki
baz rotasyonundan bağımsızdır.

## Ölçek regresyonu

Aynı generalized problem farklı K ölçekleriyle çözülür. Her ölçekte rigid mode
sayısı korunmalı ve pozitif elastic mode sabit bir Hz eşiği nedeniyle rigid
sınıfına düşmemelidir.

## Sonuç

V0.5 exit criterion, bütün eski V0.1–V0.4 testleriyle yeni modal testlerin
macOS ve Windows üzerinde birlikte geçmesidir. DSYGV'yi gerçekten çağıran
testler aynı zamanda LAPACK sağlayıcısı, Fortran ABI ve runtime link doğrulaması
görevi görür.
