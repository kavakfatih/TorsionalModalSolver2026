# Benchmark 004: İki Ataletli TVD Sistemi

Bu benchmark, tek torsional elastomer bağlantısıyla bağlanmış göbek ve atalet
halkasının fixed-hub ve serbest-serbest analitik modal sonuçlarını doğrular.

```text
Hub (J_h)
   |
   | K' — tam bağlı annüler elastomerin depolama rijitliği
   |
Inertia Ring (J_r)
```

## Doğrulanan davranışlar

- Fixed-hub sonucu mevcut tek-DOF doğal frekans yordamıyla aynıdır.
- Serbest-serbest sistem bir sıfır frekanslı rijit-cisim modu üretir.
- Elastik bağıl mod frekansı analitik iki ataletli referansla uyuşur.
- DOF sırası `[göbek, atalet halkası]` için göbek bileşeni bire ölçeklenen
  analitik mod şekilleri `[1, 1]` ve `[1, -J_h/J_r]` olur.
- `K'` dört katına çıktığında elastik frekans iki katına; iki atalet birlikte
  dört katına çıktığında elastik frekans yarıya iner.
- `J_h` çok büyüdüğünde serbest-serbest elastik mod fixed-hub sonucuna yaklaşır.

Bu benchmark yalnız `K'` kullanan frozen-property sönümsüz modal tahmindir.
`K''` sistem üstverisinde tutulabilir fakat sönüm veya kompleks özdeğer hesabına
katılmaz.

## V0.2.3 genel topoloji gösterimi

Aynı referans fizik, çözüm formülü değiştirilmeden genel sistem veri modeline
şu şekilde aktarılır:

```text
Node 1: J_h = 0.10 kg·m², serbest veya fixed-hub için constrained
Node 2: J_r = 0.20 kg·m², serbest
Element 1: node 1–2, K = K' = 1000 N·m/rad, c = 0 N·m·s/rad
```

Serbest-serbest gösterim `2`, fixed-hub gösterim `1` aktif DOF taşır. K''
kayıp rijitliği viskoz `c` ile aynı birimde olmadığından genel elemanın damping
alanına aktarılmaz. V0.2.3 topoloji regresyonu kendi başına matris assembly veya
modal çözüm içermez; aşağıdaki analitik sonuçlar sonraki solver sürümleri için
bağımsız oracle olarak korunur.

## V0.5.0 generalized modal solver eşdeğerliği

V0.5 yolu aynı modeli full K/M assembly, constraint reduction ve
`K_r phi=lambda M_r phi` çözümünden geçirir. Dense LP64 LAPACK `DSYGV`
reference backend'i analitik frekanslarla aynı eigenvalue'ları üretmelidir.

Analitik çözüm, karşılaştırmayı kolaylaştırmak için göbek genliğini `1` alır.
DSYGV ise varsayılan olarak:

```text
phi_i^T M phi_i = 1
```

mass normalization kullanır. Bu nedenle bileşenlerin sayısal büyüklükleri
farklıdır; iki sonuç ortak bir skalerle aynı fiziksel mode shape'i temsil eder.
Eigenvector işareti de arbitrary olduğundan `phi` ve `-phi` eşdeğerdir.

Doğrulama doğrudan bileşen eşitliği yerine şu invariantları kullanır:

- eigenvalue ve frekans eşitliği,
- sign-invariant modal correlation,
- `phi^T M phi=1` mass normalization,
- rigid ve elastic modlar arasında M-orthogonality,
- boyutsuz relative eigenpair residual,
- fixed-hub çözümünde constrained göbek bileşeninin physical recovery sonrası
  sıfır olması.

Bu çözüm lineer, sönümsüz ve frozen-property kapsamındadır. `G'`/K' çözüm
boyunca sabittir; damping, K'' veya öz-tutarlı frekans-malzame iterasyonu
uygulanmaz.

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)
- [Matematik modeli](../../docs/mathematics/two_inertia_modal_model.md)
- [Generalized modal matematiği](../../docs/mathematics/generalized_modal_eigenproblem.md)
- [V0.5 doğrulama kapsamı](../../docs/validation/modal_eigen_validation.md)

Sıfır olmayan analitik sonuçlar ve ölçekleme ilişkileri `1e-10` bağıl hata,
sıfır modu ile normalize mod bileşenleri `1e-12` mutlak hata sınırıyla
`test_torsional_system.f90` tarafından doğrulanır. Generalized solver
eşdeğerliği `test_generalized_eigen_solver.f90` ve `test_modal_analysis.f90`
ile residual/orthogonality toleransları altında ayrıca doğrulanır.
