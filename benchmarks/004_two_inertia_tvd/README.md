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
- DOF sırası `[göbek, atalet halkası]` için mod şekilleri `[1, 1]` ve
  `[1, -J_h/J_r]` olur.
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
alanına aktarılmaz. Bu topoloji regresyonu matris assembly veya yeni bir modal
çözüm içermez; aşağıdaki analitik sonuçlar değişmeden kalır.

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)
- [Matematik modeli](../../docs/mathematics/two_inertia_modal_model.md)

Sıfır olmayan analitik sonuçlar ve ölçekleme ilişkileri `1e-10` bağıl hata,
sıfır modu ile normalize mod bileşenleri `1e-12` mutlak hata sınırıyla
`test_torsional_system.f90` tarafından doğrulanır.
