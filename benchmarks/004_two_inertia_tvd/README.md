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

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)
- [Matematik modeli](../../docs/mathematics/two_inertia_modal_model.md)

Sıfır olmayan analitik sonuçlar ve ölçekleme ilişkileri `1e-10` bağıl hata,
sıfır modu ile normalize mod bileşenleri `1e-12` mutlak hata sınırıyla
`test_torsional_system.f90` tarafından doğrulanır.
