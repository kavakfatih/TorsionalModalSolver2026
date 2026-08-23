# Benchmark 002: Dinamik Elastomer

Bu referans veri kümesi, V0.2.0 dinamik elastomer veri modelinin tek bir EPDM
çalışma noktasını ve kayıp faktörü hesabını doğrular. Benchmark bir solver,
kompleks burulma rijitliği veya performans zamanlaması içermez.

## Veri akışı

```text
EPDM çalışma noktası
        ↓
MPa → Pa dönüşümü
        ↓
G' ve G'' saklama
        ↓
tan(delta) = G'' / G'
```

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)

Bu sabit örnek, ileride DMA veya modal test veri içe aktarma katmanları
eklenirken temel regresyon referansı olarak kullanılabilir.
