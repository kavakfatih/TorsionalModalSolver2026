# Beklenen Sonuç

## Polar alan momenti

```text
Jp = π/2 (0.05⁴ - 0.02⁴)
   = 9.56614963018092e-6 m⁴
```

## Kompleks rijitlik bileşenleri

```text
K'  = 1 000 000 × Jp / 0.01
    = 956.614963018092 N·m/rad

K'' = 100 000 × Jp / 0.01
    = 95.6614963018092 N·m/rad
```

## Kayıp faktörü

```text
G''/G' = 100 000 / 1 000 000 = 0.1
K''/K' = 95.6614963018092 / 956.614963018092 = 0.1
tan(delta) = 0.1
```

Sonuç çalışma noktası `100 Hz` ve `293.15 K` değerlerini korur. Tüm sayısal
sonuçlar `test_dynamic_torsional_stiffness.f90` tarafından yüzde `0,1`'den
küçük bağıl hata sınırıyla doğrulanır.
