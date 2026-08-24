# Beklenen Sonuç

## Annüler kauçuk burç geometri faktörü

```text
Cθ = 4π(0.01)(0.02²)(0.05²)/(0.05²-0.02²)
   = 5.98398600683770e-5 m³
```

## Kompleks rijitlik bileşenleri

```text
K'  = 1 000 000 × Cθ
    = 59.8398600683770 N·m/rad

K'' = 100 000 × Cθ
    = 5.98398600683770 N·m/rad
```

## Kayıp faktörü

```text
G''/G' = 100 000 / 1 000 000 = 0.1
K''/K' = 5.98398600683770 / 59.8398600683770 = 0.1
tan(delta) = 0.1
```

Sonuç çalışma noktası `100 Hz` ve `293.15 K` değerlerini korur. Tüm sayısal
sonuçlar `test_dynamic_torsional_stiffness.f90` tarafından `1e-10` bağıl hata
sınırıyla doğrulanır.
