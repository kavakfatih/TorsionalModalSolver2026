# Beklenen Sonuç

## SI veri değerleri

```text
G'          = 1 000 000 Pa
G''         =   100 000 Pa
frequency   =       100 Hz
temperature =       293.15 K
```

## Kayıp faktörü

```text
tan(delta) = G'' / G'
           = 100 000 / 1 000 000
           = 0.1
```

Kayıp faktörü boyutsuzdur. `test_dynamic_modulus.f90`, bu değeri çift
hassasiyetli kayan nokta toleransı içinde doğrular. Bu benchmark için kompleks
burulma rijitliği veya doğal frekans sonucu beklenmez.
