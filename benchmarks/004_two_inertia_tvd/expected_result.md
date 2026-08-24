# Beklenen Sonuç

## Fixed-hub doğal frekansı

```text
omega_fixed = sqrt(K'/J_r)
            = sqrt(1000/0.20)
            = 70.710678118654755 rad/s

f_fixed = omega_fixed/(2*pi)
        = 11.253953951963826 Hz
```

## Serbest-serbest rijit-cisim modu

```text
omega_0 = 0 rad/s
f_0 = 0 Hz
phi_0 = [1, 1]
```

İki rijit gövde birlikte döndüğü için elastomerin bağıl dönmesi ve geri
çağırıcı momenti yoktur.

## Serbest-serbest elastik mod

```text
omega_e = sqrt(K'*(1/J_h + 1/J_r))
        = sqrt(1000*(1/0.10 + 1/0.20))
        = 122.474487139158910 rad/s

f_e = omega_e/(2*pi)
    = 19.492420030841906 Hz
```

Göbek genliği bire normalize edildiğinde:

```text
phi_e = [1, -J_h/J_r]
      = [1, -0.5]
```

## Ölçekleme ve limit regresyonları

- `K' -> 4K'` olduğunda `f_e -> 2f_e`.
- `J_h -> 4J_h` ve `J_r -> 4J_r` olduğunda `f_e -> f_e/2`.
- `J_h/J_r -> sonsuz` limitinde serbest-serbest `f_e`, fixed-hub
  `f_fixed` değerine yaklaşır.

`K''`, kayıp faktörü, malzeme referans frekansı ve sıcaklığı bu sönümsüz
analitik sonuçları değiştirmez.
