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

## V0.5 generalized eigenproblem referansı

DSYGV backend'i serbest-serbest sistem için artan sırada şu özdeğerleri
üretmelidir:

```text
lambda = [0, 15000] 1/s^2
```

`phi_i^T M phi_i = 1` mass normalization uygulandığında beklenen mode
shape'ler, global `+/-` işaret belirsizliği dışında, şöyledir:

```text
phi_rigid   = [1/sqrt(0.30), 1/sqrt(0.30)]
            = [1.825741858350554, 1.825741858350554]

phi_elastic = [1/sqrt(0.15), -0.5/sqrt(0.15)]
            = [2.581988897471611, -1.290994448735806]
```

Fixed-hub reduction sonrasında tek aktif atalet `J_r` olduğundan:

```text
lambda_fixed = 5000 1/s^2
phi_reduced  = [1/sqrt(0.20)]
             = [2.236067977499790]
phi_physical = [0, 2.236067977499790]
```

Eigenvector işareti keyfidir; yukarıdaki her vektörün negatifi aynı fiziksel
mode'u temsil eder. Karşılaştırmalar bu nedenle sign-invariant modal
correlation ve mass normalization invariantlarıyla yapılır.

## Ölçekleme ve limit regresyonları

- `K' -> 4K'` olduğunda `f_e -> 2f_e`.
- `J_h -> 4J_h` ve `J_r -> 4J_r` olduğunda `f_e -> f_e/2`.
- `J_h/J_r -> sonsuz` limitinde serbest-serbest `f_e`, fixed-hub
  `f_fixed` değerine yaklaşır.

`K''`, kayıp faktörü, malzeme referans frekansı ve sıcaklığı bu sönümsüz
analitik sonuçları değiştirmez.
