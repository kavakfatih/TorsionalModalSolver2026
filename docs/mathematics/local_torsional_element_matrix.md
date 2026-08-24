# Lokal Torsional Eleman Rijitlik Matrisi

## Model ve yerel koordinatlar

İki düğümlü, lineer ve kütlesiz torsional elemanın yerel koordinat sırası:

```text
theta_e = [theta_i, theta_j]^T
```

Her iki açı aynı pozitif dönme yönünü kullanır. Elemanın bağıl dönmesi:

```text
delta_theta = theta_i - theta_j
```

Burada `theta_i` ve `theta_j` radyan, eleman rijitliği `k` ise `N·m/rad`
birimindedir.

## Enerjiden türetim

Lineer elastik elemanda depolanan potansiyel enerji:

```text
U_e = 1/2 k (theta_i - theta_j)^2
```

Radyan SI sisteminde boyutsuz kabul edildiğinden `U_e` enerji değeri joule
(`J`) birimindedir.

Yerel iç moment katkısı, enerjinin yerel koordinatlara göre gradyanıdır:

```text
tau_e = dU_e/dtheta_e = K_e theta_e
```

İkinci türev, lokal rijitlik matrisini verir:

```text
          [ 1  -1 ]
K_e = k * [       ]
          [-1   1 ]
```

Matrisin her katsayısı `N·m/rad`, `K_e theta_e` sonucu ise `N·m` birimindedir.
Hareket denklemindeki iç moment katkısı `+K_e theta_e`, elemanın düğümlere
uyguladığı fiziksel geri çağırıcı moment bunun negatifidir.

## Fiziksel invariantlar

### Simetri

```text
K_e(1,2) = K_e(2,1)
```

Simetri, lineer elastik bağlantının karşılıklılık ve konservatif enerji
özelliğini ifade eder.

### Sıfır satır toplamı

Her satırın toplamı sıfırdır. Bu nedenle ortak dönme vektörü:

```text
theta_r = [1, 1]^T
```

iç moment veya elastik enerji üretmez. Bu vektör yerel rijit-cisim null modudur.

### Pozitif yarı-tanımlılık

Her yerel açı vektörü için:

```text
theta_e^T K_e theta_e = k (theta_i - theta_j)^2 >= 0
```

`k > 0` olduğunda özdeğerler `0` ve `2k` olur. Sıfır rijit-cisim özdeğeri
nedeniyle matris pozitif tanımlı değil, pozitif yarı-tanımlıdır. Negatif
rijitlik pasif eleman modelini ve enerji koşulunu bozduğu için reddedilir.

## Uygulama ve gelecek assembly bağlantısı

`calculate_local_stiffness`, doğrulanmış `torsional_element_t` değerinden bu
2x2 katkıyı saf bir hesapla üretir. `local_matrix_2x2` yalnız katsayıları taşır;
birim, üretici yordamın fiziksel sözleşmesinden gelir.

Gelecekteki global assembly, yerel `[i,j]` uç kimliklerini aktif global DOF
indekslerine eşleyip dört katsayıyı global K matrisine ekleyecektir. V0.2.4:

- global K assembly,
- polar ataletlerden M matrisi,
- sınır koşulu eliminasyonu,
- özdeğer çözümü

uygulamaz.
