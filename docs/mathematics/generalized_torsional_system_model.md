# Genelleştirilmiş Torsional Sistem Matematik Modeli

## Ayrık koordinatlar

Her serbest torsional düğüm için bir açısal koordinat tanımlanır:

```text
theta = [theta_1, theta_2, ..., theta_n]^T
```

Sabitlenmiş düğümler topolojide korunur fakat aktif koordinat vektörüne
katılmaz. Düğüm kimlikleri topolojik etiketlerdir; V0.3.0 DOF haritası fiziksel
kimlikleri ekleme sırasındaki aktif `1..n` denklem kimliklerine dönüştürür.
Kısıtlı düğümler haritada `equation_id=0` ile korunur.

## Düğüm ataletinin enerji anlamı

`i` düğümündeki yığılmış polar atalet `J_i` için kinetik enerji:

```text
T_i = 1/2 J_i theta_i'^2
```

Global dönel atalet matrix assembly işleminde düğümün katkısı:

```text
M_ii += J_i
```

V0.3.0 aktif düğümlerin bu katkılarını diagonal global M matrisinde toplar.
Kısıtlı düğümün ataleti indirgenmiş aktif matrise eklenmez.

## Eleman rijitliği ve sönümü

`i` ve `j` düğümlerini bağlayan lineer elemanın bağıl açısı ve elastik
potansiyel enerjisi:

```text
delta_theta_e = theta_i - theta_j
U_e = 1/2 k_e delta_theta_e^2
```

Elemanın iki ucundaki geri çağırıcı momentler eşit büyüklükte ve zıt işaretlidir.
Elemanın ürettiği ve V0.3.0 global rijitlik assembly işlemine girdi olan yerel
katkı:

```text
          [ 1  -1 ]
K_e = k_e [       ]
          [-1   1 ]
```

Eşdeğer viskoz sönüm için bağıl hız ve kavramsal yerel katkı:

```text
delta_theta_e' = theta_i' - theta_j'
T_c = c_e delta_theta_e'

c_e [ 1  -1 ]
    [-1   1 ]
```

Rijitlik için lokal ve global K katkıları uygulanmıştır. Viskoz sönüm ifadesi
gelecekteki assembly arayüzünü tanımlar; V0.3.0 lokal veya global C matrisi
üretmez.

## Gelecekteki hareket denklemi

Sönümsüz, dış zorlamasız hedef formülasyon:

```text
[M] theta'' + [K] theta = 0
```

İleride eşdeğer viskoz sönüm açıkça etkinleştirilirse genel biçim:

```text
[M] theta'' + [C] theta' + [K] theta = 0
```

V0.3.0 veri topolojisini, aktif DOF eşlemesini ve dense global M/K assembly
katmanını kurar. Homojen sıfır dönme kısıtları indirgenmiş denklem matrisinde
uygulanır; sıfırdan farklı prescribed dönme ve özdeğer çözümü kapsam dışıdır.
Lokal matrisin ayrıntılı türetimi
[`local_torsional_element_matrix.md`](local_torsional_element_matrix.md)
belgesinde, global toplama ise
[`global_matrix_assembly.md`](global_matrix_assembly.md) belgesindedir.

## Kompleks kayıp rijitliği ile viskoz sönüm ayrımı

Mevcut dinamik elastomer modeli:

```text
K* = K' + iK''
```

Genel elemandaki `c`, bağıl açısal hızla çarpılan viskoz katsayıdır. Birimler:

```text
K'' : N·m/rad
c   : N·m·s/rad
```

Dolayısıyla `K''` doğrudan `c` alanına atanamaz. Frekansa bağlı eşdeğerlik
kurmak ek model seçimi gerektirir ve bu sürümde uygulanmaz. İki ataletli
topoloji köprüsü yalnız konservatif omurgayı taşır: `K = K'`, `c = 0`.

## İki ataletli özel durum

Benchmark 004 sistemi genel topolojide iki düğüm ve bir elemandır:

```text
J_1 = J_h = 0.10 kg·m²
J_2 = J_r = 0.20 kg·m²
K_1 = K'  = 1000 N·m/rad
c_1 = 0 N·m·s/rad
```

İki düğüm serbestse aktif DOF sayısı `2`; göbek düğümü sabitlenirse `1` olur.
Mevcut analitik frekans ve mod şekli çözümü
[`two_inertia_modal_model.md`](two_inertia_modal_model.md) belgesinde kalır.
Genel topoloji bu sonuçları yeniden hesaplamaz.

## SI birimleri ve geçerlilik

| Büyüklük | Sembol | SI birimi |
| --- | --- | --- |
| Açısal koordinat | `theta` | `rad` |
| Açısal hız | `theta'` | `rad/s` |
| Polar kütle ataleti | `J` | `kg·m²` |
| Torsional rijitlik | `K` | `N·m/rad` |
| Viskoz torsional sönüm | `c` | `N·m·s/rad` |

Model lineer, küçük açılı, yığılmış ataletli ve kütlesiz bağlantılı
ayrık-parametre yaklaşımıdır. Başlangıç açısı `theta(0)` veri olarak saklanır,
fakat bu sürümde denge veya zaman integrasyonu hesabına girmez.
