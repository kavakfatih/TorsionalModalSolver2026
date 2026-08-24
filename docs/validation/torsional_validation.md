# Torsional Sistem Analitik Doğrulaması

## Amaç ve doğrulama sınıfı

V0.3.1, V0.3.0 ile kurulan torsional node-element, DOF mapping ve global
mass/stiffness matrix assembly altyapısını değiştirmeden doğrulama kapsamını
güçlendirir. Amaç, üretim API'lerinin bağımsız analitik sonuçlar ve fiziksel
invariantlarla tutarlı çalıştığını göstermektir.

Bu belge ve `test_torsional_validation` bir **analytical code/equation
verification** çalışmasıdır. DMA, modal test tezgâhı veya saha ölçümleriyle
yapılan deneysel model validation çalışmasının yerini almaz.

## Matematiksel model

Sönümsüz ayrık torsional sistemin hareket denklemi:

```text
[M] theta'' + [K] theta = 0
```

İki uçlu lineer elemanın lokal rijitliği:

```text
          [ 1  -1 ]
K_e = k_e [       ]
          [-1   1 ]
```

Global rijitlik, her elemanın bağlantı matrisi `A_e` üzerinden toplanır:

```text
K = sum(A_e^T K_e A_e)
```

Düğümlerde yığılmış polar atalet modeli için:

```text
M(equation_i, equation_i) += J_i
```

| Büyüklük | Fiziksel anlam | SI birimi |
| --- | --- | --- |
| `k`, `K` | Torsional rijitlik | `N·m/rad` |
| `J`, `M` | Polar kütle ataleti | `kg·m²` |
| `theta` | Açısal yer değiştirme | `rad` |
| `omega` | Açısal frekans | `rad/s` |
| `f` | Doğal frekans | `Hz` |

## Varsayımlar ve kapsam

- Dönmeler küçüktür ve eleman davranışı lineerdir.
- Tüm düğümler ortak eksen ve pozitif dönme yönünü kullanır.
- Eleman kütlesiz, polar ataletler düğümlerde yığılmıştır.
- Düğüm ataletleri ve reel depolama rijitlikleri sonlu ve pozitiftir.
- `equation_id=0`, yalnız homojen `theta=0` kısıtının eliminasyon sentinel'idir.
- Serbest sistemlerde rijit-cisim modu korunur; indirgenmiş kısıtlı K için
  sıfır satır toplamı veya rijit-cisim modu beklenmez.
- Global damping/C matrix, K'', dış yük, nonlinear davranış, eigen solver,
  LAPACK, sparse depolama ve FEM bu doğrulamanın kapsamı dışındadır.

## Analitik referanslar

### Tek torsional eleman

`k=100 N·m/rad` için dört lokal katsayı doğrudan sınanır:

```text
K_e = [ 100  -100 ] N·m/rad
      [-100   100 ]
```

Pozitif köşegen ile negatif çapraz katsayılar işaret konvansiyonunu, matris
simetrisi ise karşılıklı lineer bağlantıyı doğrular.

### Rijit-cisim modu

Serbest iki düğüm birlikte döndüğünde elemanda bağıl dönme oluşmaz:

```text
theta_r = [1, 1]^T
K theta_r = 0
```

Bu kontrol yanlış eleman işareti, yanlış assembly konumu veya yanlış DOF
eşlemesini görünür kılar.

### Üç düğümlü global assembly

Fiziksel kimlikleri matris indislerinden farklı seçilen zincir:

```text
Node 30 -- k1=100 -- Node 10 -- k2=200 -- Node 70
```

Ekleme sırasındaki aktif denklem düzeninde analitik global rijitlik:

```text
K = [ 100  -100     0 ]
    [-100   300  -200 ] N·m/rad
    [   0  -200   200 ]
```

Merkezdeki `300` değeri iki lokal katkının toplandığını gösterir. Tamamen
serbest zincirde `K=K^T`, satır toplamları sıfır ve `[1,1,1]^T` null moddur.

### İki ataletli TVD

Serbest iki atalet ve tek torsional bağlantı için elastik mod:

```text
omega_e² = K (1/J1 + 1/J2)
f_e = omega_e/(2 pi)
phi_e = [1, -J1/J2]^T
```

Referans değerler:

```text
J1 = 0.10 kg·m²
J2 = 0.20 kg·m²
K  = 1000 N·m/rad
omega_e² = 15000 s^-2
f_e = 19.492420030841906 Hz
```

Bilinen `phi_e` üzerinde aşağıdaki generalized eigen residual ve Rayleigh
quotient kontrolleri yapılır:

```text
K phi_e - omega_e² M phi_e = 0
omega_R² = (phi_e^T K phi_e)/(phi_e^T M phi_e)
```

Test bilinmeyen bir özdeğer veya özvektör aramaz. Bu nedenle Rayleigh hesabı
yeni bir eigen solver değildir; bilinen analitik modun assembled M/K üzerinde
doğrulanmasıdır. Ayrıca `[1,1]^T` rijit mod ile `phi_e` için
`phi_r^T M phi_e=0` kütle ortogonalliği sınanır.

## DOF mapping doğrulaması

Fiziksel düğümler haritada korunur; yalnız aktif denklemler kesintisiz
numaralandırılır:

| Vaka | Fiziksel node sırası | Constraint | Beklenen equation ID | Aktif DOF |
| --- | --- | --- | --- | --- |
| A | `[30,10,70]` | Yok | `[1,2,3]` | 3 |
| B | `[30,10,70]` | Node 30 | `[0,1,2]` | 2 |

Bu kontroller `Physical Node ID -> Equation ID -> Matrix Index` ayrımını ve
kısıt sonrası denklem sürekliliğini kilitler.

## Matris kalite ölçütleri

Global K için aşağıdaki ölçütler kullanılır:

```text
Symmetry residual = ||K-K^T||_F
U = 1/2 theta^T K theta >= 0
```

Frobenius simetri normu, katsayıların karşılıklı konumlarda aynı olduğunu;
negatif olmayan `U` ise pasif lineer bağlantının yapay enerji üretmediğini
doğrular. Sıfır ve ortak dönme vektörleri sıfır enerjili olabilir; bağıl dönme
vektörlerinde enerji pozitiftir.

## Test senaryoları ve toleranslar

| Senaryo | Ölçüt | Tolerans |
| --- | --- | --- |
| Düğüm, equation ID ve boyut | Tam eşitlik | `0` |
| Lokal/global K katsayıları | Maksimum mutlak hata | `1e-10 N·m/rad` |
| Rijit-cisim kalıntısı | Maksimum mutlak hata | `1e-10 N·m` |
| K simetrisi | `||K-K^T||_F` | `1e-10 N·m/rad` |
| Elastik enerji | `U >= -tolerance` | `1e-10 N·m` |
| İki-atalet modal residual | Normalize maksimum hata | `1e-10` |
| Doğal frekans | Bağıl hata | `1e-10` |
| Kütle ortogonalliği | Mutlak iç çarpım | `1e-10 kg·m²` |

Mutlak toleranslar sıfır beklenen sonuçları güvenli biçimde sınar. Pozitif
frekans için ölçekten bağımsız bağıl tolerans kullanılır.

## Sonuç ve izlenebilirlik

Doğrulama hedefi `engine/tests/test_torsional_validation.f90` dosyasından
üretilen `tms26.torsional_validation` CTest testidir. 2026-08-25 tarihinde
GNU Fortran 16.2.0 ile yapılan temiz Debug doğrulamasında tam build ve bu hedef
başarılı olmuş, toplam `53/53` CTest geçmiştir.

İlgili teknik kaynaklar:

- [`local_torsional_element_matrix.md`](../mathematics/local_torsional_element_matrix.md)
- [`global_matrix_assembly.md`](../mathematics/global_matrix_assembly.md)
- [`two_inertia_modal_model.md`](../mathematics/two_inertia_modal_model.md)
- [`two_inertia_torsional_system.md`](../physics/two_inertia_torsional_system.md)
- [`Benchmark 004`](../../benchmarks/004_two_inertia_tvd/)

Bu kanıt, belirtilen denklemlerin ve assembly uygulamasının referans vakalarda
tutarlı olduğunu gösterir; genel eigen çözüm doğruluğu veya deneysel TVD
korelasyonu hakkında sonuç üretmez.
