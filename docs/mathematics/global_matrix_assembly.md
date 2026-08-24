# Global Torsional Matris Assembly

## Amaç

V0.3.0, genel torsional topolojiyi gelecekteki modal denkleme girdi olacak
dense global dönel atalet ve rijitlik matrislerine dönüştürür:

```text
[M] theta'' + [K] theta = 0
```

Bu sürüm matrisleri oluşturur; özdeğerleri çözmez.

## Fiziksel düğüm ve denklem kimliği

`torsional_node_t%id`, fiziksel topolojinin pozitif etiketidir. Bu değer bir
dizi indisi veya global matris satırı değildir. Örneğin:

```text
Fiziksel node ID: 10  20  30
Aktif equation ID:  1   2   3
```

`dof_map_t`, iki uzay arasındaki eşlemeyi açıkça taşır. Serbest düğümler sisteme
eklenme sırasıyla kesintisiz `1..n_active` denklem kimliği alır. Homojen sıfır
dönme ile kısıtlanmış düğüm haritada korunur fakat `equation_id=0` alır.

Sıfır kimliği yalnız kısıt işaretidir. Haritada bulunmayan fiziksel düğüm hata
üretir; böylece “kısıtlı” ile “eksik” durumları birbirine karışmaz. Tamamen
kısıtlı geçerli sistem `n_active=0` ve 0x0 global matrisler üretir.

Harita fiziksel düğüm kümesini izlenebilirlik için tam tutar, ancak yalnız
serbest düğümlerin aktif denklem uzayını `1..n_active` olarak numaralandırır.
Gelecekteki modal reduction, aktif denklem uzayı ile indirgenmiş modal baz
arasında ayrı bir dönüşüm katmanı uygulayabilir.

## Dense matris sözleşmesi

İlk gerçekleştirim `dense_matrix_t` kullanır. Katsayı depolaması private ve
allocatable'dır; boyut, indeks ve sonlu sayı kontrolleri yordamlarla yapılır.
Global stiffness/mass API'leri ham depolamaya bağlı olmadığından ileride sparse
bir gerçekleştirim eklenebilir.

Genel dense taşıyıcı fiziksel birim kodlamaz. Üst katman sözleşmeleri:

| Matris | Fiziksel anlam | Katsayı birimi |
| --- | --- | --- |
| `K` | Torsional rijitlik | `N·m/rad` |
| `M` | Polar kütle ataleti | `kg·m²` |

## Lokal rijitlik katkısı

Elemanın yerel koordinat sırası `[theta_i, theta_j]` ve bağıl açısı
`theta_i-theta_j` olduğunda:

```text
          [ 1  -1 ]
K_e = k_e [       ]
          [-1   1 ]
```

Türetim ve enerji invariantları
[`local_torsional_element_matrix.md`](local_torsional_element_matrix.md)
belgesindedir. `get_local_stiffness`, mevcut üretim hesabını standart eleman
arayüzü üzerinden döndürür.

## Lokalden globale toplama

Her elemanın uç düğümleri önce DOF haritasına gönderilir:

```text
[node_i_id, node_j_id]
          ↓ DOF map
[equation_i, equation_j]
          ↓ scatter-add
global K
```

Yerel `a,b` katsayısı aşağıdaki global konuma eklenir:

```text
K(equation_a,equation_b) += K_e(a,b)
```

Assembly toplama yapar; mevcut katsayının üzerine yazmaz. Bu nedenle aynı
düğümde buluşan veya paralel elemanların katkıları doğal olarak toplanır.

`equation_id=0` olan yerel satır veya sütun atlanır. Bu işlem homojen
`theta_constrained=0` sınır koşulu için indirgenmiş aktif K matrisini verir.
Sıfırdan farklı prescribed dönmede gerekli sağ taraf düzeltmesi bu sürümde yoktur.

## Üç düğümlü zincir

```text
Node 10 -- k1 -- Node 20 -- k2 -- Node 30
```

İki lokal katkının toplanmasıyla:

```text
    [ k1       -k1        0  ]
K = [-k1    k1 + k2     -k2 ]
    [  0       -k2       k2 ]
```

Tamamen serbest sistemde K:

- simetriktir,
- pozitif yarı-tanımlıdır,
- her satırda sıfır toplama sahiptir,
- `[1,1,1]^T` rijit-cisim dönme vektörünü null mod olarak taşır.

Bu son iki özellik homojen kısıt eliminasyonundan sonraki indirgenmiş K için
genel olarak beklenmez.

## Global dönel atalet matrisi

TMS26'nın ilk genel modeli düğümde yığılmış polar atalet kullanır:

```text
M(equation_i,equation_i) += J_i
```

Üç serbest düğüm için:

```text
    [J1  0   0 ]
M = [0   J2  0 ]
    [0   0   J3]
```

`J_i` birimi `kg·m²`'dir. Köşegen dışı katsayılar sıfırdır. Kısıtlı düğümün
ataleti indirgenmiş aktif M matrisine eklenmez. Eleman ataleti veya tutarlı mass
matrix bu sürümde uygulanmaz.

## Varsayımlar ve kapsam sınırı

- Her düğüm ortak eksen çevresinde bir torsional DOF taşır.
- Düğümler yığılmış, sonlu ve pozitif polar atalete sahiptir.
- Elemanlar lineer, kütlesiz ve aynı pozitif açı yönünü kullanır.
- Global koordinat dönüşüm matrisi gerekmez.
- Global K yalnız reel depolama rijitliği omurgasını kullanır; K'' ve viskoz c
  bu matrislere katılmaz.
- Global C assembly, dış yük/RHS, sıfırdan farklı prescribed dönme, sparse
  depolama, LAPACK ve eigen çözümü kapsam dışıdır.
