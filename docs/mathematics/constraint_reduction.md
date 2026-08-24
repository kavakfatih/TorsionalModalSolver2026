# Constraint Reduction Matematik Modeli

## Amaç

Bu belge, V0.4.0 tam torsional M/K sisteminin fixed veya prescribed Physical
DOF'lar uygulanarak aktif solver sistemine direct elimination ile
indirgenmesini tanımlar. Torsional fizik, eleman rijitliği ve atalet denklemleri
değişmez; yalnız koordinat uzayı dönüştürülür.

## Semboller ve birimler

- `n`: Physical DOF ve tam Equation ID sayısı, boyutsuz
- `n_r`: aktif Equation ID sayısı, boyutsuz
- `q`: tam torsional koordinat vektörü, `rad`
- `q_r`: aktif koordinat vektörü, `rad`
- `q_p`: prescribed bileşenleri taşıyan tam vektör, `rad`
- `K`: tam torsional rijitlik matrisi, `N·m/rad`
- `M`: tam polar kütle ataleti matrisi, `kg·m²`
- `P`: `n x n_r` boyutlu seçim matrisi, boyutsuz
- `Kr`, `Mr`: indirgenmiş rijitlik ve atalet matrisleri

Her Physical DOF `(node_id,dof_type)` ile tanımlanır. Tam Equation ID değerleri
constraint'ten bağımsız `1..n` aralığındadır. Aktif Equation ID değerleri,
constraint uygulanmamış Physical DOF'lar için `1..n_r` aralığında ayrı
oluşturulur.

## Seçim matrisi

Aktif koordinat `a`, tam denklem `i(a)` ile eşleşiyorsa seçim matrisi:

```text
P(i(a),a) = 1
P(i,a)    = 0, diğer tüm konumlarda
```

biçimindedir. Her aktif sütunda tam olarak bir bir bulunur; constraint altındaki
tam denklemlerin satırları sıfırdır. Tam durum vektörü:

```text
q = P q_r + q_p
```

ile elde edilir. `q_p`, aktif konumlarda sıfır; constraint konumlarında ilgili
prescribed değerdir.

## Direct elimination

Tam lineer hareket denklemi:

```text
M q'' + K q = f
```

olsun. Homojen constraint veya sabit prescribed değer çevresindeki homojen
perturbasyon için `q = P q_r` yazılıp denklem soldan `transpose(P)` ile
çarpılır:

```text
Mr q_r'' + Kr q_r = f_r

Kr = transpose(P) K P
Mr = transpose(P) M P
```

`P` bir seçim matrisi olduğundan eşdeğer indeks gösterimi:

```text
Kr = K(active,active)
Mr = M(active,active)
```

biçimindedir. Bu işlem yeni katsayı türetmez; tam matrislerin aktif satır ve
sütunlarını sıralı olarak seçer. K ve M simetrikse Kr ve Mr de simetriktir.

## Sıfırdan farklı prescribed değer

Genel geri kurma bağıntısı:

```text
q = P q_r + q_p
```

olarak korunur. Sıfırdan farklı `q_p`, denklem türüne göre sağ tarafta ek
terimler doğurur. Sabit bir prescribed değer için elastik düzeltme örneği:

```text
f_r,corrected = transpose(P) (f - K q_p)
```

Zamana bağlı prescribed değer ayrıca `M q_p''` ve varsa damping terimlerini
gerektirir. V0.4.0 yük vektörü veya zaman çözümü uygulamadığı için bu RHS
düzeltmelerini hesaplamaz. Prescribed değer yalnız veri doğrulaması ve tam
durum recovery altyapısında korunur.

## Tüm DOF'ların constraint altında olması

`n_r = 0` ise:

```text
P  : n x 0
Kr : 0 x 0
Mr : 0 x 0
```

geçerli sonuçtur. Bu durum hata değildir; sistemde solver'a aktarılacak aktif
koordinat olmadığını gösterir.

## Node sırasından bağımsızlık

Fiziksel node ekleme sırası değiştiğinde Equation ID değerleri ve matris satır
sırası değişebilir. Bir permütasyon matrisi `Q` için:

```text
K_new = transpose(Q) K Q
M_new = transpose(Q) M Q
```

olur. Constraint ve aktif eşleme aynı Physical DOF kimlikleriyle uygulandığında
indirgenmiş sistemler de aynı permütasyon ilişkisini taşır. Fiziksel sonuçlar
node ekleme sırasına bağlı değildir; doğrudan katsayı karşılaştırması yapılacaksa
önce ortak Physical DOF sırasına alınmalıdır.

## Result recovery ve gelecekteki modal denklem

Genel durum recovery bağıntısı `q = P q_r + q_p` biçimindedir. Serbest titreşim
özvektörü, prescribed denge durumu çevresindeki homojen perturbasyonu temsil
ettiği için:

```text
phi = P phi_r
```

olarak geri kazanılır. Gelecekteki indirgenmiş modal problem:

```text
Kr phi_r = lambda Mr phi_r
```

olacaktır. V0.4.0 yalnız Kr/Mr ve recovery eşlemesini üretir; `lambda` veya
`phi_r` çözmez.

## Matris backend sınırı

Direct elimination, yalnız matris boyutu ve katsayı erişimi gerektirir. Bu
nedenle üst seviye yordam `reduce_matrix` olarak ifade edilir; dense depolama
adını API'ye taşımaz. İlk gerçekleştirim mevcut dense matris türünü kullanır.
Gelecekte sparse/CSR backend eklenmesi matematik modelini değiştirmez.

## Geçerlilik ve kapsam

- Equation ID ve aktif Equation ID eşlemeleri benzersiz ve kesintisizdir.
- K ve M kare, aynı tam Equation ID sırasına sahip ve sonlu olmalıdır.
- Constraint hedefi mevcut bir Physical DOF olmalıdır.
- Aynı Physical DOF için birden fazla constraint geçersizdir.
- V0.4.0 yalnız direct elimination uygular.
- MPC, Lagrange multiplier, penalty, contact, RHS çözümü ve eigen solver kapsam
  dışıdır.

Mimari sorumluluklar
[`../architecture/V0.4_constraint_foundation.md`](../architecture/V0.4_constraint_foundation.md),
kalıcı tasarım kararı ise
[`../decisions/0009-constraint-reduction-architecture.md`](../decisions/0009-constraint-reduction-architecture.md)
belgelerinde açıklanır.
