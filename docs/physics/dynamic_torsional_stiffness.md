# Dinamik Burulma Rijitliği

## Amaç ve kapsam

Bu belge, V0.2.0 sırasında frekans ve sıcaklığa bağlı kompleks kayma modülünün
burulma rijitliğine aktarılması için hazırlanan modeli tanımlar. İlişkiler
V0.2.1 ile `tms_dynamic_torsional_stiffness` modülünde uygulanmıştır. Fiziksel
yorum ve doğal frekans bağlantısı
[`complex_torsional_stiffness.md`](complex_torsional_stiffness.md) belgesinde
ayrıntılandırılır.

## Kompleks burulma rijitliği

Lineer kayma idealizasyonunda kompleks burulma rijitliği şöyledir:

```text
K* = K' + iK''
K' = G'J/L
K'' = G''J/L
```

`K'`, elastik olarak enerji depolayan burulma rijitliği bileşenidir. `K''`,
çevrimsel enerji kaybıyla ilişkili burulma rijitliği bileşenidir. Her iki
bileşen de kullanılan G' ve G'' verilerinin frekans ve sıcaklığına bağlıdır.

## Fiziksel büyüklükler ve birimler

| Büyüklük | Fiziksel anlam | SI birimi |
| --- | --- | --- |
| `G'` | Depolama kayma modülü | `Pa = N/m²` |
| `G''` | Kayıp kayma modülü | `Pa = N/m²` |
| `J` | Elastomer kesitinin polar alan atalet momenti | `m⁴` |
| `L` | Elastomerin etkin burulma uzunluğu | `m` |
| `K'` | Depolama burulma rijitliği | `N·m/rad` |
| `K''` | Kayıp burulma rijitliği | `N·m/rad` |

Buradaki `J`, geometriye ait polar **alan** atalet momentidir; doğal frekans
denkleminde kullanılan `kg·m²` birimli polar **kütle** ataletiyle aynı büyüklük
değildir. Radyan SI boyut analizinde boyutsuz kabul edilse de mühendislik
anlamını açık tutmak için rijitlik birimi `N·m/rad` olarak yazılır.

## Varsayımlar ve sınırlar

- Elastomer lineer viskoelastik ve geometri küçük deformasyon bölgesindedir.
- G', G'', J ve L kullanılan çalışma noktası için geçerli kabul edilir.
- Geometrik nonlinearite, hiperelastisite ve frekanslar arası interpolasyon
  bu modelin kapsamı dışındadır.
