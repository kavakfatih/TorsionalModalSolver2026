# Dinamik Burulma Rijitliği

## Amaç ve kapsam

Bu belge, frekans ve sıcaklığa bağlı kompleks kayma modülünün gelecekte
kompleks burulma rijitliğine nasıl aktarılacağını tanımlar. V0.2.0 yalnızca
hazırlık veri modelini sağlar; aşağıdaki ilişkileri hesaplayan bir solver veya
kompleks rijitlik yordamı henüz yoktur.

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
  bu hazırlık modelinin kapsamı dışındadır.
