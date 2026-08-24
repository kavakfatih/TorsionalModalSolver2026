# Dinamik Burulma Rijitliği

## Amaç ve kapsam

Bu belge, V0.2.0 sırasında frekans ve sıcaklığa bağlı kompleks kayma modülünün
burulma rijitliğine aktarılması için hazırlanan modeli tanımlar. İlişkiler
V0.2.1 ile `tms_dynamic_torsional_stiffness` modülünde uygulanmıştır. Fiziksel
yorum ve doğal frekans bağlantısı
[`complex_torsional_stiffness.md`](complex_torsional_stiffness.md) belgesinde
ayrıntılandırılır.

## Kompleks burulma rijitliği

Tam bağlı eş merkezli silindirik elastomer tabakanın lineer kayma
idealizasyonunda kompleks burulma rijitliği şöyledir:

```text
K* = K' + iK''
Cθ = 4πLri²ro²/(ro²-ri²)
K' = G'Cθ
K'' = G''Cθ
```

`K'`, elastik olarak enerji depolayan burulma rijitliği bileşenidir. `K''`,
çevrimsel enerji kaybıyla ilişkili burulma rijitliği bileşenidir. Her iki
bileşen de kullanılan G' ve G'' verilerinin frekans ve sıcaklığına bağlıdır.

## Fiziksel büyüklükler ve birimler

| Büyüklük | Fiziksel anlam | SI birimi |
| --- | --- | --- |
| `G'` | Depolama kayma modülü | `Pa = N/m²` |
| `G''` | Kayıp kayma modülü | `Pa = N/m²` |
| `ri` | İç silindirik bağ yüzeyi yarıçapı | `m` |
| `ro` | Dış silindirik bağ yüzeyi yarıçapı | `m` |
| `L` | Bağlı elastomerin eksenel genişliği | `m` |
| `Cθ` | Annüler kauçuk burç geometri faktörü | `m³` |
| `K'` | Depolama burulma rijitliği | `N·m/rad` |
| `K''` | Kayıp burulma rijitliği | `N·m/rad` |

`Cθ`, polar alan momenti veya polar kütle ataleti değildir. Pa birimli kayma
modülüyle çarpıldığında `N·m/rad` birimli rijitlik üretir. Radyan SI boyut
analizinde boyutsuz kabul edilse de mühendislik anlamını açık tutmak için
rijitlik biriminde gösterilir.

## Varsayımlar ve sınırlar

- Elastomer lineer viskoelastik, silindirik yüzeylere tam bağlı ve geometri
  küçük deformasyon bölgesindedir; uç etkileri ihmal edilir.
- G', G'', ri, ro ve L kullanılan çalışma noktası için geçerli kabul edilir.
- Hesapta `0 < ri < ro`, `L > 0`, `G' > 0` ve `G'' >= 0` koşulları
  doğrulanır; geçersiz girdiler reddedilir.
- `ri = 0`, bağlı iç silindirik yüzey bulunmadığı için kapsam dışıdır.
- Annüler milin eksen boyunca burulmasına ait `GJp/ℓ` modeli ayrıdır ve
  TVD solver'larında kullanılmaz.
- Geometrik nonlinearite, hiperelastisite ve frekanslar arası interpolasyon
  bu modelin kapsamı dışındadır.
