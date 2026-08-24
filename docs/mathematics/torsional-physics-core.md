# Torsional Physics Core

## Amaç ve kapsam

Bu belge, eksenel simetrik annüler bir TVD için analitik hesap zincirini
tanımlar. Model; homojen atalet halkası, lineer elastik elastomer ve tek
serbestlik dereceli sönümsüz doğal frekans kabulüne dayanır. V0.2.1.2 fizik
düzeltmesi, TVD elastomerini tam bağlı silindirik kauçuk burç olarak ele alır.
Tüm girdiler hesap yordamlarına SI birimleriyle verilmelidir.

## Annüler halkanın kütlesi ve polar ataleti

İç yarıçapı `ri`, dış yarıçapı `ro` ve eksenel genişliği `b` olan homojen
halkanın hacmi ve kütlesi:

```text
V = π (ro² - ri²) b
m = ρ V
```

Halkanın dönme eksenine göre polar kütle atalet momenti:

```text
J = 1/2 m (ro² + ri²)
```

`ri`, `ro` ve `b` metre; yoğunluk `ρ` kg/m³; kütle `m` kg; polar kütle atalet
momenti `J` kg·m² cinsindendir. Model homojen yoğunluk ve eksenel simetri kabul
eder. Geçerli girdiler `0 <= ri < ro`, `b > 0` ve `ρ > 0` koşullarını sağlar.

## Annüler burulma modellerinin ayrımı

Annüler bir geometri için kullanılacak rijitlik denklemi, momentin hangi
yüzeylerden aktarıldığına bağlıdır. TMS26 iki farklı modeli açıkça ayırır.

### A) Annüler mil veya uç yüzey burulması

Prizmatik bir annüler milin ekseni boyunca Saint-Venant burulmasında moment
uç yüzeylerden aktarılır. Kesitin polar alan momenti ve rijitliği:

```text
Jp = π/2 (ro⁴ - ri⁴)
Kshaft = G Jp / ℓ
```

Burada `Jp` birimi `m⁴`, eksenel burulma boyu `ℓ` birimi `m` ve `Kshaft`
birimi `N·m/rad` değeridir. `ri = 0` katı dairesel mil sınırı olarak bu
modelde geçerlidir. `calculate_rubber_polar_area_moment` yalnızca bu geometrik
`Jp` değerini hesaplar; TMS26 V0.2.x TVD rijitlik solver'ları bu denklemi
kullanmaz.

### B) Tam bağlı annüler TVD kauçuk burcu

TMS26 ana TVD modelinde elastomer, rijit iç göbeğin dış silindirik yüzeyi
ile rijit atalet halkasının iç silindirik yüzeyine tam bağlıdır. Moment,
silindirik ara yüzlerden aktarılır ve göbek ile halka birbirine göre döner.

Eksenel simetrik lineer elastisite çözümünde çevresel yer değiştirme alanı
`uθ(r) = Ar + B/r`, kayma gerilmesi ise
`τrθ = G(duθ/dr - uθ/r)` biçimindedir. İç ve dış yüzeylerdeki tam
bağlılık koşulları uygulandığında bağıl dönme başına moment:

```text
Cθ = 4π L ri² ro² / (ro² - ri²)
Kθ = G Cθ
   = 4π G L ri² ro² / (ro² - ri²)
```

`Cθ` geometri faktörünün birimi `m³`, `G` birimi `Pa` ve `Kθ` birimi
`N·m/rad` değeridir. `L`, eksenel bir burulma yolu değil, bağlı silindirik
yüzeyin eksenel genişliğidir; bu nedenle rijitlik `L` ile doğru orantılıdır.

Model; eş merkezli rijit silindirler, homojen ve izotrop lineer elastomer,
küçük bağıl dönme, ara yüzlerde kaymasız tam bağ ve ihmal edilen uç
etkileri varsayımlarına dayanır. Geçerlilik koşulları `0 < ri < ro`,
`L > 0` ve `G > 0` değerleridir.

`ri = 0`, moment aktaracak iç silindirik bağ yüzeyi bulunmadığı için bu
burç modelinin fiziksel kapsamı dışındadır ve solver tarafından reddedilir.
`ri`, `ro` değerine yaklaştığında payda sıfıra yaklaşır ve ideal lineer
model rijitliği kuvvetli biçimde artar. Tam çakışma `ro <= ri` geçersizdir.
Dinamik bileşenler
[`dynamic_elastomer_model.md`](dynamic_elastomer_model.md) belgesinde verilir.

## Doğal frekans

Atalet halkası ile elastomer rijitliğinden oluşan tek serbestlik dereceli,
sönümsüz burulma sistemi için açısal ve çevrimsel doğal frekanslar:

```text
ωn = sqrt(kθ / J)
fn = ωn / (2π)
```

`kθ` N·m/rad, `J` kg·m², `ωn` rad/s ve `fn` Hz cinsindendir. Model `kθ > 0`
ve `J > 0` koşullarını gerektirir. Sönüm, nonlinear elastomer davranışı, bağlı
çoklu ataletler ve eigen çözümü kapsam dışındadır.

## Sayısal doğrulama

Her fizik yordamı bağımsız analitik değerle sınanır. Kabul ölçütü, hesaplanan
değer ile referans değer arasındaki bağıl hatanın `0,001` değerinden, yani yüzde
`0,1`'den küçük olması temel mühendislik kabulüdür. Dinamik burulma
rijitliği analitik doğrulamasında daha sıkı `1e-10` bağıl hata sınırı
uygulanır. Ortak referans model `benchmarks/001_simple_annular_tvd/` altında
tanımlanmıştır.
