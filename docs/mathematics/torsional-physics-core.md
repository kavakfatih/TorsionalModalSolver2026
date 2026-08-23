# Torsional Physics Core

## Amaç ve kapsam

V0.1.2 fizik çekirdeği, eksenel simetrik annüler bir TVD için ilk analitik
hesap zincirini tanımlar. Model; homojen atalet halkası, lineer elastik
elastomer ve tek serbestlik dereceli sönümsüz doğal frekans kabulüne dayanır.
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

## Elastomer polar alan momenti ve burulma rijitliği

Annüler elastomer kesitin polar alan atalet momenti:

```text
Jp = π/2 (ro⁴ - ri⁴)
```

Lineer elastik kayma modeli için burulma rijitliği:

```text
kθ = G' Jp / L
```

`Jp` m⁴, storage shear modulus `G'` Pa, efektif uzunluk `L` m ve burulma
rijitliği `kθ` N·m/rad cinsindendir. Elastomer homojen ve izotrop kabul edilir;
şekil değiştirmeler küçük, davranış lineer elastiktir. Geçerli girdiler
`0 <= ri < ro`, `L > 0` ve `G' > 0` koşullarını sağlar. Loss shear modulus `G''`
statik sonuçta kullanılmaz. V0.2.1 dinamik solver'ı aynı ortak `Jp` hesabını
kullanarak K' ve K'' bileşenlerini ayrı üretir; ayrıntılı bağıntılar
[`dynamic_elastomer_model.md`](dynamic_elastomer_model.md) altında verilir.
Ortak geometri yordamı negatif yarıçapı ve `ro <= ri` durumunu reddeder.
Dinamik solver ayrıca `L > 0`, `G' > 0` ve `G'' >= 0` koşullarını çalışma
zamanında zorunlu kılar.

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
