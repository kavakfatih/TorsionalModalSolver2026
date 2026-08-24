# Karar 0005: Tam Bağlı Annüler Kauçuk Burç Modeli

- Durum: Kabul edildi
- Tarih: 2026-08-24
- Geçersiz kıldığı karar: Karar 0004 içindeki ortak `Jp` kullanımı

## Bağlam

V0.2.1 statik ve dinamik TVD rijitlik solver'ları, annüler kesitin eksen
boyunca Saint-Venant burulmasına ait `K = GJp/ℓ` denklemini kullanıyordu. Bu
denklem momentin kesitin uç yüzeylerinden aktarıldığı prizmatik bir mil içindir.
TMS26 TVD geometrisinde ise elastomer, rijit iç göbek ve rijit dış halkanın
silindirik yüzeylerine tam bağlıdır ve moment bu yüzeylerden aktarılır.

## Karar

TMS26 V0.2.x ana TVD modeli, tam bağlı eş merkezli silindirik kauçuk burç
denklemini kullanacaktır:

```text
Cθ = 4πLri²ro²/(ro²-ri²)
Kθ = G Cθ
K' = G' Cθ
K'' = G'' Cθ
```

`calculate_annular_bush_torsion_geometry_factor`, `m³` birimli `Cθ` değerini
hesaplayan ortak saf geometri yordamı olacaktır. Statik ve dinamik TVD
solver'ları bu yordamı kullanacaktır.

`calculate_rubber_polar_area_moment` silinmeyecektir. Yordam, annüler mil veya
uç yüzey burulması gibi farklı sınır koşulları ve ilerideki doğrulama
senaryoları için `Jp` geometri değerini hesaplamayı sürdürecektir.

## Fiziksel varsayımlar

- İç göbek ve dış halka rijit, eş merkezli silindirlerdir.
- Elastomer homojen, izotrop ve lineerdir; deformasyonlar küçüktür.
- Elastomer iki silindirik yüzeye kaymasız ve ayrılmasız tam bağlıdır.
- Eksenel uç etkileri ihmal edilir.
- Geçerli geometri `0 < ri < ro` ve `L > 0` koşullarını sağlar.
- `ri = 0`, moment aktaracak bir iç silindirik bağ yüzeyi bulunmadığından
  model kapsamı dışındadır.

## Sonuçlar

- Eksenel bağlı genişlik `L` arttıkça rijitlik doğrusal artar.
- Radyal boşluk azaldıkça ideal model rijitliği kuvvetli biçimde artar.
- `K''/K' = G''/G'` eşitliği ortak `Cθ` faktörü nedeniyle korunur.
- Eski TVD benchmark ve doğal frekans referansları yeni denklemle
  güncellenmelidir.
- Sonuç veri türü ve statik/dinamik solver API imzaları değişmez.
