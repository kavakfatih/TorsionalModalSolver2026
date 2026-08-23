# Çekirdek Veri Modeli

## Amaç

V0.2.0 çekirdek veri modeli, torsional vibration damper (TVD) bileşenlerinin
geometrik ve malzeme özelliklerini fizik yordamlarına taşır. Veri türleri
doğrulama veya hesap yapmaz; fizik davranışı `engine/src/solver/` altındaki
modüllerde tutulur.

## Modül bağımlılıkları

- `tms_kinds`, proje genelindeki çift hassasiyetli `dp` türünü tanımlar.
- `tms_constants`, `tms_kinds` türüyle matematik ve dönüşüm sabitlerini tanımlar.
- `tms_units`, `tms_constants` çarpanlarını kullanan birim dönüşümlerini sunar.
- `tms_geometry`, `tms_kinds` tabanlı TVD geometri türlerini tanımlar.
- `tms_dynamic_modulus`, G', G'', frekans ve sıcaklık alanlarıyla kompleks
  kayma modülü veri türünü ve kayıp faktörü hesabını tanımlar.
- `tms_material_frequency`, dinamik modül türünü bir malzeme çalışma noktası
  olarak genişletir.
- `tms_material`, geriye uyumlu elastomer türü içinde sıfır veya daha fazla
  frekans-sıcaklık çalışma noktası saklar.
- `tms_inertia`, `inertia_ring_geometry_t` ile halka kütlesi ve ataletini
  hesaplar.
- `tms_torsional_stiffness`, `rubber_geometry_t` ve
  `dynamic_rubber_material_t` ile lineer burulma rijitliğini hesaplar.
- `tms_frequency_solver`, hesaplanan rijitlik ve ataletten doğal frekansı bulur.

Mevcut torsional fizik akışı aşağıdaki sırayı izler:

1. Mühendislik birimleri `tms_units` ile SI birimlerine dönüştürülür.
2. SI değerleri mevcut geometri ve malzeme türlerine yazılır.
3. Atalet ve burulma rijitliği birbirinden bağımsız hesaplanır.
4. Bu iki skaler sonuç doğal frekans yordamına verilir.

Fizik yordamları için yeni bir derived type oluşturulmamıştır.

## Geometri türleri

İlk geometri modeli eksenel simetrik, eş merkezli bileşenleri temsil eder:

- `rubber_geometry_t`: iç yarıçap, dış yarıçap ve eksenel uzunluk
- `inertia_ring_geometry_t`: iç yarıçap, dış yarıçap ve eksenel uzunluk
- `hub_geometry_t`: delik yarıçapı, dış yarıçap ve eksenel uzunluk
- `tvd_geometry_t`: yukarıdaki üç bileşeni bir arada tutan bileşik tür

Tüm geometrik alanlar metre cinsindedir. Sıfır başlangıç değerleri yalnızca
deterministik ilk durumu sağlar; fiziksel olarak geçerli bir geometri anlamına
gelmez.

## Dinamik elastomer türleri

`dynamic_rubber_material_t`, tek bir sıcaklık ve frekans çalışma noktasındaki
aşağıdaki verileri taşır:

- malzeme adı
- yoğunluk (`kg/m³`)
- storage shear modulus G' (`Pa`)
- loss shear modulus G'' (`Pa`)
- sıcaklık (`K`)
- frekans (`Hz`)

Bu alanlar V0.1.x istemcileriyle kaynak uyumluluğunu korumak için değişmeden
bırakılmıştır. V0.2.0 ile eklenen `frequency_points` dizisi,
`material_frequency_point` elemanları üzerinden birden fazla çalışma noktası
saklar. Eski tek noktalı alanlarla yeni dizi arasında otomatik eşzamanlama
yoktur; çağıran kod hangi temsili kullandığını açıkça seçer.

`material_frequency_point`, `dynamic_shear_modulus` türünü genişletir ve şu SI
alanlarını devralır:

- `storage_modulus`: G' (`Pa`)
- `loss_modulus`: G'' (`Pa`)
- `frequency`: (`Hz`)
- `temperature`: (`K`)

Bu ayrım DMA ile Dewesoft veya Brüel & Kjær modal test verilerinin ortak veri
alanlarına alınmasına hazırlık sağlar. Mevcut burulma rijitliği hesabı, geriye
uyumlu tek çalışma noktasındaki G' değerini kullanmayı sürdürür; yeni dizi ve
G'' henüz solver akışına bağlanmamıştır.
