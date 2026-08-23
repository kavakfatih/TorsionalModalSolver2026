# Çekirdek Veri Modeli

## Amaç

V0.1.2 çekirdek veri modeli, torsional vibration damper (TVD) bileşenlerinin
geometrik ve malzeme özelliklerini fizik yordamlarına taşır. Veri türleri
doğrulama veya hesap yapmaz; fizik davranışı `engine/src/solver/` altındaki
modüllerde tutulur.

## Modül bağımlılıkları

- `tms_kinds`, proje genelindeki çift hassasiyetli `dp` türünü tanımlar.
- `tms_constants`, `tms_kinds` türüyle matematik ve dönüşüm sabitlerini tanımlar.
- `tms_units`, `tms_constants` çarpanlarını kullanan birim dönüşümlerini sunar.
- `tms_geometry`, `tms_kinds` tabanlı TVD geometri türlerini tanımlar.
- `tms_material`, `tms_kinds` tabanlı dinamik elastomer türünü tanımlar.
- `tms_inertia`, `inertia_ring_geometry_t` ile halka kütlesi ve ataletini
  hesaplar.
- `tms_torsional_stiffness`, `rubber_geometry_t` ve
  `dynamic_rubber_material_t` ile lineer burulma rijitliğini hesaplar.
- `tms_frequency_solver`, hesaplanan rijitlik ve ataletten doğal frekansı bulur.

V0.1.2 fizik akışı aşağıdaki sırayı izler:

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

## Dinamik elastomer türü

`dynamic_rubber_material_t`, tek bir sıcaklık ve frekans çalışma noktasındaki
aşağıdaki verileri taşır:

- malzeme adı
- yoğunluk (`kg/m³`)
- storage shear modulus G' (`Pa`)
- loss shear modulus G'' (`Pa`)
- sıcaklık (`K`)
- frekans (`Hz`)

Gelecek sürümlerde doğrulama ve frekans-sıcaklık bağımlı malzeme tabloları ayrı
API katmanlarında ele alınacaktır. V0.1.2 burulma rijitliği hesabı, seçili
çalışma noktasındaki G' değerini lineer kayma modülü olarak kullanır; G'' henüz
çözümde yer almaz.
