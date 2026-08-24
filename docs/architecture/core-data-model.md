# Çekirdek Veri Modeli

## Amaç

V0.2.4 çekirdek veri modeli, torsional vibration damper (TVD) bileşenlerinin
geometrik ve malzeme özelliklerini fizik yordamlarına, genel düğüm-eleman
topolojisini ise gelecekteki sistem analiz katmanlarına taşır. Geometri ve
malzeme türleri yalnız veri taşır; genel sistem yönetim yordamları topolojik ve
temel fiziksel önkoşulları doğrular.

## Modül bağımlılıkları

- `tms_kinds`, proje genelindeki çift hassasiyetli `dp` türünü tanımlar.
- `tms_constants`, `tms_kinds` türüyle matematik ve dönüşüm sabitlerini tanımlar.
- `tms_units`, `tms_constants` çarpanlarını kullanan birim dönüşümlerini sunar.
- `tms_geometry`, TVD geometri türlerini, annüler kesit polar alan momentini
  ve tam bağlı annüler kauçuk burç geometri faktörünü tanımlar.
- `tms_dynamic_modulus`, G', G'', frekans ve sıcaklık alanlarıyla kompleks
  kayma modülü veri türünü ve kayıp faktörü hesabını tanımlar.
- `tms_material_frequency`, dinamik modül türünü bir malzeme çalışma noktası
  olarak genişletir.
- `tms_material`, geriye uyumlu elastomer türü içinde sıfır veya daha fazla
  frekans-sıcaklık çalışma noktası saklar.
- `tms_inertia`, `hub_geometry_t` ve `inertia_ring_geometry_t` ile göbek ve
  halka kütle özelliklerini ortak annüler rijit gövde hesabından üretir.
- `tms_torsional_stiffness`, `rubber_geometry_t` ve
  `dynamic_rubber_material_t` ile lineer burulma rijitliğini hesaplar.
- `tms_dynamic_torsional_stiffness`, aynı girdilerden K', K'', kayıp faktörü
  ve çalışma noktası üstverisini üretir.
- `tms_frequency_solver`, hesaplanan rijitlik ve ataletten doğal frekansı bulur.
- `tms_local_matrix`, iki uçlu bir elemanın 2x2 lokal matris katsayılarını
  sistem katmanından bağımsız ve sabit boyutlu bir veri türünde taşır.
- `tms_torsional_node`, yığılmış polar atalet, başlangıç açısı ve dönel sınır
  koşulu taşıyan genel düğüm türünü tanımlar.
- `tms_torsional_element`, iki düğüm arasındaki lineer rijitlik ve eşdeğer
  viskoz sönüm bağlantısını tanımlar ve 2x2 lokal rijitlik katkısını üretir.
- `tms_generalized_torsional_system`, private düğüm/eleman koleksiyonlarını,
  ekleme-okuma yordamlarını, aktif DOF sayımını ve sistem doğrulamasını sağlar.
- `tms_torsional_system`, bileşen sonuçlarını iki ataletli sistem türünde
  birleştirir; fixed-hub ve serbest-serbest analitik modal sonuçları üretir ve
  bu özel modeli genel topolojiye dönüştürür.

Mevcut torsional fizik akışı aşağıdaki sırayı izler:

1. Mühendislik birimleri `tms_units` ile SI birimlerine dönüştürülür.
2. SI değerleri mevcut geometri ve malzeme türlerine yazılır.
3. Göbek ve halka ataletleri ile kompleks burulma rijitliği mevcut bileşen
   yordamlarında birbirinden bağımsız hesaplanır.
4. Sistem builder'ı bu sonuçları `two_inertia_tvd_system_t` içinde birleştirir.
5. Fixed-hub ve serbest-serbest modal yordamlar yalnız K' depolama rijitliğini
   mevcut doğal frekans yordamına aktarır; K'' sistem üstverisinde korunur.
6. İstenirse iki ataletli özel sistem, iki düğüm ve bir elemandan oluşan genel
   torsional topolojiye dönüştürülür; bu adım frekansı yeniden çözmez.

Kompleks sonuç, reel ve sanal bileşenleri açıkça adlandıran
`complex_torsional_stiffness_t` derived type değeriyle taşınır.

## Geometri türleri

İlk geometri modeli eksenel simetrik, eş merkezli bileşenleri temsil eder:

- `rubber_geometry_t`: iç yarıçap, dış yarıçap ve eksenel uzunluk
- `inertia_ring_geometry_t`: iç yarıçap, dış yarıçap ve eksenel uzunluk
- `hub_geometry_t`: delik yarıçapı, dış yarıçap ve eksenel uzunluk
- `tvd_geometry_t`: yukarıdaki üç bileşeni bir arada tutan bileşik tür

Tüm geometrik alanlar metre cinsindedir. Sıfır başlangıç değerleri yalnızca
deterministik ilk durumu sağlar; fiziksel olarak geçerli bir geometri anlamına
gelmez. Geometri modülü iki ayrı fizik modelini karıştırmayan yardımcı
yordamlar sunar:

- `calculate_rubber_polar_area_moment`, eksen boyunca Saint-Venant burulması
  uygulanacak annüler kesit için `Jp = π/2 (ro⁴-ri⁴)` değerini (`m⁴`)
  hesaplar. TVD rijitlik solver'ları bu yordamı kullanmaz.
- `calculate_annular_bush_torsion_geometry_factor`, rijit iç göbek ve dış
  halkaya tam bağlı elastomer için
  `Cθ = 4πLri²ro²/(ro²-ri²)` değerini (`m³`) hesaplar. Statik ve
  dinamik TVD rijitlik solver'ları bu faktörü kullanır.

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
uyumlu tek çalışma noktasındaki G' değerini kullanmayı sürdürür. V0.2.1 dinamik
solver'ı aynı noktadaki G' ve G'' değerlerini kullanır; yeni dizi üzerinde seçim
veya interpolasyon henüz yapılmaz.

## Kompleks rijitlik sonuç türü

`complex_torsional_stiffness_t` aşağıdaki alanları taşır:

- `storage_stiffness`: K' (`N·m/rad`)
- `loss_stiffness`: K'' (`N·m/rad`)
- `loss_factor`: K''/K' (`boyutsuz`)
- `frequency`: kullanılan çalışma noktası (`Hz`)
- `temperature`: kullanılan çalışma noktası (`K`)

Bu sonuç türü solver katmanında tanımlıdır. Geometri ve malzeme veri türleri
solver sonuçlarına bağımlı değildir; bağımlılık yönü çekirdek veriden solver'a
doğrudur.

## Torsional sistem ve modal sonuç türleri

`two_inertia_tvd_system_t`, aşağıdaki SI alanlarıyla bileşen hesaplarını tek bir
modal sistem durumunda birleştirir:

- göbek ve atalet halkası polar kütle ataletleri (`kg·m²`)
- K' depolama ve K'' kayıp burulma rijitlikleri (`N·m/rad`)
- boyutsuz kayıp faktörü
- dinamik malzemenin referans frekansı (`Hz`) ve sıcaklığı (`K`)

`two_inertia_modal_result_t`, DOF sırası `[göbek, atalet halkası]` olacak
biçimde rijit-cisim ve elastik mod frekanslarını (`Hz`) ve normalize mod
şekillerini taşır. Rijit-cisim modu `[1,1]`, elastik mod ise
`[1,-J_h/J_r]` biçimindedir.

Sistem builder'ı atalet ve dinamik rijitlik denklemlerini yinelemez; mevcut
üretim yordamlarını çağırıp sonuçlarını aktarır. Modal katman K'' değerinden
sönüm veya kompleks özdeğer türetmez. G' ve K', malzeme referans çalışma
noktasında sabitlenir; hesaplanan doğal frekansa göre interpolasyon veya
iterasyon yapılmaz.

## Genel torsional topoloji türleri

`torsional_node_t`, pozitif benzersiz kimlik, polar kütle ataleti (`kg·m²`),
başlangıç açısı (`rad`) ve sabitlenmişlik bilgisi taşır. Her sabitlenmemiş düğüm
bir aktif torsional DOF oluşturur. Düğüm kimliği doğrudan DOF indeksi değildir.

`torsional_element_t`, pozitif benzersiz kimlik, iki farklı uç düğüm kimliği,
lineer rijitlik (`N·m/rad`) ve eşdeğer viskoz sönüm (`N·m·s/rad`) taşır.
Paralel elemanlar farklı kimliklerle temsil edilebilir. Elemanın saf
`calculate_local_stiffness` yordamı, `[theta_i, theta_j]` yerel sırası için
`Ke = k[[1,-1],[-1,1]]` matrisini `local_matrix_2x2` değeri olarak döndürür.

`torsional_system_t` içindeki koleksiyonlar private'dır. Public saf yordamlar:

- düğüm ve eleman ekler,
- düğüm/eleman sayısını ve ekleme sırasındaki eleman kopyasını döndürür,
- sabitlenmemiş düğümlerden aktif DOF sayısını belirler,
- kimlik benzersizliğini ve bağlantı uçlarının varlığını doğrular.

Eleman katmanı yalnız lokal K katkısını üretir. Genel sistem katmanı bu katkıyı
global K matrisine birleştirmez; M veya C matrisi oluşturmaz ve modal çözüm
yapmaz. Bağımlılık yönü `tms_kinds -> tms_local_matrix ->
tms_torsional_element -> tms_generalized_torsional_system` biçimindedir.
Mevcut iki ataletli dönüşüm `J_h`, `J_r` ve K' değerlerini taşır. K''
`N·m/rad`, viskoz `c` ise `N·m·s/rad` olduğundan damping alanına doğrudan
aktarılmaz.
