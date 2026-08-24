# Çekirdek Veri Modeli

## Amaç

V0.4.0 çekirdek veri modeli, torsional vibration damper (TVD) bileşenlerinin
geometrik ve malzeme özelliklerini fizik yordamlarına, genel düğüm-eleman
topolojisini ise tam M/K assembly, constraint yönetimi ve indirgenmiş Kr/Mr
katmanına taşır. Geometri ve malzeme türleri yalnız veri taşır; sistem, matris
ve constraint yordamları topolojik, boyutsal ve temel fiziksel önkoşulları
doğrular.

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
- `tms_matrix_types`, private allocatable depolamalı genel dense matris türünü,
  boyut ve katsayı erişim yordamlarını sağlar.
- `tms_dof_types`, anlamlı DOF türlerini ve `(node_id,dof_type)` Physical DOF
  kimliğini tanımlar. İlk desteklenen tür `TORSIONAL_ROTATION` değeridir.
- `tms_torsional_node`, yığılmış polar atalet, başlangıç açısı ve dönel sınır
  koşulu taşıyan genel düğüm türünü tanımlar.
- `tms_torsional_element`, iki düğüm arasındaki lineer rijitlik ve eşdeğer
  viskoz sönüm bağlantısını tanımlar ve 2x2 lokal rijitlik katkısını üretir.
- `tms_generalized_torsional_system`, private düğüm/eleman koleksiyonlarını,
  ekleme-okuma yordamlarını, aktif DOF sayımını ve sistem doğrulamasını sağlar.
- `tms_dof_map`, Physical DOF'ları constraint durumundan bağımsız tam Equation
  ID değerlerine eşler ve sistem-harita uyumunu doğrular.
- `tms_constraint_types`, fixed/prescribed constraint kayıtlarını, hedef
  Physical DOF'u ve radyan cinsindeki prescribed değeri taşır.
- `tms_constraint_manager`, constraint kayıtlarını doğrular ve tam Equation ID
  değerlerinden ayrı Active Equation ID haritasını üretir.
- `tms_stiffness_matrix`, lokal 2x2 katkıları global K matrisine toplar.
- `tms_mass_matrix`, düğüm polar ataletlerini global diagonal M matrisine ekler.
- `tms_matrix_assembly`, sistem topolojisini ve DOF haritasını tüketerek global
  tam torsional M/K matrislerini üretir.
- `tms_matrix_reduction`, tam matrisi aktif denklem haritasıyla storage-bağımsız
  direct elimination üzerinden indirger.
- `tms_reduced_system`, Kr/Mr çiftini ve tam fiziksel sonuç recovery bilgisini
  tek bir doğrulanabilir solver girdisinde birleştirir.
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
7. `tms_dof_map`, her Physical DOF'a constraint'ten bağımsız tam Equation ID
   verir.
8. Elemanların lokal K katkıları tam denklem kimlikleri üzerinden global K
   matrisine; tüm düğümlerin polar ataletleri tam diagonal M matrisine eklenir.
9. `tms_constraint_manager`, fixed veya prescribed constraint kayıtlarını
   doğrular ve Active Equation ID haritasını üretir.
10. `tms_matrix_reduction`, tam K/M matrislerinden Kr/Mr aktif alt sistemini
    çıkarır.
11. `tms_reduced_system`, Kr/Mr ile `q=Pq_r+q_p` ve gelecekteki
    `phi=Pphi_r` recovery bilgisini birlikte taşır.

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
başlangıç açısı (`rad`) ve geriye uyumlu sabitlenmişlik bilgisi taşır. Her düğüm
bir `TORSIONAL_ROTATION` Physical DOF oluşturur; aktif olup olmadığı constraint
katmanında belirlenir. Düğüm kimliği doğrudan DOF veya matris indeksi değildir.

`torsional_element_t`, pozitif benzersiz kimlik, iki farklı uç düğüm kimliği,
lineer rijitlik (`N·m/rad`) ve eşdeğer viskoz sönüm (`N·m·s/rad`) taşır.
Paralel elemanlar farklı kimliklerle temsil edilebilir. Elemanın saf
`calculate_local_stiffness` ve geriye uyumlu standart
`get_local_stiffness` arayüzü, `[theta_i, theta_j]` yerel sırası için
`Ke = k[[1,-1],[-1,1]]` matrisini `local_matrix_2x2` değeri olarak döndürür.

`torsional_system_t` içindeki koleksiyonlar private'dır. Public saf yordamlar:

- düğüm ve eleman ekler,
- düğüm/eleman sayısını ve ekleme sırasındaki eleman kopyasını döndürür,
- düğüm ve geriye uyumlu sabitlenmişlik bilgisini sorgular,
- kimlik benzersizliğini ve bağlantı uçlarının varlığını doğrular.

Eleman katmanı yalnız lokal K katkısını üretir. `tms_dof_map`, ekleme sırasındaki
tüm Physical DOF'lara kesintisiz tam Equation ID verir. `tms_matrix_assembly`,
lokal K katkılarını tam global K matrisine ve tüm düğüm J değerlerini diagonal
tam M matrisine toplar. Constraint manager bundan sonra ayrı Active Equation ID
haritasını kurar; reduction katmanı Kr/Mr matrislerini üretir. Genel C assembly
veya modal çözüm yapılmaz.

Temel bağımlılık katmanında `tms_dof_types` bağımsız fiziksel kimlikleri,
`tms_kinds` ise `tms_local_matrix` ve `tms_matrix_types` sayısal türlerini
besler. Lokal matris dalı
`tms_torsional_element -> tms_generalized_torsional_system -> tms_dof_map`,
dense dalı ise `tms_matrix_types -> {tms_stiffness_matrix,tms_mass_matrix}`
biçimindedir. Sistem, eleman, tam DOF haritası ve iki global matris türü
`tms_matrix_assembly` içinde birleşir. Constraint türleri ve manager aktif
haritayı oluşturur; `tms_matrix_reduction` ile `tms_reduced_system` bu haritayı
tam matrislere bağlar. Private depolama sayesinde sparse gerçekleştirim ileride
eklenebilir.

V0.3.0'da `equation_id=0` ile doğrudan aktif matris assembly yaklaşımı
kullanılıyordu. V0.4.0 tam Equation ID bilgisini korur ve constraint
eliminasyonunu assembly sonrasına taşır. Tarihsel karar
[`../decisions/0008-global-matrix-assembly-design.md`](../decisions/0008-global-matrix-assembly-design.md),
yeni sorumluluk ayrımı ise
[`V0.4_constraint_foundation.md`](V0.4_constraint_foundation.md) ile
[`../decisions/0009-constraint-reduction-architecture.md`](../decisions/0009-constraint-reduction-architecture.md)
belgelerinde açıklanır.
Mevcut iki ataletli dönüşüm `J_h`, `J_r` ve K' değerlerini taşır. K''
`N·m/rad`, viskoz `c` ise `N·m·s/rad` olduğundan damping alanına doğrudan
aktarılmaz.
