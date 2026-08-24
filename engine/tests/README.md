# TMS26 Hesap Motoru Testleri

Bu dizindeki Fortran test programları CMake ile derlenir ve CTest üzerinden
çalıştırılır. Yeni veya değiştirilen her davranış uygun test kategorisiyle
kapsanmalıdır.

## Unit test

Tek modül, tür veya yordam izole olarak doğrulanır. `test_kinds`,
`test_constants`, `test_units`, `test_geometry`, `test_material`,
`test_dynamic_modulus`, `test_torsional_node` ve `test_torsional_element` bu
kategoridedir. `test_dynamic_modulus`, G' ile G''
değerlerinin SI biriminde saklanmasını, çoklu çalışma noktası altyapısını ve
`tan(delta) = G'' / G'` hesabını doğrular. `test_geometry`, geometri veri
alanlarının yanında annüler kesit polar alan momentini ve tam bağlı kauçuk
burç `Cθ` faktörünü ayrı bağımsız analitik değerlerle sınar.

## Physics validation

Fiziksel veya matematiksel yordam, bağımsız analitik sonuç ve açık bir hata
sınırıyla doğrulanır. `test_inertia`, `test_hub_inertia`,
`test_torsional_stiffness`, `test_dynamic_torsional_stiffness`,
`test_frequency_solver` ve `test_torsional_system` bu kategoridedir.
Genel fizik doğrulamalarında bağıl hata yüzde `0,1`'den küçük olmalıdır;
dinamik burulma rijitliği testi `1e-10` bağıl hata sınırı kullanır. Bu test
üretim `tms_dynamic_torsional_stiffness` modülünü doğrudan çağırarak Cθ,
K', K'', frekans, sıcaklık ve
`G''/G' = K''/K' = tan(delta)` eşitliğini birlikte doğrular. Ayrıca K'
bileşeninin statik solver sonucuyla aynı kaldığını, `G'' = 0` için kayıp
bileşenlerinin sıfırlandığını ve `ri`, `ro` değerine yaklaştıkça ideal
model rijitliğinin kuvvetli biçimde arttığını sınar. Dinamik test ayrıca
`K' ∝ G'`, `K' ∝ L` ve `K'' ∝ L` ölçekleme regresyonlarını üretim
fonksiyonu sonuçlarını karşılaştırarak doğrular. Statik test, G ve L iki katına
çıkarıldığında K değerinin de iki katına çıktığını doğrular.

Negatif veya sıfır iç yarıçap, sırasız/eşit yarıçaplar, pozitif olmayan
eksenel genişlik, pozitif olmayan depolama modülü ve negatif kayıp modülü
ayrı CTest vakalarıdır.
Bu vakalarda üretim yordamının `error stop` ile sonlanması beklenir ve
`WILL_FAIL` özelliği beklenen reddi test başarısına dönüştürür.

`test_hub_inertia`, homojen annüler göbeğin hacim, kütle ve polar kütle
ataletini bağımsız analitik sabitlerle doğrular. `test_torsional_system`,
builder entegrasyonunu, fixed-hub eşdeğerliğini, serbest-serbest sıfır ve
elastik modları, normalize mod şekillerini, K'/atalet ölçeklemelerini ve büyük
göbek ataleti limitini sınar. K'' ve kayıp faktörü değişiminin sönümsüz modal
sonucu değiştirmediği de ayrı bir regresyonla korunur. Pozitif olmayan atalet,
K' ve yoğunluklar ile geçersiz göbek/halka geometrileri ayrı `WILL_FAIL`
regresyonlarıdır.

`test_torsional_node`, kimlik, polar atalet, başlangıç açısı ve sınır koşulu
alanlarını; `test_torsional_element`, uç düğüm kimlikleri, K rijitliği ve
eşdeğer viskoz c alanını doğrular. Pozitif olmayan veya sonlu olmayan
büyüklükler, self-connection ve negatif sönüm ayrı `WILL_FAIL` vakalarıdır.

`test_local_stiffness_matrix`, `k = 100 N·m/rad` için üretim yordamının
`[[100,-100],[-100,100]]` lokal katkısını verdiğini doğrular. Simetri, sıfır
satır toplamı, bağıl dönmede pozitif enerji ve ortak dönmede sıfır enerji ayrı
assertion'larla korunur. Negatif rijitlik aynı üretim yordamını çağıran ayrı bir
`WILL_FAIL` regresyonudur.

`test_generalized_torsional_system`, private koleksiyonların public yönetim
yordamlarını, aktif DOF sayımını ve Benchmark 004 iki-ataletli sisteminin genel
iki-düğüm/bir-eleman gösterimini doğrular. K'' kayıp rijitliğinin viskoz sönüm
alanına aktarılmadığı, serbest-serbest gösterimin iki ve fixed-hub gösteriminin
bir aktif DOF taşıdığı sınanır. Boş sistem, yinelenen kimlik ve tanımsız eleman
ucu hata regresyonlarıdır.

`test_matrix_assembly`, fiziksel node ID ile equation ID ayrımını, tek eleman
ve üç düğümlü zincir için global torsional K matrisini, K simetrisini, serbest
sistem sıfır satır toplamını ve rijit-cisim null modunu doğrular. Global M
matrisinde düğüm polar ataletlerinin `kg·m²` birimiyle diagonal saklandığını ve
köşegen dışı katsayıların sıfır kaldığını sınar. Kısıtlı uç için indirgenmiş
`K=[k]`, `M=[J_free]` ve tamamen kısıtlı sistem için 0x0 matrisler korunur.
Geçersiz boyut, eksik node lookup, aralık dışı denklem, negatif atalet ve
uyumsuz DOF haritası ayrı `WILL_FAIL` regresyonlarıdır.

`test_torsional_validation`, V0.3.0 foundation katmanlarını değiştirmeden
uçtan uca analitik doğrulama sağlar. Tek eleman işaret konvansiyonu, serbest
rijit-cisim modu, üç düğümlü global katkı toplamı, iki constraint durumunda
DOF mapping sürekliliği, Frobenius simetri normu ve
`U=1/2 theta^T K theta>=0` enerji koşulu birlikte sınanır. İki ataletli TVD
referansı ayrıca assembled M/K modal residual, bilinen analitik modun Rayleigh
quotient değeri ve kütle ortogonalliği üzerinden mevcut analitik solver ile
çapraz doğrulanır. Bu test genel eigen çözümü yapmaz ve deneysel model
validation yerine analytical code verification sunar.

`test_numerical_hardening`, V0.3.2 girdi ve sayı aralığı sağlamlaştırmasını
doğrular. Nominal test akışı; mutlak `1e-12 m` ve bağıl `1e-9` toleranslı TVD
geometri ara yüzlerini, uç atalet/rijitlik ölçeklerinde sonlu eşdeğer atalet ile
doğal frekansı ve fixed-DOF için indirgenmiş `K=[k]`, `M=[J_free]` matrislerini
sınar.

Atalet, yoğunluk, yarıçap, uzunluk, rijitlik, modül, frekans ve sıcaklık
alanlarındaki IEEE `NaN`, pozitif sonsuz ve negatif sonsuz değerleri ayrı
`tms26.numerical_hardening.rejects_*` CTest süreçleridir. Bu süreçler yalnız
üretim doğrulayıcısının kontrollü `error stop` sonucunu başarı kabul eder.
Assertion yardımcıları da sonlu olmayan actual, expected veya tolerance
değerlerinin testi yanlış biçimde geçirmesine izin vermez. Ayrıntılar
[`../../docs/validation/numerical_robustness_validation.md`](../../docs/validation/numerical_robustness_validation.md)
belgesindedir.

## Benchmark regression

Birden fazla fizik adımını temsil eden sabit referans modelin sonuçları zaman
içinde korunur. `benchmarks/001_simple_annular_tvd/` girdileri; kütle, polar
atalet, rijitlik ve doğal frekans testlerinin ortak regresyon temelidir. Model
veya kabul edilen formül değişirse benchmark girdileri, beklenen sonuçlar ve
ilgili testler aynı commit içinde güncellenir.

`benchmarks/002_dynamic_elastomer/` dinamik modül veri noktasını,
`benchmarks/003_dynamic_torsional_stiffness/` ise bu noktanın annüler geometri
üzerindeki kompleks rijitlik sonucunu tanımlar.
`benchmarks/004_two_inertia_tvd/`, fixed-hub ve serbest-serbest iki ataletli
sistemin analitik frekansları ile mod şekillerini ve V0.2.3 genel topoloji
eşlemesini tanımlar.

## Test ekleme

Yeni Fortran testi `test_<konu>.f90` biçiminde adlandırılır ve
`engine/CMakeLists.txt` içindeki `tms26_add_fortran_test` yordamıyla CTest'e
kaydedilir. Tüm testler aşağıdaki komutla çalıştırılır:

```sh
ctest --test-dir build --output-on-failure
```
