# V0.3.2 Numerik Güvenilirlik Doğrulaması

## Amaç ve doğrulama sınıfı

V0.3.2, mevcut torsional fizik denklemlerini ve node-element assembly
mimarisini değiştirmeden girdi doğrulamasını ve sayısal değerlendirmeyi
sağlamlaştırır. Bu belgenin kapsamı:

- IEEE sonlu olmayan girdilerin deterministik olarak reddedilmesi,
- TVD bileşen ara yüzlerinin açık bir toleransla eşleşmesi,
- büyük atalet veya rijitlik oranlarında matematiksel olarak eşdeğer, daha
  kararlı frekans değerlendirmesi,
- homojen sabitlenmiş DOF için indirgenmiş global matrislerin korunmasıdır.

Bu çalışma bir **analytical and numerical code verification** çalışmasıdır.
Deneysel DMA, modal test tezgâhı veya saha verisiyle fiziksel model validation
çalışmasının yerini almaz.

## Sonlu girdi sözleşmesi

Fiziksel hesaplara giren gerçek sayılar, işaret veya sıralama önkoşulları
uygulanmadan önce sonlu olmalıdır:

```text
is_finite(x) = true
```

IEEE `NaN`, pozitif sonsuz ve negatif sonsuz değerleri geçerli mühendislik
girdisi değildir. Karşılaştırmalarda `NaN` sonucu her iki yönde de false
olabildiği için yalnız `x > 0` veya `x >= 0` kontrolü yeterli değildir.

| Girdi sınıfı | SI birimi | Sonluluk sonrası korunan fiziksel koşul |
| --- | --- | --- |
| Polar kütle ataleti | `kg·m²` | `J > 0` |
| Yoğunluk | `kg/m³` | `rho > 0` |
| Yarıçap ve uzunluk | `m` | İlgili geometri modelinin işaret ve sıralama koşulu |
| Torsional rijitlik | `N·m/rad` | `K > 0` |
| Depolama modülü | `Pa` | `G` veya `G' > 0` |
| Kayıp modülü | `Pa` | `G'' >= 0` |
| Frekans | `Hz` | İlgili çalışma noktası sözleşmesi |
| Mutlak sıcaklık | `K` | İlgili çalışma noktası sözleşmesi |

V0.3.2 sonluluk katmanı, mevcut sonlu değer alanlarını veya fiziksel denklemleri
genişletmez. Her alan için `NaN`, pozitif sonsuz ve negatif sonsuz ayrı süreçte
üretim doğrulayıcısına gönderilir. Beklenen sonuç kontrollü `error stop` ile
reddedilmedir.

## Toleranslı TVD geometri ara yüzleri

Bileşik TVD geometrisinde bağlı yüzeylerin nominal olarak çakışması gerekir:

```text
hub.outer_radius       ~= rubber.inner_radius
rubber.outer_radius    ~= inertia_ring.inner_radius
```

Ondalık girdi dönüşümü ve kayan nokta yuvarlaması nedeniyle doğrudan bit eşitliği
kullanılmaz. İki uzunluk `a` ve `b` aşağıdaki koşulla eş kabul edilir:

```text
abs(a-b) <= absolute_tolerance_m
            + relative_tolerance * max(abs(a),abs(b))

absolute_tolerance_m = 1e-12 m
relative_tolerance   = 1e-9
```

Mutlak terim sıfıra yakın değerlerde metre cinsinden bir alt sınır, bağıl terim
ise geometri ölçeği büyüdükçe orantılı pay sağlar. Her iki yüzey çifti için:

- tam eşleşme kabul edilir,
- birleşik toleransın içindeki küçük fark kabul edilir,
- toleransı aşan fark reddedilir,
- sonlu olmayan yarıçap tolerans karşılaştırmasına girmeden reddedilir.

Bu kontrol temas, mesh veya deformasyon çözümü değildir; yalnız tanımlanan
eş merkezli bileşenlerin veri sürekliliği önkoşuludur.

## Kararlı doğal frekans değerlendirmesi

Tek serbestlik dereceli sönümsüz doğal frekansın fizik denklemi değişmez:

```text
f = 1/(2*pi) * sqrt(K/J)
```

Pozitif ve sonlu `K` ile `J` için aynı sonuç, ara `K/J` oranını oluşturmadan
şöyle değerlendirilir:

```text
f = sqrt(K) / (2*pi*sqrt(J))
```

Bu biçim, son frekans temsil edilebilir olduğu halde `K/J` ara değerinin taşma
veya alt-taşma riskini azaltır. Matematiksel model ve SI birimleri değişmez.

Serbest-serbest iki ataletli sistemde eşdeğer atalet:

```text
J_eq = J_h*J_r/(J_h+J_r)
```

doğrudan çarpım yerine pozitif ataletlerin büyüklük sırasından yararlanılarak
değerlendirilir:

```text
J_min = min(J_h,J_r)
J_max = max(J_h,J_r)
J_eq  = J_min/(1 + J_min/J_max)

f_e = sqrt(K')/(2*pi*sqrt(J_eq))
```

İki ifade cebirsel olarak eşdeğerdir. İkinci biçim `J_h*J_r` ara çarpımını ve
büyük/küçük atalet oranlarında gereksiz sayı aralığı kaybını önler.

## Uç ölçek regresyonları

`test_numerical_hardening`, seçilen çok farklı fakat `real(dp)` aralığında
temsil edilebilir atalet ve rijitlik değerleri için aşağıdaki koşulları korur:

- eşdeğer atalet sonlu ve pozitiftir,
- tek-DOF ve iki-ataletli doğal frekanslar sonlu ve pozitiftir,
- kararlı değerlendirme bağımsız cebirsel referansla belirtilen bağıl tolerans
  içinde uyuşur,
- fiziksel ölçekleme ve limit davranışı değişmez.

Test değerleri sonucun kendisinin temsil edilemediği taşma bölgesini geçerli
sonuç olarak kabul etmez. Amaç, temsil edilebilir sonuçtan önce oluşabilecek
gereksiz ara işlem taşmasını önlemektir.

## Fixed-DOF indirgenmiş matris regresyonu

İki düğümlü tek elemanda ilk düğüm homojen olarak sabit, ikinci düğüm serbestse:

```text
equation_id = [0,1]

K_reduced = [K]
M_reduced = [J_free]
```

Bu sonuç V0.3.0 ile eklenen homojen kısıt eliminasyonunun regresyonudur; yeni bir
boundary-condition solver değildir. Sıfırdan farklı prescribed dönme ve sağ
taraf düzeltmesi kapsam dışı kalır.

## Assertion ve `WILL_FAIL` güvenliği

Test assertion yordamları yalnız beklenen ve hesaplanan değerleri değil,
toleransın kendisini de sonlu ve negatif olmayan bir sayı olarak doğrular.
`NaN` veya sonsuz tolerans başarı koşuluna dönüşemez.

Üretim yordamının geçersiz girdiyi `error stop` ile reddettiği her senaryo ayrı
CTest sürecidir ve `WILL_FAIL` özelliği taşır. Her test adı doğrulanan alanı ve
IEEE sınıfını açıkça belirtir. Bilinmeyen komut seçicisi beklenen fiziksel red
olarak yorumlanmaz; böylece CMake ile test programı arasındaki adlandırma hatası
yanlış başarı üretmez.

## Başarı ölçütleri ve izlenebilirlik

V0.3.2 kabulü için:

- tüm nominal ve uç ölçek regresyonları geçmeli,
- kayıtlı her `NaN` ve sonsuz girdi vakası üretim katmanında reddedilmeli,
- temiz Debug build sonrasında tüm CTest takımı başarılı olmalıdır.

Ana doğrulama hedefi `engine/tests/test_numerical_hardening.f90` dosyasından
üretilen `tms26.numerical_hardening` CTest testidir. Ayrı geçersiz girdi
çağrıları aynı yürütülebilir dosyanın kayıtlı `WILL_FAIL` vakalarıdır. Kesin test
sayısı, CMake kaydı tamamlandıktan sonra `ctest -N` ile belirlenir ve görev
raporunda verilir.

İlgili teknik kaynaklar:

- [`torsional_validation.md`](torsional_validation.md)
- [`torsional_frequency_model.md`](../mathematics/torsional_frequency_model.md)
- [`two_inertia_modal_model.md`](../mathematics/two_inertia_modal_model.md)
- [`global_matrix_assembly.md`](../mathematics/global_matrix_assembly.md)
- [`two_inertia_torsional_system.md`](../physics/two_inertia_torsional_system.md)

Bu doğrulama; genel eigen çözümü, deneysel elastomer korelasyonu, nonlinear
malzeme, sparse solver veya FEM doğruluğu hakkında sonuç üretmez.
