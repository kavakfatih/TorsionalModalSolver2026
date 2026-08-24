# İki Ataletli TVD Torsional Sistem Modeli

## Fiziksel düzen

V0.2.2 modeli üç eş merkezli bileşeni tek bir ayrık-parametreli sistemde
birleştirir:

```text
Hub
 |
 | tam bağlı annüler elastomer
 |
Inertia Ring
```

Göbek ile halka homojen rijit gövdeler, elastomer ise bağıl dönmeye karşı
moment taşıyan lineer viskoelastik bağlantı olarak idealize edilir. Sistem
builder'ı `tvd_geometry_t`, iki rijit gövde yoğunluğu ve
`dynamic_rubber_material_t` verilerini mevcut atalet ve dinamik rijitlik
yordamlarına aktarır; bu denklemleri sistem katmanında yeniden uygulamaz.
Elastomerin dağıtılmış kütlesi ve polar ataleti iki rijit gövdeli kütle
matrisinde ihmal edilir.

## Fixed-hub sınır koşulu

Göbek dış yapıya ideal biçimde sabitlenir. Tek aktif koordinat atalet halkası
dönmesi `theta_r` olur:

```text
J_r theta_r'' + K' theta_r = 0
```

Bu model, elastomerin depoladığı enerji ile halkanın kinetik enerjisi
arasındaki sönümsüz salınımı temsil eder. Doğal frekans mevcut tek-DOF
frekans yordamıyla hesaplanır.

## Serbest-serbest sınır koşulu

Göbek ve halka dış dünyaya bağlanmadığında iki dönel serbestlik derecesi
vardır. Kütle ve rijitlik matrisleri:

```text
M = [ J_h   0   ]       K = [ K'  -K' ]
    [  0   J_r  ]           [-K'   K' ]
```

Sistem iki fiziksel mod üretir:

1. Rijit-cisim modu: Göbek ve halka `[1, 1]` şekliyle birlikte döner.
   Elastomer şekil değiştirmediği için geri çağırıcı moment ve frekans sıfırdır.
2. Elastik bağıl mod: Göbek genliği `1` seçildiğinde halka genliği
   `-J_h/J_r` olur. Zıt yönlü hareket elastomeri burar ve pozitif bir doğal
   frekans oluşturur.

V0.2.2 bu iki modu analitik çözer. Genel matris özdeğer çözümü, seyrek matris
altyapısı veya FEM kullanılmaz.

## K' ve K'' kullanım sınırı

Kompleks burulma rijitliği:

```text
K* = K' + i K''
```

Burada `K'` çevrim sırasında depolanan elastik enerjiyi, `K''` ise kaybedilen
enerjiyi temsil eder. V0.2.2 doğal frekans hesabı reel ve sönümsüzdür;
bu nedenle hareket denkleminin rijitlik teriminde yalnız `K'` kullanılır.

`K''` ve `tan(delta) = K''/K'` sistem veri tipinde korunur, fakat bunlardan
sönüm oranı, sönümlü doğal frekans veya kompleks özdeğer türetilmez. `K''`
değerini reel bir özdeğer denklemine doğrudan eklemek, bu sürümde bulunmayan
kompleks modal formülasyonu örtük biçimde varsayacağı için fiziksel olarak
uygun değildir.

## Frozen-property modal tahmin

Elastomerin `G'` ve `G''` değerleri belirli bir malzeme referans frekansı ve
sıcaklığında tanımlıdır. Sistem builder'ı bu çalışma noktasını sonuç verisine
aktarır. Modal çözüm sırasında:

- `G'` ve dolayısıyla `K'` sabit tutulur.
- Hesaplanan doğal frekans ile malzeme referans frekansı karşılaştırılarak
  interpolasyon yapılmaz.
- Yeni bir `G'(f)` seçmek için iterasyon veya öz-tutarlı bağlaşım yapılmaz.

Sonuç bu nedenle **frozen-property undamped modal estimate** niteliğindedir.
Frekans bağımlı malzeme interpolasyonu ve iteratif bağlaşım ileriki bir sürümün
konusudur.

## Geçerlilik koşulları ve SI birimleri

- `J_h` ve `J_r` pozitif olmalı ve `kg·m²` cinsinden verilmelidir.
- `K'` pozitif olmalı ve `N·m/rad` cinsinden verilmelidir.
- Göbek ve halka geometrileri `0 <= ri < ro`, `L > 0`; yoğunluklar
  `kg/m³` cinsinden pozitif olmalıdır.
- Elastomer geometrisi, tam bağlı annüler burç modelinin koşullarını sağlamalıdır.
- Elastomer kütlesi ve polar ataleti ihmal edilir; kütle matrisi yalnız
  göbek ve atalet halkası ataletlerinden oluşur.
- Frekans `Hz`, sıcaklık `K`, açısal genlikler `rad` cinsindendir; normalize
  mod şekli bileşenleri boyutsuzdur.
- Model lineer, küçük deformasyonlu, eksenel simetrik ve dış zorlamasızdır.

Analitik türetim
[`../mathematics/two_inertia_modal_model.md`](../mathematics/two_inertia_modal_model.md),
sayısal doğrulama ise
[`../../benchmarks/004_two_inertia_tvd/`](../../benchmarks/004_two_inertia_tvd/)
altındadır. Model, standart ayrık-parametreli torsional vibration teorisine
dayanır.
