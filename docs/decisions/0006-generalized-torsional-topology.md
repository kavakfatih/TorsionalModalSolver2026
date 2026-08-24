# Karar 0006: Genelleştirilmiş Torsional Düğüm-Eleman Topolojisi

- Durum: Kabul edildi
- Tarih: 2026-08-24

## Bağlam

V0.2.2 yalnız göbek ve atalet halkasına özel iki ataletli bir sistem türü ve
analitik solver sağlıyordu. Gelecekteki M/K matrix assembly ve genel modal
çözüm için düğüm sayısından bağımsız, doğrulanabilir bir topoloji sözleşmesine
ihtiyaç vardır. Mevcut `tms_torsional_system` modül adı ve public iki-ataletli
API'si geriye uyumluluk için korunmalıdır.

## Karar

- Genel veri modeli `torsional_node_t`, `torsional_element_t` ve
  `torsional_system_t` türlerinden oluşacaktır.
- Yeni ana veri modülü, mevcut modül adıyla çakışmayı önlemek için
  `tms_generalized_torsional_system` adını kullanacak; dosyası istenen
  `engine/src/system/tms_torsional_system.f90` yolunda bulunacaktır.
- Sistem koleksiyonları private olacak; ekleme, sayma, okuma ve doğrulama public
  saf yordamlarla yapılacaktır.
- Aktif DOF sayısı sabitlenmemiş düğüm sayısıdır. Kimlikler DOF indeksi değildir.
- Paralel elemanlar, ayrık alt sistemler, tam sabitlenmiş topolojiler ve
  serbest-serbest rijit-cisim davranışı veri modeli tarafından yasaklanmaz.
- Mevcut `tms_torsional_system` iki-ataletli API'si korunacak ve genel topolojiye
  tek yönlü bir dönüşüm yordamı ekleyecektir.
- İki-atalet dönüşümü `J_h`, `J_r` ve `K'` değerlerini taşır. `K''` kayıp
  rijitliği viskoz sönüm katsayısı değildir; genel eleman sönümü sıfır bırakılır.
- V0.2.3 matris assembly veya özdeğer çözümü içermeyecektir.

## Sonuçlar

- Gelecekteki matris ve modal katmanlar bileşene özel TVD türlerinden bağımsız
  bir node-element topolojisini tüketebilir.
- Mevcut fixed-hub ve serbest-serbest analitik solver API'leri değişmeden kalır.
- Topolojik invariantlar dışarıdan doğrudan koleksiyon mutasyonuyla atlanamaz.
- Kayıp rijitliğinden viskoz sönüme dönüşüm için ileride frekans ve model
  seçimini açıkça tanımlayan ayrı bir karar gerekir.
- Benchmark 004 sonuçları değişmez; yalnız aynı fiziksel sistemin genel
  topolojideki gösterimi eklenir.
