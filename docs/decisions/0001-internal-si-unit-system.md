# Karar 0001: İç SI Birim Sistemi

- Durum: Kabul edildi
- Tarih: 2026-08-23

## Bağlam

TVD verileri uygulama ve kullanıcı kaynaklarında mm, MPa ve derece gibi farklı
mühendislik birimleriyle gelebilir. Birim bilgisi alan adlarında veya hesap
adımlarında belirsiz kalırsa ölçek hataları fiziksel olarak yanlış sonuçlara yol
açabilir.

## Karar

Hesap motorunun çekirdek veri türleri yalnızca SI değerleri saklayacaktır:

- uzunluk: metre (`m`)
- yoğunluk: kilogram/metreküp (`kg/m³`)
- gerilme ve modül: paskal (`Pa`)
- sıcaklık: kelvin (`K`)
- frekans: hertz (`Hz`)
- açı: radyan (`rad`)

Fiziksel alan adları uygun olduğunda birim son eki taşıyacaktır. Mühendislik
birimleri çekirdek türlere yazılmadan önce `tms_units` ile dönüştürülecektir.

## Sonuçlar

- Çekirdek hesaplarda örtük birim dönüşümü yapılmayacaktır.
- Birim kaynaklı ölçek hatalarının incelenmesi kolaylaşacaktır.
- Girdi ve kullanıcı arayüzü katmanları birim dönüşümünden sorumlu olacaktır.
- V0.1.1 türleri fiziksel geçerlilik kontrolü yapmayacaktır; doğrulama daha sonra
  ayrı bir davranış olarak eklenecektir.
