# Karar 0002: İlk Torsional Fizik Modeli

- Durum: Kabul edildi
- Tarih: 2026-08-23

## Bağlam

TMS26 hesap motorunun ileride çok serbestlik dereceli ve daha gelişmiş elastomer
modellerini desteklemesi hedeflenir. İlk fizik sürümünün, daha karmaşık çözüm
yöntemleri eklenmeden önce birim sözleşmesini ve modül sınırlarını doğrulayacak
analitik bir temel sağlaması gerekir.

## Karar

V0.1.2 aşağıdaki en küçük fizik modelini kullanacaktır:

- Atalet halkası homojen, eksenel simetrik annüler katıdır.
- Elastomer homojen, izotrop ve lineer elastiktir; G' kayma modülü kullanılır.
- Doğal frekans modeli sönümsüz ve tek serbestlik derecelidir.
- Hesap yordamları mevcut geometri ve malzeme türlerini kullanır; yeni derived
  type tanımlamaz.
- Girdiler çekirdek yordamlara SI birimleriyle verilir.

## Sonuçlar

- Formüller bağımsız analitik testlerle doğrulanabilir.
- Geometri, malzeme ve fizik sorumlulukları ayrı modüllerde kalır.
- Girdi doğrulama bu sürümde yordamların önkoşuludur; ayrı doğrulama API'si
  henüz yoktur.
- Sönüm, G'', nonlinear davranış, FEM ve eigen çözümü sonraki sürümlere bırakılır.
