# Değişiklik Günlüğü

Bu dosyada TMS26 projesindeki kullanıcıya dönük önemli değişiklikler kaydedilir.

Biçim, [Keep a Changelog](https://keepachangelog.com/tr-TR/1.1.0/) yaklaşımını
temel alır ve proje anlamsal sürümleme ilkelerini izlemeyi hedefler.

## [Yayımlanmamış]

### Düzeltildi

- Windows CI iş akışında MinGW64 GNU Fortran paketi ve araç zinciri PATH
  doğrulaması eklendi.

## [0.1.2] - 2026-08-23

### Eklendi

- macOS ve Windows üzerinde CMake, Ninja ve GNU Fortran ile derleme/test yapan
  GitHub Actions CI iş akışları.
- Birim, entegrasyon, benchmark ve CI doğrulama kapsamını tanımlayan test
  stratejisi belgesi.
- Homojen annüler halka için kütle ve polar kütle atalet momenti hesabı.
- Lineer elastik annüler elastomer için burulma rijitliği hesabı.
- Tek serbestlik dereceli, sönümsüz sistem için doğal frekans hesabı.
- Üç torsional fizik yordamı için yüzde 0,1 hata sınırına sahip analitik testler.
- Basit annüler TVD referans benchmark tanımı ve beklenen sonuçları.

## [0.1.1] - 2026-08-23

### Eklendi

- Pi ve mühendislik birim dönüşüm sabitlerini sağlayan `tms_constants` modülü.
- mm → m, MPa → Pa ve derece → radyan dönüşümlerini sağlayan saf ve eleman bazlı
  `tms_units` yordamları.
- Elastomer, atalet halkası, göbek ve bileşik TVD geometrisi veri türleri.
- Dinamik elastomer için yoğunluk, G', G'', sıcaklık ve frekans veri türü.
- Sabitler, birimler, geometri ve malzeme modülleri için CTest testleri.
- SI tabanlı iç birim sözleşmesi ile çekirdek veri modeli belgeleri.

## [0.1.0] - 2026-08-23

### Eklendi

- Fortran 2018, CMake, Ninja ve CTest tabanlı başlangıç altyapısı.
- İlk tür tanımları modülü ve derleyici doğrulama testi.
- Temel proje dokümantasyonu ve geliştirme kuralları.
