# Değişiklik Günlüğü

Bu dosyada TMS26 projesindeki kullanıcıya dönük önemli değişiklikler kaydedilir.

Biçim, [Keep a Changelog](https://keepachangelog.com/tr-TR/1.1.0/) yaklaşımını
temel alır ve proje anlamsal sürümleme ilkelerini izlemeyi hedefler.

## [Yayımlanmamış]

### Değiştirildi

- Geliştirme görevleri için temiz Debug build, tam derleme, CTest ve commit
  sonrası GitHub Actions kontrolünü zorunlu kılan Definition of Done eklendi.
- Yeni Fortran modülleri için CMake kaydı, bağımlılık sırası, `pure` matematik
  yordamı ve fizik testi kuralları standartlaştırıldı.
- Annüler elastomer polar alan momenti geometri testinde doğrudan analitik
  referansla doğrulanacak şekilde test kapsamı güçlendirildi.
- Dinamik burulma rijitliği analitik testinin bağıl hata sınırı `1e-10`
  seviyesine sıkılaştırıldı ve fiziksel girdi sınırları CTest'e eklendi.
- V0.2.1.3 bakımında dinamik rijitlik regresyon kapsamı, `K' ∝ G'`,
  `K' ∝ L` ve `K'' ∝ L` ölçekleme kontrolleriyle tamamlandı.

### Düzeltildi

- TVD elastomer rijitliğinde eksen boyunca Saint-Venant burulmasına ait
  `GJp/ℓ` denklemi yerine, rijit göbek ve dış halkaya tam bağlı annüler
  kauçuk burç için `4πGLri²ro²/(ro²-ri²)` denklemi uygulanmaya başlandı.
- Statik ve dinamik solver'lar `m³` birimli ortak `Cθ` geometri faktörüne
  bağlandı; Benchmark 001/003 ve doğal frekans referansları düzeltildi.
- `ri = 0` değerinin tam bağlı silindirik burç modeli için fiziksel olarak
  geçersiz olduğu belgelenerek girdi doğrulamasına eklendi.
- `calculate_rubber_polar_area_moment` yordamının public PURE arayüzü dış ve iç
  yarıçapı metre cinsinden alan iki skaler argümanla uyumlu hale getirildi.
- Dinamik rijitlik testinin üretim yordamı yerine yerel bir hesap kopyasını
  sınaması giderildi; test doğrudan üretim modülüne bağlandı.
- Negatif veya sırasız yarıçaplar, pozitif olmayan etkin uzunluk ve depolama
  modülü ile negatif kayıp modülünün geçersiz sonuç üretmesi engellendi.

## [0.2.1] - 2026-08-23

### Eklendi

- K', K'', kayıp faktörü, frekans ve sıcaklığı taşıyan
  `complex_torsional_stiffness_t` veri türü.
- Dinamik elastomer malzemesi ile annüler kauçuk geometrisinden kompleks
  burulma rijitliği hesaplayan `calculate_dynamic_torsional_stiffness` yordamı.
- `rubber_geometry_t` için yeniden kullanılabilir polar alan momenti hesabı.
- G''/G' ile K''/K' eşitliğini yüzde 0,1'den küçük hata sınırında doğrulayan
  CTest testi ve EPDM Benchmark 003 referansı.
- Kompleks rijitliğin fiziksel etkileri, matematiksel bağlantısı, mimari veri
  akışı ve API kararı için dokümantasyon.

### Değiştirildi

- Mevcut statik burulma rijitliği solver'ı davranışı korunarak ortak polar
  alan momenti yordamını kullanacak şekilde düzenlendi.
- Proje sürümü `0.2.1` olarak güncellendi ve CTest kapsamı on teste çıkarıldı.

## [0.2.0] - 2026-08-23

### Eklendi

- Kompleks dinamik kayma modülünün G', G'', frekans ve sıcaklık bileşenlerini
  saklayan `dynamic_shear_modulus` veri türü.
- DMA ve modal test verilerine hazırlanmak için `material_frequency_point`
  veri türü ve malzeme içinde çoklu çalışma noktası altyapısı.
- Boyutsuz `tan(delta) = G'' / G'` kayıp faktörü hesabı ve CTest birim testi.
- Dinamik elastomer matematik modeli, kompleks burulma rijitliği hazırlık
  belgesi, mimari karar kaydı ve EPDM referans benchmark'ı.

### Değiştirildi

- `dynamic_rubber_material_t`, mevcut tek noktalı alanları korunarak dinamik
  frekans-sıcaklık veri noktalarını saklayacak şekilde genişletildi.
- Proje sürümü `0.2.0` olarak güncellendi ve CTest kapsamı dokuz teste çıkarıldı.

## [0.1.3] - 2026-08-23

### Eklendi

- Mimari, matematik, fizik, geliştirme ve karar dokümantasyonu dizin indeksleri.
- Tek serbestlik dereceli torsional frekans modeli matematik belgesi.
- Unit test, physics validation ve benchmark regression kategorilerini açıklayan
  hesap motoru test belgesi.

### Değiştirildi

- GitHub otomasyon kuralları zorunlu build, test ve commit sonrası CI kalite
  kapılarıyla netleştirildi.
- Basit annüler TVD benchmark belgesine fiziksel model ve hesap zinciri eklendi.

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
