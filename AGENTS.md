# TMS26 Geliştirme Kuralları

Bu kurallar deponun tamamı için geçerlidir.

## Dil ve standart

- Tüm Fortran kaynakları Fortran 2018 standardına uygun olmalıdır.
- Her program birimi `implicit none` kullanmalıdır.
- Taşınabilir tür tanımları merkezi modüllerden alınmalı; derleyiciye özgü tür
  numaraları doğrudan kullanılmamalıdır.
- Kod açıklamaları, dokümantasyon yorumları ve geliştiriciye dönük açıklamalar
  Türkçe yazılmalıdır. Tanımlayıcı adları için tutarlı ve anlaşılır İngilizce
  adlar kullanılabilir.

## Matematik ve fizik açıklamaları

- Fiziksel veya matematiksel hesap yapan her fonksiyon ve alt yordam için
  açıklama zorunludur.
- Açıklama; kullanılan modeli veya denklemi, varsayımları, girdileri, çıktıları,
  fiziksel birimleri ve varsa geçerlilik sınırlarını belirtmelidir.
- Kaynak gösterilmesi gereken yöntemler `docs/mathematics/` altında ayrıntılı
  olarak belgelenmeli ve kod açıklamasından bu belgeye yönlendirme yapılmalıdır.

## Testler

- Yeni veya değiştirilen her davranış uygun bir otomatik testle kapsanmalıdır.
- Hata düzeltmeleri, hatayı yeniden üreten bir regresyon testi içermelidir.
- Testler CTest üzerinden çalıştırılmalı ve desteklenen platformlarda aynı sonucu
  vermelidir.
- Bir değişiklik tamamlanmadan önce yapılandırma, derleme ve test komutları
  çalıştırılmalıdır.

## Dokümantasyon

- Mimari değişikliklerde `docs/architecture/`, matematiksel değişikliklerde
  `docs/mathematics/`, geliştirme süreci değişikliklerinde `docs/development/`
  güncellenmelidir.
- Önemli ve kalıcı teknik kararlar `docs/decisions/` altında karar kaydı olarak
  belgelenmelidir.
- Kullanıcıya dönük önemli değişiklikler `CHANGELOG.md` dosyasına eklenmelidir.

## Derleme sistemi

- Derleme ve test tanımları CMake ile yönetilmelidir; Ninja önerilen üreticidir.
- Kaynak dosyaları hedef tabanlı CMake komutlarıyla ilgili hedefe eklenmelidir.
- Derleyiciye veya platforma özgü seçenekler koşullu ve mümkün olduğunca dar
  kapsamlı tutulmalıdır.
- Kaynak ağacında derleme çıktısı üretilmemeli; ayrı bir `build/` dizini
  kullanılmalıdır.

## Definition of Done

Bir geliştirme görevi, aşağıdaki kalite kapılarının tamamı başarıyla geçmeden
tamamlanmış kabul edilmez.

1. Proje kökündeki önceki build dizini silinerek temiz başlangıç hazırlanır:

   ```sh
   rm -rf build
   ```

2. Debug yapılandırması CMake ve Ninja ile oluşturulur:

   ```sh
   cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
   ```

3. Tüm hedefler tam olarak derlenir:

   ```sh
   cmake --build build
   ```

4. Tüm testler hata ayrıntılarıyla çalıştırılır:

   ```sh
   ctest --test-dir build --output-on-failure
   ```

5. Yerel doğrulama sonrasında commit ve push yapılır; push ile tetiklenen
   GitHub Actions işlerinin sonucu kontrol edilir.

Yerel yapılandırma, build veya test adımlarından biri başarısız olursa:

- Commit yapılmaz.
- Push yapılmaz.
- Görev tamamlandı olarak raporlanmaz.

Push sonrasında GitHub Actions başarısız olursa görev tamamlandı olarak
raporlanmaz; hata giderilir ve Definition of Done akışı temiz build ile baştan
uygulanır.

## Fortran Module Development Rules

Yeni bir Fortran modülü eklendiğinde aşağıdaki kurallar uygulanır:

1. Kaynak dosyası ilgili hedefin `CMakeLists.txt` içindeki `target_sources`
   listesine eklenir.
2. Modül bağımlılık sırası sağlayıcı modülden tüketici modüle doğru kontrol
   edilir.
3. `use` ile alınan tüm modüllerin tüketici derlenmeden önce üretildiği, CMake
   bağımlılık taraması ve temiz Ninja build ile doğrulanır.
4. Yan etkisiz matematiksel yordamlar `pure`, eleman bazında uygulanabilenler
   uygun olduğunda `elemental` olarak tanımlanır. Saf yordamlar I/O yapmaz ve
   global durum değiştirmez.
5. Her yeni fiziksel veya matematiksel hesap için analitik referanslı otomatik
   test eklenir ve CTest'e kaydedilir.

## GitHub Workflow Automation

Tüm geliştirme görevlerinde aşağıdaki sıra uygulanacaktır:

1. Kod değişikliklerini gerçekleştir.
2. CMake configure çalıştır.
3. Projeyi derle.
4. CTest çalıştır.
5. Testler başarılı ise `git status` kontrol et.
6. Conventional Commit formatında commit oluştur.
7. GitHub remote kontrol et.
8. GitHub'a push yap.
9. GitHub Actions sonuçlarını kontrol et.
10. Sonucu raporla.

### Kurallar

- CMake configure ve build adımları her geliştirme görevinde zorunludur.
- CTest çalıştırılması ve tüm testlerin geçmesi zorunludur.
- Build başarısızsa commit yapılmaz.
- Test başarısızsa push yapılmaz.
- Commit ve push sonrasında GitHub Actions sonuçları kontrol edilmelidir.
- Fiziksel veya matematiksel hesap değişikliklerinde `docs/mathematics/`
  altındaki ilgili dokümantasyonun güncellenmesi zorunludur.
- Her commit tek bir amacı temsil eder.

### Commit örnekleri

- `feat: yeni özellik`
- `fix: hata düzeltme`
- `test: test ekleme`
- `docs: dokümantasyon`
- `ci: CI/CD değişikliği`

### Görev sonunda rapor

Görev sonunda aşağıdaki bilgiler raporlanmalıdır:

- Değişen dosyalar
- Build sonucu
- Test sonucu
- Commit hash
- Push sonucu
- GitHub Actions sonucu
