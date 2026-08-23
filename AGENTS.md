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

- Build başarısızsa commit yapılmaz.
- Test başarısızsa push yapılmaz.
- Fiziksel veya matematiksel hesap değişikliklerinde dokümantasyon güncellenir.
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
