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
