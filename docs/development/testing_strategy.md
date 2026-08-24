# TMS26 Test Stratejisi

## Amaç

TMS26 test yaklaşımı, hesap motorunun fiziksel veri sözleşmelerini, sayısal
yordamlarını ve platformlar arası derlenebilirliğini erken aşamada doğrulamayı
amaçlar. Her davranış değişikliği, uygun test seviyesi ve CI doğrulaması ile
birlikte teslim edilmelidir.

## Unit test

Birim testleri, tek bir Fortran modülünün veya yordamının beklenen davranışını
izole biçimde doğrular. Dönüşüm, veri türü, sınır durum ve hata düzeltmesi gibi
davranışlar CTest üzerinden ayrı çalıştırılabilir test hedefleri olarak eklenir.
Fiziksel veya matematiksel hesap yapan yordamlar, denklemin bilinen sonuçlarını
ve birim sözleşmesini kapsayan testlere sahip olmalıdır.

## Fizik regresyon testi

Fizik regresyon testleri, analitik referans değerlere ek olarak modelin
ölçekleme davranışlarını korur. Bu kontroller üretim denklemini test kodunda
yeniden kurmaz; yalnızca farklı girdilerle alınan üretim fonksiyonu sonuçlarını
karşılaştırır.

Tam bağlı annüler TVD kauçuk burç modeli için aşağıdaki ölçekleme
değişmezliği kontrolleri `1e-10` bağıl hata sınırıyla uygulanır:

- Modül ölçekleme: Geometri ile diğer malzeme alanları sabitken
  `K'(2G') = 2K'(G')` olmalıdır.
- Eksenel genişlik ölçekleme: Malzeme ve yarıçaplar sabitken
  `K'(2L) = 2K'(L)` ve `K''(2L) = 2K''(L)` olmalıdır.

## Integration test

Entegrasyon testleri, birden fazla modülün birlikte kullanımını doğrular.
`test_generalized_torsional_system`, node/element koleksiyon yönetimini ve
mevcut iki ataletli TVD'nin konservatif alt modelinin genel topolojiye kayıpsız
dönüşümünü sınar. Benchmark 004 için `J_h`, `J_r` ve K' aktarılır; boyutsal
olarak farklı K'' değerinin viskoz `c` alanına aktarılmadığı doğrulanır.

Genel sistem testleri ayrıca sabitlenmemiş düğümlerin aktif DOF sayısını,
fixed-hub dönüşümünde göbek kısıtını, yinelenen kimlikleri ve tanımsız eleman
uçlarını kapsar. Bu testler birim testlerinin yerini almaz ve matris assembly
veya eigen çözümü uygulamaz.

## Benchmark test

Benchmark testleri, temsilî TVD modellerinde yürütme süresi ve bellek kullanım
eğilimlerinin yanı sıra analitik referans sonuçların korunmasını sağlar.
Başlangıç senaryoları `benchmarks/` altında tekrarlanabilir girdiler ve beklenen
sonuçlar olarak tutulur. Donanım farklılıklarına bağlı performans değerleri
henüz ana CI geçiş koşulu değildir; analitik referanslar ilgili CTest fizik
doğrulamalarında sınanır.

## CI validation

GitHub Actions, `main` ve `develop` dallarına yapılan push işlemlerinde ve tüm
pull request'lerde aşağıdaki işlemleri yapar:

1. Kaynak kodunu alır.
2. macOS üzerinde Homebrew GNU Fortran; Windows üzerinde MinGW64 GNU Fortran
   araç zincirini kurar.
3. CMake ve Ninja ile ayrı `build/` dizininde yapılandırır ve derler.
4. CTest test takımını `--output-on-failure` seçeneğiyle çalıştırır.

Bir CI işi başarısız olursa değişiklik birleştirilmeden önce hata yerelde
tekrar üretilmeli, düzeltilmeli ve ilgili test eklenmelidir.

## Definition of Done doğrulaması

Her geliştirme görevi, önceden üretilmiş modül ve nesne dosyalarının sonucu
gizlemesini önlemek için temiz bir Debug build ile doğrulanır:

```sh
rm -rf build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

Yerel kalite kapıları geçmeden commit veya push yapılmaz. Push sonrasında macOS
ve Windows GitHub Actions sonuçları commit SHA ile eşleştirilerek kontrol edilir;
CI başarısızsa görev tamamlanmış sayılmaz.

Yeni Fortran modüllerinde kaynak kaydı, modül bağımlılık sırası ve üretilen
`.mod` dosyalarının tüketici hedeflere erişimi temiz build çıktısıyla
doğrulanır. Saf matematik yordamlarının `pure` niteliği korunur ve yeni fizik
hesapları analitik CTest kapsamına alınır.
