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

## Integration test

Entegrasyon testleri, birden fazla modülün birlikte kullanımını doğrular.
Örneğin gelecekte geometri, malzeme ve modal çözüm bileşenlerinin aynı modelde
etkileşimini test eder. Bu testler yalnızca modüller arası davranış oluştuğunda
eklenir; birim testlerinin yerini almaz.

## Benchmark test

Benchmark testleri, temsilî TVD modellerinde yürütme süresi ve bellek kullanım
eğilimlerini izler. Başlangıçta `benchmarks/` altında tekrarlanabilir senaryolar
olarak tutulur. Donanım farklılıkları nedeniyle benchmark sonuçları henüz ana
CI geçiş koşulu değildir; performans regresyonu izleme amacı taşır.

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
