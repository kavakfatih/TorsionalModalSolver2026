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

Lokal torsional eleman matrisi regresyonu, `k = 100 N·m/rad` için bilinen 2x2
katsayıları doğrudan doğrular. Ek fizik invariantları şunlardır:

- karşılıklı bağlantı için matris simetrisi,
- ortak rijit-cisim dönmesi için sıfır satır toplamı ve sıfır enerji,
- bağıl dönme için `theta^T Ke theta > 0`,
- negatif rijitliğin üretim matris yordamı tarafından reddedilmesi.

Bu test global assembly veya özdeğer çözümü yapmaz; yalnız elemanın lokal
fiziksel katkısını sınar.

Global matrix assembly regresyonu, fiziksel node ID değerlerini matris indisi
olarak kullanmadan aşağıdaki davranışları `1e-10` mutlak toleransla doğrular:

- `k=100 N·m/rad` tek elemanın 2x2 global K matrisi,
- `k1=100`, `k2=200 N·m/rad` üç düğümlü zincirin analitik global K matrisi,
- serbest sistemde K simetrisi, sıfır satır toplamı ve rijit-cisim null modu,
- `J=[0.1,0.2,0.3] kg·m²` için diagonal global dönel M matrisi,
- kısıtlı düğümün equation ID sıfır ile homojen eliminasyonu,
- tamamen kısıtlı geçerli sistemde 0x0 M/K matrisleri.

Negatif matris boyutu, bilinmeyen node lookup, aralık dışı equation ID, negatif
atalet ve başka sisteme ait DOF haritası ayrı `WILL_FAIL` vakalarıdır.

V0.3.1 foundation doğrulaması, mevcut unit ve assembly regresyonlarının yerine
geçmeden üretim katmanlarını tek bir analitik kanıt zincirinde birleştirir.
`test_torsional_validation` aşağıdaki kontrolleri `1e-10` seviyesinde uygular:

- tek elemanın `K_e=k[[1,-1],[-1,1]]` katsayıları, işaretleri ve simetrisi,
- serbest iki düğüm için `K[1,1]^T=0` rijit-cisim kalıntısı,
- farklı fiziksel node ID değerleriyle üç düğümlü global K katkı toplamı,
- tamamen serbest ve ilk düğümü kısıtlı DOF mapping sürekliliği,
- Frobenius normuyla K simetrisi ve `U=1/2 theta^T K theta>=0` enerji koşulu,
- iki ataletli analitik frekansın mevcut solver, assembled M/K modal residual,
  Rayleigh quotient ve kütle ortogonalliği ile çapraz doğrulanması.

Bu test bilinmeyen özdeğer aramaz; bilinen analitik modun Rayleigh değerini
hesaplar. Dolayısıyla yeni eigen solver davranışı oluşturmaz. Kapsam,
analytical code verification niteliğindedir ve deneysel model validation
çalışmasının yerine geçmez. Ayrıntılar
[`../validation/torsional_validation.md`](../validation/torsional_validation.md)
belgesindedir.

## Integration test

Entegrasyon testleri, birden fazla modülün birlikte kullanımını doğrular.
`test_generalized_torsional_system`, node/element koleksiyon yönetimini ve
mevcut iki ataletli TVD'nin konservatif alt modelinin genel topolojiye kayıpsız
dönüşümünü sınar. Benchmark 004 için `J_h`, `J_r` ve K' aktarılır; boyutsal
olarak farklı K'' değerinin viskoz `c` alanına aktarılmadığı doğrulanır.

Genel sistem testleri ayrıca sabitlenmemiş düğümlerin aktif DOF sayısını,
fixed-hub dönüşümünde göbek kısıtını, yinelenen kimlikleri ve tanımsız eleman
uçlarını kapsar. Bu testler birim testlerinin yerini almaz. Global assembly
`test_matrix_assembly` içinde ayrı doğrulanır; eigen çözümü uygulanmaz.

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
