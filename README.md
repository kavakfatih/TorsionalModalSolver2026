# TorsionalModalSolver2026 (TMS26)

TMS26, elastomer esaslı burulma titreşimi sistemleri için geliştirilen bir
mühendislik hesaplama yazılımıdır. Projenin hesap motoru modern Fortran 2018
ile geliştirilecek; derleme ve test süreçleri CMake ile yönetilecektir.

Güncel geliştirme sürümü `0.2.0`, elastomerlerin frekans ve sıcaklığa bağlı
kompleks kayma modülü verilerini temsil edecek mühendislik temelini sağlar.
Bu sürüm veri saklama ve kayıp faktörü hesabıyla sınırlıdır; interpolasyon,
eğri uydurma ve kompleks burulma rijitliği çözümü henüz uygulanmamıştır.

## V0.2.0 çekirdek kapsamı

- `tms_kinds`: taşınabilir çift hassasiyetli `dp` türü
- `tms_constants`: pi ve temel mühendislik birim dönüşüm sabitleri
- `tms_units`: mm → m, MPa → Pa ve derece → radyan dönüşümleri
- `tms_geometry`: elastomer, atalet halkası, göbek ve bileşik TVD geometrisi
- `tms_dynamic_modulus`: G', G'', frekans ve sıcaklık ile kayıp faktörü hesabı
- `tms_material_frequency`: frekans-sıcaklık bağımlı malzeme veri noktası
- `tms_material`: tek çalışma noktası alanları ve dinamik veri noktaları
- `tms_inertia`: homojen annüler halkanın kütlesi ve polar kütle ataleti
- `tms_torsional_stiffness`: lineer elastomer bölgenin burulma rijitliği
- `tms_frequency_solver`: tek serbestlik dereceli doğal frekans hesabı

Hesap motorunun iç veri sözleşmesi SI birimlerini kullanır. Uzunluk metre,
yoğunluk `kg/m³`, kayma modülleri Pa, sıcaklık K ve frekans Hz cinsinden
saklanır. Dışarıdan alınan mühendislik birimleri, veri yapılarına yazılmadan
önce `tms_units` yordamlarıyla dönüştürülmelidir.

## Geliştirme Durumu

TMS26 şu anda V0.2.0 aşamasındadır. Hesap motorunun temel veri sözleşmesi, ilk
analitik torsional hesap zinciri ve dinamik elastomer veri modeli kullanılabilir.
Kompleks modül yaklaşımı bu sürümde çözüm algoritmasına bağlanmamıştır.

### Tamamlananlar

- Fortran 2018, CMake, Ninja ve CTest tabanlı derleme/test altyapısı
- SI birim dönüşümleri ile geometri ve dinamik elastomer veri türleri
- Frekans ve sıcaklık çalışma noktalarında G' ve G'' saklama altyapısı
- `tan(delta) = G'' / G'` kayıp faktörü hesabı ve analitik testi
- Homojen annüler atalet halkası için kütle ve polar kütle ataleti hesabı
- Lineer elastomer bölgesi için burulma rijitliği hesabı
- Tek serbestlik dereceli, sönümsüz doğal frekans hesabı
- Analitik referans testleri ve basit annüler TVD benchmark'ı
- macOS ve Windows için GitHub Actions derleme/test iş akışları
- Mimari, matematik, fizik, geliştirme ve karar belgeleri için dizin indeksleri

### Henüz kapsam dışında

- İnterpolasyon, eğri uydurma, Prony serisi ve nonlinear elastomer modeli
- Kompleks burulma rijitliği ve G'' kullanımına dayalı dinamik çözüm
- Çok serbestlik dereceli modal/eigen çözümü
- FEM, DXF ve Qt tabanlı kullanıcı arayüzü

Her fiziksel hesap genişletmesi, ilgili matematik belgesi, otomatik test ve
benchmark güncellemesiyle birlikte eklenir.

## Gereksinimler

- Fortran 2018 destekli bir derleyici
  - macOS: güncel GNU Fortran (`gfortran`) veya Intel `ifx`
  - Windows: güncel GNU Fortran (`gfortran`) veya Intel `ifx`
- CMake 3.25 veya üzeri
- Ninja
- VS Code (önerilen geliştirme ortamı)

## Derleme ve test

Tek yapılandırmalı Ninja üreticisiyle macOS, Linux, Git Bash veya PowerShell
üzerinde aşağıdaki komutlar kullanılabilir:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

### Commit öncesi yerel kontrol

Bir değişikliği commit etmeden önce, temiz bir `build/` dizininde aşağıdaki
komutlar çalıştırılmalıdır:

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

Birden fazla Fortran derleyicisi kuruluysa yapılandırma sırasında derleyici
açıkça seçilebilir:

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_COMPILER=gfortran
```

Windows PowerShell üzerinde Intel oneAPI derleyicisi kullanılacaksa komutlar,
oneAPI ortamı etkinleştirildikten sonra `-DCMAKE_Fortran_COMPILER=ifx` seçeneği
ile çalıştırılabilir.

## Continuous Integration

GitHub Actions, `main` ve `develop` dallarına yapılan her push ile tüm pull
request'lerde derleme ve test doğrulaması yapar.

- macOS iş akışı Homebrew ile GNU Fortran, CMake ve Ninja kurar; projeyi derler
  ve CTest testlerini çalıştırır.
- Windows iş akışı MSYS2 MinGW64 ortamında GNU Fortran, CMake ve Ninja kurar;
  projeyi derler ve CTest testlerini çalıştırır.
- Otomatik testler, iki platformda da `ctest --output-on-failure` ile raporlanır.

İş akışı tanımları [macOS CI](.github/workflows/macos-build.yml) ve
[Windows CI](.github/workflows/windows-build.yml) dosyalarındadır.

## Dizin yapısı

- `docs/`: mimari, matematik, fizik, geliştirme ve karar kayıtları
- `engine/src/`: Fortran hesap motoru kaynakları
- `engine/tests/`: hesap motoru testleri
- `benchmarks/`: performans ölçümleri
- `examples/`: örnek kullanım senaryoları
- `third_party/`: üçüncü taraf bileşenler
- `tools/`: geliştirme yardımcı araçları

## Geliştirme ilkeleri

Kaynak kod Fortran 2018 standardına uygun yazılır. Kod açıklamaları Türkçedir;
fiziksel veya matematiksel hesap yapan her yordam, dayandığı modeli,
varsayımları, girdileri, çıktıları ve birimleri açıklamak zorundadır. Her yeni
davranış test ve ilgili dokümantasyon güncellemesiyle birlikte eklenir. Ayrıntılı
kurallar için [`AGENTS.md`](AGENTS.md) dosyasına bakın.

Birim sözleşmesinin gerekçesi `docs/decisions/0001-internal-si-unit-system.md`,
dönüşüm denklemleri ise `docs/mathematics/unit-conversions.md` altında
belgelenmiştir.

İlk torsional fizik modelinin denklemleri
[`docs/mathematics/torsional-physics-core.md`](docs/mathematics/torsional-physics-core.md),
tek serbestlik dereceli frekans modeli
[`docs/mathematics/torsional_frequency_model.md`](docs/mathematics/torsional_frequency_model.md),
analitik referans problemi ise
[`benchmarks/001_simple_annular_tvd/`](benchmarks/001_simple_annular_tvd/)
altında açıklanmıştır. Dinamik elastomer veri örneği
[`benchmarks/002_dynamic_elastomer/`](benchmarks/002_dynamic_elastomer/),
kompleks modül modeli ise
[`docs/mathematics/dynamic_elastomer_model.md`](docs/mathematics/dynamic_elastomer_model.md)
altında açıklanmıştır.

## Lisans

Projenin lisans koşulları henüz belirlenmemiştir.
