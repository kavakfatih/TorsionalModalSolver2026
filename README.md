# TorsionalModalSolver2026 (TMS26)

TMS26, elastomer esaslı burulma titreşimi sistemleri için geliştirilen bir
mühendislik hesaplama yazılımıdır. Projenin hesap motoru modern Fortran 2018
ile geliştirilecek; derleme ve test süreçleri CMake ile yönetilecektir.

Güncel geliştirme sürümü `0.1.1` temel matematik ve fizik veri altyapısını
içerir. Fiziksel çözüm algoritmaları, FEM, kullanıcı arayüzü ve DXF desteği
henüz uygulanmamıştır.

## V0.1.1 çekirdek kapsamı

- `tms_kinds`: taşınabilir çift hassasiyetli `dp` türü
- `tms_constants`: pi ve temel mühendislik birim dönüşüm sabitleri
- `tms_units`: mm → m, MPa → Pa ve derece → radyan dönüşümleri
- `tms_geometry`: elastomer, atalet halkası, göbek ve bileşik TVD geometrisi
- `tms_material`: dinamik elastomer malzemenin çalışma noktası verileri

Hesap motorunun iç veri sözleşmesi SI birimlerini kullanır. Uzunluk metre,
yoğunluk `kg/m³`, kayma modülleri Pa, sıcaklık K ve frekans Hz cinsinden
saklanır. Dışarıdan alınan mühendislik birimleri, veri yapılarına yazılmadan
önce `tms_units` yordamlarıyla dönüştürülmelidir.

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

- `docs/`: mimari, matematik, geliştirme ve karar kayıtları
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

## Lisans

Projenin lisans koşulları henüz belirlenmemiştir.
