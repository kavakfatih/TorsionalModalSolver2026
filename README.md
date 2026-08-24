# TorsionalModalSolver2026 (TMS26)

TMS26, elastomer esaslı burulma titreşimi sistemleri için geliştirilen bir
mühendislik hesaplama yazılımıdır. Projenin hesap motoru modern Fortran 2018
ile geliştirilecek; derleme ve test süreçleri CMake ile yönetilecektir.

Güncel geliştirme sürümü `0.2.4`, mevcut fixed-hub ve serbest-serbest iki
ataletli TVD fiziğini korurken genel torsional elemanların 2x2 lokal rijitlik
matrisi katkısını üretmesini sağlar. Yığılmış polar atalet düğümleri, lineer
elemanlar ve lokal K katkısı gelecekteki global M/K assembly altyapısına
hazırdır. Bu sürümde global assembly veya genel özdeğer çözümü uygulanmamıştır.

## V0.2.4 çekirdek kapsamı

- `tms_kinds`: taşınabilir çift hassasiyetli `dp` türü
- `tms_constants`: pi ve temel mühendislik birim dönüşüm sabitleri
- `tms_units`: mm → m, MPa → Pa ve derece → radyan dönüşümleri
- `tms_geometry`: elastomer, atalet halkası, göbek ve bileşik TVD geometrisi
- `tms_dynamic_modulus`: G', G'', frekans ve sıcaklık ile kayıp faktörü hesabı
- `tms_material_frequency`: frekans-sıcaklık bağımlı malzeme veri noktası
- `tms_material`: tek çalışma noktası alanları ve dinamik veri noktaları
- `tms_inertia`: homojen annüler göbek ve halkanın kütle özellikleri
- `tms_torsional_stiffness`: lineer elastomer bölgenin burulma rijitliği
- `tms_dynamic_torsional_stiffness`: K', K'' ve kayıp faktörü hesabı
- `tms_frequency_solver`: tek serbestlik dereceli doğal frekans hesabı
- `tms_local_matrix`: iki uçlu elemanlar için 2x2 lokal matris veri taşıyıcısı
- `tms_torsional_node`: genel düğüm kimliği, polar atalet ve sınır koşulu
- `tms_torsional_element`: lineer K/c bağlantısı ve 2x2 lokal rijitlik hesabı
- `tms_generalized_torsional_system`: koleksiyon yönetimi, aktif DOF sayımı ve
  topoloji doğrulaması
- `tms_torsional_system`: TVD builder'ı ile fixed-hub ve serbest-serbest
  analitik modal çözüm ve genel topolojiye geriye uyumlu dönüşüm

Hesap motorunun iç veri sözleşmesi SI birimlerini kullanır. Uzunluk metre,
yoğunluk `kg/m³`, kayma modülleri Pa, sıcaklık K ve frekans Hz cinsinden
saklanır. Dışarıdan alınan mühendislik birimleri, veri yapılarına yazılmadan
önce `tms_units` yordamlarıyla dönüştürülmelidir.

## Geliştirme Durumu

TMS26 şu anda V0.2.4 aşamasındadır. Dinamik elastomer ve kompleks burulma
rijitliği hesabı, annüler rijit gövde ataletleri ve iki sınır koşullu analitik
TVD modeli kullanılabilir. Aynı iki ataletli sistem artık iki düğüm ve bir
eleman olarak genel topolojide temsil edilebilir. K'' sistem verisinde korunur;
boyutsal olarak farklı viskoz sönüm katsayısına doğrudan dönüştürülmez. Her
lineer eleman, global sisteme henüz birleştirilmeyen lokal K matrisini üretir.

### Tamamlananlar

- Fortran 2018, CMake, Ninja ve CTest tabanlı derleme/test altyapısı
- SI birim dönüşümleri ile geometri ve dinamik elastomer veri türleri
- Frekans ve sıcaklık çalışma noktalarında G' ve G'' saklama altyapısı
- `tan(delta) = G'' / G'` kayıp faktörü hesabı ve analitik testi
- Annüler mil/end-face burulması için ayrı polar alan momenti hesabı
- Tam bağlı annüler TVD kauçuk burcu için `Cθ` geometri faktörü
- `K' = G'Cθ` ve `K'' = G''Cθ` kompleks rijitlik bileşenleri
- `G''/G' = K''/K'` eşitliğinin analitik doğrulaması
- Homojen annüler göbek ve atalet halkası için kütle ve polar kütle ataleti
  hesabı
- Lineer elastomer bölgesi için burulma rijitliği hesabı
- Tek serbestlik dereceli, sönümsüz doğal frekans hesabı
- Geometri, yoğunluk ve dinamik malzemeyi birleştiren iki ataletli TVD builder'ı
- Fixed-hub 1-DOF doğal frekansı ve serbest-serbest iki analitik mod
- `[1,1]` rijit-cisim ve `[1,-J_h/J_r]` elastik mod şekilleri
- K', atalet ve büyük göbek ataleti limit regresyonları
- Genel torsional node/element veri türleri ve private sistem koleksiyonları
- Sabitlenmemiş düğümlerden aktif DOF sayımı ve topoloji doğrulaması
- İki ataletli TVD için serbest-serbest ve fixed-hub genel topoloji köprüsü
- `Ke = k[[1,-1],[-1,1]]` biçiminde 2x2 lokal eleman rijitlik katkısı
- Lokal matris simetrisi, sıfır satır toplamı ve pozitif yarı-tanımlılık testleri
- Analitik referans testleri ve annüler TVD benchmark'ları
- macOS ve Windows için GitHub Actions derleme/test iş akışları
- Mimari, matematik, fizik, geliştirme ve karar belgeleri için dizin indeksleri

### Henüz kapsam dışında

- İnterpolasyon, eğri uydurma, Prony serisi ve nonlinear elastomer modeli
- Kompleks rijitlik kullanan sönümlü doğal frekans veya frekans cevabı çözümü
- Global M/K/C matrix assembly ve sınır koşulu eliminasyonu
- Genel çok serbestlik dereceli modal/eigen çözümü
- Frekansa bağlı G'(f) interpolasyonu ve öz-tutarlı modal iterasyon
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
rm -rf build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
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
- `engine/src/matrix/`: küçük lokal matris veri türleri
- `engine/src/system/`: genel torsional node, element ve sistem topolojisi
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
altında açıklanmıştır. Kompleks rijitlik fiziksel modeli
[`docs/physics/complex_torsional_stiffness.md`](docs/physics/complex_torsional_stiffness.md),
EPDM analitik örneği ise
[`benchmarks/003_dynamic_torsional_stiffness/`](benchmarks/003_dynamic_torsional_stiffness/)
altında bulunur. İki ataletli analitik modal türetim
[`docs/mathematics/two_inertia_modal_model.md`](docs/mathematics/two_inertia_modal_model.md),
fiziksel sistem sınırları
[`docs/physics/two_inertia_torsional_system.md`](docs/physics/two_inertia_torsional_system.md)
ve sayısal referans
[`benchmarks/004_two_inertia_tvd/`](benchmarks/004_two_inertia_tvd/)
altında açıklanmıştır. Genel torsional sistemin fiziksel sözleşmesi
[`docs/physics/generalized_torsional_system.md`](docs/physics/generalized_torsional_system.md),
gelecekteki M/K bağlantısı
[`docs/mathematics/generalized_torsional_system_model.md`](docs/mathematics/generalized_torsional_system_model.md)
ve kalıcı mimari karar
[`docs/decisions/0006-generalized-torsional-topology.md`](docs/decisions/0006-generalized-torsional-topology.md)
altında bulunur. Lokal eleman rijitlik türetimi
[`docs/mathematics/local_torsional_element_matrix.md`](docs/mathematics/local_torsional_element_matrix.md),
matris taşıyıcısı ve bağımlılık kararı ise
[`docs/decisions/0007-local-element-matrix-design.md`](docs/decisions/0007-local-element-matrix-design.md)
altında bulunur.

## Lisans

Projenin lisans koşulları henüz belirlenmemiştir.
