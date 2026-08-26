# TorsionalModalSolver2026 (TMS26)

TMS26, elastomer esaslı burulma titreşimi sistemleri için geliştirilen bir
mühendislik hesaplama yazılımıdır. Projenin hesap motoru modern Fortran 2018
ile geliştirilecek; derleme ve test süreçleri CMake ile yönetilecektir.

Güncel geliştirme sürümü `0.5.0`, V0.4 constraint foundation tarafından
üretilen indirgenmiş `K_r/M_r` sistemini gerçek simetrik genelleştirilmiş
özdeğer problemi olarak çözer. Dense LAPACK `DSYGV`, küçük modeller için
doğrulanabilir reference backend olarak kullanılır; modal API backend
ayrıntılarını açmaz ve gelecekteki sparse/Lanczos-family çözücülere hazırlanır.

## V0.5.0 generalized modal eigen solver kapsamı

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
- `tms_matrix_types`: private allocatable depolamalı genel dense matris türü
- `tms_dof_types`: anlamlı DOF türü ve `(node_id,dof_type)` Physical DOF kimliği
- `tms_dof_map`: Physical DOF ile constraint'ten bağımsız tam Equation ID
- `tms_stiffness_matrix`: global torsional K matrisi ve lokal katkı toplama
- `tms_mass_matrix`: düğüm polar ataletlerinden global diagonal M matrisi
- `tms_matrix_assembly`: topoloji, tam DOF haritası ve full M/K assembly
- `tms_matrix_reduction`: aktif denklem haritasıyla storage-bağımsız indirgeme
- `tms_constraint_types`: fixed ve prescribed value constraint veri modeli
- `tms_constraint_manager`: constraint doğrulaması ve Active Equation ID haritası
- `tms_reduced_system`: Kr/Mr ile tam fiziksel sonuç recovery bilgisi
- `tms_generalized_eigen_problem`: sonlu ve simetrik K/M problem sözleşmesi
- `tms_eigen_solution`: backend-neutral eigenvalue/eigenvector taşıyıcısı
- `tms_lapack_dsygv_backend`: LP64 LAPACK DSYGV reference backend'i
- `tms_generalized_eigen_solver`: backend ayrıntısını gizleyen solver facade'ı
- `tms_modal_validation`: rijit mod, residual ve kütle ortogonalliği kontrolü
- `tms_modal_result`: frekans, sınıflandırma, tanı ve mod şekli sonuçları
- `tms_modal_analysis`: reduced sistemden çözüme ve fiziksel mode recovery'ye
  uzanan modal analiz orkestrasyonu
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

TMS26 şu anda V0.5.0 aşamasındadır. Dinamik elastomer ve kompleks burulma
rijitliği hesabı, annüler rijit gövde ataletleri ve iki sınır koşullu analitik
TVD modeli kullanılabilir. Aynı iki ataletli sistem artık iki düğüm ve bir
eleman olarak genel topolojide temsil edilebilir. K'' sistem verisinde korunur;
boyutsal olarak farklı viskoz sönüm katsayısına doğrudan dönüştürülmez.

Lineer elemanların lokal K katkıları ve düğüm polar ataletleri önce
constraint'ten bağımsız tam Equation ID uzayında full K/M matrislerine
birleştirilir. Constraint manager bundan sonra ayrı Active Equation ID
haritasını kurar; direct elimination Kr/Mr matrislerini üretir. V0.5.0 modal
analiz yolu bu matrisleri `K_r phi = lambda M_r phi` biçiminde çözer, modları
ölçeğe duyarlı biçimde sınıflandırır ve `phi=P phi_r` ile Physical DOF uzayına
geri açar. Prescribed statik offset mode shape'e eklenmez.

Modal sonuç **linear**, **undamped** ve **frozen-property** kapsamındadır.
Elastomer rijitliği hangi frekans-sıcaklık çalışma noktasından üretildiyse o
değer çözüm boyunca sabit kabul edilir; `G'(f) -> K(f) -> eigenfrequency`
öz-tutarlı iterasyonu yapılmaz.

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
- Fiziksel düğüm kimliklerinden bağımsız, kısıtları destekleyen DOF haritası
- Private dense depolama ile boyut ve sonlu katsayı doğrulaması
- Açık lokal→denklem→global dönüşümüyle dense global torsional K assembly
- Düğümde yığılmış polar ataletlerden diagonal global M assembly
- Homojen sıfır dönme kısıtı için indirgenmiş aktif M/K matrisleri
- Tek eleman işareti, rijit-cisim modu, DOF mapping, üç düğümlü assembly,
  K simetrisi/enerjisi ve iki ataletli analitik frekansı birlikte doğrulayan
  V0.3.1 foundation regresyonu
- Atalet, yoğunluk, yarıçap, uzunluk, rijitlik, modül, frekans ve sıcaklık
  girdileri için IEEE sonluluk doğrulaması
- Göbek-elastomer ve elastomer-atalet halkası ara yüzlerinde toleranslı geometri
  sürekliliği doğrulaması
- Uç atalet/rijitlik oranları için taşmaya dirençli eşdeğer atalet ve doğal
  frekans değerlendirmesi
- Sabitlenmiş tek uçta indirgenmiş `K=[k]` ve `M=[J_free]` regresyonu
- `(node_id,dof_type)` Physical DOF kimliği ve `TORSIONAL_ROTATION` türü
- Constraint'ten bağımsız tam Equation ID ile ayrı Active Equation ID haritası
- Fixed ve prescribed value kayıtlarını doğrulayan constraint manager
- Constraint'ten bağımsız full K/M assembly ve storage-bağımsız direct
  elimination ile Kr/Mr üretimi
- Tüm DOF'lar fixed olduğunda geçerli `0x0` indirgenmiş sistem
- `q=Pq_r+q_p` durum recovery ve `phi=Pphi_r` modal recovery eşlemeleri
- Node sırası değişimine karşı Physical DOF tabanlı indirgeme regresyonu
- Constraint-aware `K_r/M_r` matrislerini doğrudan alan generalized
  eigenproblem sözleşmesi
- `ITYPE=1`, `JOBZ='V'`, `UPLO='U'` kullanan LP64 LAPACK `DSYGV` dense
  reference backend'i ve platformdan bağımsız CMake bağlantısı
- Original K/M matrislerini değiştirmeyen workspace-query ve working-copy
  çözüm akışı
- Serbest sistemlerde singular/positive-semidefinite K'yi ve birden çok
  rigid-body mode'u geçerli kabul eden ölçeğe duyarlı sınıflandırma
- `phi^T M phi=1` kütle normalizasyonu, boyutsuz relative residual ve
  `Phi^T M Phi` ortogonallik tanıları
- Mode shape işaret belirsizliği ile repeated eigenvalue eigenspace
  eşdeğerliğini dikkate alan regresyonlar
- Fixed 1-DOF, free-free iki atalet, üç düğümlü zincir, constrained recovery,
  repeated mode, invalid K/M ve tamamen constrained sistem testleri
- `phi=Pphi_r` ile constrained bileşenleri sıfır bırakan fiziksel mode recovery
- Analitik referans testleri ve annüler TVD benchmark'ları
- macOS LP64 LAPACK ve Windows LP64 OpenBLAS sağlayıcılarıyla GitHub Actions
  derleme/test iş akışları
- Mimari, matematik, fizik, geliştirme ve karar belgeleri için dizin indeksleri

### Henüz kapsam dışında

- İnterpolasyon, eğri uydurma, Prony serisi ve nonlinear elastomer modeli
- Kompleks rijitlik kullanan sönümlü doğal frekans veya frekans cevabı çözümü
- Global C/damping assembly ve sıfırdan farklı prescribed dönme için RHS
  düzeltmesi veya yük çözümü
- Sparse matris depolaması, CSR ve production iterative eigen solver
- Lanczos, Block Lanczos, Krylov–Schur, LOBPCG, ARPACK ve SLEPc backend'leri
- MPC, constraint equation, Lagrange multiplier, penalty ve contact
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
- LP64 (32-bit Fortran `INTEGER`) arayüzlü LAPACK ve BLAS
- VS Code (önerilen geliştirme ortamı)

LAPACK, `K_r phi = lambda M_r phi` gerçek simetrik genelleştirilmiş özdeğer
problemini `DSYGV` ile çözmek için kullanılır. Kaynak kod veya platforma özgü
`-llapack/-lblas` yolları repoya gömülmez; CMake sağlayıcıyı
`find_package(LAPACK REQUIRED)` ile bulur ve `LAPACK::LAPACK` hedefini bağlar.
ILP64/OpenBLAS64 bu sürümün ABI sözleşmesiyle uyumlu değildir.

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

### macOS LAPACK kurulumu

Homebrew reference LAPACK paketi macOS sistem `Accelerate.framework` ile
çakışmaması için keg-only kurulur. Cellar sürümü veya Apple Silicon/Intel yolu
elle yazılmamalı; güncel prefix `brew --prefix` ile alınmalıdır:

```sh
brew install gcc cmake ninja lapack
lapack_prefix="$(brew --prefix lapack)"
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DBLA_SIZEOF_INTEGER=4 \
  -DBLA_VENDOR=Generic \
  -DCMAKE_PREFIX_PATH="$lapack_prefix"
cmake --build build
ctest --test-dir build --output-on-failure
```

CMake, açık Homebrew prefix'i verilmediğinde uyumlu bir sistem LP64 sağlayıcısı
bulabilir. Yukarıdaki komut CI ile aynı reference LAPACK seçimini tekrarlar.

### Windows MSYS2 MinGW64 LAPACK kurulumu

MSYS2 `MINGW64` shell içinde standard LP64 OpenBLAS paketi kullanılmalıdır:

```sh
pacman -S --needed mingw-w64-x86_64-gcc \
  mingw-w64-x86_64-gcc-fortran \
  mingw-w64-x86_64-cmake \
  mingw-w64-x86_64-ninja \
  mingw-w64-x86_64-openblas
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_Fortran_COMPILER=gfortran \
  -DBLA_SIZEOF_INTEGER=4 \
  -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH=/mingw64
cmake --build build
ctest --test-dir build --output-on-failure
```

`mingw-w64-x86_64-openblas64` paketi ILP64 olduğu için kullanılmamalıdır.

## Continuous Integration

GitHub Actions, `main` ve `develop` dallarına yapılan her push ile tüm pull
request'lerde derleme ve test doğrulaması yapar.

- macOS iş akışı Homebrew ile GNU Fortran, CMake, Ninja ve keg-only LP64
  LAPACK kurar; provider prefix'ini dinamik belirler.
- Windows iş akışı MSYS2 MinGW64 ortamında GNU Fortran, CMake, Ninja ve LP64
  OpenBLAS/LAPACK kurar.
- İki iş akışı da LAPACK provider'ını configure sırasında CMake ile doğrular,
  projeyi derler ve DSYGV kullanan otomatik testleri çalıştırır.
- Otomatik testler, iki platformda da `ctest --output-on-failure` ile raporlanır.

İş akışı tanımları [macOS CI](.github/workflows/macos-build.yml) ve
[Windows CI](.github/workflows/windows-build.yml) dosyalarındadır.

## Dizin yapısı

- `docs/`: mimari, matematik, fizik, doğrulama, geliştirme ve karar kayıtları
- `engine/src/`: Fortran hesap motoru kaynakları
- `engine/src/constraints/`: constraint veri modeli, yönetimi ve reduced system
- `engine/src/eigen/`: backend-neutral eigenproblem/solution ve DSYGV backend
- `engine/src/modal/`: modal analiz, doğrulama, sonuç ve fiziksel recovery
- `engine/src/matrix/`: lokal/global matris, DOF haritası, assembly ve reduction
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
genel M/K bağlantısının matematiksel temeli
[`docs/mathematics/generalized_torsional_system_model.md`](docs/mathematics/generalized_torsional_system_model.md)
ve kalıcı mimari karar
[`docs/decisions/0006-generalized-torsional-topology.md`](docs/decisions/0006-generalized-torsional-topology.md)
altında bulunur. Lokal eleman rijitlik türetimi
[`docs/mathematics/local_torsional_element_matrix.md`](docs/mathematics/local_torsional_element_matrix.md),
matris taşıyıcısı ve bağımlılık kararı ise
[`docs/decisions/0007-local-element-matrix-design.md`](docs/decisions/0007-local-element-matrix-design.md)
altında bulunur. DOF eşlemesi ile global M/K assembly matematiği
[`docs/mathematics/global_matrix_assembly.md`](docs/mathematics/global_matrix_assembly.md),
kalıcı mimari karar ise
[`docs/decisions/0008-global-matrix-assembly-design.md`](docs/decisions/0008-global-matrix-assembly-design.md)
altında bulunur. V0.4.0 constraint ve reduced-system sorumlulukları
[`docs/architecture/V0.4_constraint_foundation.md`](docs/architecture/V0.4_constraint_foundation.md),
direct-elimination matematiği
[`docs/mathematics/constraint_reduction.md`](docs/mathematics/constraint_reduction.md)
ve bu mimarinin kalıcı kararı
[`docs/decisions/0009-constraint-reduction-architecture.md`](docs/decisions/0009-constraint-reduction-architecture.md)
altında açıklanır. V0.3.1 analitik foundation doğrulama modeli, toleransları ve
sonuç sözleşmesi
[`docs/validation/torsional_validation.md`](docs/validation/torsional_validation.md)
belgesinde açıklanır. V0.3.2 sonlu girdi, geometri ara yüzü ve uç ölçek
doğrulamaları
[`docs/validation/numerical_robustness_validation.md`](docs/validation/numerical_robustness_validation.md)
belgesinde açıklanır.

V0.5.0 modal katman akışı ve backend sınırı
[`docs/architecture/V0.5_modal_eigen_solver.md`](docs/architecture/V0.5_modal_eigen_solver.md),
genelleştirilmiş özdeğer denklemi, kütle normalizasyonu ve residual tanımları
[`docs/mathematics/generalized_modal_eigenproblem.md`](docs/mathematics/generalized_modal_eigenproblem.md),
analitik ve sayısal doğrulama kapsamı
[`docs/validation/modal_eigen_validation.md`](docs/validation/modal_eigen_validation.md),
DSYGV reference backend ile gelecek sparse/Lanczos-family facade kararı ise
[`docs/decisions/0010-generalized-eigen-solver-backend.md`](docs/decisions/0010-generalized-eigen-solver-backend.md)
belgesinde açıklanır.

## Lisans

Projenin lisans koşulları henüz belirlenmemiştir.
