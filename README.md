# TorsionalModalSolver2026 (TMS26)

TMS26, elastomer esaslı burulma titreşimi sistemleri için geliştirilen bir
mühendislik hesaplama yazılımıdır. Projenin hesap motoru modern Fortran 2018
ile geliştirilecek; derleme ve test süreçleri CMake ile yönetilecektir.

Güncel geliştirme sürümü `0.8.0`, doğrulanmış reference master curve'ü
tabulated `log10(a_T)`, WLF veya Arrhenius temperature-shift provider'ıyla
bileştirir. Canonical koordinat dönüşümü `f_r=a_T(T)f` log-space'te
doğrulanır; physical `f,T` ile reduced lookup `f_r,T_ref` trace içinde ayrı
tutulur. Material-aware dynamic stiffness
`Z_r(f)=K'_r(f)-omega^2 M_r+i(K''_r(f)+omega C_r)` biçimindedir. Dense LAPACK
`ZSYSVX` mevcut complex-symmetric reference backend olarak yeniden kullanılır;
V0.5 `DSYGV` modal ve V0.6 frozen harmonic yolları değişmeden korunur.

## V0.8.0 thermorheological runtime kapsamı

- `tms_kinds`: taşınabilir çift hassasiyetli `dp` türü
- `tms_constants`: pi ve temel mühendislik birim dönüşüm sabitleri
- `tms_units`: mm → m, MPa → Pa ve derece → radyan dönüşümleri
- `tms_geometry`: elastomer, atalet halkası, göbek ve bileşik TVD geometrisi
- `tms_dynamic_modulus`: G', G'', frekans ve sıcaklık ile kayıp faktörü hesabı
- `tms_material_frequency`: frekans-sıcaklık bağımlı malzeme veri noktası
- `tms_material`: tek çalışma noktası alanları ve dinamik veri noktaları
- `tms_dynamic_material_metadata`: dataset operating-state ve test traceability
- `tms_dynamic_modulus_provider`: genişletilebilir constitutive sorgu sınırı
- `tms_tabulated_dynamic_modulus_provider`: linear/log-frequency interpolation,
  measured-isotherm ve no-extrapolation sözleşmesi
- `tms_temperature_shift_provider`: ortak shift evaluation, validation ve
  explicit sıcaklık-domain sözleşmesi
- `tms_tabulated_temperature_shift`: sıcaklık ekseninde linear `log10(a_T)`
  interpolation ve no-extrapolation
- `tms_wlf_temperature_shift`: positive `C1/C2` ve pole-korumalı WLF modeli
- `tms_arrhenius_temperature_shift`: SI birimli activation-energy modeli
- `tms_thermorheological_dynamic_modulus_provider`: master curve ile shift
  modelini mevcut dynamic modulus provider API'sinde bileştiren runtime katmanı
- `tms_inertia`: homojen annüler göbek ve halkanın kütle özellikleri
- `tms_torsional_stiffness`: lineer elastomer bölgenin burulma rijitliği
- `tms_dynamic_torsional_stiffness`: K', K'' ve kayıp faktörü hesabı
- `tms_dynamic_torsional_property_binding`: eleman ID, provider ve önceden
  hesaplanmış bonded-annular geometri katsayısı bağlantısı
- `tms_frequency_solver`: tek serbestlik dereceli doğal frekans hesabı
- `tms_local_matrix`: iki uçlu elemanlar için 2x2 lokal matris veri taşıyıcısı
- `tms_matrix_types`: private allocatable depolamalı genel dense matris türü
- `tms_dof_types`: anlamlı DOF türü ve `(node_id,dof_type)` Physical DOF kimliği
- `tms_dof_map`: Physical DOF ile constraint'ten bağımsız tam Equation ID
- `tms_stiffness_matrix`: global torsional storage K' matrisi
- `tms_loss_stiffness_matrix`: ayrı global torsional loss K'' matrisi
- `tms_damping_matrix`: ayrı global viskoz torsional C matrisi
- `tms_mass_matrix`: düğüm polar ataletlerinden global diagonal M matrisi
- `tms_matrix_assembly`: topoloji ve full K'/K''/C/M assembly
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
- `tms_reduced_dynamic_system`: reduced K'/K''/C/M ve complex recovery bağlamı
- `tms_dynamic_stiffness`: `Z=K'-omega^2M+i(K''+omega C)` oluşturucusu
- `tms_harmonic_excitation`: Physical DOF tabanlı complex peak torque assembly
- `tms_complex_linear_problem` / `tms_complex_linear_solution`: backend-neutral
  çoklu-RHS doğrusal problem, status ve diagnostics sözleşmesi
- `tms_lapack_zsysvx_backend`: LP64 complex-symmetric dense reference backend
- `tms_complex_linear_solver`: ZSYSVX ayrıntısını gizleyen solver facade'ı
- `tms_harmonic_response`: sweep status, response, diagnostics ve TVD çıktıları
- `tms_harmonic_analysis`: direct frequency sweep çözüm orkestrasyonu
- `tms_material_aware_harmonic_analysis`: prevalidated frequency-dependent
  K'/K'' override, mixed/multiple provider assembly ve ZSYSVX orkestrasyonu
- `tms_material_aware_harmonic_response` / `tms_material_state_trace`: V0.6
  harmonic sonuç bileşimi ile dataset, G*/K* ve interpolation izleri
- `tms_frf`: tek tanımlı torque input channel için rotational FRF yardımcıları
- `tms_torsional_node`: genel düğüm kimliği, polar atalet ve sınır koşulu
- `tms_torsional_element`: ayrı K', K'' ve c bağlantı kanalları
- `tms_generalized_torsional_system`: koleksiyon yönetimi, aktif DOF sayımı ve
  topoloji doğrulaması
- `tms_torsional_system`: TVD builder'ı ile fixed-hub ve serbest-serbest
  analitik modal çözüm ve genel topolojiye geriye uyumlu dönüşüm

Hesap motorunun iç veri sözleşmesi SI birimlerini kullanır. Uzunluk metre,
yoğunluk `kg/m³`, kayma modülleri Pa, sıcaklık K ve frekans Hz cinsinden
saklanır. Dışarıdan alınan mühendislik birimleri, veri yapılarına yazılmadan
önce `tms_units` yordamlarıyla dönüştürülmelidir.

## Geliştirme Durumu

TMS26 şu anda V0.8.0 aşamasındadır. Dinamik elastomer master curve'ü,
temperature-shift provider'ları, bonded-annular dynamic element binding'i,
generalized node-element topolojisi, constraint reduction, sönümsüz modal
analiz ve complex harmonic response aynı çekirdekte kullanılabilir. K'' ile
boyutsal olarak farklı viskoz `c` ayrı eleman alanları, global matrisler ve
doğrulama kanalları olarak korunur.

Lineer elemanların lokal K', K'' ve C katkıları ile düğüm polar ataletleri önce
constraint'ten bağımsız Full Equation ID uzayında birleştirilir. Constraint
manager aynı retained indekslerle `K'_r/K''_r/C_r/M_r` matrislerini üretir.
V0.5 modal yolu `K_r phi=lambda M_r phi` denklemini çözmeye devam eder. V0.6
frozen harmonic yolu her explicit pozitif frekansta `Z_r theta_hat=T_hat`
denklemini doğrudan çözer. V0.7 material-aware yol aynı ZSYSVX backend'ini
kullanır; yalnız bound elemanların `K'(f),K''(f)` katkılarını provider'dan
günceller. V0.8 bu provider sınırını değiştirmeden, externally prescribed
operating sıcaklıkta reduced-frequency master-curve lookup'u ekler. Her iki yol
`theta_hat=P theta_hat_r` ile Physical DOF uzayına geri açılır; prescribed
statik offset harmonic phasor'a eklenmez.

Harmonic sonuç **direct**, **full-order**, **linear**, **small-amplitude** ve
**frequency-domain** kapsamındadır. Complex genlikler `exp(+i*omega*t)`
konvansiyonunda peak amplitude değerleridir; RMS değildir. V0.6 API'sinde
özellikler frozen kalır. Ayrı V0.7 API'si `LINEAR_FREQUENCY` veya seçimlik
`LINEAR_LOG_FREQUENCY` ile tek measured isotherm üzerinde G'/G'' interpolate
eder. V0.8 provider'ı validated `T_ref` master curve ve shift modelini
horizontal-only TTS ile bileştirir. Temperature ve reduced-frequency
extrapolation, vertical shift, self-heating ve master-curve fitting yapmaz.

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
- K', K'' ve viskoz C için ayrı lokal/global matris türleri, full assembly ve
  aynı retained indekslerle direct reduction
- Explicit, sonlu, pozitif ve kesin artan frequency grid ile complex peak
  nodal torque excitation assembly
- `Z_r=K'_r-omega^2M_r+i(K''_r+omega C_r)` complex-symmetric dynamic stiffness
- `FACT='N'`, `UPLO='U'` ve workspace query kullanan LP64 LAPACK `ZSYSVX`
  backend'i; çoklu RHS'ye açık backend-neutral complex solver facade'ı
- `SOLVED`, `SOLVED_ILL_CONDITIONED` ve `SINGULAR` frequency-point durumları;
  singular noktada uydurma response olmadan sweep'e devam etme
- RCOND, RHS başına FERR/BERR ve backend-independent scaled relative residual
- Reduced/physical complex response, magnitude, phase, velocity ve acceleration
- Oriented TVD relative angle, complex element torque, transmitted magnitude,
  average dissipated power ve dissipated energy-per-cycle yardımcıları
- Tek tanımlı torque input channel için rotational receptance, mobility ve
  accelerance; genel forced response ile FRF'nin açık ayrımı
- Viscous-only, K''-only, combined damping, fixed/free two-inertia, passivity,
  exact singular, ill-conditioned ve V0.5 modal cross-validation regresyonları
- Canonical SI `G'(f),G''(f)` tablosu için private-copy provider abstraction,
  explicit dataset/test metadata ve SHEAR-only direct torsional contract
- Default linear-frequency ve seçimlik linear-log-frequency axis interpolation;
  exact-point machine tolerance, measured-isotherm ve no-extrapolation kuralları
- Bound eleman nominal K'/K'' değerlerini frozen yol için koruyan, material-aware
  yolda double-count etmeden dynamic K'/K'' ile override eden ayrı binding
- Mixed constant/dynamic sistemler, birden çok provider, full-sweep
  prevalidation ve frequency-independent M/C/DOF/constraint preparation
- Solved veya singular bütün requested frequency noktalarında dataset kimliği,
  G'/G'', tan(delta), K'/K'' ve interpolation bracket/alpha material trace'i
- Exact/interpolated provider, invalid data, annular mapping, 1-DOF full chain,
  complex 2x2 mixed/multiple system ve singular trace analitik regresyonları
- `a_T=tau(T)/tau(T_ref)` ve `f_r=a_Tf` convention'ıyla tabulated
  `log10(a_T)`, WLF ve Arrhenius temperature-shift provider'ları
- Shift sıcaklık domain'i ile reduced-frequency master-curve domain'ini ayrı
  doğrulayan, çarpım taşması yaratmayan log-space runtime evaluation
- Physical frequency/operating temperature ile reduced lookup
  frequency/reference temperature'ı ayıran geriye uyumlu material trace
- Mevcut binding ve `analyze_material_aware_harmonic_response()` API'sini
  yeniden kullanan thermorheological 1-DOF, multiple-provider, passivity ve
  singular-trace analitik regresyonları
- Analitik referans testleri ve annüler TVD benchmark'ları
- macOS LP64 LAPACK ve Windows LP64 OpenBLAS sağlayıcılarıyla GitHub Actions
  derleme/test iş akışları
- Mimari, matematik, fizik, geliştirme ve karar belgeleri için dizin indeksleri

### Henüz kapsam dışında

- Spline/PCHIP, curve fitting, Prony serisi ve nonlinear elastomer modeli
- Sönümlü kompleks özdeğer problemi ve mode-superposition harmonic çözüm
- Master-curve identification, empirical isotherm alignment, TRS quality
  metrics, vertical shift, amplitude/prestrain interpolation ve self-heating
- Sıfırdan farklı dynamic prescribed dönme için RHS correction ve reaction
  torque recovery
- Sparse matris depolaması, CSR ve production iterative eigen solver
- Sparse complex direct solver, GMRES/BiCGSTAB ve complex Krylov yöntemleri
- Lanczos, Block Lanczos, Krylov–Schur, LOBPCG, ARPACK ve SLEPc backend'leri
- MPC, constraint equation, Lagrange multiplier, penalty ve contact
- Frekansa bağlı malzemeyle nonlinear/öz-tutarlı modal eigenproblem iterasyonu
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

LAPACK, `K_r phi=lambda M_r phi` gerçek simetrik genelleştirilmiş özdeğer
problemini `DSYGV`, complex-symmetric `Z_r theta_hat=T_hat` problemini ise
`ZSYSVX` ile çözmek için kullanılır. Kaynak kod veya platforma özgü
`-llapack/-lblas` yolları repoya gömülmez; CMake sağlayıcıyı
`find_package(LAPACK REQUIRED)` ile bulur ve `LAPACK::LAPACK` hedefini bağlar.
Her iki backend default 32-bit Fortran `INTEGER` kullanan LP64 ABI gerektirir;
ILP64/OpenBLAS64 bu sürümün sözleşmesiyle uyumlu değildir.

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
  projeyi derler ve DSYGV ile ZSYSVX kullanan otomatik testleri çalıştırır.
- Otomatik testler, iki platformda da `ctest --output-on-failure` ile raporlanır.

İş akışı tanımları [macOS CI](.github/workflows/macos-build.yml) ve
[Windows CI](.github/workflows/windows-build.yml) dosyalarındadır.

## Dizin yapısı

- `docs/`: mimari, malzeme, matematik, fizik, doğrulama, geliştirme ve karar kayıtları
- `engine/src/`: Fortran hesap motoru kaynakları
- `engine/src/constraints/`: constraint veri modeli, yönetimi ve reduced system
- `engine/src/eigen/`: backend-neutral eigenproblem/solution ve DSYGV backend
- `engine/src/harmonic/`: dynamic stiffness, complex solver, frozen/material-aware
  orchestration ve harmonic sonuç/trace katmanları
- `engine/src/materials/`: legacy malzeme türleri, dataset metadata'sı ve
  dynamic modulus provider'ları
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

V0.6.0 direct harmonic katman akışı
[`docs/architecture/V0.6_frequency_domain_response.md`](docs/architecture/V0.6_frequency_domain_response.md),
phasor denklemleri ve enerji/FRF bağıntıları
[`docs/mathematics/harmonic_torsional_response.md`](docs/mathematics/harmonic_torsional_response.md),
K'/K''/C fiziksel ayrımı
[`docs/physics/torsional_damping_models.md`](docs/physics/torsional_damping_models.md),
analitik ve solver-status doğrulamaları
[`docs/validation/harmonic_response_validation.md`](docs/validation/harmonic_response_validation.md),
complex-symmetric ZSYSVX backend kararı ise
[`docs/decisions/0011-frequency-domain-complex-solver.md`](docs/decisions/0011-frequency-domain-complex-solver.md)
belgelerinde açıklanır.

V0.7.0 provider/binding ve material-aware assembly akışı
[`docs/architecture/V0.7_dynamic_material_provider.md`](docs/architecture/V0.7_dynamic_material_provider.md),
dataset/test metadata ve DMA→TVD aktarım sınırları
[`docs/materials/tabulated_dynamic_elastomer_material.md`](docs/materials/tabulated_dynamic_elastomer_material.md),
interpolation, isotherm, passivity ve causality matematiği
[`docs/mathematics/dynamic_modulus_interpolation.md`](docs/mathematics/dynamic_modulus_interpolation.md),
analitik doğrulama kapıları
[`docs/validation/dynamic_material_provider_validation.md`](docs/validation/dynamic_material_provider_validation.md),
kalıcı provider kararı ise
[`docs/decisions/0012-tabulated-dynamic-material-provider.md`](docs/decisions/0012-tabulated-dynamic-material-provider.md)
belgelerinde açıklanır.

V0.8.0 thermorheological provider bileşimi ve trace akışı
[`docs/architecture/V0.8_thermorheological_runtime.md`](docs/architecture/V0.8_thermorheological_runtime.md),
dynamic elastomer varsayımları
[`docs/materials/thermorheological_dynamic_elastomer.md`](docs/materials/thermorheological_dynamic_elastomer.md),
tabulated/WLF/Arrhenius denklemleri
[`docs/mathematics/temperature_shift_functions.md`](docs/mathematics/temperature_shift_functions.md),
analitik doğrulama kapıları
[`docs/validation/thermorheological_runtime_validation.md`](docs/validation/thermorheological_runtime_validation.md),
canonical convention kararı ise
[`docs/decisions/0013-thermorheological-runtime-convention.md`](docs/decisions/0013-thermorheological-runtime-convention.md)
belgelerinde açıklanır. Uygulama ve doğrulama özeti
[`docs/development/V0.8.0_development_report.md`](docs/development/V0.8.0_development_report.md)
altında tutulur.

## Lisans

Projenin lisans koşulları henüz belirlenmemiştir.
