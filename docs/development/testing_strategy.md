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

V0.3.2 numerik güvenilirlik regresyonu, fizik denklemlerini değiştirmeden girdi
ve sayı aralığı sınırlarını doğrular. `test_numerical_hardening` aşağıdaki
başlıkları kapsar:

- atalet, yoğunluk, yarıçap, uzunluk, rijitlik, modül, frekans ve sıcaklık için
  IEEE `NaN`, pozitif sonsuz ve negatif sonsuz girdilerin reddi,
- göbek-elastomer ve elastomer-atalet halkası ara yüzlerinde mutlak `1e-12 m`
  ile bağıl `1e-9` birleşik toleransının iç/dış sınır vakaları,
- uç atalet ve rijitlik oranlarında sonlu eşdeğer atalet ile doğal frekans,
- ilk düğümü sabit iki düğümlü sistemde `K=[k]` ve `M=[J_free]` indirgemesi.

`error stop` beklenen her geçersiz girdi ayrı CTest sürecinde çalışır ve
`WILL_FAIL` olarak kaydedilir. Test seçicileri alan ile IEEE sınıfını açıkça
adlandırır. Bilinmeyen seçici beklenen red sayılmaz; böylece CMake kaydı ile
test programı arasındaki yazım hatası yanlış başarıya dönüşmez.

Assertion yardımcıları; hesaplanan değer, beklenen değer ve toleransın sonlu
olmasını zorunlu tutar. Tolerans ayrıca negatif olamaz. Bu sözleşme, özellikle
`NaN` karşılaştırmalarının false dönmesi nedeniyle oluşabilecek sessiz test
başarılarını önler. Ayrıntılı doğrulama kapsamı
[`../validation/numerical_robustness_validation.md`](../validation/numerical_robustness_validation.md)
belgesindedir.

V0.4.0 constraint foundation regresyonu, V0.3 tam fizik ve assembly
invariantlarını korurken constraint'ten bağımsız tam denklem uzayı ile aktif
solver uzayını ayrı doğrular. `test_constraint_foundation` aşağıdaki zinciri
üretim yordamları üzerinden kapsar:

- Physical DOF `(node_id,dof_type)` ile tam Equation ID eşlemesinin constraint
  değiştiğinde korunması,
- iki düğüm ve tek elemanda ilk düğüm fixed iken `Kr=[k]` ve
  `Mr=[J_free]`,
- üç düğümlü zincirde tam sistemden iki aktif denkleme direct elimination,
- tüm DOF'lar fixed iken geçerli `0x0` Kr/Mr,
- node ekleme sırası değiştiğinde ortak Physical DOF sırasındaki eşdeğer K/M
  sonucu ve fiziksel kimlik tabanlı result recovery,
- `q=Pq_r+q_p` ile prescribed bileşenleri içeren tam vektör recovery ve
  `phi=Pphi_r` homojen modal recovery altyapısı,
- sıfırdan farklı prescribed değerin recovery metadata'sında korunurken
  principal Kr/Mr katsayılarına eklenmemesi,
- test assertion ve recovery girdilerinde IEEE `NaN` değerinin reddedilmesi.

Aynı Physical DOF kümesini farklı sırada taşıyan başka bir sisteme ait full
DOF haritası, K/M satırlarının yanlış elenmesini önlemek için map-system
bütünlük doğrulamasında reddedilir.

Olmayan node, aynı Physical DOF için yinelenen constraint, geçersiz DOF türü
ve geçersiz prescribed değerler ayrı
`tms26.constraint_foundation.rejects_*` CTest süreçlerinde `WILL_FAIL` olarak
kaydedilir. Bilinmeyen test seçicisi beklenen fiziksel red sayılmaz.

Bu regresyon Kr/Mr üretimini doğrular; özdeğer aramaz, LAPACK çağırmaz ve
sıfırdan farklı prescribed değer için RHS düzeltmesi çözmez. Direct-elimination
matematiği
[`../mathematics/constraint_reduction.md`](../mathematics/constraint_reduction.md),
modül sorumlulukları ise
[`../architecture/V0.4_constraint_foundation.md`](../architecture/V0.4_constraint_foundation.md)
belgelerindedir.

V0.5.0 modal doğrulaması, V0.4 `K_r/M_r` çıktısını production generalized
eigen ve modal-analysis yordamları üzerinden çözer. Ana CTest programları
`test_generalized_eigen_solver` ve `test_modal_analysis` olup aşağıdaki
birbirini tamamlayan kanıtları sağlar:

- fixed 1-DOF sistemde `lambda=k/J`, frekans, mass normalization ve residual,
- free-free iki atalette bir rigid ile bir elastic mode ve mevcut analitik
  solver ile sign-invariant çapraz doğrulama,
- eş üç düğümlü zincirde `[0,k/J,3k/J]` eigenvalue sırası,
- constrained zincirde reduced çözüm ve `phi=Pphi_r` physical recovery,
- repeated eigenvalue için tekil eigenvector yerine multiplicity ve eigenspace
  eşdeğerliği,
- ayrık serbest alt sistemlerde birden çok rigid-body mode,
- `Phi^T M Phi~=I` mass orthogonality ve her eigenpair için boyutsuz relative
  residual,
- original K/M matrislerinin DSYGV çağrısından sonra değişmemesi,
- tamamen constrained `0x0` sistemde LAPACK çağrılmadan temiz tanı.

Mode shape işareti fiziksel olarak arbitrary olduğundan `phi` ile `-phi`
eşdeğer kabul edilir. Repeated eigenspace içinde seçilen baz da unique değildir;
test exact sütun karşılaştırması yapmaz. Nonsymmetric/non-finite K/M, boyut
uyuşmazlığı, non-SPD M ve anlamlı negatif eigenvalue durumları ayrı beklenen-hata
CTest süreçleridir. Non-SPD, fully constrained ve significant-negative yolları
yalnız exit code ile değil, beklenen TMS26 tanı metniyle de doğrulanır. Singular
symmetric positive-semidefinite K ise free-free rigid mode taşıyabildiği için
geçerlidir. Modal baz için residual ve `Phi^T M Phi~=I` yanında
`Phi^T K Phi~=diag(lambda)` bağıntısı da bağımsız test hesabıyla sınanır.

Bu kapsam lineer, sönümsüz ve frozen-property çözümü doğrular. Toleranslar ve
analitik model tablosu
[`../validation/modal_eigen_validation.md`](../validation/modal_eigen_validation.md),
matematiksel tanımlar ise
[`../mathematics/generalized_modal_eigenproblem.md`](../mathematics/generalized_modal_eigenproblem.md)
belgelerindedir.

V0.6.0 harmonic doğrulaması, mevcut modal regresyonları değiştirmeden ayrı bir
direct/full-order complex response kanıt zinciri ekler. Test katmanları:

- element `K''` ve `C` lokal katkıları ile passive-property doğrulaması,
- üç düğümlü full/reduced `K'`, `K''`, `C`, `M` exact assembly regresyonu,
- dynamic stiffness için reel/sanal bileşenler ve complex-symmetry invariantı,
- LP64 LAPACK ZSYSVX facade için bilinen çözüm, multiple RHS, input
  immutability, exact singular ve working-precision ill-conditioned status,
- fixed 1-DOF viscous-only, loss-only ve combined damping analitik cevapları,
- fixed-hub ve free-free iki ataletli TVD harmonic referansları,
- excitation scatter-add, constraint tanıları ve explicit frequency-array
  sözleşmesi,
- physical complex recovery, magnitude/phase, velocity/acceleration, FRF,
  relative angle, element torque ve dissipated power/energy.

Exact singular nokta beklenen program hatası değildir; `SINGULAR` analysis
status ile normal test akışında doğrulanır. `SOLVED_ILL_CONDITIONED` durumda
hesaplanan cevap ile RCOND/FERR/BERR ve backend-independent residual korunur.
Exact `RCOND`, `FERR` veya `BERR` değerleri platformlar arasında zorlanmaz.
Invalid input vakaları ise ayrı CTest selector süreçlerinde ve gerektiğinde
beklenen diagnostic regex'iyle sınanır.

Harmonic test yardımcıları kompleks değerlerin reel ve sanal bileşenlerini
IEEE finite kontrolünden geçirir. Relative residual, matris/vector normları ve
passivity toleransları problem ölçeğine duyarlı olmalıdır. Ayrıntılı vaka
listesi ve tolerans ilkeleri
[`../validation/harmonic_response_validation.md`](../validation/harmonic_response_validation.md),
denklem sözleşmesi
[`../mathematics/harmonic_torsional_response.md`](../mathematics/harmonic_torsional_response.md)
belgelerindedir.

V0.7.0 material-aware harmonic doğrulaması üç yeni CTest programı ve ayrı
expected-failure süreçleriyle aşağıdaki kapıları kapsar:

- tabulated provider exact point, one-ULP numerical match,
  `LINEAR_FREQUENCY` ve `LINEAR_LOG_FREQUENCY` analitik interpolation,
- iki policy için alt/üst no-extrapolation, measured-isotherm temperature ve
  strictly increasing/passive dataset validation,
- negatif G'' için clipping olmadan red ve input/private-copy immutability,
- bonded-annular `G'/G'' -> K'/K''` mapping, loss-factor eşitliği ve SHEAR-only
  direct torsional contract,
- dynamic override/no-double-counting ile V0.6 frozen nominal yolun yan yana
  regresyonu; frequency-independent viscous c'nin korunması,
- constant+dynamic mixed zincir ve iki farklı provider/element binding'i için
  bağımsız complex 2x2 analitik çözüm,
- tüm provider domain'inin ZSYSVX öncesinde prevalidation'ı,
- solved ve exact-singular frequency noktalarında dataset/material state trace.

Interpolation toleransları fiziksel deney uncertainty'si değildir. Exact
eşleme machine epsilon ölçeğinde, harmonic analitik karşılaştırmalar ise
double-precision dense solver roundoff'una uygun scale-aware toleransla yapılır.
Ayrıntılar
[`../validation/dynamic_material_provider_validation.md`](../validation/dynamic_material_provider_validation.md)
belgesindedir.

V0.8.0 thermorheological runtime doğrulaması, V0.7 provider abstraction'ını
değiştirmeden temperature-shift katmanını aşağıdaki ayrı kapılarla sınar:

- WLF için `T_ref` identity, hand-calculated `log10(a_T)`, physical sıcaklık
  yönü, positive `C1/C2`, explicit domain ve WLF-pole rejection,
- Arrhenius için hand-calculated `Ea/R(1/T-1/T_ref)`, `T_ref` identity,
  physical yön ve invalid activation-energy/domain rejection,
- tabulated shift için exact/reference noktalar, `log10(a_T)` linear
  interpolation, machine-equivalent temperature, strictly increasing T,
  NaN/Inf/duplicate rejection ve no extrapolation,
- abstract provider sınırında operating/reference T, model kind, bracket ve
  `log10(a_T)`/`a_T` tutarlılığı; bozuk test-double evaluation rejection,
- `log10(f_r)=log10(f)+log10(a_T)` analytical reduced-frequency hesabı,
  log-space numerical safety ve ayrı shift-temperature/master-frequency
  domain kontrolleri,
- master curve ile shift provider reference-temperature consistency ve
  `T=T_ref` durumunda aynı V0.7 `G'/G''` sonucu,
- physical frequency/operating temperature ile reduced lookup
  frequency/reference temperature'ın trace içinde ayrılması,
- bütün geçerli sıcaklıklarda `G'>0`, `G''>=0`, `K'>0`, `K''>=0` passivity,
- mevcut dynamic binding ve `analyze_material_aware_harmonic_response()` ile
  bağımsız fixed-hub 1-DOF analytical chain,
- aynı modelde farklı WLF/Arrhenius provider'ları ve singular response
  noktasında material trace retention.

Tabulated shift ordinatlarına monotonicity zorlanmaz; yalnız temperature axis
strictly increasing olmalıdır. Testler temperature/master-curve extrapolation,
endpoint clamp veya nearest-value fallback'i kabul etmez. Büyük shift
senaryoları invalid point'in ara overflow üretmeden kontrollü reddedildiğini
doğrular. Ayrıntılar
[`../validation/thermorheological_runtime_validation.md`](../validation/thermorheological_runtime_validation.md)
belgesindedir.

V0.8.1 experimental identification doğrulaması dokuz ayrı CTest ailesiyle
aşağıdaki kapıları ekler:

- family-level authoritative common state ve explicit point quality,
- quality gap üzerinden interpolation yapılmayan contiguous channel segments,
- measured-domain feasible shift ve exact piecewise-linear L2 integral,
- deterministic coarse scan ve yalnız valid interior bracket ile Brent,
- joint G'/G'' shift, storage-only fallback, overlap ve curvature evidence,
- explicit reference zero shift, colder/hotter adjacent chain ve broken-chain,
- provenance-preserving master cloud ile low/high/both/no-extension stitching,
- runtime-valid adjacent interval union, internal-hole rejection, edge-domain
  shrink ve cross-isotherm coverage bridge,
- VGP/Cole-Cole analytical points, exact TRS, deterministic non-TRS ve weak
  identifiability,
- independent generalized-Maxwell equations ile curved exact-TRS recovery,
- existing V0.8.0 provider nesneleri üzerinden exact ve intermediate/off-
  reference runtime round-trip.

Exact synthetic shift toleransı `1.5e-6` boyutsuz düzeydedir ve
`8*sqrt(epsilon(dp))` numerical stopping ölçeğini kapsar; measurement
uncertainty değildir. Exact integral/identity testleri `1e-12`–`1e-14`
ölçeğindedir. Universal TRS acceptance threshold kullanılmaz. Bütün yeni
testler `-fcheck=all` build'inde de çalıştırılır. Ayrıntılar
[`../validation/V0.8.1_tts_validation.md`](../validation/V0.8.1_tts_validation.md)
belgesindedir.

V0.8.2 parametric shift-law doğrulaması dört normal test ailesi ve iki
expected-failure runtime-domain kaydı ekler:

- analytical Arrhenius beta/Ea_app recovery ve iki measured reference altında
  apparent activation-energy invariance,
- profiled WLF C1/C2 recovery, exact reference reparameterization ve physical
  relative-prediction invariance,
- large-C2 linear-limit fixture'ında düşük residual ile parameter
  identifiability'nin ayrılması,
- pole-boundary/no-interior, insufficient, duplicate ve nonfinite input
  status yolları,
- V0.8.1 pair results extraction, cumulative empirical shifts'in fit
  observation olmadığının ve input result'un immutable kaldığının kontrolü,
- en az beş temperature'da Leave-One-Temperature-Out predictive diagnostics,
- existing V0.8.0 Arrhenius/WLF providers ile measured-domain parametric
  runtime round-trip ve domain dışı sorgu rejection.

Arrhenius/WLF exact recovery toleransları numerical stopping scale'ine göre
ayarlanır; universal rubber acceptance threshold'u değildir. Pair
curvature/overlap statistical variance olarak kullanılmaz. Ayrıntılar
[`../validation/V0.8.2_shift_law_validation.md`](../validation/V0.8.2_shift_law_validation.md)
belgesindedir.

V0.8.3 repeatability doğrulaması beş CTest ailesi ekler:

- `tms26.tts_sample_statistics`: mean, `n-1` sample SD, SE, median, MAD,
  scaled MAD, single/empty/nonfinite availability semantics,
- `tms26.tts_bootstrap`: portable RNG sequence, seed reproducibility,
  replacement draw, Type-7 quantile ve valid/unavailable draw accounting,
- `tms26.tts_repeatability`: complete V0.8.1 campaign, physical-temperature
  mapping, reference normalization, immutability ve incompatible-state yolları,
- `tms26.tts_parametric_repeatability`: valid Arrhenius ile identifiable WLF
  cohort statistics ve placeholder içermeyen partial availability,
- `tms26.tts_repeatability_bootstrap`: independent campaign ile same-specimen
  rerun population ayrımı ve deterministic interval reproducibility.

Statistical sample unit complete independent campaign'dir. Frequency point,
isotherm ve adjacent pair bağımsız replicate sayılmaz. Aynı specimen rerun'ı
default descriptive evidence'a alınabilir, fakat independent bootstrap
population'ına alınmaz. Same-specimen-only study bağımsız confidence interval
üretemez.

Tek cluster draw planı bütün adjacent shift, common-reference shift,
Arrhenius ve WLF niceliklerinde paylaşılır. `A=[1,2,3]`, `B=10A` campaign
fixture'ı her draw'da `mean(B)=10 mean(A)` koşulunu doğrular; böylece
quantity-level ayrı resampling regresyonu yakalanır. Derived fit availability
önceden population filtrelemez: whole-campaign draw'dan sonra usability
uygulanır ve ikiden az usable değer kalan draw explicit unavailable olur.

Common-reference noktasında normalization gereği mean/SD ve bootstrap interval
sıfır olabilir. Test, bunu `is_reference_anchor=true` ve
`uncertainty_informative=false` ile yapısal ankraj olarak doğrular; zero
measurement uncertainty yorumu yapmaz. Bootstrap confidence interval da
engineering acceptance tolerance değildir. Ayrıntılar
[`../validation/V0.8.3_repeatability_validation.md`](../validation/V0.8.3_repeatability_validation.md)
belgesindedir.

## Integration test

Entegrasyon testleri, birden fazla modülün birlikte kullanımını doğrular.
`test_generalized_torsional_system`, node/element koleksiyon yönetimini ve
mevcut iki ataletli TVD'nin genel topolojiye kayıpsız dönüşümünü sınar.
Benchmark 004 için `J_h`, `J_r`, K' ve K'' aktarılır; K'' ayrı kayıp rijitliği
alanında korunurken boyutsal olarak farklı viskoz `c` değeri sıfır kalır.

Genel sistem testleri ayrıca sabitlenmemiş düğümlerin aktif DOF sayısını,
fixed-hub dönüşümünde göbek kısıtını, yinelenen kimlikleri ve tanımsız eleman
uçlarını kapsar. Bu testler birim testlerinin yerini almaz. Global assembly
`test_matrix_assembly` içinde ayrı doğrulanır; eigen çözümü uygulanmaz.

V0.4.0 için `test_constraint_foundation`, tam M/K assembly, constraint manager,
aktif Equation ID haritası, storage-bağımsız reduction ve result recovery
katmanlarını birlikte doğrulayan entegrasyon testidir. Önceki V0.3 assembly ve
validation testleri tarihsel regresyon koruması olarak çalışmayı sürdürür.

V0.5.0 için `test_modal_analysis`, reduced-system builder, generalized solver
facade, DSYGV reference backend, modal validation/result ve mevcut
`recover_mode_shape` katmanlarını uçtan uca birleştirir. LAPACK backend yalnız
sayısal K/M problemine bağımlıdır; test Geometry, Material, Physical DOF ve
constraint sorumluluklarının backend'e sızmadığını da katman kullanımıyla
korur.

V0.6.0 için complex linear solver testleri facade ile ZSYSVX reference
backend'ini, harmonic integration testleri ise full dynamic assembly'den
status-aware frequency sweep ve physical response recovery'ye kadar olan yolu
birleştirir. Low-level solver birden çok RHS'yi sınarken public harmonic analiz
tek load case'i koruyabilir. Singular bir frekans noktası bütün sweep'i
sonlandırmamalıdır.

V0.8.0 için thermorheological provider entegrasyonu constitutive layer ile
mevcut V0.7 binding/harmonic orchestration arasındaki public sınırı sınar.
Entegrasyon testi yeni bir harmonic API çağırmaz; aşağıdaki zinciri production
provider, binding ve
`analyze_material_aware_harmonic_response()` üzerinden yürütür. Multiple shift
provider ve singular-trace senaryoları solver reuse ile state trace'in birlikte
korunduğunu kanıtlar.

```text
physical f,T -> a_T -> f_r -> G'/G'' -> K'/K'' -> Z -> theta
```

## Benchmark test

Benchmark testleri, temsilî TVD modellerinde yürütme süresi ve bellek kullanım
eğilimlerinin yanı sıra analitik referans sonuçların korunmasını sağlar.
Başlangıç senaryoları `benchmarks/` altında tekrarlanabilir girdiler ve beklenen
sonuçlar olarak tutulur. Donanım farklılıklarına bağlı performans değerleri
henüz ana CI geçiş koşulu değildir; analitik referanslar ilgili CTest fizik
doğrulamalarında sınanır.

Benchmark 004, hub-normalized analitik iki-atalet modlarını DSYGV'nin
mass-normalized ve sign-arbitrary modlarıyla eigenvalue, residual, mass inner
product ve physical recovery üzerinden eşleştirir. Böylece normalize bileşen
değerlerinin farklı olması fizik hatası olarak değerlendirilmez.

V0.6 harmonic regression referansları aynı iki-atalet fiziğinin fixed-hub ve
free-free finite-frequency cevaplarını kullanır. Harmonic benchmark verisi;
peak complex torque, açık frequency array, frozen K'/K''/c/M ve response/phase/
relative-angle/energy sonuçlarını modal oracle'dan kavramsal olarak ayrı
tutmalıdır.

## CI validation

GitHub Actions, `main` ve `develop` dallarına yapılan push işlemlerinde ve tüm
pull request'lerde aşağıdaki işlemleri yapar:

1. Kaynak kodunu alır.
2. macOS üzerinde Homebrew GNU Fortran ile keg-only LP64 LAPACK; Windows
   üzerinde MinGW64 GNU Fortran ile LP64 OpenBLAS/LAPACK araç zincirini kurar.
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
