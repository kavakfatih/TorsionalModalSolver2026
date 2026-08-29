# TMS26 Hesap Motoru Testleri

Bu dizindeki Fortran test programları CMake ile derlenir ve CTest üzerinden
çalıştırılır. Yeni veya değiştirilen her davranış uygun test kategorisiyle
kapsanmalıdır.

## Unit test

Tek modül, tür veya yordam izole olarak doğrulanır. `test_kinds`,
`test_constants`, `test_units`, `test_geometry`, `test_material`,
`test_dynamic_modulus`, `test_torsional_node` ve `test_torsional_element` bu
kategoridedir. `test_dynamic_modulus`, G' ile G''
değerlerinin SI biriminde saklanmasını, çoklu çalışma noktası altyapısını ve
`tan(delta) = G'' / G'` hesabını doğrular. `test_geometry`, geometri veri
alanlarının yanında annüler kesit polar alan momentini ve tam bağlı kauçuk
burç `Cθ` faktörünü ayrı bağımsız analitik değerlerle sınar.

## Physics validation

Fiziksel veya matematiksel yordam, bağımsız analitik sonuç ve açık bir hata
sınırıyla doğrulanır. `test_inertia`, `test_hub_inertia`,
`test_torsional_stiffness`, `test_dynamic_torsional_stiffness`,
`test_frequency_solver` ve `test_torsional_system` bu kategoridedir.
Genel fizik doğrulamalarında bağıl hata yüzde `0,1`'den küçük olmalıdır;
dinamik burulma rijitliği testi `1e-10` bağıl hata sınırı kullanır. Bu test
üretim `tms_dynamic_torsional_stiffness` modülünü doğrudan çağırarak Cθ,
K', K'', frekans, sıcaklık ve
`G''/G' = K''/K' = tan(delta)` eşitliğini birlikte doğrular. Ayrıca K'
bileşeninin statik solver sonucuyla aynı kaldığını, `G'' = 0` için kayıp
bileşenlerinin sıfırlandığını ve `ri`, `ro` değerine yaklaştıkça ideal
model rijitliğinin kuvvetli biçimde arttığını sınar. Dinamik test ayrıca
`K' ∝ G'`, `K' ∝ L` ve `K'' ∝ L` ölçekleme regresyonlarını üretim
fonksiyonu sonuçlarını karşılaştırarak doğrular. Statik test, G ve L iki katına
çıkarıldığında K değerinin de iki katına çıktığını doğrular.

Negatif veya sıfır iç yarıçap, sırasız/eşit yarıçaplar, pozitif olmayan
eksenel genişlik, pozitif olmayan depolama modülü ve negatif kayıp modülü
ayrı CTest vakalarıdır.
Bu vakalarda üretim yordamının `error stop` ile sonlanması beklenir ve
`WILL_FAIL` özelliği beklenen reddi test başarısına dönüştürür.

`test_hub_inertia`, homojen annüler göbeğin hacim, kütle ve polar kütle
ataletini bağımsız analitik sabitlerle doğrular. `test_torsional_system`,
builder entegrasyonunu, fixed-hub eşdeğerliğini, serbest-serbest sıfır ve
elastik modları, normalize mod şekillerini, K'/atalet ölçeklemelerini ve büyük
göbek ataleti limitini sınar. K'' ve kayıp faktörü değişiminin sönümsüz modal
sonucu değiştirmediği de ayrı bir regresyonla korunur. Pozitif olmayan atalet,
K' ve yoğunluklar ile geçersiz göbek/halka geometrileri ayrı `WILL_FAIL`
regresyonlarıdır.

`test_torsional_node`, kimlik, polar atalet, başlangıç açısı ve sınır koşulu
alanlarını; `test_torsional_element`, uç düğüm kimlikleri ile ayrı K', K'' ve
viskoz c alanlarını doğrular. Pozitif olmayan K', negatif K''/c veya sonlu
olmayan büyüklükler ve self-connection ayrı `WILL_FAIL` vakalarıdır.

`test_local_stiffness_matrix`, `k = 100 N·m/rad` için üretim yordamının
`[[100,-100],[-100,100]]` lokal katkısını verdiğini doğrular. Simetri, sıfır
satır toplamı, bağıl dönmede pozitif enerji ve ortak dönmede sıfır enerji ayrı
assertion'larla korunur. Negatif rijitlik aynı üretim yordamını çağıran ayrı bir
`WILL_FAIL` regresyonudur.

`test_generalized_torsional_system`, private koleksiyonların public yönetim
yordamlarını, aktif DOF sayımını ve Benchmark 004 iki-ataletli sisteminin genel
iki-düğüm/bir-eleman gösterimini doğrular. K'' kayıp rijitliğinin ayrı element
alanına aktarıldığı, viskoz c'nin sıfır kaldığı, serbest-serbest gösterimin iki
ve fixed-hub gösteriminin bir aktif DOF taşıdığı sınanır. Boş sistem, yinelenen
kimlik ve tanımsız eleman ucu hata regresyonlarıdır.

`test_matrix_assembly`, fiziksel node ID ile equation ID ayrımını, tek eleman
ve üç düğümlü zincir için global torsional K matrisini, K simetrisini, serbest
sistem sıfır satır toplamını ve rijit-cisim null modunu doğrular. Global M
matrisinde düğüm polar ataletlerinin `kg·m²` birimiyle diagonal saklandığını ve
köşegen dışı katsayıların sıfır kaldığını sınar. Kısıtlı uç için indirgenmiş
`K=[k]`, `M=[J_free]` ve tamamen kısıtlı sistem için 0x0 matrisler korunur.
Geçersiz boyut, eksik node lookup, aralık dışı denklem, negatif atalet ve
uyumsuz DOF haritası ayrı `WILL_FAIL` regresyonlarıdır.

`test_torsional_validation`, V0.3.0 foundation katmanlarını değiştirmeden
uçtan uca analitik doğrulama sağlar. Tek eleman işaret konvansiyonu, serbest
rijit-cisim modu, üç düğümlü global katkı toplamı, iki constraint durumunda
DOF mapping sürekliliği, Frobenius simetri normu ve
`U=1/2 theta^T K theta>=0` enerji koşulu birlikte sınanır. İki ataletli TVD
referansı ayrıca assembled M/K modal residual, bilinen analitik modun Rayleigh
quotient değeri ve kütle ortogonalliği üzerinden mevcut analitik solver ile
çapraz doğrulanır. Bu test genel eigen çözümü yapmaz ve deneysel model
validation yerine analytical code verification sunar.

`test_numerical_hardening`, V0.3.2 girdi ve sayı aralığı sağlamlaştırmasını
doğrular. Nominal test akışı; mutlak `1e-12 m` ve bağıl `1e-9` toleranslı TVD
geometri ara yüzlerini, uç atalet/rijitlik ölçeklerinde sonlu eşdeğer atalet ile
doğal frekansı ve fixed-DOF için indirgenmiş `K=[k]`, `M=[J_free]` matrislerini
sınar.

Atalet, yoğunluk, yarıçap, uzunluk, rijitlik, modül, frekans ve sıcaklık
alanlarındaki IEEE `NaN`, pozitif sonsuz ve negatif sonsuz değerleri ayrı
`tms26.numerical_hardening.rejects_*` CTest süreçleridir. Bu süreçler yalnız
üretim doğrulayıcısının kontrollü `error stop` sonucunu başarı kabul eder.
Assertion yardımcıları da sonlu olmayan actual, expected veya tolerance
değerlerinin testi yanlış biçimde geçirmesine izin vermez. Ayrıntılar
[`../../docs/validation/numerical_robustness_validation.md`](../../docs/validation/numerical_robustness_validation.md)
belgesindedir.

`test_constraint_foundation`, V0.4.0 tam sistemden constraint-aware aktif
sisteme geçişini doğrular. Test, iki düğümlü tek elemanda ilk düğüm fixed iken
`Kr=[k]` ve `Mr=[J_free]` sonucunu; üç düğümlü zincirde iki aktif DOF'u ve tüm
DOF'lar fixed olduğunda geçerli `0x0` Kr/Mr matrislerini sınar. Physical DOF,
tam Equation ID ve Active Equation ID değerleri ayrı kontrol edilir.

Node ekleme sırası değiştirilmiş eşdeğer model, katsayıları ortak Physical DOF
sırasına aldıktan sonra aynı indirgenmiş fizik ve fiziksel kimlik tabanlı
recovery sonucunu vermelidir. Recovery testi genel durum için `q=Pq_r+q_p`,
homojen modal vektör için `phi=Pphi_r` eşlemesini doğrular; prescribed değer
Kr/Mr katsayılarına eklenmez ve eigen çözümü yapılmaz. Assertion yardımcıları
ile recovery girdileri IEEE `NaN` değerini kabul etmez.

Olmayan node, aynı Physical DOF için yinelenen constraint ve geçersiz DOF türü
ayrı `tms26.constraint_foundation.rejects_*` CTest süreçleridir. Bu süreçlerde
yalnız üretim doğrulayıcısının kontrollü `error stop` sonucu başarı kabul
edilir. Farklı node sırasındaki başka bir sisteme ait full DOF haritası da
yanlış K/M satır eliminasyonunu önlemek için reddedilir. Ayrıntılı mimari ve
matematik sözleşmeleri
[`../../docs/architecture/V0.4_constraint_foundation.md`](../../docs/architecture/V0.4_constraint_foundation.md)
ve
[`../../docs/mathematics/constraint_reduction.md`](../../docs/mathematics/constraint_reduction.md)
belgelerindedir.

`test_generalized_eigen_solver`, backend-neutral problem ve solution
sözleşmeleriyle LP64 LAPACK DSYGV reference backend'ini doğrudan doğrular.
K/M karelik, ortak boyut, sonluluk ve simetri önkoşulları; M positive definite
zorunluluğu; singular positive-semidefinite K kabulü; input immutability;
ascending eigenpair eşleşmesi ve anlamlı LAPACK diagnostic'leri bu test
grubundadır. Invalid selector vakaları ayrı
`tms26.generalized_eigen_solver.rejects_*` CTest süreçleri olarak çalışır.
Non-SPD M vakalarında CMake sarmalayıcısı hem nonzero çıkışı hem de beklenen
`M positive definite` TMS26 tanısını doğrular. Backend-neutral solution
sözleşmesi `1<=mode_count<=DOF_count` partial-spectrum aralığını korur;
DSYGV V0.5'te yine bütün spectrum'u çözer.

`test_modal_analysis`, V0.4 reduced system ile V0.5 modal katmanı arasındaki
uçtan uca bağlantıyı sınar. Fixed 1-DOF, free-free iki atalet, üç düğümlü
zincir, ilk düğümü constrained zincir, repeated eigenvalue, birden çok rigid
mode ve tamamen constrained sistem senaryolarını kapsar. Kontroller şunlardır:

- eigenvalue ve Hz frekansının analitik referansla uyuşması,
- `phi^T M phi=1` mass normalization,
- boyutsuz relative eigenpair residual,
- `Phi^T M Phi~=I` mass orthogonality,
- `Phi^T K Phi~=diag(lambda)` modal stiffness diagonalization,
- mode shape işaretinden bağımsız correlation,
- repeated eigenvalue için basis-independent eigenspace eşdeğerliği,
- `phi=Pphi_r` physical recovery ve constrained bileşenlerin sıfır kalması.

Anlamlı negatif eigenvalue ve aktif DOF bulunmaması temiz modal tanılarla
reddedilir; bu iki yolun tanı metni CMake expected-failure sarmalayıcısıyla
kilitlenir. Küçük roundoff kaynaklı negatif değerler ölçeğe duyarlı rigid-mode
toleransı içinde ele alınır. Requested mode count seçimi ve gelecekteki
partial-spectrum backend'lere açık `1<=m<=n` ortak tolerans sözleşmesi de
regresyon kapsamındadır. Bu testler lineer, sönümsüz ve frozen-property
çözümü doğrular; frequency-dependent elastomer iterasyonu veya damping çözümü
yapmaz. Ayrıntılar
[`../../docs/validation/modal_eigen_validation.md`](../../docs/validation/modal_eigen_validation.md)
belgesindedir.

V0.6 harmonic test grubu, modal testlerden ayrı direct/full-order complex
response yolunu doğrular. `test_complex_linear_solver`, backend-neutral
problem/solution contract'ı ile LP64 LAPACK ZSYSVX reference backend'ini
gerçekten çağırır. Kapsamı:

- complex symmetric fakat Hermitian olmayan bilinen çözüm,
- multiple RHS ve per-RHS `FERR/BERR/residual`,
- original A/B input immutability,
- exact singular `SINGULAR` status ve fabricated response bulunmaması,
- working-precision `SOLVED_ILL_CONDITIONED` status ile çözüm ve diagnostics'in
  korunması,
- karelik, boyut, transpose-symmetry ve kompleks IEEE finite önkoşullarıdır.

Harmonic physics/integration testleri ayrıca şu zinciri üretim API'leri
üzerinden sınar:

- lokal ve global `K''/C` assembly ile aynı-index constraint reduction,
- `Z=K'-omega^2 M+i(K''+omega C)` ve `Z^T=Z`, `Z^H/=Z` ayrımı,
- fixed 1-DOF viscous/loss/combined analitik cevapları,
- fixed-hub ve free-free iki ataletli TVD referansları,
- finite, pozitif ve strictly increasing explicit frequency array,
- complex nodal torque scatter-add ile unknown/unsupported/constrained target
  tanıları,
- status-aware sweep, reduced/physical response recovery ve input
  immutability,
- magnitude, phase, velocity, acceleration, rotational FRF, relative angle,
  dynamic element torque, transmitted-torque magnitude ve passive dissipated
  power/energy.

Singular bir frequency point normal analysis-state sonucudur ve CTest
expected-failure vakası değildir. Geçersiz frequency, excitation, negatif veya
nonfinite `K''/c` ise ayrı selector süreçlerinde reddedilir. Platformlar arası
testler exact `RCOND/FERR/BERR` değerini değil, status, sonluluk, işaret ve
backend-independent relative residual sözleşmesini kullanır. Ayrıntılar
[`../../docs/validation/harmonic_response_validation.md`](../../docs/validation/harmonic_response_validation.md)
belgesindedir.

V0.7 test grubu `test_tabulated_dynamic_modulus_provider`,
`test_dynamic_torsional_property_binding` ve
`test_material_aware_harmonic_analysis` programlarından oluşur. Exact ve
machine-equivalent frequency noktaları, iki interpolation policy, strict
no-extrapolation, measured-isotherm temperature, passive data validation,
private-copy immutability ve SHEAR-only binding ayrı doğrulanır.

Uçtan uca test; interpolate edilen G'/G'' değerlerini bonded-annular K'/K''
ile 1-DOF dynamic stiffness'e bağlar ve ZSYSVX cevabını bağımsız analitik
formülle karşılaştırır. Deliberately farklı nominal element K'/K'' değerleri
no-double-counting'i, aynı modelin V0.6 frozen analizi geriye uyumluluğu
kanıtlar. Mixed constant/dynamic ve iki-provider seri sistemleri bağımsız 2x2
complex inverse ile sınanır. Exact singular dynamic stiffness noktasında
harmonic response unavailable iken material trace'in korunduğu doğrulanır.
Duplicate/unknown binding, full-sweep domain eksikliği ve temperature mismatch
solver çağrısından önce ayrı expected-failure CTest süreçlerinde reddedilir.
Ayrıntılar
[`../../docs/validation/dynamic_material_provider_validation.md`](../../docs/validation/dynamic_material_provider_validation.md)
belgesindedir.

V0.8 thermorheological test grubu temperature-shift provider'larını,
thermorheological dynamic-modulus provider bileşimini ve mevcut
material-aware harmonic zincirini kapsar. Ana doğrulamalar:

- WLF `T_ref` identity, sıcaklık yönü, hand-calculated shift ve pole/domain
  hataları,
- Arrhenius identity, analytical `Ea/R` sonucu, sıcaklık yönü ve invalid
  input/domain hataları,
- tabulated `log10(a_T)` exact/midpoint interpolation, reference zero shift,
  machine-equivalent temperature, input doğrulama ve no extrapolation,
- common provider boundary'de model kind, sıcaklık/bracket ve
  `log10(a_T)`/`a_T` evaluation invariant'ları,
- master curve ile shift provider `T_ref` consistency,
- `physical_frequency_hz` ile `lookup_frequency_hz=f_r` ayrımı ve returned
  modulus'un physical `f,T` durumunu koruması,
- reduced-frequency domain'in önce log-space'te doğrulanması ve extreme shift
  girdilerinin ara overflow/underflow olmadan reddedilmesi,
- bütün geçerli sorgularda G'/G'' ve K'/K'' passivity,
- fixed-hub 1-DOF independent analytical response, aynı modelde multiple shift
  provider ve singular response'ta shift/material trace retention,
- thermorheological provider'ın mevcut
  `analyze_material_aware_harmonic_response()` API'siyle çözülmesi.

V0.7 unshifted provider testleri `log10(a_T)=0`, `a_T=1` ve
`lookup_frequency=physical_frequency` defaults'unu geriye uyumluluk olarak
korur. V0.5 modal ve V0.6 frozen harmonic test grupları değişmeden çalışır.
Ayrıntılar
[`../../docs/validation/thermorheological_runtime_validation.md`](../../docs/validation/thermorheological_runtime_validation.md)
belgesindedir.

V0.8.1 test grubu `tts_data_model`, `tts_scalar_minimizer`, `tts_pair_shift`,
`tts_shift_chain`, `tts_master_curve`, `tts_diagnostics`, `tts_identification`
`tts_generalized_maxwell` ve `tts_runtime_roundtrip` CTest ailelerinden oluşur.
Ortak synthetic helper
yalnız test target'larına derlenir; production material modeli değildir.

Pair testleri fixed grid sampling yerine exact linear residual integralini,
irregular frequency spacing ve farklı point density altında doğrular. Loss
quality gap ayrı segmentlere bölünür; `VALID,G''=0` runtime için korunurken
log-objective'e alınmaz. G'' support yoksa storage-only, iki channel varsa
joint shift ve ayrı diagnostic shifts kontrol edilir. No-support ve
no-interior-minimum clean status ile döner.

Shift-chain ve master testleri explicit reference, hot/cold cumulative shifts,
broken chain, provenance, low/high/both/no extension, duplicate priority,
strict ordering, boundary diagnostics, runtime interval-union coverage,
internal/edge/cross-isotherm gap davranışı ve checked exponentiation
kapsamındadır. Top-level test exact
TRS ile deterministic non-TRS residual/discrepancy ayrımını ve plateau
curvature evidence'ını sınar. Generalized-Maxwell testi bağımsız eğrisel
viscoelastic denklemlerle known shifts'i doğrular. Round-trip testi V0.8.1
output'larını mevcut V0.8.0 providers ile exact point yanında intermediate ve
off-reference physical `(f,T)` sorgusuna taşır. Ayrıntılar
[`../../docs/validation/V0.8.1_tts_validation.md`](../../docs/validation/V0.8.1_tts_validation.md)
belgesindedir.

## Benchmark regression

Birden fazla fizik adımını temsil eden sabit referans modelin sonuçları zaman
içinde korunur. `benchmarks/001_simple_annular_tvd/` girdileri; kütle, polar
atalet, rijitlik ve doğal frekans testlerinin ortak regresyon temelidir. Model
veya kabul edilen formül değişirse benchmark girdileri, beklenen sonuçlar ve
ilgili testler aynı commit içinde güncellenir.

`benchmarks/002_dynamic_elastomer/` dinamik modül veri noktasını,
`benchmarks/003_dynamic_torsional_stiffness/` ise bu noktanın annüler geometri
üzerindeki kompleks rijitlik sonucunu tanımlar.
`benchmarks/004_two_inertia_tvd/`, fixed-hub ve serbest-serbest iki ataletli
sistemin analitik frekansları ile hub-normalized mod şekillerini, V0.2.3 genel
topoloji eşlemesini ve V0.5 DSYGV mass-normalized/sign-arbitrary modal sonucunun
aynı fiziksel eigenspace'i verdiğini tanımlar.

`benchmarks/005_tts_identification/`, exact horizontal-collapse ve farklı
storage/loss shifts taşıyan deterministic non-TRS family'leri tanımlar.
Pair/chain, experimental cloud, TRS evidence, stitching ve V0.8.0 provider
round-trip testleri aynı known-truth girdileri kullanır.

## Test ekleme

Yeni Fortran testi `test_<konu>.f90` biçiminde adlandırılır ve
`engine/CMakeLists.txt` içindeki `tms26_add_fortran_test` yordamıyla CTest'e
kaydedilir. Tüm testler aşağıdaki komutla çalıştırılır:

```sh
ctest --test-dir build --output-on-failure
```
