# Thermorheological Dynamic Elastomer

## Fiziksel amaç

V0.8.0, daha önce doğrulanmış bir dinamik elastomer reference master curve'ünü
operating temperature'a yatay olarak taşır. Constitutive model:

\[
G^*(f,T)=G'(f,T)+iG''(f,T),
\qquad G'>0,\quad G''\geq0
\]

olup canonical birimler `G'`, `G''` için Pa; physical ve reduced frequency
için Hz; sıcaklık için K'dir. Horizontal shift yeni bir modulus değeri
uydurmaz; reference curve'ü `f_r=a_T(T)f` noktasında değerlendirir.

## Validated reference-state sözleşmesi

Thermorheological provider iki doğrulanmış bileşimi taşır:

1. `T_ref` sıcaklığında bir V0.7 tabulated `G'(f_r),G''(f_r)` master curve,
2. açık `T_min..T_max` domain'li bir temperature-shift provider.

Master-curve metadata sıcaklığı ile shift modelinin reference temperature'ı
machine-equivalent olmalıdır. `dynamic_material_metadata_t` içindeki
`dataset_temperature_k`, thermorheological provider'da operating temperature
değil reference/master-curve temperature anlamına gelir.

Temperature shifting yalnız frequency/temperature koordinatını değiştirir.
Aşağıdaki material state alanları reference dataset ile aynı kalır:

- material ve compound kimliği,
- specimen ve source,
- dynamic strain amplitude,
- static prestrain,
- deformation mode,
- conditioning ve material state.

V0.8.0 amplitude normalization, prestrain normalization veya Payne correction
yapmaz. Doğrudan torsional binding SHEAR verisi kullanır; tensile-to-shear
conversion uygulanmaz.

## Physical ve lookup koordinatları

Provider sorgusu physical excitation frequency `f` ve externally prescribed
operating temperature `T` ile yapılır. Shift modelinden
`s=log10(a_T)` alınır ve reduced coordinate:

\[
\log_{10}(f_r)=\log_{10}(f)+s
\]

olarak hesaplanır. Master curve `f_r,T_ref` ile sorgulanır; dönen material
state ise physical `f,T` taşır. Reduced coordinate yalnız evaluation/trace
bağlamındadır.

Örneğin `f=100 Hz` ve `a_T=0.1` için master-curve lookup `f_r=10 Hz`
çevresindeki bracket'ta yapılır. Returned modulus yine `100 Hz` physical
frequency ve gerçek operating temperature'ı temsil eder.

## Passivity ve horizontal-only kabulü

Validated master curve `G'>0`, `G''>=0` ise yatay kaydırma bu işaretleri
korur. Pozitif bonded-annular geometri faktörü ile:

\[
K'=C_\theta G'>0,\qquad K''=C_\theta G''\geq0
\]

olur. V0.8.0'da vertical factor `b_T=1`'dir. Density correction, thermal
expansion, temperature-dependent inertia/density ve automatic vertical
alignment uygulanmaz.

## Deney verisi ve TRS sorumluluğu

Runtime provider, girdinin önceden thermorheologically simple olduğunun ve
master curve'ün kabul edilebilir doğrulukta oluşturulduğunun doğrulandığını
varsayar. V0.8.0 aşağıdakileri yapmaz:

- measured isotherm shifting ve empirical `a_T` identification,
- master-curve construction veya global shift optimization,
- overlap-error, Van Gurp–Palmen veya Cole–Cole diagnostics,
- WLF/Arrhenius model seçimi ve TRS acceptance,
- noisy DMA robustness değerlendirmesi.

Bu çalışmalar V0.8.1 Master-Curve Identification & TRS Validation kapsamıdır.
Specimen DMA durumunun gerçek TVD compound, cure, strain amplitude, prestrain
ve conditioning koşullarına aktarılabilirliği ayrıca test correlation ile
kanıtlanmalıdır.

## Thermal coupling sınırı

Operating temperature dışarıdan verilen, uniform ve sweep boyunca prescribed
bir büyüklüktür. V0.6 dissipated energy sonucu temperature rise hesabına
bağlanmaz; `G'' -> heat -> T -> G''` self-heating feedback loop'u yoktur.

## Commercial solver ilkeleriyle karşılaştırma

Karşılaştırma yalnız kamuya açık model ilkeleri düzeyindedir; proprietary
internal implementation eşdeğerliği iddia edilmez.

- [ANSYS Material Reference — Viscoelasticity](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_mat/evis.html),
  thermorheologically-simple frequency–temperature superposition ile
  experimental complex modulus/temperature shifting ve harmonic bağlamında
  WLF ve Tool–Narayanaswamy türü shift modellerini belgeler. ANSYS shift
  notation'ı TMS26 `a_T=tau(T)/tau(T_ref)` convention'ıyla aynı isim veya
  işarete sahip
  varsayılmaz; veri aktarımında dönüşüm açıkça yapılmalıdır.
- [Hexagon Marc 2024.2, Volume C — Program Input](https://documentation-be.hexagon.com/bundle/Marc_2024.2-Volume_C_Program_Input/raw/resource/enus/Marc_2024.2-Volume_C_Program_Input.pdf),
  frequency-domain storage/loss property ile thermo-rheologically simple
  `SHIFT FUNCTION` katmanını; WLF, Arrhenius, combined modelleri ve
  temperature-dependent tabulated logarithmic shift-factor seçeneklerini
  belgeler. TMS26 V0.8.0 yalnız tabulated, WLF ve Arrhenius doğruluk-first
  subset'ini kapsar; combined model kapsam dışıdır.

Bu karşılaştırma dosya formatı, parameter sign convention veya iki ticari
çözücüyle bire bir sayısal eşdeğerlik anlamına gelmez.
