# Tabulated Dynamic Elastomer Material

## Amaç

V0.7 provider'ı deneysel dinamik **shear** verisini doğrudan constitutive
girdi olarak kullanır:

\[
G^*(f)=G'(f)+iG''(f), \qquad G'>0,\quad G''\geq 0.
\]

`G'` depolanan elastik enerjiyle, `G''` bir çevrimde kaybedilen enerjiyle
ilişkilidir. Birimler Pa, frekans Hz ve mutlak sıcaklık K'dir. Loss factor
primary tablo kolonu değildir; `tan(delta)=G''/G'` olarak türetilir.

Provider yalnız dynamic shear modulus kabul eder. Tensile `E'`/`E''`, bulk
modulus veya complex Poisson ratio direct torsional girdiye dönüştürülmez;
özellikle `G=E/3` gibi örtük bir kabul yapılmaz.

## Tek operating-state sözleşmesi

Bir provider tek bir doğrulanmış deney durumunu temsil eder. Aynı frequency
curve içindeki bütün noktalar şu dataset-level durumda ortak olmalıdır:

- dataset temperature,
- dynamic shear strain amplitude,
- static shear prestrain,
- deformation mode,
- conditioning/precycling,
- specimen/material state.

Farklı genlik, prestrain, sıcaklık veya conditioning ölçümleri ayrı dataset ve
provider olarak tutulur. V0.7 bu eksenlerde interpolation yapmaz.

## Metadata ve availability

`dynamic_material_metadata_t` dataset/material kimliğini ve sıcaklığı zorunlu
tutar. Specimen kimliği, strain amplitude, prestrain, conditioning, material
state, test source/method, standard/reference ve notes alanları açık `has_*`
availability bayrakları taşır. Bilinmeyen veri `-1` veya `-999` gibi magic
değerlerle ifade edilmez.

Metadata ASTM D5992 ve ISO 4664-1 gibi dinamik kauçuk test çerçevelerine
izlenebilirlik sağlayabilir; bu standartların frekans aralıkları yazılım hard
limit'i değildir. Geçerli analiz aralığını gerçek measured dataset belirler.

## Veri kalitesi ve passivity

Provider en az iki nokta, sonlu ve strictly increasing `f>0`, sonlu `G'>0`,
sonlu `G''>=0` ve sonlu pozitif dataset sıcaklığı gerektirir. Duplicate veya
unordered frekanslar reddedilir. Negatif deneysel `G''` sıfıra clamp edilmez;
data-quality hatasıdır.

Input noktaları ve metadata provider private storage'ına bağımsız kopyalanır.
Caller dizisi veya getter kopyası daha sonra değiştirilse de authoritative
provider durumu değişmez. Legacy `dynamic_rubber_material_t` alanları ve
`frequency_points(:)` korunur; ADR 0003 uyarınca aralarında otomatik
eşzamanlama yoktur.

## DMA → TVD transfer sınırı

DMA specimen verisinin bonded annular TVD'ye aktarılması, specimen material
state ile gerçek component operating state'in yeterince benzer olduğu
varsayımını içerir. Strain amplitude, prestrain, temperature, conditioning,
cure state, orientation ve specimen preparation farkları sonucu etkileyebilir.
Provider bu eşdeğerliği otomatik doğrulamaz; component/test correlation gerekir.

CSV/Excel veya vendor parser, automatic column/unit guessing, Hz↔rad/s ve
Pa↔MPa tahmini V0.7 kapsamı dışındadır.

## ANSYS / Marc ile prensip karşılaştırması

Karşılaştırma yalnız doğrulanabilir public prensip düzeyindedir; proprietary
solver internals hakkında tahmin yapılmaz.

- ANSYS tarafındaki experimental complex shear modulus yaklaşımıyla ortak
  kavramlar `G'(f)`, `G''(f)`, frequency-domain viscoelastic property ve
  piecewise-linear frequency interpolation'dır. Temperature shifting/TTS ayrı
  bir model katmanıdır ve V0.7 provider'ına örtük olarak dahil edilmez.
- Marc-style table yaklaşımıyla ortak kavramlar storage/loss property'nin
  frequency ve açık operating-condition metadata'sına bağlı tabulated veri
  olması ile deterministic linear table interpolation'dır. Frequency,
  amplitude, pre-deformation ve temperature bağımlılık eksenleri kavramsal
  olarak ayrı tutulur; V0.7 yalnız tek amplitude/prestrain/isotherm curve çözer.

Bu benzerlikler dosya formatı veya iki ticari çözümleyiciyle bire bir sayısal
eşdeğerlik iddiası değildir. Karşılaştırma için aynı canonical SI dataset,
interpolation policy, boundary conditions ve harmonic sign convention açıkça
eşleştirilmelidir.
