# Karar 0013: Thermorheological Runtime Convention

- Durum: Kabul edildi
- Tarih: 2026-08-28
- Sürüm: V0.8.0

## Bağlam

V0.7, tek measured isotherm üzerinde tabulated `G'(f),G''(f)` sağlayan
`dynamic_modulus_provider_t` sınırını ve mevcut material-aware harmonic
analysis akışını kurdu. Operating temperature değiştiğinde validated reference
master curve'ün hangi yönde ve hangi koordinatta sorgulanacağı için tek,
izlenebilir ve ticari solver notation'larından bağımsız bir core convention
gerekiyordu.

## Karar

TMS26 canonical shift convention'ı:

\[
a_T(T)=\frac{\tau(T)}{\tau(T_{ref})},\qquad a_T(T_{ref})=1,
\]

\[
t_r=\frac{t}{a_T},\qquad f_r=a_Tf.
\]

Authoritative storage ve hesap büyüklüğü `s(T)=log10(a_T)`'dir. `a_T=10^s`
türetilmiş trace değeridir; reduced-frequency domain öncelikle
`log10(f_r)=log10(f)+s` ile değerlendirilir.

V0.8.0 üç temperature-shift provider destekler:

- tabulated `log10(a_T)` ve temperature-axis linear interpolation,
- WLF,
- Arrhenius.

Her provider explicit validated `T_min,T_max` domain'i taşır. Temperature
extrapolation yapılmaz. WLF domain'i paydanın strictly positive kaldığı pole
öncesi bölgeyle sınırlıdır. Shift değerleri için genel monotonicity zorlanmaz.

Master curve, `T_ref` metadata'sı taşıyan V0.7 tabulated dynamic modulus
provider'dır. Master-curve reference temperature ile shift-provider reference
temperature machine-equivalent olmalıdır. Thermorheological provider mevcut
`dynamic_modulus_provider_t` interface'ini uygular, mevcut dynamic torsional
binding'e bağlanır ve mevcut `analyze_material_aware_harmonic_response()`
tarafından çözülür. Yeni harmonic API oluşturulmaz.

Returned modulus physical `f,T` taşır; lookup trace reduced `f_r,T_ref`
taşır. V0.7 unshifted evaluation için `s=0`, `a_T=1` ve `f_lookup=f_physical`
default'ları korunur.

## Sonuçlar

- Physical operating state ile master-curve lookup koordinatı ayrılır.
- Büyük shift değerlerinde reduced-frequency domain kontrolü log-space'te
  güvenli yapılır.
- Temperature ve reduced-frequency domain'lerinin ikisi de zorunludur;
  extrapolation, endpoint clamp ve nearest-value selection yoktur.
- Horizontal shift validated master curve'ün `G'>0`, `G''>=0` passivity'sini
  ve mevcut positive `C_theta` mapping'inin `K'>0`, `K''>=0` sonucunu korur.
- V0.5 modal, V0.6 frozen harmonic ve V0.7 unshifted material semantics geriye
  uyumlu kalır.

## Açıkça kapsam dışı kararlar

V0.8.0 aşağıdakileri içermez:

- master-curve fitting veya TRS identification,
- vertical shift, density correction ya da automatic modulus alignment,
- self-heating veya thermal-mechanical feedback,
- Payne effect, amplitude veya prestrain interpolation,
- combined WLF–Arrhenius,
- transient viscoelasticity,
- nonlinear harmonic iteration.

Master-curve identification ve TRS validation V0.8.1 kapsamıdır.

## Commercial convention notu

[ANSYS Material Reference — Viscoelasticity](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/ans_mat/evis.html)
ve [Hexagon Marc 2024.2 Volume C](https://documentation-be.hexagon.com/bundle/Marc_2024.2-Volume_C_Program_Input/raw/resource/enus/Marc_2024.2-Volume_C_Program_Input.pdf)
temperature-shift/superposition modellerini belgeler. Bu ürünlerin parameter
adları veya işaret convention'ları TMS26 core değişkenlerine doğrudan
eşitlenmez. Dış veri dönüşümü açık formül karşılaştırmasıyla yapılır;
proprietary internal implementation eşdeğerliği varsayılmaz.
