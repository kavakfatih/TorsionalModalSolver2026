# Malzeme Modelleri

Bu dizin, TMS26 içinde kullanılan constitutive malzeme verilerinin fiziksel
anlamını, canonical SI sözleşmesini, deney metadata'sını ve solver'a aktarım
sınırlarını açıklar.

## Belgeler

- [`tabulated_dynamic_elastomer_material.md`](tabulated_dynamic_elastomer_material.md):
  tek measured isotherm üzerinde tabulated dinamik shear-modulus provider'ı,
  operating-state metadata'sı ve DMA → TVD aktarım sınırları
- [`thermorheological_dynamic_elastomer.md`](thermorheological_dynamic_elastomer.md):
  validated reference master curve, temperature shift, physical/reduced
  material-state ayrımı, horizontal-only TTS ve V0.8.1 TRS sınırı
- [`tts_experimental_data_model.md`](tts_experimental_data_model.md):
  family-level common state, isotherm/specimen provenance, explicit point
  quality, contiguous segment ve zero-loss semantics
- [`parametric_temperature_shift_laws.md`](parametric_temperature_shift_laws.md):
  empirical/parametric ayrımı, `Ea_app`, WLF identifiability, material-state
  applicability ve measured-domain runtime sınırı
- [`experimental_repeatability_and_uncertainty.md`](experimental_repeatability_and_uncertainty.md):
  intralaboratory repeatability, campaign independence/provenance, standards
  context'i ve DMA evidence ile bonded-TVD validation sınırı

Malzeme dosyası içe aktarma, vendor formatı yorumlama ve birim tahmini core
material modelinin parçası değildir. Core API frekansı Hz, modülleri Pa ve
sıcaklığı K cinsinden bekler.
