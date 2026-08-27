# Malzeme Modelleri

Bu dizin, TMS26 içinde kullanılan constitutive malzeme verilerinin fiziksel
anlamını, canonical SI sözleşmesini, deney metadata'sını ve solver'a aktarım
sınırlarını açıklar.

## Belgeler

- [`tabulated_dynamic_elastomer_material.md`](tabulated_dynamic_elastomer_material.md):
  tek measured isotherm üzerinde tabulated dinamik shear-modulus provider'ı,
  operating-state metadata'sı ve DMA → TVD aktarım sınırları

Malzeme dosyası içe aktarma, vendor formatı yorumlama ve birim tahmini core
material modelinin parçası değildir. Core API frekansı Hz, modülleri Pa ve
sıcaklığı K cinsinden bekler.
