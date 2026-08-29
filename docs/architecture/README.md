# Mimari Dokümantasyon

Bu dizin, TMS26 modüllerinin sorumluluklarını, bağımlılıklarını, veri akışını ve
katman sınırlarını açıklar.

## Mevcut belgeler

- [`core-data-model.md`](core-data-model.md): çekirdek geometri, malzeme,
  torsional fizik, lokal/global matris, DOF eşleme, genel düğüm-eleman
  topolojisi ve iki ataletli sistem modüllerinin veri modeli
- [`V0.4_constraint_foundation.md`](V0.4_constraint_foundation.md): Physical
  DOF, tam ve aktif denklem kimlikleri, constraint yönetimi, direct elimination,
  indirgenmiş M/K sistemi ve result recovery mimarisi
- [`V0.5_modal_eigen_solver.md`](V0.5_modal_eigen_solver.md): reduced Kr/Mr,
  backend-neutral generalized solver facade, dense DSYGV reference backend,
  modal doğrulama/sonuç katmanı ve fiziksel mode recovery akışı
- [`V0.6_frequency_domain_response.md`](V0.6_frequency_domain_response.md):
  ayrı K'/K''/C/M assembly ve reduction, complex dynamic stiffness, ZSYSVX
  facade, status-aware sweep, physical response recovery ve derived TVD sonuçları
- [`V0.7_dynamic_material_provider.md`](V0.7_dynamic_material_provider.md):
  tabulated G'(f)/G''(f) provider sınırı, ayrı eleman binding'i, material-aware
  K'/K'' preparation/assembly, frozen API uyumluluğu ve material trace
- [`V0.8_thermorheological_runtime.md`](V0.8_thermorheological_runtime.md):
  tabulated/WLF/Arrhenius shift provider sınırı, validated master-curve
  bileşimi, physical/reduced coordinate semantics, trace evolution ve mevcut
  V0.7 harmonic API reuse
- [`V0.8.1_master_curve_identification.md`](V0.8.1_master_curve_identification.md):
  experimental/runtime boundary, adjacent-pair identification, master cloud,
  deterministic stitching, result/status modeli ve V0.8.0 provider export
- [`V0.8.2_parametric_shift_law_identification.md`](V0.8.2_parametric_shift_law_identification.md):
  empirical/parametric ayrımı, pair-space Arrhenius/WLF fit, LOTO,
  identifiability ve explicit measured-domain runtime export

Yeni modül veya katman eklendiğinde sorumluluk, bağımlılık yönü ve mevcut API'ye
etkisi burada belgelenmelidir.
