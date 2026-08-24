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

Yeni modül veya katman eklendiğinde sorumluluk, bağımlılık yönü ve mevcut API'ye
etkisi burada belgelenmelidir.
