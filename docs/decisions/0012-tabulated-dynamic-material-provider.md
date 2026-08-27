# Karar 0012: Tabulated Dynamic Material Provider

- Durum: Kabul edildi
- Tarih: 2026-08-28
- Sürüm: V0.7.0

## Bağlam

V0.6 harmonic solver frozen K'/K''/C/M ile çalışıyordu. TVD test correlation
için deneysel dinamik elastomer `G'(f),G''(f)` verisinin explicit operating
state ve izlenebilir interpolation ile harmonic assembly'ye bağlanması
gerekiyordu. Bu davranışın V0.5 symmetric generalized eigenproblem'ini veya
V0.6 frozen API'sini değiştirmemesi gerekiyordu.

## Karar

- Primary constitutive büyüklükler direct dynamic shear `G'` ve `G''`'dir.
- Harmonic convention `exp(+i*omega*t)`; passive işaret `G''>=0`, `K''>=0`.
- Default policy `LINEAR_FREQUENCY`, optional policy
  `LINEAR_LOG_FREQUENCY`'dir; yalnız frequency axis logaritmiktir.
- Extrapolation, temperature interpolation, amplitude/prestrain interpolation,
  WLF, Arrhenius, TTS, Prony ve spline yoktur.
- Her provider tek measured isotherm ve tek validated material operating state
  temsil eder.
- Dataset/test mode explicit'tir; direct torsional binding yalnız SHEAR kabul
  eder ve tensile-to-shear conversion yapmaz.
- Provider abstraction harmonic solver'dan interpolation ayrıntısını gizler.
- Element tablosu topoloji nesnesine gömülmez; provider + element ID + bir kez
  hesaplanan annular C_theta ayrı binding'de tutulur.
- Material-aware harmonic orchestration ayrıdır. Bound elemanın stored nominal
  K'/K'' değerleri frozen yol için korunur, dynamic yolda override edilir ve
  double-count edilmez; c ayrı frequency-independent kanaldır.
- Tüm provider domain'i solver çağrısından önce değerlendirilir. M, C, DOF,
  constraint ve topology sweep boyunca yeniden oluşturulmaz.
- Material-aware result, mevcut harmonic result'i bileşimle kullanır ve
  singular çözüm noktalarında dahi material state trace'i korur.
- V0.5 modal solver frozen-property ve DSYGV tabanlı kalır.

## Sonuçlar

Bu sınır V0.8'de master-curve veya temperature-shift provider eklenirken
harmonic solver API'sinin yeniden yazılmasını gerektirmez. Buna karşılık V0.7
tablosu causal time-domain model veya Kramers–Kronig consistency garanti etmez.
Import/parsing ve unit inference ayrı I/O katmanının sorumluluğudur.

ANSYS ve Marc karşılaştırması yalnız doğrulanabilir prensip düzeyindedir:
experimental storage/loss properties, tabulated dependency ve ayrı shift
katmanı kavramlarıyla uyum hedeflenir; proprietary solver internals hakkında
varsayım yapılmaz.
