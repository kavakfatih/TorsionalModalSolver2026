# Fizik Dokümantasyonu

Bu dizin, TMS26 içindeki matematiksel denklemlerin temsil ettiği fiziksel
sistemleri, idealizasyonları ve mühendislik geçerlilik sınırlarını açıklar.

TMS26 fizik kapsamı; homojen annüler göbek ve atalet halkası, lineer elastik
elastomer, kompleks dinamik elastomer rijitliği, fixed-hub/serbest-serbest
analitik modeller, genel torsional düğüm-eleman topolojisi ve 2x2 lokal eleman
rijitlik katkısını içerir. DOF eşleme ile dense global dönel M/K assembly
uygulanmıştır. V0.5 modal tahmin frozen-property ve sönümsüzdür. V0.6 direct
harmonic yolunda K', K'' ve viskoz C ayrı tutulur; çözüm lineer,
frequency-domain ve frozen-property kapsamındadır.
Ayrıntılı denklemler [`../mathematics/`](../mathematics/) altında; tek-DOF ve
iki ataletli referans modeller sırasıyla
[`Benchmark 001`](../../benchmarks/001_simple_annular_tvd/) ve
[`Benchmark 004`](../../benchmarks/004_two_inertia_tvd/) altında tutulur.

## Mevcut belgeler

- [`dynamic_torsional_stiffness.md`](dynamic_torsional_stiffness.md): kompleks
  kayma modülünden kompleks burulma rijitliğine geçiş için V0.2.0 hazırlığı
- [`complex_torsional_stiffness.md`](complex_torsional_stiffness.md): kompleks
  rijitlik, enerji depolama/kaybı ve TVD doğal frekansına etkisi
- [`two_inertia_torsional_system.md`](two_inertia_torsional_system.md):
  fixed-hub ve serbest-serbest TVD sınır koşulları, iki analitik mod ve
  frozen-property geçerlilik sınırı
- [`generalized_torsional_system.md`](generalized_torsional_system.md):
  genel torsional düğüm, eleman, aktif DOF, sistem doğrulaması ve iki-ataletli
  TVD topoloji köprüsü
- [`torsional_damping_models.md`](torsional_damping_models.md): K' depolama,
  K'' kayıp rijitliği ve viskoz c ayrımı; passive energy/power bağıntıları ve
  frozen-property sönüm kapsamı

Yeni fizik modeli eklenirken modelin fiziksel anlamı, kabulleri, geçerlilik
aralığı ve kapsam dışı davranışları bu dizinde belgelenmelidir.
