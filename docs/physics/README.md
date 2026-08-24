# Fizik Dokümantasyonu

Bu dizin, TMS26 içindeki matematiksel denklemlerin temsil ettiği fiziksel
sistemleri, idealizasyonları ve mühendislik geçerlilik sınırlarını açıklar.

V0.2.2 fizik kapsamı; homojen annüler göbek ve atalet halkası, lineer elastik
elastomer, kompleks dinamik elastomer rijitliği ile fixed-hub ve serbest-serbest
analitik torsional sistem modellerini içerir. Modal tahmin frozen-property ve
sönümsüzdür; yalnız K' kullanılır, K'' kompleks özdeğer hesabına bağlanmaz.
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

Yeni fizik modeli eklenirken modelin fiziksel anlamı, kabulleri, geçerlilik
aralığı ve kapsam dışı davranışları bu dizinde belgelenmelidir.
