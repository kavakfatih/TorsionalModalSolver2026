# Fizik Dokümantasyonu

Bu dizin, TMS26 içindeki matematiksel denklemlerin temsil ettiği fiziksel
sistemleri, idealizasyonları ve mühendislik geçerlilik sınırlarını açıklar.

V0.2.1 fizik kapsamı; homojen annüler atalet halkası, lineer elastik elastomer,
tek serbestlik dereceli sönümsüz torsional titreşim modeli ve kompleks dinamik
elastomer rijitliğiyle sınırlıdır. Kompleks rijitlik henüz doğal frekans
çözümüne eklenmemiştir.
Ayrıntılı denklemler [`../mathematics/`](../mathematics/) altında, referans model
ise [`../../benchmarks/001_simple_annular_tvd/`](../../benchmarks/001_simple_annular_tvd/)
altında tutulur.

## Mevcut belgeler

- [`dynamic_torsional_stiffness.md`](dynamic_torsional_stiffness.md): kompleks
  kayma modülünden kompleks burulma rijitliğine geçiş için V0.2.0 hazırlığı
- [`complex_torsional_stiffness.md`](complex_torsional_stiffness.md): kompleks
  rijitlik, enerji depolama/kaybı ve TVD doğal frekansına etkisi

Yeni fizik modeli eklenirken modelin fiziksel anlamı, kabulleri, geçerlilik
aralığı ve kapsam dışı davranışları bu dizinde belgelenmelidir.
