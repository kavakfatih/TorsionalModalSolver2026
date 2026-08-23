# Benchmark 001: Basit Annüler TVD

Bu benchmark, V0.1.2 torsional fizik zincirini tek bir analitik örnekte
doğrular. Homojen atalet halkasının kütlesi ve polar ataleti, lineer elastomer
burulma rijitliği ve tek serbestlik dereceli doğal frekans ardışık hesaplanır.

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)

Kabul ölçütü, her hesap sonucundaki bağıl hatanın yüzde `0,1`'den küçük
olmasıdır. Bu tanım performans zamanlaması içermez; sayısal referans benchmark
olarak kullanılacaktır.
