# Benchmark 001: Basit Annüler TVD

Bu benchmark, düzeltilmiş TMS26 torsional fizik zincirini tek bir analitik
örnekte doğrular. Homojen atalet halkasının kütlesi ve polar ataleti, tam
bağlı annüler elastomer burulma rijitliği ve tek serbestlik dereceli doğal
frekans ardışık hesaplanır.

## Model

```text
Hub
 |
Elastomer
 |
Inertia Ring
```

Göbek referans dönme bileşenidir. Elastomer, rijit göbeğin dış silindirik
yüzeyi ile rijit atalet halkasının iç silindirik yüzeyine tam bağlı annüler
kauçuk burç olarak lineer burulma rijitliği sağlar. Inertia ring ise sistemin
polar kütle ataletini temsil eder.

## Hesap zinciri

```text
Geometry
   ↓
Material
   ↓
Stiffness
   ↓
Inertia
   ↓
Natural Frequency
```

Geometri ve malzeme verileri SI birimlerine dönüştürüldükten sonra burulma
rijitliği ile polar kütle ataleti hesaplanır. Bu iki sonuç, sönümsüz tek
serbestlik dereceli doğal frekans modelinin girdileridir.

TVD elastomer rijitliği, eksen boyunca annüler mil burulmasından farklı olarak
aşağıdaki tam bağlı silindirik burç denklemini kullanır:

```text
Cθ = 4πLri²ro²/(ro²-ri²)
Kθ = G Cθ
```

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)

Kabul ölçütü, her hesap sonucundaki bağıl hatanın yüzde `0,1`'den küçük
olmasıdır. Bu tanım performans zamanlaması içermez; sayısal referans benchmark
olarak kullanılacaktır.
