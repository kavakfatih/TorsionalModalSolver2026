# Kompleks Burulma Rijitliği

## Kompleks rijitlik kavramı

Harmonik burulma altında lineer viskoelastik elastomer, uygulanan dönmeye hem
faz içi hem de faz dışı moment tepkisi verir. Bu davranış kompleks burulma
rijitliği ile gösterilir:

```text
K* = K' + iK''
```

Annüler elastomer için geometrik dönüşüm şöyledir:

```text
Jp = π/2 (ro⁴ - ri⁴)
K' = G' Jp / L
K'' = G'' Jp / L
```

Burada `Jp` kauçuk kesitin polar alan momenti (`m⁴`), `L` etkin elastomer
uzunluğu (`m`), G' ve G'' ise sırasıyla depolama ve kayıp kayma modülleridir
(`Pa`). K' ve K'' birimi `N·m/rad` değeridir.

## Elastomer enerji depolama davranışı

Depolama rijitliği K', moment tepkisinin dönme ile faz içi bileşenidir.
Elastomerin çevrim sırasında elastik enerji depolayıp geri verme kabiliyetini
ve TVD'nin geri çağırıcı momentini temsil eder.

## Elastomer enerji kaybı

Kayıp rijitliği K'', moment tepkisinin dönmeye göre faz dışı bileşenidir.
Elastomer içindeki viskoz kayıplar nedeniyle mekanik enerjinin bir bölümünün
her çevrimde ısıya dönüşmesiyle ilişkilidir.

## Kayıp faktörü

Aynı çalışma noktasında G' ve G'' bileşenlerine aynı `Jp/L` katsayısı
uygulandığından iki kayıp faktörü tanımı eşittir:

```text
tan(delta) = G''/G' = K''/K'
```

Kayıp faktörü boyutsuzdur ve lineer viskoelastik modelde malzemenin göreli
sönüm davranışını gösterir. Geçerlilik için G' ve K' sıfırdan büyük olmalıdır.

## TVD doğal frekansına etkisi

Depolama rijitliği K', mevcut tek serbestlik dereceli sönümsüz yaklaşımda
doğal frekansın rijitlik girdisine karşılık gelir:

```text
fn ≈ 1/(2π) sqrt(K'/Jmass)
```

K'' ise sönüm ve faz davranışını etkiler; frekans cevabındaki rezonans genliği,
fazı ve sönümlü rezonans konumu değerlendirilirken dikkate alınır. V0.2.1
`tms_frequency_solver` yordamını değiştirmez: bu yordam gerçek skaler rijitlikle
sönümsüz doğal frekans hesaplamayı sürdürür. Kompleks eigen veya frekans cevabı
çözümü bu sürümün kapsamı dışındadır.

## Girdi doğrulaması

Hesabın sonlu ve fiziksel anlamlı bir kompleks rijitlik üretmesi için aşağıdaki
önkoşullar zorunludur:

- İç ve dış yarıçaplar negatif olamaz.
- Annüler kesitte `ro > ri` olmalıdır.
- Etkin elastomer uzunluğu için `L > 0` olmalıdır.
- Depolama modülü için `G' > 0` olmalıdır.
- Pasif kayıp modeli için `G'' >= 0` olmalıdır.

`calculate_dynamic_torsional_stiffness` ve ortak polar alan momenti yordamı bu
koşulları `pure` niteliklerini koruyarak denetler; geçersiz girdide hesap
`error stop` ile sonlanır. Frekans ve sıcaklık bu modelde hesap girdisi olarak
dönüştürülmez, malzeme çalışma noktasından sonuca aynen aktarılır.

## Uygulama sınırları

- Malzeme homojen ve lineer viskoelastik kabul edilir.
- Geometri eksenel simetrik annüler kesit ve küçük deformasyon kabulüne dayanır.
- Solver, `dynamic_rubber_material_t` içindeki tek çalışma noktası alanlarını
  kullanır; `frequency_points` dizisinde seçim veya interpolasyon yapmaz.
- FEM, nonlinear hiperelastisite ve Prony serisi uygulanmaz.
