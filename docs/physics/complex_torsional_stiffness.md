# Kompleks Burulma Rijitliği

## Kompleks rijitlik kavramı

Harmonik burulma altında lineer viskoelastik elastomer, uygulanan dönmeye hem
faz içi hem de faz dışı moment tepkisi verir. Bu davranış kompleks burulma
rijitliği ile gösterilir:

```text
K* = K' + iK''
```

Rijit iç göbek ve rijit dış atalet halkasına tam bağlı annüler elastomer
burç için geometrik dönüşüm şöyledir:

```text
Cθ = 4π L ri² ro² / (ro² - ri²)
K' = G' Cθ
K'' = G'' Cθ
```

Burada `ri` ve `ro` silindirik bağ yüzeylerinin yarıçapları (`m`), `L`
bağlı yüzeyin eksenel genişliği (`m`) ve `Cθ` burç geometri faktörüdür
(`m³`). G' ve G'' sırasıyla depolama ve kayıp kayma modülleridir (`Pa`).
K' ve K'' birimi `N·m/rad` değeridir.

Bu modelde moment, elastomerin uç yüzeylerinden değil iç ve dış silindirik
bağ yüzeylerinden aktarılır. Bu nedenle annüler milin eksen boyunca
burulmasına ait `GJp/ℓ` denklemi TVD solver'ında kullanılmaz.

## Elastomer enerji depolama davranışı

Depolama rijitliği K', moment tepkisinin dönme ile faz içi bileşenidir.
Elastomerin çevrim sırasında elastik enerji depolayıp geri verme kabiliyetini
ve TVD'nin geri çağırıcı momentini temsil eder.

## Elastomer enerji kaybı

Kayıp rijitliği K'', moment tepkisinin dönmeye göre faz dışı bileşenidir.
Elastomer içindeki viskoz kayıplar nedeniyle mekanik enerjinin bir bölümünün
her çevrimde ısıya dönüşmesiyle ilişkilidir.

## Kayıp faktörü

Aynı çalışma noktasında G' ve G'' bileşenlerine aynı `Cθ` faktörü
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
fazı ve sönümlü rezonans konumu değerlendirilirken dikkate alınır. V0.2.2
`tms_torsional_system`, fixed-hub ve serbest-serbest analitik modal tahminlerde
kompleks sonuçtan yalnız K' bileşenini mevcut `tms_frequency_solver` yordamına
aktarır. K'' sistem verisinde korunur ancak sönümsüz özdeğer denklemine dahil
edilmez. Kompleks özdeğer veya frekans cevabı çözümü bu sürümün kapsamı
dışındadır; frozen-property sınırı
[`two_inertia_torsional_system.md`](two_inertia_torsional_system.md) belgesinde
açıklanır.

## Girdi doğrulaması

Hesabın sonlu ve fiziksel anlamlı bir kompleks rijitlik üretmesi için aşağıdaki
önkoşullar zorunludur:

- İç yarıçap pozitif olmalıdır; `ri = 0` bağlı iç silindirik yüzey
  oluşturmadığından bu modelin kapsamı dışındadır.
- Annüler kesitte `ro > ri` olmalıdır.
- Bağlı eksenel genişlik için `L > 0` olmalıdır.
- Depolama modülü için `G' > 0` olmalıdır.
- Pasif kayıp modeli için `G'' >= 0` olmalıdır.

`calculate_dynamic_torsional_stiffness` ve
`calculate_annular_bush_torsion_geometry_factor` bu koşulları `pure`
niteliklerini koruyarak denetler; geçersiz girdide hesap `error stop` ile
sonlanır. Frekans ve sıcaklık bu modelde hesap girdisi olarak dönüştürülmez,
malzeme çalışma noktasından sonuca aynen aktarılır.

## Uygulama sınırları

- Malzeme homojen ve lineer viskoelastik kabul edilir.
- Geometri eş merkezli rijit silindirler, tam bağlı annüler elastomer,
  küçük deformasyon ve ihmal edilen uç etkileri kabulüne dayanır.
- Solver, `dynamic_rubber_material_t` içindeki tek çalışma noktası alanlarını
  kullanır; `frequency_points` dizisinde seçim veya interpolasyon yapmaz.
- FEM, nonlinear hiperelastisite ve Prony serisi uygulanmaz.
