# Dynamic Elastomer Model

## Amaç ve kapsam

Elastomerin harmonik uyarım altındaki kayma davranışı, frekans ve sıcaklığa
bağlı kompleks dinamik kayma modülüyle temsil edilir. V0.2.0 bu büyüklükleri
saklayan veri yapısını ve kayıp faktörü hesabını tanımlar; interpolasyon, eğri
uydurma, Prony serisi ve çözüm algoritması içermez. V0.7 bu legacy veri türünü
değiştirmeden, ayrı tabulated provider ve harmonic orchestration ekler;
interpolation sözleşmesi
[`dynamic_modulus_interpolation.md`](dynamic_modulus_interpolation.md)
belgesindedir.

## Kompleks kayma modülü

TMS26, kompleks kayma modülü için aşağıdaki işaret gösterimini kullanır:

```text
G* = G' + iG''
```

Burada `i` sanal birimdir. `G'` kompleks modülün reel, `G''` ise sanal
bileşenidir. Her iki bileşen de ölçüm frekansı `f` ve mutlak sıcaklık `T` ile
birlikte değerlendirilir:

```text
G* = G*(f, T)
```

## Storage modulus

Depolama modülü `G'`, harmonik bir çevrim sırasında elastik olarak depolanan ve
geri kazanılabilen enerjiyle ilişkili bileşendir. Malzemenin dinamik kayma
rijitliği davranışını temsil eder.

## Loss modulus

Kayıp modülü `G''`, harmonik bir çevrim sırasında iç sürtünme nedeniyle
dağıtılan enerjiyle ilişkili bileşendir. Malzemenin viskoz kayıp davranışını
temsil eder.

## Loss factor

Kayıp faktörü, kayıp modülünün depolama modülüne oranıdır:

```text
tan(delta) = G'' / G'
```

Bu boyutsuz oran sönüm davranışının göstergesidir. Hesabın geçerli olması için
`G' > 0` önkoşulu uygulanır. Genel amaçlı `calculate_loss_factor` yardımcı
yordamı bu önkoşulun çağıran tarafından sağlandığını varsayar. V0.2.1
dinamik burulma rijitliği solver'ı ise kendi hesabından önce `G' > 0` ve
`G'' >= 0` koşullarını doğrular.

## Kompleks burulma rijitliğine bağlantı

V0.2.1 ana TVD modeli, aynı çalışma noktasındaki kompleks kayma modülünü
rijit iç göbek ve rijit dış halkaya tam bağlı annüler elastomer burca
aşağıdaki şekilde uygular:

```text
Cθ = 4π L ri² ro² / (ro² - ri²)
K* = K' + iK''
K' = G' Cθ
K'' = G'' Cθ
```

`ri` ve `ro` metre birimli iç/dış yarıçap, `L` bağlı silindirik
yüzeyin metre birimli eksenel genişliğidir. `Cθ` birimi `m³`, K' ve K''
birimleri `N·m/rad` değeridir. Aynı geometri faktörü iki modül bileşenine de
uygulandığı için:

```text
tan(delta) = G''/G' = K''/K'
```

Bu eşitlik `G' > 0`, `G'' >= 0`, `0 < ri < ro` ve `L > 0` koşullarında
geçerlidir. `calculate_dynamic_torsional_stiffness`, malzemenin mevcut tek
çalışma noktası alanlarını kullanır; frekans noktası seçimi veya interpolasyon
yapmaz. Hesap, `ri <= 0`, `ro <= ri`, `L <= 0`, `G' <= 0` ve `G'' < 0`
girdilerini kabul etmez.

V0.2.2 sistem builder'ı K', K'', kayıp faktörü, referans frekansı ve sıcaklığı
iki ataletli TVD sistemine birlikte aktarır. Sönümsüz modal tahmin yalnız K'
kullanır; K'' kompleks özdeğer hesabına katılmaz ve hesaplanan frekansa göre
yeni bir G' seçimi veya iterasyon yapılmaz. Bu frozen-property sınırı
[`two_inertia_modal_model.md`](two_inertia_modal_model.md) belgesinde
ayrıntılandırılır.

Annüler bir milin eksen boyunca Saint-Venant burulmasına ait
`K = GJp/ℓ` denklemi farklı sınır koşullarına dayanır ve TVD burç solver'larında
kullanılmaz. İki modelin ayrıntılı ayrımı
[`torsional-physics-core.md`](torsional-physics-core.md) belgesindedir.

## Birimler

| Büyüklük | Sembol | SI birimi |
| --- | --- | --- |
| Depolama modülü | `G'` | `Pa` |
| Kayıp modülü | `G''` | `Pa` |
| Burç geometri faktörü | `Cθ` | `m³` |
| Frekans | `f` | `Hz` |
| Mutlak sıcaklık | `T` | `K` |
| Kayıp faktörü | `tan(delta)` | boyutsuz |

Mühendislik verisi MPa cinsindeyse çekirdek veri yapısına yazılmadan önce
`tms_units` modülündeki `mpa_to_pa` yordamıyla Pa birimine dönüştürülür.

## Varsayımlar ve geçerlilik sınırları

- Her veri noktası yalnızca kendi frekans ve sıcaklık koşulunu temsil eder.
- G' ve G'' aynı frekans, sıcaklık ve kayma koşulunda elde edilmiş kabul edilir.
- Genlik bağımlılığı ve nonlinear hiperelastik davranış modellenmez.
- Legacy single-point yordam noktalar arasında interpolation yapmaz. V0.7
  provider yalnız aynı measured isotherm üzerinde açık linear-frequency veya
  linear-log-frequency-axis interpolation yapar; sıcaklık-frekans kaydırması,
  WLF/Arrhenius ve TTS uygulamaz.
