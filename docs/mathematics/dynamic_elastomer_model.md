# Dynamic Elastomer Model

## Amaç ve kapsam

Elastomerin harmonik uyarım altındaki kayma davranışı, frekans ve sıcaklığa
bağlı kompleks dinamik kayma modülüyle temsil edilir. V0.2.0 bu büyüklükleri
saklayan veri yapısını ve kayıp faktörü hesabını tanımlar; interpolasyon, eğri
uydurma, Prony serisi ve çözüm algoritması içermez.

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

V0.2.1, aynı çalışma noktasındaki kompleks kayma modülünü annüler elastomer
geometrisine aşağıdaki şekilde uygular:

```text
Jp = π/2 (ro⁴ - ri⁴)
K* = K' + iK''
K' = G' Jp / L
K'' = G'' Jp / L
```

`Jp` annüler kauçuk kesitinin `m⁴` birimli polar alan momenti, `L` ise `m`
birimli etkin burulma uzunluğudur. K' ve K'' değerleri `N·m/rad` birimindedir.
Aynı geometri katsayısı iki modül bileşenine de uygulandığı için:

```text
tan(delta) = G''/G' = K''/K'
```

Bu eşitlik `G' > 0`, geçerli annüler geometri ve `L > 0` koşullarında
geçerlidir. `calculate_dynamic_torsional_stiffness`, malzemenin mevcut tek
çalışma noktası alanlarını kullanır; frekans noktası seçimi veya interpolasyon
yapmaz. Hesap, negatif yarıçapı, `ro <= ri`, `L <= 0`, `G' <= 0` ve
`G'' < 0` girdilerini kabul etmez.

## Birimler

| Büyüklük | Sembol | SI birimi |
| --- | --- | --- |
| Depolama modülü | `G'` | `Pa` |
| Kayıp modülü | `G''` | `Pa` |
| Frekans | `f` | `Hz` |
| Mutlak sıcaklık | `T` | `K` |
| Kayıp faktörü | `tan(delta)` | boyutsuz |

Mühendislik verisi MPa cinsindeyse çekirdek veri yapısına yazılmadan önce
`tms_units` modülündeki `mpa_to_pa` yordamıyla Pa birimine dönüştürülür.

## Varsayımlar ve geçerlilik sınırları

- Her veri noktası yalnızca kendi frekans ve sıcaklık koşulunu temsil eder.
- G' ve G'' aynı frekans, sıcaklık ve kayma koşulunda elde edilmiş kabul edilir.
- Genlik bağımlılığı ve nonlinear hiperelastik davranış modellenmez.
- Noktalar arasında interpolasyon veya sıcaklık-frekans kaydırması yapılmaz.
