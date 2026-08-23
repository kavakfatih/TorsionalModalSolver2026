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
`G' > 0` önkoşulu uygulanır. V0.2.0 veri modeli negatif veya fiziksel olmayan
girdileri kendiliğinden doğrulamaz.

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
