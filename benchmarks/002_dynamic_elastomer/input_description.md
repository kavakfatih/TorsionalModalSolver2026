# Girdi Açıklaması

## Malzeme çalışma noktası

| Alan | Değer | Birim |
| --- | ---: | --- |
| Malzeme | EPDM örneği | - |
| Sıcaklık | 293,15 | K |
| Frekans | 100 | Hz |
| Depolama modülü G' | 1 | MPa |
| Kayıp modülü G'' | 0,1 | MPa |

Modül değerleri çekirdek veri yapısına aktarılmadan önce `tms_units` içindeki
`mpa_to_pa` yordamıyla Pa birimine dönüştürülür. G' ve G'' değerlerinin aynı
frekans ve sıcaklıkta elde edildiği kabul edilir.

Bu benchmark kayma genliği, geometri, yoğunluk veya deney cihazına özgü üstveri
tanımlamaz; interpolasyon, eğri uydurma ve Prony serisi uygulanmaz.
