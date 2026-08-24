# Girdi Açıklaması

## Analitik referans sistem

| Alan | Sembol | Değer | SI birimi |
| --- | --- | ---: | --- |
| Göbek polar kütle ataleti | `J_h` | 0.10 | `kg·m²` |
| Atalet halkası polar kütle ataleti | `J_r` | 0.20 | `kg·m²` |
| Depolama burulma rijitliği | `K'` | 1000 | `N·m/rad` |

Modal model için gerekli fiziksel girdiler bu üç değerdir. Regresyon testindeki
sistem üstverisi, dinamik rijitlik sözleşmesini temsil etmek üzere
`K'' = 100 N·m/rad`, `tan(delta) = 0.1`, referans frekansı `100 Hz` ve
sıcaklık `293.15 K` kullanır. Bu dört üstveri alanı V0.2.2 sönümsüz
frekanslarına veya mod şekillerine etki etmez.

## Sınır koşulları

- Fixed-hub durumda göbek açısal yer değiştirmesi sıfırdır.
- Serbest-serbest durumda göbek ve halka için dış torsional mesnet yoktur.
- Elastomer lineer ve küçük deformasyon bölgesindedir.
- Elastomerin dağıtılmış kütlesi ve polar ataleti ihmal edilir.
- `K'`, belirtilen malzeme çalışma noktasında sabitlenmiş kabul edilir;
  frekansa bağlı interpolasyon veya iterasyon yapılmaz.
