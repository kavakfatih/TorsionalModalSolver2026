# Temel Birim Dönüşümleri

## Amaç ve kapsam

TMS26 hesap motoru fiziksel verileri SI birimleriyle saklar. V0.1.1 sürümünde
desteklenen dönüşümler uzunluk, basınç veya modül ve düzlem açısı içindir. Bu
dönüşümler yalnızca ölçek değiştirir; fiziksel model veya malzeme davranışı
uygulamaz.

## Milimetreden metreye

Milimetre cinsinden uzunluk `L_mm`, metre cinsinden `L_m` değerine aşağıdaki
eşitlikle dönüştürülür:

```text
L_m = L_mm × 10⁻³
```

Girdi birimi mm, çıktı birimi m'dir. Dönüşüm işareti korur ve tüm sonlu gerçel
değerler için doğrusaldır.

## Megapaskaldan paskala

MPa cinsinden basınç, gerilme veya modül `S_MPa`, Pa cinsinden `S_Pa` değerine
aşağıdaki eşitlikle dönüştürülür:

```text
S_Pa = S_MPa × 10⁶
```

Girdi birimi MPa, çıktı birimi `Pa = N/m²`'dir. Dönüşüm, değerin sıcaklık veya
frekans bağımlılığını değiştirmez.

## Dereceden radyana

Derece cinsinden düzlem açısı `θ_deg`, radyan cinsinden `θ_rad` değerine
aşağıdaki eşitlikle dönüştürülür:

```text
θ_rad = θ_deg × π / 180
```

Radyan boyutsuz bir SI türetilmiş birimidir. Dönüşüm işareti ve devir sayısını
korur; sonucu belirli bir açı aralığına indirgemez.
