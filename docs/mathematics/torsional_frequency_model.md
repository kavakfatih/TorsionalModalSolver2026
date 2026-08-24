# Tek Serbestlik Dereceli Torsional Frekans Modeli

## Hareket denklemi

Tek serbestlik dereceli lineer torsional sistemin serbest titreşim denklemi:

```text
J * theta'' + c * theta' + k * theta = 0
```

Burada:

- `J`: dönme eksenine göre polar kütle atalet momenti (`kg·m²`)
- `c`: torsional viskoz sönüm katsayısı (`N·m·s/rad`)
- `k`: torsional stiffness (`N·m/rad`)
- `theta`: açısal yer değiştirme (`rad`)
- `theta'`: açısal hız (`rad/s`)
- `theta''`: açısal ivme (`rad/s²`)

Bu denklem, atalet momenti, sönüm momenti ve elastik geri çağırıcı moment
toplamının dış moment bulunmadığında sıfır olmasını ifade eder.

## Sönümsüz doğal frekans

V0.1.3 hesap çekirdeğinde doğal frekans için `c = 0` kabul edilir. Açısal doğal
frekans ve hertz cinsinden doğal frekans:

```text
omega_n = sqrt(k / J)
f_n = 1 / (2 * pi) * sqrt(k / J)
```

Pozitif ve sonlu `k` ile `J` için V0.3.2 uygulaması aynı matematiksel sonucu
aşağıdaki cebirsel olarak eşdeğer biçimde değerlendirir:

```text
omega_n = sqrt(k) / sqrt(J)
f_n = sqrt(k) / (2 * pi * sqrt(J))
```

Bu değerlendirme, sonuç temsil edilebilir olduğu halde `k/J` ara oranının sayı
aralığını aşma veya alt-taşma riskini azaltır. Fizik denklemi, birimler ve model
varsayımları değişmez.

Birimler:

- `k`: `N·m/rad`
- `J`: `kg·m²`
- `omega_n`: `rad/s`
- `f_n`: `Hz`

## Varsayımlar ve geçerlilik

- Malzeme ve sistem davranışı lineer elastiktir.
- Açısal yer değiştirmeler ve deformasyonlar küçüktür.
- Sistem tek serbestlik derecelidir.
- Sönümsüz doğal frekans hesabında `c = 0` kabul edilir.
- Burulma rijitliği ile polar kütle atalet momenti pozitif olmalıdır.
- Burulma rijitliği ile polar kütle atalet momenti sonlu olmalıdır.

Sönümlü doğal frekans, nonlinear malzeme davranışı, çok serbestlik dereceli
sistem ve eigen çözümü bu sürümün kapsamı dışındadır. Mevcut uygulama
`engine/src/solver/tms_frequency_solver.f90`, analitik doğrulama ise
`engine/tests/test_frequency_solver.f90` dosyasındadır. Uç ölçek sayısal
doğrulaması
[`../validation/numerical_robustness_validation.md`](../validation/numerical_robustness_validation.md)
belgesinde açıklanır.
