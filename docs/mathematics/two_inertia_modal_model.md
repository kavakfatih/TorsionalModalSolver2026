# İki Ataletli Analitik Modal Model

## Kapsam

Bu belge, TMS26 V0.2.2 içindeki fixed-hub tek serbestlik dereceli ve
serbest-serbest iki serbestlik dereceli torsional vibration damper (TVD)
modellerinin analitik temelini tanımlar. Genel bir özdeğer çözücü kullanılmaz.

Serbest-serbest modelin genelleştirilmiş koordinat sırası şöyledir:

```text
q = [theta_h, theta_r]^T
```

Burada `theta_h` göbek, `theta_r` atalet halkası açısal yer değiştirmesidir.

## Fixed-hub model

Göbek sabitlendiğinde yalnız atalet halkası hareket eder:

```text
J_r theta_r'' + K' theta_r = 0
```

Sönümsüz doğal açısal frekans ve doğal frekans:

```text
omega_n = sqrt(K' / J_r)
f_n = omega_n / (2 pi)
```

Üretim kodu bu sonucu yeniden uygulamak yerine mevcut
`calculate_natural_frequency` yordamına `K'` ve `J_r` değerlerini verir.

## Serbest-serbest iki ataletli model

Göbek ile atalet halkası serbest olduğunda hareket denklemleri:

```text
J_h theta_h'' + K'(theta_h - theta_r) = 0
J_r theta_r'' + K'(theta_r - theta_h) = 0
```

Matris biçimi:

```text
M = [ J_h   0   ]       K = [ K'  -K' ]
    [  0   J_r  ]           [-K'   K' ]
```

Harmonik hareket `q = phi exp(i omega t)` kabulüyle karakteristik denklem:

```text
det(K - omega^2 M) = 0

omega^2 [omega^2 J_h J_r - K'(J_h + J_r)] = 0
```

İki analitik kök elde edilir:

```text
omega_0 = 0

omega_e^2 = K' (1/J_h + 1/J_r)
          = K' (J_h + J_r)/(J_h J_r)

f_e = omega_e/(2 pi)
```

Elastik mod için eşdeğer atalet şu biçimde de yazılabilir:

```text
J_eq = 1/(1/J_h + 1/J_r)
f_e = 1/(2 pi) sqrt(K'/J_eq)
```

Bu eşdeğerlik, üretim kodunun mevcut tek-DOF doğal frekans yordamını güvenli
biçimde yeniden kullanmasını sağlar.

## Mod şekilleri

DOF sırası `[göbek, atalet halkası]` ve göbek genliği `1` olacak biçimde
normalizasyon kullanılır.

Rijit-cisim modunda elastomer bağıl dönme yapmaz:

```text
phi_0 = [1, 1]^T
```

Elastik modda iki rijit gövde zıt yönde döner. Analitik özvektör bağıntısı:

```text
phi_e = [1, -J_h/J_r]^T
```

Bu şekil aynı zamanda rijit-cisim modu ile kütle ortogonalliğini sağlar:

```text
phi_0^T M phi_e = J_h - J_h = 0
```

## Semboller ve SI birimleri

| Sembol | Fiziksel anlam | SI birimi |
| --- | --- | --- |
| `J_h` | Göbek polar kütle ataleti | `kg·m²` |
| `J_r` | Atalet halkası polar kütle ataleti | `kg·m²` |
| `K'` | Depolama burulma rijitliği | `N·m/rad` |
| `K''` | Kayıp burulma rijitliği | `N·m/rad` |
| `theta_h`, `theta_r` | Açısal yer değiştirme | `rad` |
| `omega` | Doğal açısal frekans | `rad/s` |
| `f` | Doğal frekans | `Hz` |

## Varsayımlar ve frozen-property sınırı

- Göbek ve atalet halkası rijit, eksenel simetrik gövdelerdir.
- Elastomer lineer, küçük deformasyon bölgesinde ve iki rijit yüzeye tam
  bağlıdır.
- Elastomerin dağıtılmış kütlesi ve polar ataleti M matrisinde ihmal edilir;
  yalnız göbek ve atalet halkası rijit gövde ataletleri kullanılır.
- Modal tahmin sönümsüzdür ve yalnız depolama rijitliği `K'` kullanılır.
- `K''` sistem verisinde enerji kaybı üstverisi olarak korunur; sönüm oranı,
  sönümlü frekans veya kompleks özdeğer hesabı yapılmaz.
- `G'`, malzemenin belirtilen referans frekansı ve sıcaklığındaki değeriyle
  modal hesap boyunca sabit tutulur. Hesaplanan frekans referans frekanstan
  farklı olsa bile interpolasyon veya öz-tutarlı iterasyon yapılmaz.
- `J_h > 0`, `J_r > 0` ve `K' > 0` olmalıdır.

Bu nedenle V0.2.2 sonucu bir **frozen-property undamped modal estimate** olarak
yorumlanmalıdır.

## Kaynak notu ve doğrulama

İki rijit rotorun tek torsional yayla bağlandığı bu ayrık-parametreli model,
standart torsional vibration ve çok serbestlik dereceli titreşim teorisine
dayanır. Genel kuramsal başvuru olarak S. S. Rao'nun *Mechanical Vibrations*
eserindeki çok serbestlik dereceli sistemler ve torsional vibration bölümleri
kullanılabilir.

Sayısal referanslar
[`benchmarks/004_two_inertia_tvd/`](../../benchmarks/004_two_inertia_tvd/),
fiziksel yorum ve geçerlilik sınırları ise
[`../physics/two_inertia_torsional_system.md`](../physics/two_inertia_torsional_system.md)
altında verilir.
