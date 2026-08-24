# Genelleştirilmiş Torsional Sistem

## Amaç ve kapsam

V0.2.3, özel iki ataletli TVD modelini değiştirmeden gelecekteki çok düğümlü
torsional sistemler için genel bir topoloji veri modeli sağlar. Bu katman fiziksel
düğümleri ve aralarındaki lineer bağlantıları saklar, aktif serbestlik derecesi
sayısını belirler ve topolojik bütünlüğü doğrular.

Bu sürüm matris assembly, özdeğer çözümü veya yeni bir modal fizik algoritması
içermez.

## Torsional düğüm

`torsional_node_t`, ortak bir dönme ekseni çevresinde tek skaler açısal
koordinat taşıyan yığılmış atalet noktasını temsil eder.

| Alan | Fiziksel ve matematiksel anlam | SI birimi |
| --- | --- | --- |
| `id` | Benzersiz topolojik düğüm etiketi | boyutsuz |
| `polar_inertia_kg_m2` | Açısal ivmeye karşı polar kütle ataleti `J` | `kg·m²` |
| `initial_angle_rad` | Başlangıç açısı `theta(0)` | `rad` |
| `constrained` | İdeal dönel kinematik sınır koşulu | boyutsuz |

Düğüm kimliği bir dizi indeksi veya doğrudan DOF sıra numarası değildir.
`initial_angle_rad` yalnız başlangıç koşulu üstverisidir; V0.2.3 bunu ön-gerilme,
denge çözümü veya modal rijitlik değişimine dönüştürmez.

## Torsional eleman

`torsional_element_t`, sistemde bulunan iki farklı düğüm arasındaki kütlesiz ve
karşılıklı lineer torsional bağlantıyı temsil eder.

| Alan | Fiziksel ve matematiksel anlam | SI birimi |
| --- | --- | --- |
| `id` | Benzersiz topolojik eleman etiketi | boyutsuz |
| `node_i_id`, `node_j_id` | Bağlantının uç düğüm kimlikleri | boyutsuz |
| `stiffness_nm_per_rad` | Bağıl açıya karşı rijitlik `K` | `N·m/rad` |
| `damping_nms_per_rad` | Bağıl açısal hıza karşı viskoz katsayı `c` | `N·m·s/rad` |

Fiziksel moment bağıntıları:

```text
T_K = K (theta_i - theta_j)
T_c = c (theta_i' - theta_j')
```

Farklı kimlikli paralel elemanlar fiziksel olarak geçerlidir ve yasaklanmaz.

## Sistem koleksiyonu ve DOF

`torsional_system_t`, düğüm ve eleman koleksiyonlarını private tutar. Ekleme,
okuma, sayma ve doğrulama yalnız public yönetim yordamları üzerinden yapılır.
Bu sınır, yinelenen kimlik veya tanımsız bağlantı ucu eklenmesini önler.

Her sabitlenmemiş düğüm bir aktif torsional DOF oluşturur:

```text
n_dof = sabitlenmemiş düğüm sayısı
```

Tamamen sabitlenmiş bir sistem `0` aktif DOF, iki serbest düğümlü sistem `2`
aktif DOF taşır. Serbest-serbest bir sistemin rijit-cisim modu veya ayrık alt
sistemler veri modeli açısından hata değildir. Bağlantılılık ve rijit-cisim
modu sayısı ilerideki analiz katmanının sorumluluğudur.

## İki ataletli TVD entegrasyonu

Mevcut `two_inertia_tvd_system_t`, yeni topolojide aşağıdaki gibi temsil edilir:

```text
Node 1: hub,          J = J_h
Node 2: inertia ring, J = J_r
Element 1: rubber,    K = K', c = 0
```

Serbest-serbest durumda iki düğüm serbesttir. Fixed-hub gösteriminde yalnız
Node 1 sabitlenir. `build_generalized_two_inertia_system` bu dönüşümü yapar;
mevcut fixed-hub ve serbest-serbest analitik solver'ları değiştirmez veya
yeniden uygulamaz.

Kompleks modeldeki `K''` kayıp rijitliği `N·m/rad`, genel elemandaki `c` ise
`N·m·s/rad` birimindedir. Bu iki büyüklük doğrudan eşdeğer değildir. Frekansa
bağlı bir dönüşüm veya viskoz model tanımlanmadığı için V0.2.3 köprüsü K''yi
sönüm alanına aktarmaz; `c = 0` kullanır. K'', kayıp faktörü, referans frekansı
ve sıcaklık mevcut iki ataletli veri türünde korunur.

## Fiziksel varsayımlar

- Her düğüm ortak eksen çevresinde tek torsional açı taşır.
- Polar atalet düğümde yığılmış, pozitif, sonlu ve zamandan bağımsızdır.
- Eleman kütlesiz, lineer ve pasiftir; `K > 0`, `c >= 0` olmalıdır.
- Açı ve şekil değiştirmeler küçüktür.
- Gyroskopik etkiler, dış zorlama ve nonlinear bağlantılar modellenmez.
- Düğüm ve eleman kimlikleri pozitif ve kendi kümelerinde benzersizdir.
- Her eleman sistemde bulunan iki farklı düğümü bağlar.

Gelecekteki M/K bağlantısının matematiksel çerçevesi
[`../mathematics/generalized_torsional_system_model.md`](../mathematics/generalized_torsional_system_model.md),
mimari karar ise
[`../decisions/0006-generalized-torsional-topology.md`](../decisions/0006-generalized-torsional-topology.md)
altında açıklanır.
