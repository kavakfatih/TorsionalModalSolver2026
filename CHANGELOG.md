# Değişiklik Günlüğü

Bu dosyada TMS26 projesindeki kullanıcıya dönük önemli değişiklikler kaydedilir.

Biçim, [Keep a Changelog](https://keepachangelog.com/tr-TR/1.1.0/) yaklaşımını
temel alır ve proje anlamsal sürümleme ilkelerini izlemeyi hedefler.

## [Yayımlanmamış]

## [0.4.0] - 2026-08-25

### Eklendi

- `(node_id,dof_type)` ile tanımlanan Physical DOF ve anlamlı
  `TORSIONAL_ROTATION` türü.
- Constraint'ten bağımsız tam Equation ID ile constraint sonrası ayrı Active
  Equation ID eşlemesi.
- Fixed ve prescribed value kayıtları için constraint veri modeli, hedef/DOF
  doğrulaması ve aktif denklem yönetimi.
- Full K/M sisteminden `Kr=K(active,active)` ve `Mr=M(active,active)` üreten
  storage-bağımsız direct-elimination katmanı.
- Kr/Mr, aktif Physical DOF listesi ve `q=Pq_r+q_p` / `phi=Pphi_r` result
  recovery bilgisini taşıyan reduced-system altyapısı.
- İki ve üç düğümlü sistem, tamamen constrained `0x0` sistem, node
  permütasyonu, recovery ve geçersiz constraint durumları için otomatik
  regresyon kapsamı.
- Constraint foundation mimari belgesi, reduction matematik belgesi ve Karar
  0009.

### Değiştirildi

- Proje sürümü `0.4.0` olarak güncellendi.
- Global M/K assembly constraint uygulamasından ayrıldı; tam Equation ID
  uzayında full matrisler üretilirken aktif Kr/Mr üretimi ayrı reduction
  katmanına taşındı.
- V0.3.0 Karar 0008'deki `equation_id=0` tabanlı active-only assembly yaklaşımı
  Karar 0009 ile değiştirildi; Karar 0008 tarihsel kayıt olarak korundu.
- Mevcut fixed-node giriş davranışı sıfır değerli constraint olarak korunurken
  tam ve aktif denklem kimlikleri public modelde açıkça ayrıldı.

### Sınırlamalar

- Prescribed değerler veri ve recovery altyapısında korunur; sıfırdan farklı
  değerler için RHS düzeltmesi veya yük çözümü yapılmaz.
- Eigen solver, LAPACK, sparse/CSR backend, MPC, Lagrange multiplier, penalty
  yöntemi ve contact bu sürümün kapsamında değildir.

## [0.3.2] - 2026-08-25

### Eklendi

- Atalet, yoğunluk, yarıçap, uzunluk, rijitlik, modül, frekans ve sıcaklık
  alanlarında IEEE `NaN`, pozitif sonsuz ve negatif sonsuz girdileri kapsayan
  sonluluk regresyonları.
- Göbek dış yarıçapı ile elastomer iç yarıçapı ve elastomer dış yarıçapı ile
  atalet halkası iç yarıçapı için mutlak `1e-12 m` ve bağıl `1e-9` toleranslı
  ara yüz sürekliliği doğrulaması.
- Uç atalet ve rijitlik ölçeklerinde sonlu sonuç, kararlı eşdeğer atalet,
  doğal frekans ve fixed-DOF indirgenmiş M/K regresyonları.
- Numerik güvenilirlik modeli, test sınırları ve izlenebilirlik tablosu.

### Değiştirildi

- Proje sürümü `0.3.2` olarak güncellendi; fizik denklemleri ve node-element
  assembly mimarisi değiştirilmeden girdi doğrulaması ile sayısal değerlendirme
  sağlamlaştırıldı.
- Test assertion yordamları, sonlu olmayan sonuç veya toleransların yanlış
  başarı üretmesini önleyecek şekilde sıkılaştırıldı.

## [0.3.1] - 2026-08-25

### Eklendi

- Tek torsional eleman işaret konvansiyonu, serbest rijit-cisim modu, üç
  düğümlü global assembly, DOF mapping sürekliliği ve matris kalite
  invariantlarını birlikte doğrulayan V0.3.1 regresyon testi.
- İki ataletli analitik frekansı mevcut solver ile assembled M/K modal
  residual, Rayleigh quotient ve kütle ortogonalliği üzerinden çapraz
  doğrulayan foundation benchmark kapsamı.
- Matematiksel modelleri, varsayımları, toleransları ve doğrulama sınırlarını
  açıklayan torsional validation belgesi ve doğrulama dizini indeksi.

### Değiştirildi

- Proje sürümü `0.3.1` olarak güncellendi; production solver API'leri ve
  V0.3.0 node-element/assembly mimarisi değiştirilmeden CTest kapsamı 53 teste
  çıkarıldı.

## [0.3.0] - 2026-08-24

### Eklendi

- Fiziksel torsional düğüm kimliklerini aktif solver denklem kimliklerinden
  ayıran, kısıtlı düğümlerde sıfır sentinel kullanan DOF haritası.
- Private allocatable depolama, boyut denetimi ve güvenli katsayı erişimi
  sağlayan genel dense matris türü.
- Lokal eleman katkılarını toplayan global torsional rijitlik matrisi ile düğüm
  polar ataletlerini diagonal yerleştiren global dönel atalet matrisi.
- Açık `node_id -> equation_id -> local -> global` dönüşümü kullanan saf
  stiffness ve inertia assembly yordamları.
- Tek eleman, üç düğümlü zincir, serbest sistem invariantları, kısıt eliminasyonu,
  0x0 tam kısıtlı sistem ve geçersiz girdi regresyonları.
- Global M/K assembly matematik belgesi ve Karar 0008 tasarım kaydı.

### Değiştirildi

- Proje sürümü `0.3.0` olarak güncellendi; CTest kapsamı 52 teste çıkarıldı.
- Mevcut `calculate_local_stiffness` korunarak standart eleman arayüzü için saf
  `get_local_stiffness` sarmalayıcısı eklendi.

## [0.2.4] - 2026-08-24

### Eklendi

- İki torsional uç için sabit boyutlu `local_matrix_2x2` veri taşıyıcısı.
- Lineer elemandan `Ke = k[[1,-1],[-1,1]]` lokal rijitlik katkısını üreten
  saf `calculate_local_stiffness` yordamı.
- Bilinen matris katsayılarını, simetriyi, sıfır satır toplamını, pozitif
  yarı-tanımlı enerjiyi, rijit-cisim null modunu ve negatif rijitlik reddini
  doğrulayan CTest kapsamı.
- Lokal torsional eleman matrisi türetimi ve Karar 0007 tasarım kaydı.

### Değiştirildi

- Proje sürümü `0.2.4` olarak güncellendi; CTest kapsamı 46 teste çıkarıldı.
- Genel torsional sistem belgeleri, lokal eleman katkısı ile henüz uygulanmayan
  global matrix assembly arasındaki sınırı açıklayacak şekilde genişletildi.

## [0.2.3] - 2026-08-24

### Eklendi

- V0.2.3 — Generalized Torsional System Foundation kapsamında polar atalet,
  başlangıç açısı ve sınır koşulu taşıyan genel torsional düğüm türü.
- İki düğüm arasında K rijitliği ve eşdeğer viskoz c katsayısı taşıyan genel
  torsional eleman türü.
- Private düğüm/eleman koleksiyonları, güvenli ekleme-okuma yordamları, aktif
  DOF sayımı ve sistem bütünlüğü doğrulaması sağlayan genel sistem modülü.
- Mevcut iki ataletli TVD sistemini serbest-serbest veya fixed-hub genel
  topolojiye dönüştüren geriye uyumlu köprü.
- Node, element, topoloji, geçersiz girdi ve Benchmark 004 dönüşüm regresyonları
  ile genel sistem fizik/matematik belgeleri ve Karar 0006.

### Değiştirildi

- Proje sürümü `0.2.3` olarak güncellendi; CTest kapsamı 44 teste çıkarıldı.
- Benchmark 004, sayısal sonuçları değişmeden genel iki-düğüm/bir-eleman
  gösterimini de belgeleyecek şekilde genişletildi.

## [0.2.2] - 2026-08-24

### Eklendi

- V0.2.2 — Torsional System Model Foundation kapsamında geometri, rijit gövde
  ataletleri ve kompleks elastomer rijitliğini birleştiren iki ataletli TVD
  sistem türü ile builder yordamı.
- Mevcut doğal frekans yordamını yeniden kullanan fixed-hub çözümü ve sıfır
  frekanslı rijit-cisim modu ile elastik bağıl modu analitik çözen
  serbest-serbest iki ataletli sistem yordamı.
- Homojen annüler göbek için hacim, kütle ve polar kütle ataleti hesabı.
- Fixed-hub ve serbest-serbest frekansları, `[1,1]` ile `[1,-J_h/J_r]` mod
  şekillerini, ölçekleme davranışlarını ve büyük göbek ataleti limitini
  doğrulayan testler ile Benchmark 004.
- İki ataletli modal denklemleri ve frozen-property sönümsüz yaklaşımın fiziksel
  sınırlarını açıklayan matematik, fizik ve mimari dokümantasyonu.

### Değiştirildi

- Proje sürümü `0.2.2` olarak güncellendi; CTest kapsamı 28 teste çıkarıldı.
- Annüler halka ve göbek kütle özellikleri ortak özel `pure` yardımcı yordamda
  birleştirildi; pozitif yoğunluk ve geçerli rijit gövde geometrisi zorunlu
  hale getirildi.
- Geliştirme görevleri için temiz Debug build, tam derleme, CTest ve commit
  sonrası GitHub Actions kontrolünü zorunlu kılan Definition of Done eklendi.
- Yeni Fortran modülleri için CMake kaydı, bağımlılık sırası, `pure` matematik
  yordamı ve fizik testi kuralları standartlaştırıldı.
- Annüler elastomer polar alan momenti geometri testinde doğrudan analitik
  referansla doğrulanacak şekilde test kapsamı güçlendirildi.
- Dinamik burulma rijitliği analitik testinin bağıl hata sınırı `1e-10`
  seviyesine sıkılaştırıldı ve fiziksel girdi sınırları CTest'e eklendi.
- V0.2.1.3 bakımında dinamik rijitlik regresyon kapsamı, `K' ∝ G'`,
  `K' ∝ L` ve `K'' ∝ L` ölçekleme kontrolleriyle tamamlandı.

### Düzeltildi

- TVD elastomer rijitliğinde eksen boyunca Saint-Venant burulmasına ait
  `GJp/ℓ` denklemi yerine, rijit göbek ve dış halkaya tam bağlı annüler
  kauçuk burç için `4πGLri²ro²/(ro²-ri²)` denklemi uygulanmaya başlandı.
- Statik ve dinamik solver'lar `m³` birimli ortak `Cθ` geometri faktörüne
  bağlandı; Benchmark 001/003 ve doğal frekans referansları düzeltildi.
- `ri = 0` değerinin tam bağlı silindirik burç modeli için fiziksel olarak
  geçersiz olduğu belgelenerek girdi doğrulamasına eklendi.
- `calculate_rubber_polar_area_moment` yordamının public PURE arayüzü dış ve iç
  yarıçapı metre cinsinden alan iki skaler argümanla uyumlu hale getirildi.
- Dinamik rijitlik testinin üretim yordamı yerine yerel bir hesap kopyasını
  sınaması giderildi; test doğrudan üretim modülüne bağlandı.
- Negatif veya sırasız yarıçaplar, pozitif olmayan etkin uzunluk ve depolama
  modülü ile negatif kayıp modülünün geçersiz sonuç üretmesi engellendi.

## [0.2.1] - 2026-08-23

### Eklendi

- K', K'', kayıp faktörü, frekans ve sıcaklığı taşıyan
  `complex_torsional_stiffness_t` veri türü.
- Dinamik elastomer malzemesi ile annüler kauçuk geometrisinden kompleks
  burulma rijitliği hesaplayan `calculate_dynamic_torsional_stiffness` yordamı.
- `rubber_geometry_t` için yeniden kullanılabilir polar alan momenti hesabı.
- G''/G' ile K''/K' eşitliğini yüzde 0,1'den küçük hata sınırında doğrulayan
  CTest testi ve EPDM Benchmark 003 referansı.
- Kompleks rijitliğin fiziksel etkileri, matematiksel bağlantısı, mimari veri
  akışı ve API kararı için dokümantasyon.

### Değiştirildi

- Mevcut statik burulma rijitliği solver'ı davranışı korunarak ortak polar
  alan momenti yordamını kullanacak şekilde düzenlendi.
- Proje sürümü `0.2.1` olarak güncellendi ve CTest kapsamı on teste çıkarıldı.

## [0.2.0] - 2026-08-23

### Eklendi

- Kompleks dinamik kayma modülünün G', G'', frekans ve sıcaklık bileşenlerini
  saklayan `dynamic_shear_modulus` veri türü.
- DMA ve modal test verilerine hazırlanmak için `material_frequency_point`
  veri türü ve malzeme içinde çoklu çalışma noktası altyapısı.
- Boyutsuz `tan(delta) = G'' / G'` kayıp faktörü hesabı ve CTest birim testi.
- Dinamik elastomer matematik modeli, kompleks burulma rijitliği hazırlık
  belgesi, mimari karar kaydı ve EPDM referans benchmark'ı.

### Değiştirildi

- `dynamic_rubber_material_t`, mevcut tek noktalı alanları korunarak dinamik
  frekans-sıcaklık veri noktalarını saklayacak şekilde genişletildi.
- Proje sürümü `0.2.0` olarak güncellendi ve CTest kapsamı dokuz teste çıkarıldı.

## [0.1.3] - 2026-08-23

### Eklendi

- Mimari, matematik, fizik, geliştirme ve karar dokümantasyonu dizin indeksleri.
- Tek serbestlik dereceli torsional frekans modeli matematik belgesi.
- Unit test, physics validation ve benchmark regression kategorilerini açıklayan
  hesap motoru test belgesi.

### Değiştirildi

- GitHub otomasyon kuralları zorunlu build, test ve commit sonrası CI kalite
  kapılarıyla netleştirildi.
- Basit annüler TVD benchmark belgesine fiziksel model ve hesap zinciri eklendi.

### Düzeltildi

- Windows CI iş akışında MinGW64 GNU Fortran paketi ve araç zinciri PATH
  doğrulaması eklendi.

## [0.1.2] - 2026-08-23

### Eklendi

- macOS ve Windows üzerinde CMake, Ninja ve GNU Fortran ile derleme/test yapan
  GitHub Actions CI iş akışları.
- Birim, entegrasyon, benchmark ve CI doğrulama kapsamını tanımlayan test
  stratejisi belgesi.
- Homojen annüler halka için kütle ve polar kütle atalet momenti hesabı.
- Lineer elastik annüler elastomer için burulma rijitliği hesabı.
- Tek serbestlik dereceli, sönümsüz sistem için doğal frekans hesabı.
- Üç torsional fizik yordamı için yüzde 0,1 hata sınırına sahip analitik testler.
- Basit annüler TVD referans benchmark tanımı ve beklenen sonuçları.

## [0.1.1] - 2026-08-23

### Eklendi

- Pi ve mühendislik birim dönüşüm sabitlerini sağlayan `tms_constants` modülü.
- mm → m, MPa → Pa ve derece → radyan dönüşümlerini sağlayan saf ve eleman bazlı
  `tms_units` yordamları.
- Elastomer, atalet halkası, göbek ve bileşik TVD geometrisi veri türleri.
- Dinamik elastomer için yoğunluk, G', G'', sıcaklık ve frekans veri türü.
- Sabitler, birimler, geometri ve malzeme modülleri için CTest testleri.
- SI tabanlı iç birim sözleşmesi ile çekirdek veri modeli belgeleri.

## [0.1.0] - 2026-08-23

### Eklendi

- Fortran 2018, CMake, Ninja ve CTest tabanlı başlangıç altyapısı.
- İlk tür tanımları modülü ve derleyici doğrulama testi.
- Temel proje dokümantasyonu ve geliştirme kuralları.
