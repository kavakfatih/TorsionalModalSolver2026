# Karar 0009: Constraint ve İndirgenmiş Sistem Mimarisi

- Durum: Kabul edildi
- Tarih: 2026-08-25

## Bağlam

Karar 0008, V0.3.0 için kısıtlı düğümlere `equation_id=0` vererek aktif
denklem uzayını assembly öncesinde oluşturdu. Global K ve M yordamları sıfır
kimlikli satır/sütunları atlıyor ve doğrudan indirgenmiş matris üretiyordu. Bu
yaklaşım homojen fixed-node kapsamı için yeterliydi; ancak tam fiziksel denklem
kimliğini korumuyor, constraint politikasını assembly davranışına bağlıyor ve
recovery bilgisini ayrı bir solver girdisi olarak sunmuyordu.

V0.4.0 sonrasında modal, harmonik ve daha genel torsional network çözümlerine
hazırlanmak için Physical DOF, tam Equation ID ve Active Equation ID
kavramlarının ayrı yaşam döngülerine sahip olması gerekir.

## Karar

- Physical DOF `(node_id,dof_type)` çiftiyle tanımlanacaktır.
- İlk anlamlı DOF türü `TORSIONAL_ROTATION` olacaktır; çıplak tamsayı
  kullanımı public istemci sözleşmesine yayılmayacaktır.
- Tam Equation ID değerleri constraint durumundan bağımsız ve kesintisiz
  `1..n_physical` olarak atanacaktır.
- Constraint sonrasında serbest koordinatlar için ayrı, kesintisiz
  `1..n_active` Active Equation ID haritası üretilecektir.
- Fixed constraint, prescribed value veri modelinin sıfır değerli özel durumu
  olacaktır. Constraint kimliği, hedef node, DOF türü, tür ve değer açıkça
  saklanacaktır.
- Global K ve M, tam Equation ID uzayında constraint'ten bağımsız assemble
  edilecektir.
- Constraint reduction assembly içine gömülmeyecek; tam K/M ve aktif map
  tüketen ayrı bir direct-elimination katmanı olacaktır.
- V0.3 public `assemble_stiffness` / `assemble_inertia` yordamları kaynak
  uyumluluğu için full assembly ile ayrı reduction katmanını çağıran legacy
  wrapper'lar olarak korunacaktır; canonical V0.4 yolu `assemble_full_*` ile
  başlayacaktır.
- İndirgenmiş matrisler `Kr=K(active,active)` ve `Mr=M(active,active)` ile
  üretilecektir. Tüm DOF'ların constraint altında olduğu `0x0` sonuç geçerlidir.
- Üst seviye reduction API'si matris depolamasını adında kodlamayacaktır.
  Mevcut dense backend ilk gerçekleştirim olarak kalacaktır.
- Reduced-system sonucu Kr, Mr, aktif Physical DOF listesi, tam/aktif denklem
  eşlemesi ve result-recovery bilgisini birlikte taşıyacaktır.
- Genel durum recovery bağıntısı `q=P q_r+q_p`, gelecekteki homojen modal
  recovery bağıntısı `phi=P phi_r` olacaktır.
- Sıfırdan farklı prescribed değer saklanacak ve recovery sırasında
  kullanılabilecektir; RHS düzeltmesi veya yük çözümü V0.4.0 kapsamına
  alınmayacaktır.
- MPC, constraint equation, Lagrange multiplier, penalty, contact, sparse
  gerçekleştirim ve eigen solver bu kararın uygulama kapsamına alınmayacaktır.

Bu karar, Karar 0008'in fiziksel node ile denklem kimliğini ayırma, private
matris depolama ve lokal-global assembly ilkelerini korur. Karar 0008'deki
`equation_id=0` ile **aktif-only numaralandırma** ve **assembly içinde homojen
kısıt eliminasyonu** kararlarının yerini V0.4.0 için bu kayıt alır. Karar 0008
tarihsel V0.3.0 tasarım kaydı olarak değiştirilmeden korunur.

## Sonuçlar

- Aynı tam sistem farklı constraint kümeleriyle assembly yeniden yapılmadan
  indirgenebilir.
- Fiziksel Equation ID, constraint sonrası aktif numaralandırmadan bağımsız
  olarak izlenebilir kalır.
- Assembly, constraint ve solver sorumlulukları ayrı test edilebilir.
- Result recovery, gelecekteki eigen solver'dan tam fiziksel mod şekline dönüş
  için gerekli eşlemeyi hazırlar.
- Dense matris ilk backend olarak kalırken reduction sözleşmesi gelecekteki
  sparse backend'i engellemez.
- V0.3 istemcilerinin `equation_id=0` anlamına dayanan davranışı yeni tam
  numaralandırmada aynı semantiği taşımaz; geriye uyumlu fixed-node giriş yolu
  korunurken tam ve aktif kimliklerin yeni API üzerinden açıkça seçilmesi
  gerekir.
- Prescribed değer kabul edilmesi, non-homogeneous sistemin çözüldüğü anlamına
  gelmez; gereken RHS düzeltmesi sonraki bir kararın konusudur.
