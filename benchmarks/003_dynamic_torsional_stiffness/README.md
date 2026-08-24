# Benchmark 003: Dinamik Burulma Rijitliği

Bu benchmark, EPDM örnek elastomerin tek frekans-sıcaklık çalışma noktasını
rijit iç göbek ile rijit dış halkaya tam bağlı annüler kauçuk burç
geometrisine uygulayarak kompleks burulma rijitliği bileşenlerini doğrular.

## Hesap zinciri

```text
G', G'', f, T
      ↓
Annüler geometri: ri, ro, L
      ↓
Cθ = 4πLri²ro²/(ro²-ri²)
      ↓
K' = G'Cθ     K'' = G''Cθ
      ↓
tan(delta) = K''/K' = G''/G'
```

## Analitik modelin fiziksel temeli

Formül, eş merkezli iki rijit silindir arasındaki elastomerin çevresel yer
değiştirme alanı `uθ(r) = Ar + B/r` olan eksenel simetrik lineer elastisite
çözümünden gelir. İç ve dış silindirik yüzeylerde kaymasız tam bağ
uygulanır; moment bu silindirik yüzeylerden aktarılır.

Varsayımlar:

- İç göbek ve dış halka rijit ve eş merkezlidir.
- Elastomer homojen, lineer viskoelastik ve küçük deformasyon bölgesindedir.
- Silindirik ara yüzlerde tam bağ vardır; kayma veya ayrılma yoktur.
- Eksenel uç etkileri ve geometrik nonlinearite ihmal edilir.
- `0 < ri < ro` ve `L > 0` koşulları geçerlidir.

Annüler milin uç yüzeylerden yüklenen Saint-Venant burulmasına ait
`K = GJp/ℓ` denklemi farklı sınır koşullarına sahiptir ve bu benchmarkta
kullanılmaz.

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)

Kabul ölçütü, Cθ, K', K'' ve kayıp faktörü sonuçlarında bağıl hatanın
`1e-10` değerinden küçük olmasıdır. Benchmark interpolasyon, FEM veya
kompleks doğal frekans çözümü içermez.
