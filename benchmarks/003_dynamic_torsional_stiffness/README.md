# Benchmark 003: Dinamik Burulma Rijitliği

Bu benchmark, EPDM örnek elastomerin tek frekans-sıcaklık çalışma noktasını
annüler geometriye uygulayarak kompleks burulma rijitliği bileşenlerini doğrular.

## Hesap zinciri

```text
G', G'', f, T
      ↓
Annüler geometri: ri, ro, L
      ↓
Jp = π/2 (ro⁴ - ri⁴)
      ↓
K' = G'Jp/L     K'' = G''Jp/L
      ↓
tan(delta) = K''/K' = G''/G'
```

- [Girdi açıklaması](input_description.md)
- [Beklenen sonuç](expected_result.md)

Kabul ölçütü, Jp, K', K'' ve kayıp faktörü sonuçlarında bağıl hatanın yüzde
`0,1`'den küçük olmasıdır. Benchmark interpolasyon, FEM veya kompleks doğal
frekans çözümü içermez.
