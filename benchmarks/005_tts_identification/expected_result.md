# Beklenen Sonuç

## Exact TRS

- Adjacent relative shifts exact synthetic farkları verir.
- Absolute shift chain `[2,1,0,-1,-2]` değerlerini numerical tolerance içinde
  geri kazanır.
- Reference satırı tam `s=0`, `a_T=1` olur.
- Joint residual sıfıra machine precision ölçeğinde yaklaşır.
- Storage/loss diagnostic shifts aynı olur.
- Bütün 35 original point experimental cloud'da provenance ile korunur.
- Runtime table strictly increasing ve single-valued olur.
- Existing V0.8.0 provider round-trip `f_r=a_Tf` ile independent truth
  G'/G'' değerlerini verir.

## Non-TRS

- Numerical identification sonucu üretilebilir.
- Joint residual exact TRS'den yüksek kalır.
- `|delta_s_storage-delta_s_loss|` sıfırdan belirgin biçimde farklıdır.
- VGP/Cole-Cole temperature/source clouds ayrışmayı korur.
- Universal TRS PASS/FAIL veya artificial smoothing uygulanmaz.

Recovered shift tolerance `1.5e-6` boyutsuz; exact integral/identity kontrolleri
`1e-12`–`1e-14` ölçeğindedir.
