# Benchmark 005 — Synthetic TTS Identification

Bu benchmark V0.8.1 experimental master-curve identification zincirini iki
deterministic synthetic family ile tanımlar:

1. tek horizontal shift ile exact collapse olan TRS family,
2. storage ve loss branch'leri farklı shift taşıyan non-TRS family.

Amaç gerçek EPDM/NR parametresi temsil etmek değil; pair objective, shift chain,
TRS evidence, provenance, runtime stitching ve V0.8.0 provider round-trip için
tekrarlanabilir known-truth sağlamaktır.

Girdi ayrıntıları [`input_description.md`](input_description.md), beklenen
sonuçlar [`expected_result.md`](expected_result.md) içindedir. Otomatik
regresyonlar `tms26.tts_identification` ve `tms26.tts_runtime_roundtrip`
CTest ailelerinde çalışır.
