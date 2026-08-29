# Girdi Tanımı

## Common State

- Material: `SYNTHETIC-ELASTOMER`
- Batch/state: `BATCH-EXACT-01`
- Dynamic strain amplitude: `0.005` [-]
- Static prestrain: `0.02` [-]
- Mode: dynamic shear
- Birimler: frequency Hz, temperature K, modulus Pa

## Isotherm Seti

- Temperature: `[253.15, 273.15, 293.15, 313.15, 333.15] K`
- Explicit reference: `293.15 K`
- Exact TRS `s=log10(a_T)`: `[2, 1, 0, -1, -2]`
- Log-frequency grid:
  `[-2.0, -1.35, -0.65, 0.0, 0.55, 1.25, 2.0]`

Synthetic truth:

```text
log10(G'/Pa)  = 6.2 + 0.25 (log10(f/Hz) + s)
log10(G''/Pa) = 5.4 + 0.15 (log10(f/Hz) + s)
```

Non-TRS variant storage shifts'i korur; loss shifts
`[2.6,1.25,0,-1.25,-2.6]` kullanır. Bunlar production elastomer constants
değildir.
