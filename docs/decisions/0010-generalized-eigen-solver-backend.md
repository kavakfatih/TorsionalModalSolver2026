# Karar 0010: Generalized Eigen Solver Backend

- Durum: Kabul edildi
- Tarih: 2026-08-26

## Bağlam

V0.4.0, constraint-aware reduced `K_r/M_r` sistemi ile fiziksel mode recovery
eşlemesini üretmektedir. V0.5.0'ın bu sistemi doğrulanabilir biçimde çözmesi,
fakat public modal API'yi belirli bir dense veya sparse algoritmaya kilitlememesi
gerekir.

## Karar

- V0.5.0 reference dense backend olarak LAPACK `DSYGV` kullanacaktır.
- Problem `ITYPE=1`, `JOBZ='V'` ve bütün backend boyunca `UPLO='U'` ile
  çözülecektir.
- LAPACK integer ABI LP64 olacaktır; CMake `BLA_SIZEOF_INTEGER=4` kullanacaktır.
- Public modal-analysis ve modal-result API'leri DSYGV, LAPACK ve raw dense
  storage kavramlarını açmayacaktır.
- Generalized problem, solver facade, DSYGV backend, eigen solution, modal
  validation ve modal result ayrı sorumluluklar olacaktır.
- Original K/M korunacak; backend yalnız çalışma kopyalarını değiştirecektir.
- K symmetric positive semidefinite olabilir. Singular K ve rigid-body mode
  hata değildir.
- M symmetric positive definite olmalıdır; authoritative SPD tanısı DSYGV
  `INFO` semantiğinden üretilecektir. İkinci bir Cholesky yazılmayacaktır.
- Mode shape'ler mass-normalized olacak; residual, M-orthogonality, repeated
  eigenspace ve sign ambiguity backend bağımsız katmanda doğrulanacaktır.
- Sonuç `linear`, `undamped`, `frozen-property` olarak etiketlenecektir.

## DSYGV seçiminin gerekçesi

DSYGV, küçük dense benchmarklarda kararlı ve yaygın bir generalized symmetric-
definite solver sağlar. V0.5.0'ın amacı büyük sparse production performansı
değil; doğru problem contract'ı, analytical verification ve gelecekteki
backend'ler için reference sonuç üretmektir.

## Gelecek sparse backend

Sparse matrix capability sonrasında iterative eigen extraction adayları ayrı
benchmark edilecektir:

- Lanczos / Block Lanczos
- Krylov–Schur
- LOBPCG
- ARPACK veya eşdeğeri
- SLEPc

Block Lanczos şimdiden zorunlu algoritma olarak kilitlenmez. Seçim;
convergence robustness, bellek, requested-mode count, clustered eigenvalue,
rigid modes, factorization ihtiyacı, platform desteği, lisans ve bakım maliyeti
ölçütleriyle yapılacaktır.

## Sonuçlar

- DSYGV kalıcı small-model reference backend olarak korunur.
- Yeni sparse/Lanczos-family backend eklemek modal-analysis veya result recovery
  kodunu yeniden yazmayı gerektirmez.
- V0.5.0 dense matris kopyası maliyetini kabul eder; sparse storage, CSR ve
  iterative eigensolver bu kararın uygulama kapsamı dışındadır.
- Frequency-dependent elastomer self-consistent iteration çözülmez.
