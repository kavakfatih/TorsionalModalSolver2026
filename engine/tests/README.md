# TMS26 Hesap Motoru Testleri

Bu dizindeki Fortran test programları CMake ile derlenir ve CTest üzerinden
çalıştırılır. Yeni veya değiştirilen her davranış uygun test kategorisiyle
kapsanmalıdır.

## Unit test

Tek modül, tür veya yordam izole olarak doğrulanır. `test_kinds`,
`test_constants`, `test_units`, `test_geometry`, `test_material` ve
`test_dynamic_modulus` bu kategoridedir. `test_dynamic_modulus`, G' ile G''
değerlerinin SI biriminde saklanmasını, çoklu çalışma noktası altyapısını ve
`tan(delta) = G'' / G'` hesabını doğrular.

## Physics validation

Fiziksel veya matematiksel yordam, bağımsız analitik sonuç ve açık bir hata
sınırıyla doğrulanır. `test_inertia`, `test_torsional_stiffness` ve
`test_frequency_solver` bu kategoridedir. V0.1.3 fizik doğrulamalarında bağıl
hata yüzde `0,1`'den küçük olmalıdır.

## Benchmark regression

Birden fazla fizik adımını temsil eden sabit referans modelin sonuçları zaman
içinde korunur. `benchmarks/001_simple_annular_tvd/` girdileri; kütle, polar
atalet, rijitlik ve doğal frekans testlerinin ortak regresyon temelidir. Model
veya kabul edilen formül değişirse benchmark girdileri, beklenen sonuçlar ve
ilgili testler aynı commit içinde güncellenir.

## Test ekleme

Yeni Fortran testi `test_<konu>.f90` biçiminde adlandırılır ve
`engine/CMakeLists.txt` içindeki `tms26_add_fortran_test` yordamıyla CTest'e
kaydedilir. Tüm testler aşağıdaki komutla çalıştırılır:

```sh
ctest --test-dir build --output-on-failure
```
