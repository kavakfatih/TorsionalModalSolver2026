# Karar 0003: Dinamik Elastomer Veri Modeli

- Durum: Kabul edildi
- Tarih: 2026-08-23

## Bağlam

V0.1.x malzeme türü G', G'', frekans ve sıcaklığı tek bir çalışma noktası
olarak saklar. DMA ve modal test verileri ise aynı malzeme için birden fazla
frekans-sıcaklık noktası gerektirir. Veri modeli genişlerken mevcut solver ve
istemci kaynaklarının bozulmaması gerekir.

## Karar

- Kompleks dinamik kayma modülünün dört temel alanı
  `dynamic_shear_modulus` türünde tanımlanacaktır.
- `material_frequency_point`, bu türü genişleterek gelecekte ölçüm kaynağına
  özgü üstverilerin eklenebileceği bir malzeme noktası sınırı oluşturacaktır.
- `dynamic_rubber_material_t` içindeki mevcut alanlar korunacak ve çoklu
  çalışma noktaları allocatable `frequency_points` dizisinde saklanacaktır.
- Eski alanlarla yeni dizi otomatik eşzamanlanmayacaktır.
- Bu sürüm yalnızca depolama ve kayıp modüllerinden kayıp faktörünü hesaplayacak;
  interpolasyon, eğri uydurma ve solver bağlantısı eklemeyecektir.

## Sonuçlar

- Mevcut tek noktalı kullanım kaynak uyumluluğunu korur.
- Bir malzeme sıfır veya daha fazla dinamik çalışma noktası taşıyabilir.
- Aynı fiziksel değerler için iki temsil bulunduğundan çağıran kodun seçimi
  açık olmalıdır; sessiz eşzamanlama kaynaklı belirsizlik oluşmaz.
- Gelecekteki veri doğrulama, interpolasyon ve deney üstverisi kararları bu
  temel türlerin üzerinde ayrı olarak ele alınabilir.
