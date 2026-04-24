# SmartEpitope: Viral Epitope Discovery Pipeline (Prototip)

SmartEpitope, viral proteinler üzerinde potansiyel terapötik hedefleri (epitopları) tespit etmek amacıyla biyoinformatik araçları ve makine öğrenmesi modellerini birleştiren melez bir analiz sistemidir.

> **⚠️ Önemli Not:** 
> Bu proje şu an bir **prototip** aşamasındadır ve sunulan web arayüzü bir **demo** niteliği taşımaktadır. Dizilim hizalamaları, model çıkarımları (inference) ve optimizasyon gibi ağır veri işleme adımları tamamen arka planda (backend/lokal makine) gerçekleştirilmiş olup; web arayüzü yalnızca önceden hesaplanmış bu karmaşık verilerin araştırmacılar tarafından kolayca incelenebilmesi amacıyla host edilmektedir.

## Proje Hakkında

SmartEpitope'un temel amacı, bir virüsün yapısal bütünlüğünden ödün vermeden mutasyona uğratamayacağı "evrimsel çıkmaz sokakları" tespit etmektir. Prototip, şu an için aşağıdaki virüslerin analizini içermektedir:
- **SARS-CoV-2** (Spike Proteini)
- **SARS-CoV-1** 
- **Influenza A** (H1N1)

## Analiz Boru Hattı (Methodology)

Sistem, protein sekanslarını arka planda dört aşamalı bir doğrulama sürecinden geçirir:

### 1. Evrimsel Korunmuşluk (NCBI & MAFFT)
Virüse ait varyant sekansları indirilir ve **MAFFT** algoritması kullanılarak hizalanır. Hizalama sonucunda, her amino asit pozisyonu için Python `math` kütüphanesi ile **Shannon Entropisi** hesaplanır. Düşük entropiye sahip pozisyonlar, tarihsel olarak mutasyona direnmiş bölgeleri temsil eder.

### 2. Yapısal Yaşamsallık (Meta ESM-2)
Evrimsel korunmuşluk, bir bölgenin tek başına iyi bir hedef olduğunu kanıtlamaz. Bu aşamada Meta'nın **650 Milyon parametreli ESM-2** (`esm2_t33_650M_UR50D`) protein dil modeli devreye girer. Python üzerinden `Transformers (Hugging Face)` ve `PyTorch` kullanılarak çalıştırılan bu model, her pozisyon için in-silico (bilgisayar ortamında) mutasyonlar uygulayarak amino asidin yapısal açıdan değiştirilemezliğini tahmin eder.

### 3. Yüzey Erişilebilirliği (SASA)
Tahmin edilen bölgenin fiziksel olarak ilaçlara veya antikorlara açık olduğundan emin olmak için 3 boyutlu yapı üzerinden **SASA (Solvent Accessible Surface Area)** hesaplaması yapılır. Bu işlem, **Biopython (`Bio.PDB.SASA`)** kütüphanesindeki Shrake-Rupley algoritması ile gerçekleştirilir.

### 4. Hibrit Skorlama ve Optimizasyon
İstatistiksel korunmuşluk ve ESM-2 yapay zeka skorları dinamik bir ağırlıklandırma sistemiyle birleştirilir. Modelin doğruluğunu (Accuracy) maksimize etmek amacıyla, Python'un `itertools` modülü kullanılarak bir **Grid Search** (Izgara Arama) algoritması çalıştırılır. Sistem, IEDB ve PDB kaynaklı "Altın Standart" (Ground Truth) verilere karşı 120 farklı parametre kombinasyonunu test ederek en uygun konfigürasyonu kendi kendine belirler.

## Kullanılan Teknolojiler

- **Backend / Veri İşleme:** Python 3.11, PyTorch, Hugging Face Transformers, Biopython, MAFFT
- **Frontend / Demo Arayüzü:** HTML5, Vanilla CSS, JavaScript, Molstar Plugin (3D Görselleştirme)

## Dosya Yapısı ve Kullanım

Sistem ağır hesaplamalar gerektirdiği için kodlar doğrudan tarayıcıda çalışmaz. Arka plan betikleri şunlardır:
- `fetch_data.py`: NCBI ve UniProt'tan sekansları çeker.
- `ai_analysis.py`: ESM-2 dil modelini çalıştırarak yapısal tahminleri yapar.
- `sasa_analysis.py`: PDB verileri üzerinden 3D yüzey erişilebilirlik metriklerini hesaplar.
- `optimize_pipeline.py`: Kapalı döngü hiperparametre optimizasyonunu (Grid Search) gerçekleştirir.

### Demo Arayüzünü İnceleme
Analiz sonuçlarını görselleştirmek için ek bir web sunucusuna ihtiyaç yoktur. Bilgisayarınızda (veya GitHub Pages üzerinde) `index.html` dosyasını açarak, işlenmiş verilerin 3 boyutlu sunumunu interaktif olarak inceleyebilirsiniz.
