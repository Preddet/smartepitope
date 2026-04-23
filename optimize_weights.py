import json
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression

def optimize_sars_weights(json_file):
    # 1. Veriyi Yükle
    with open(json_file, 'r') as f:
        data = json.load(f)
    df = pd.DataFrame(data)

    # 2. IEDB "Ground Truth" (Gerçek Doğrular) Etiketlemesi
    # Literatürde doğrulanmış başarılı epitop bölgeleri
    verified_ranges = [
        (177, 196), (331, 356), (370, 395), (403, 427), 
        (437, 461), (483, 493), (525, 546), (653, 666), 
        (1028, 1049), (1086, 1105)
    ]

    def check_if_verified(pos):
        for start, end in verified_ranges:
            if start <= pos <= end:
                return 1 # Başarılı Epitop
        return 0 # Başarısız/Bilinmeyen

    # AI ve Stats değerlerini çekelim (varsayımsal olarak verinde olduklarını kabul ediyoruz)
    # Eğer json'da sadece 'final_score' varsa, önceki kodlarından 'stats_val' ve 'ai_val'ı da eklemelisin.
    df['is_epitope'] = df['pos'].apply(check_if_verified)

    # 3. Lojistik Regresyon Eğitimi
    # Bağımsız değişkenler: Stats ve AI skorları (Burada item içindeki değerleri kullanıyoruz)
    # Not: final_targets.json oluştururken stats_val ve ai_val değerlerini de kaydettiysen burası çalışır.
    X = df[['stats_val', 'ai_val']] 
    y = df['is_epitope']

    model = LogisticRegression()
    model.fit(X, y)

    # 4. Katsayıları (Ağırlıkları) Al
    weights = model.coef_[0]
    # Katsayıları normalize edelim (toplamları 1 olacak şekilde)
    total = sum(abs(weights))
    w_stats = round(abs(weights[0]) / total, 2)
    w_ai = round(abs(weights[1]) / total, 2)

    print(f"\n--- SARS-CoV-2 İÇİN OPTİMİZE EDİLMİŞ AĞIRLIKLAR ---")
    print(f"Stats (Geçmiş) Ağırlığı: {w_stats}")
    print(f"AI (Gelecek) Ağırlığı: {w_ai}")
    print(f"\nYeni Formülün: S_final = {w_stats} * Stats + {w_ai} * AI")

    return w_stats, w_ai

# Not: final_targets.json dosyan stats_val ve ai_val içermelidir.
# optimize_sars_weights("final_targets.json")