import json
import itertools
import os

def load_json(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)

def calculate_overlap(predicted_pos, window_size, ground_truth_ranges):
    # Calculate the range of the predicted snippet
    half_window = window_size // 2
    pred_start = predicted_pos - half_window
    pred_end = predicted_pos + half_window
    
    # Check if this range overlaps with any ground truth range
    for gt in ground_truth_ranges:
        gt_start, gt_end = gt["range"]
        # Overlap condition
        if max(pred_start, gt_start) <= min(pred_end, gt_end):
            return True
    return False

def run_grid_search():
    print("--- SmartEpitope Kapalı Döngü Optimizasyon (Grid Search) ---")
    
    ground_truth = load_json("data/ground_truth.json")
    
    # Grid Search Parametreleri
    sasa_thresholds = [15.0, 20.0, 25.0, 30.0, 35.0]
    ai_biases = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
    window_sizes = [3, 5, 7, 9]
    
    viruses = {
        "sars2": "data/sars2_final.json",
        "sars1": "data/sars1_final.json",
        "flu": "data/flu_final.json"
    }
    
    optimization_results = {}
    
    for v_key, json_path in viruses.items():
        print(f"\n[{v_key.upper()}] Optimizasyon başlıyor...")
        
        if not os.path.exists(json_path):
            print(f"Hata: {json_path} bulunamadı.")
            continue
            
        pool = load_json(json_path)
        gt_ranges = ground_truth.get(v_key, [])
        
        if not gt_ranges:
            print(f"Uyarı: {v_key} için ground truth bulunamadı.")
            continue
            
        best_acc = -1.0
        best_params = {}
        best_top3 = []
        
        # Orijinal ağırlıkları havuzun ilk elemanından al
        base_w_ai = pool[0].get("w_ai", 0.5)
        base_w_stats = pool[0].get("w_stats", 0.5)
        
        total_combinations = len(sasa_thresholds) * len(ai_biases) * len(window_sizes)
        print(f"Toplam test edilecek kombinasyon: {total_combinations}")
        
        for sasa_thr, ai_bias, win_size in itertools.product(sasa_thresholds, ai_biases, window_sizes):
            
            # 1. Ağırlıkları sapma (bias) ile tekrar hesapla
            biased_w_ai = base_w_ai * ai_bias
            total_w = biased_w_ai + base_w_stats
            new_w_ai = biased_w_ai / total_w
            new_w_stats = base_w_stats / total_w
            
            # 2. S_final skorlarını güncelle ve filtrele
            valid_candidates = []
            for item in pool:
                # SASA filtresi
                if item.get("sasa", 0) <= sasa_thr:
                    continue
                    
                s_stats = item["stats_score"]
                s_ai = item["ai_score"]
                
                new_final = (new_w_stats * s_stats) + (new_w_ai * s_ai)
                
                valid_candidates.append({
                    "pos": item["pos"],
                    "sasa": item["sasa"],
                    "final_score": new_final,
                    "w_ai": new_w_ai,
                    "w_stats": new_w_stats
                })
                
            # 3. Sırala ve Top 3 al
            valid_candidates.sort(key=lambda x: x["final_score"], reverse=True)
            top3 = valid_candidates[:3]
            
            if not top3:
                continue
                
            # 4. Accuracy (Doğruluk) hesapla
            hits = 0
            for cand in top3:
                if calculate_overlap(cand["pos"], win_size, gt_ranges):
                    hits += 1
                    
            accuracy = (hits / len(top3)) * 100
            
            # En iyi parametreleri güncelle
            if accuracy > best_acc:
                best_acc = accuracy
                best_params = {
                    "sasa_threshold": sasa_thr,
                    "ai_bias": ai_bias,
                    "window_size": win_size,
                    "w_ai_final": new_w_ai,
                    "w_stats_final": new_w_stats
                }
                best_top3 = top3
                
                # Eğer %100 bulduysak daha aramaya gerek yok (isteğe bağlı, ama dursun)
                # if accuracy == 100.0: break
                
        print(f"--- {v_key.upper()} Optimizasyon Sonucu ---")
        print(f"Maksimum Doğruluk: %{best_acc:.1f}")
        print(f"En İyi Parametreler: SASA > {best_params['sasa_threshold']}, AI_Bias = {best_params['ai_bias']}x, Pencere = {best_params['window_size']}")
        print(f"Sonuç Ağırlıkları: Stats={best_params['w_stats_final']:.2f}, AI={best_params['w_ai_final']:.2f}")
        
        optimization_results[v_key] = {
            "accuracy": best_acc,
            "best_parameters": best_params,
            "top3_predicted_positions": [c["pos"] for c in best_top3]
        }
        
    with open("data/optimization_results.json", "w") as f:
        json.dump(optimization_results, f, indent=4)
        
    print("\n[OK] Optimizasyon tamamlandı. Sonuçlar data/optimization_results.json dosyasına kaydedildi.")

if __name__ == "__main__":
    run_grid_search()
