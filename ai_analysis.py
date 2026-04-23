import torch
from transformers import EsmTokenizer, EsmForMaskedLM
import json
import os
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

model_name = "facebook/esm2_t30_150M_UR50D"
tokenizer = EsmTokenizer.from_pretrained(model_name)
model = EsmForMaskedLM.from_pretrained(model_name)

def get_pdb_sequence(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    try:
        chain = structure[0]['A']
    except:
        chain = list(structure[0].get_chains())[0]
        
    seq = ""
    res_ids = []
    for res in chain:
        if res.id[0] == ' ':
            res_name = res.get_resname()
            try:
                aa = protein_letters_3to1.get(res_name, '')
                if aa:
                    seq += aa
                    res_ids.append(res.id[1])
            except:
                pass
    return seq, res_ids

def get_ai_scores(sequence):
    print("   -> ESM-2 Analizi Başlatılıyor...")
    inputs = tokenizer(sequence, return_tensors="pt")
    
    with torch.no_grad():
        results = model(**inputs).logits
    
    log_probs = torch.log_softmax(results, dim=-1)
    
    ai_scores = []
    for i, token_id in enumerate(inputs['input_ids'][0][1:-1]):
        score = log_probs[0, i+1, token_id].item()
        normalized_ai = 1 / (1 + abs(score)) 
        ai_scores.append(round(normalized_ai, 4))
        
    return ai_scores

def process_ai_for_virus(pdb_file, stats_json, output_json):
    if not os.path.exists(stats_json) or not os.path.exists(pdb_file):
        print(f"Hata: Gerekli dosyalar bulunamadı.")
        return

    pdb_seq, res_ids = get_pdb_sequence(pdb_file)
    ai_scores = get_ai_scores(pdb_seq)
    
    ai_dict = {}
    for i, res_id in enumerate(res_ids):
        if i < len(ai_scores):
            ai_dict[res_id] = ai_scores[i]
    
    with open(stats_json, 'r') as f:
        stats_data = json.load(f)
    
    hybrid_results = []
    for item in stats_data:
        res_id = item['pos']
        stats_val = item['score']
        ai_val = ai_dict.get(res_id, 0.0)
        
        hybrid_score = (stats_val * 0.5) + (ai_val * 0.5)
        
        hybrid_results.append({
            "pos": res_id,
            "stats_score": stats_val,
            "ai_score": ai_val,
            "final_score": round(hybrid_score, 4)
        })
    
    with open(output_json, "w") as f:
        json.dump(hybrid_results, f)
    print(f"[OK] Hibrit analiz tamamlandı: {output_json}")

if __name__ == "__main__":
    tasks = [
        ("data/7K8S.pdb", "data/sars2_conservation.json", "data/sars2_hybrid.json"),
        ("data/5X58.pdb", "data/sars1_conservation.json", "data/sars1_hybrid.json"),
        ("data/1RVX.pdb", "data/flu_conservation.json", "data/flu_hybrid.json")
    ]
    
    for pdb, stats, out in tasks:
        print(f"\n--- {out} İçin AI Analizi ---")
        process_ai_for_virus(pdb, stats, out)