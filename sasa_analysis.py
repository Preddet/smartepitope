from Bio.PDB import PDBParser, SASA
from Bio.PDB.Polypeptide import protein_letters_3to1
import json
import os

def calculate_sasa_and_merge(pdb_file, hybrid_json, output_json, virus_name="Virus"):
    if not os.path.exists(pdb_file) or not os.path.exists(hybrid_json):
        print(f"Uyarı: Dosyalar bulunamadı.")
        return

    print(f"\n--- {virus_name} İçin SASA ve Altın Epitop Analizi ---")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(virus_name, pdb_file)
    
    sr = SASA.ShrakeRupley()
    sr.compute(structure, level="R")
    
    with open(hybrid_json, 'r') as f:
        targets = json.load(f)
    
    try:
        chain_a = structure[0]['A']
    except:
        chain_a = list(structure[0].get_chains())[0]

    res_list = []
    for res in chain_a:
        if res.id[0] == ' ':
            try:
                aa = protein_letters_3to1.get(res.get_resname(), '')
                if aa:
                    res_list.append((res.id[1], aa))
            except:
                pass
                
    def get_seq_snippet(center_id):
        idx = -1
        for i, (rid, _) in enumerate(res_list):
            if rid == center_id:
                idx = i
                break
        if idx == -1: return "N-A"
        
        snippet = []
        for i in range(idx - 2, idx + 3):
            if 0 <= i < len(res_list):
                snippet.append(res_list[i][1])
            else:
                snippet.append("-")
        return "-".join(snippet)

    all_results = []
    golden_epitopes = []
    
    for item in targets:
        pos = item['pos']
        try:
            res = chain_a[pos]
            sasa_val = float(res.sasa)
            
            item['sasa'] = round(sasa_val, 2)
            item['is_surface'] = bool(sasa_val > 25.0)
            item['seq_snippet'] = get_seq_snippet(pos)
            
            all_results.append(item)
            
            if item['is_surface'] and float(item['final_score']) > 0.80:
                golden_epitopes.append(item)
                
        except KeyError:
            continue

    golden_epitopes.sort(key=lambda x: x['final_score'], reverse=True)
    top_3 = golden_epitopes[:3]

    with open(output_json, "w") as f:
        json.dump(all_results, f) 
        
    summary_file = output_json.replace(".json", "_top3.json")
    with open(summary_file, "w") as f:
        json.dump(top_3, f)

    print(f"[OK] En iyi 3 aday {summary_file} dosyasına yazıldı.")

if __name__ == "__main__":
    tasks = [
        ("data/7K8S.pdb", "data/sars2_hybrid.json", "data/sars2_final.json", "SARS-CoV-2"),
        ("data/5X58.pdb", "data/sars1_hybrid.json", "data/sars1_final.json", "SARS-CoV-1"),
        ("data/1RVX.pdb", "data/flu_hybrid.json", "data/flu_final.json", "Influenza")
    ]
    
    for pdb, hyb, out, name in tasks:
        calculate_sasa_and_merge(pdb, hyb, out, name)