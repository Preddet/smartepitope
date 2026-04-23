import subprocess
import os
import numpy as np
import json
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import protein_letters_3to1

mafft_exe = os.path.join(os.getcwd(), "mafft.bat")

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
        if res.id[0] == ' ': # Standard amino acid
            res_name = res.get_resname()
            try:
                aa = protein_letters_3to1.get(res_name, '')
                if aa:
                    seq += aa
                    res_ids.append(res.id[1])
            except:
                pass
    return seq, res_ids

def analyze_virus(pdb_file, input_fasta, aligned_fasta, output_json):
    if not os.path.exists(input_fasta):
        print(f"Hata: {input_fasta} bulunamadı.")
        return

    print(f"\n--- {input_fasta} PDB ile Hizalanıyor ---")
    
    pdb_seq, res_ids = get_pdb_sequence(pdb_file)
    
    pdb_rec = SeqRecord(Seq(pdb_seq), id="PDB_REF", description="Reference from 3D Model")
    records = list(SeqIO.parse(input_fasta, "fasta"))
    records.insert(0, pdb_rec)
    
    combined_fasta = input_fasta.replace(".fasta", "_combined.fasta")
    SeqIO.write(records, combined_fasta, "fasta")

    try:
        with open(aligned_fasta, "w") as f_out:
            subprocess.run([mafft_exe, "--auto", combined_fasta], stdout=f_out, check=True, shell=True)
        print(f"[OK] Hizalama tamamlandı: {aligned_fasta}")
    except Exception as e:
        print(f"Hizalama hatası: {e}")
        return

    print(f"--- Korunmuşluk Analizi Yapılıyor ---")
    try:
        alignment = AlignIO.read(aligned_fasta, "fasta")
        
        pdb_aln = None
        for rec in alignment:
            if rec.id == "PDB_REF":
                pdb_aln = str(rec.seq)
                break
                
        if not pdb_aln:
            print("Hata: PDB_REF hizalamada bulunamadı.")
            return

        seq_length = alignment.get_alignment_length()
        results = []
        pdb_seq_idx = 0 
        
        for i in range(seq_length):
            pdb_char = pdb_aln[i]
            if pdb_char == '-':
                continue
                
            res_id = res_ids[pdb_seq_idx]
            pdb_seq_idx += 1
            
            column = alignment[:, i].upper()
            residues = [res for res in column if res in "ACDEFGHIKLMNPQRSTVWY"]
            
            if len(residues) == 0:
                entropy = 4.322
            else:
                _, counts = np.unique(residues, return_counts=True)
                probs = counts / len(residues)
                entropy = -np.sum(probs * np.log2(probs))
            
            score = 1 - (entropy / 4.322)
            
            results.append({
                "pos": res_id,
                "score": round(float(score), 4)
            })

        with open(output_json, "w") as f:
            json.dump(results, f)
        print(f"[OK] Sonuçlar kaydedildi: {output_json}")

    except Exception as e:
        print(f"Analiz hatası: {e}")

if __name__ == "__main__":
    tasks = [
        ("data/7K8S.pdb", "data/sars2_sequences.fasta", "data/sars2_aligned.fasta", "data/sars2_conservation.json"),
        ("data/5X58.pdb", "data/sars1_sequences.fasta", "data/sars1_aligned.fasta", "data/sars1_conservation.json"),
        ("data/1RVX.pdb", "data/flu_sequences.fasta", "data/flu_aligned.fasta", "data/flu_conservation.json")
    ]
    
    for pdb, inp, aln, out in tasks:
        analyze_virus(pdb, inp, aln, out)