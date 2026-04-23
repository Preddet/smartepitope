from Bio import Entrez, SeqIO
import time
import os

# NCBI İletişim Bilgisi
Entrez.email = "gkcelebi69@gmail.com" 

# Virüs Konfigürasyonları
VIRUS_CONFIGS = {
    "sars2": {
        "query": 'SARS-CoV-2[Organism] AND (Spike OR "surface glycoprotein" OR "S protein")',
        "len_range": (1200, 1300),
        "filename": "data/sars2_sequences.fasta"
    },
    "sars1": {
        "query": 'SARS-CoV[Organism] AND (Spike OR "surface glycoprotein" OR "S protein")',
        "len_range": (1200, 1300),
        "filename": "data/sars1_sequences.fasta"
    },
    "flu": {
        "query": 'Influenza A virus (H1N1)[Organism] AND (Hemagglutinin OR "HA protein")',
        "len_range": (500, 600), # Hemagglutinin genellikle 560 aa civarıdır
        "filename": "data/flu_sequences.fasta"
    }
}

def fetch_sequences(virus_key, limit=300):
    config = VIRUS_CONFIGS.get(virus_key)
    if not config:
        print(f"Hata: {virus_key} için konfigürasyon bulunamadı.")
        return

    print(f"\n--- {virus_key.upper()} İçin Sekanslar İndiriliyor ---")
    try:
        # 1. NCBI Araması
        handle = Entrez.esearch(db="protein", term=config["query"], retmax=limit)
        record = Entrez.read(handle)
        ids = record["IdList"]
        print(f"Toplam {len(ids)} potansiyel eşleşme bulundu.")

        if not ids:
            return

        # 2. Sekansları İndir
        fetch_handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
        all_sequences = list(SeqIO.parse(fetch_handle, "fasta"))
        
        # 3. Filtreleme
        cleaned = []
        for rec in all_sequences:
            min_len, max_len = config["len_range"]
            if (min_len <= len(rec.seq) <= max_len) and ("X" not in rec.seq):
                cleaned.append(rec)
        
        print(f"Filtreleme sonrası: {len(cleaned)} / {len(all_sequences)} sekans hazır.")

        # 4. Kaydet
        if not os.path.exists('data'):
            os.makedirs('data')
            
        SeqIO.write(cleaned, config["filename"], "fasta")
        print(f"[OK] {config['filename']} dosyasina kaydedildi.")
        return config["filename"]

    except Exception as e:
        print(f"Hata: {e}")
        return None

if __name__ == "__main__":
    # Demo amaçlı tüm virüsler için indir
    for v_key in VIRUS_CONFIGS.keys():
        fetch_sequences(v_key, limit=100) # Test için limiti düşük tutuyoruz