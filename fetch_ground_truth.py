"""
IEDB API'sinden doğrulanmış epitop verilerini çekip ground_truth.json dosyasını günceller.

Filtreler:
- Yalnızca Spike / Hemagglutinin proteinini hedefleyen epitoplar
- Pozisyon bilgisi olan kayıtlar (starting_position + ending_position dolu)
- En az bir "Positive" ölçüm sonucu olan kayıtlar
"""

import urllib.request
import json
import time
import os

IEDB_API = "https://query-api.iedb.org/epitope_search"

# organism_name: IEDB'nin kullandığı tam isim
# spike_keywords: Spike/HA proteinine filtre için anahtar kelimeler
VIRUS_CONFIGS = {
    "sars2": {
        "org_encoded": "SARS-CoV2",
        "spike_keywords": ["P0DTC2", "spike", "surface glycoprotein"],
        "limit": 2000,
    },
    "sars1": {
        "org_encoded": "SARS-CoV1",
        "spike_keywords": ["P59594", "spike", "surface glycoprotein"],
        "limit": 1000,
    },
    "flu": {
        "org_encoded": "Influenza%20A%20virus",
        "spike_keywords": ["hemagglutinin", "HA protein"],
        "limit": 1000,
    },
}

def is_spike_antigen(source_antigens, keywords):
    for ag in source_antigens:
        name = (ag.get("name") or "").lower()
        iri = (ag.get("iri") or "").lower()
        for kw in keywords:
            if kw.lower() in name or kw.lower() in iri:
                return True
    return False

def fetch_all(org_encoded, limit, timeout=120):
    url = (
        f"{IEDB_API}"
        f"?source_organism_names=cs.%7B{org_encoded}%7D"
        f"&structure_type=eq.Linear%20peptide"
        f"&limit={limit}"
    )
    for attempt in range(2):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except Exception as e:
            if attempt == 0:
                print(f"  Hata, tekrar deneniyor... ({e})")
                time.sleep(2)
                continue
            raise e

def process_epitopes(all_data, keywords):
    results = []
    seen_ranges = set()

    for entry in all_data:
        source_antigens = entry.get("curated_source_antigens") or []

        if not is_spike_antigen(source_antigens, keywords):
            continue

        measures = entry.get("qualitative_measures") or []
        if measures and all(m == "Negative" for m in measures):
            continue

        for ag in source_antigens:
            start = ag.get("starting_position")
            end = ag.get("ending_position")

            if start is None or end is None:
                continue
            if (end - start) < 4:
                continue

            range_key = (start, end)
            if range_key in seen_ranges:
                continue
            seen_ranges.add(range_key)

            seq = entry.get("linear_sequence", "")
            results.append({
                "name": seq if seq else f"IEDB_{entry.get('structure_id', '')}",
                "range": [start, end],
                "source": "IEDB",
                "iedb_id": entry.get("structure_id"),
            })

    results.sort(key=lambda x: x["range"][0])
    return results


def fetch_iedb_epitopes(virus_key, config):
    print(f"\n--- {virus_key.upper()} için IEDB sorgulanıyor ---")
    try:
        all_data = fetch_all(config["org_encoded"], config["limit"])
        print(f"  Toplam ham epitop: {len(all_data)}")
        results = process_epitopes(all_data, config["spike_keywords"])
        print(f"  Spike epitopları (pozisyonlu, pozitif): {len(results)}")
        return results
    except Exception as e:
        print(f"  Hata: {e}")
        return []


def run():
    # Mevcut ground_truth.json'u oku (var olan verileri koru)
    output_path = os.path.join("data", "ground_truth.json")
    existing = {}
    if os.path.exists(output_path):
        with open(output_path, "r", encoding="utf-8") as f:
            existing = json.load(f)

    output = dict(existing)  # mevcut verileri başlangıç noktası al

    for virus_key, config in VIRUS_CONFIGS.items():
        epitopes = fetch_iedb_epitopes(virus_key, config)
        if epitopes:  # sadece başarılı fetch'lerde güncelle
            output[virus_key] = epitopes
        else:
            print(f"  [{virus_key}] fetch başarısız, önceki veri korundu.")
        time.sleep(1)

    os.makedirs("data", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=4, ensure_ascii=False)

    print(f"\n[OK] {output_path} güncellendi.")
    for k, v in output.items():
        print(f"  {k}: {len(v)} epitop")


if __name__ == "__main__":
    run()
