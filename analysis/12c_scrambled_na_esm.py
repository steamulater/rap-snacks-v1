"""Scramble native_ala sequences and fold with ESMFold API. N_SCRAMBLES=3 per bar."""
import random, time
import numpy as np
import pandas as pd
import requests
from pathlib import Path

ROOT           = Path(__file__).parent.parent
ENRICHED       = ROOT / "data/aggregated_lines_v2_enriched.csv"
CANDIDATES_CSV = ROOT / "data/phase2_candidates.csv"
OUT_CSV        = ROOT / "outputs/bioreason/scrambled_na_esm.csv"  # reuse outputs dir
VALID_AA       = set("ACDEFGHIKLMNPQRSTVWY")
ESM_URL        = "https://api.esmatlas.com/foldSequence/v1/pdb/"
N_SCRAMBLES    = 3
SEED           = 42
DELAY          = 1.5

def shuffle_seq(seq, seed):
    rng = random.Random(seed)
    s = list(seq); rng.shuffle(s)
    return "".join(s)

def fold_esm(seq):
    try:
        r = requests.post(ESM_URL, data=seq, timeout=60,
                          headers={"Content-Type": "application/x-www-form-urlencoded"})
        if r.status_code != 200: return None
        bf = [float(l[60:66]) for l in r.text.splitlines()
              if l.startswith("ATOM") and len(l) >= 66]
        return float(np.mean(bf)) if bf else None
    except Exception as e:
        print(f"    ESM error: {e}"); return None

enriched   = pd.read_csv(ENRICHED)
candidates = pd.read_csv(CANDIDATES_CSV)
bar_ids    = list(candidates["bar_id"])
sub        = enriched[enriched["bar_id"].isin(bar_ids)].set_index("bar_id")

rows = []
for bar_id in bar_ids:
    na_seq = str(sub.loc[bar_id, "fasta_seq_native_alanine"]).upper().strip()
    if not set(na_seq) <= VALID_AA:
        print(f"  [SKIP] {bar_id}"); continue
    for sc_i in range(N_SCRAMBLES):
        sc_seq = shuffle_seq(na_seq, seed=SEED + sc_i)
        name   = f"{bar_id}_sc_na_{sc_i}"
        print(f"  ESMFold {name}", end="  ", flush=True)
        time.sleep(DELAY)
        plddt = fold_esm(sc_seq)
        print(f"pLDDT={'N/A' if plddt is None else f'{plddt:.3f}'}")
        rows.append({"name": name, "bar_id": bar_id, "bucket": "scrambled_na",
                     "sequence": sc_seq, "esm_plddt": plddt})

df = pd.DataFrame(rows)
OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT_CSV, index=False)
print(f"\nSaved {len(df)} rows → {OUT_CSV}")
print(f"mean pLDDT: {df['esm_plddt'].mean():.3f}")
print("\n--- Per-bar ---")
for bar_id, g in df.groupby("bar_id"):
    print(f"  {bar_id:8s}  mean={g['esm_plddt'].mean():.3f}")
