"""
prep_boltz_fasta_native_ala.py
-------------------------------
Generates Boltz-2-compatible input for the native_alanine condition.

native_alanine: letters that are canonical AAs pass through unchanged;
BJOZUX (non-standard, non-mappable) are replaced with Alanine (A).
This is the deterministic "what if we just Ala-sub the invalids" baseline.

Outputs:
  outputs/boltz_inputs_native_ala/na0.fasta  na1.fasta  ...
  data/boltz_id_map_native_ala.json              <- na0 -> bar_0 mapping

Colab usage:
  !boltz predict boltz_inputs_native_ala/ \\
      --use_msa_server \\
      --output_format pdb \\
      --diffusion_samples 5 \\
      --out_dir boltz_outputs_native_ala/

Usage:
  python pipeline/prep_boltz_fasta_native_ala.py
"""

import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import convert_native_alanine, clean as clean_bar

SNAPSHOT    = Path("data/bar_index_snapshot.json")
OUT_DIR     = Path("outputs/boltz_inputs_native_ala")
ID_MAP      = Path("data/boltz_id_map_native_ala.json")
FASTA_WIDTH = 60
LAMBDA_VAL  = 2.0
SEED        = 42

with open(SNAPSHOT) as f:
    snapshot = json.load(f)

bars = sorted(snapshot.items(), key=lambda x: x[1]["bar_n"])

OUT_DIR.mkdir(parents=True, exist_ok=True)

id_map   = {}
written  = 0
skipped  = 0

for i, (bar_id, data) in enumerate(bars):
    boltz_id = f"na{i}"
    id_map[boltz_id] = bar_id

    lyric_cleaned = data.get("lyric_cleaned", "")
    if not lyric_cleaned:
        print(f"  [SKIP] {bar_id} — no lyric_cleaned")
        skipped += 1
        continue

    rng = random.Random(SEED)
    seq, mask = convert_native_alanine(lyric_cleaned, LAMBDA_VAL, rng)

    if not seq:
        print(f"  [SKIP] {bar_id} — empty sequence after conversion")
        skipped += 1
        continue

    wrapped = "\n".join(seq[j:j+FASTA_WIDTH] for j in range(0, len(seq), FASTA_WIDTH))
    fasta_text = f">{boltz_id}|protein\n{wrapped}\n"

    out_path = OUT_DIR / f"{boltz_id}.fasta"
    out_path.write_text(fasta_text)
    written += 1

ID_MAP.write_text(json.dumps(id_map, indent=2))

print(f"Wrote {written} individual FASTA files to {OUT_DIR}/")
if skipped:
    print(f"Skipped {skipped} bars")
print(f"Wrote {ID_MAP}  ({len(id_map)} mappings)")

# Spot-check: show a few examples
print("\nSpot-check (concordance vs native_alanine):")
spot_bars = ["bar_27", "bar_11", "bar_38"]
for bar_id in spot_bars:
    conc = snapshot[bar_id].get("fasta_seq_concordance", "")
    boltz_id = next(k for k, v in id_map.items() if v == bar_id)
    na_seq_lines = (OUT_DIR / f"{boltz_id}.fasta").read_text().splitlines()[1:]
    na_seq = "".join(na_seq_lines)
    n_diff = sum(a != b for a, b in zip(conc, na_seq))
    print(f"  {bar_id} ({boltz_id}):  conc={conc[:40]}...")
    print(f"  {'':>10}  natv={na_seq[:40]}...")
    print(f"  {'':>10}  len={len(na_seq)}  positions differing: {n_diff}/{len(conc)}")
    print()

print("Upload boltz_inputs_native_ala/ to Colab, then run:")
print("  !boltz predict boltz_inputs_native_ala/ \\")
print("      --use_msa_server \\")
print("      --output_format pdb \\")
print("      --diffusion_samples 5 \\")
print("      --out_dir boltz_outputs_native_ala/")
print(f"\nExample mapping:  na0 -> {id_map.get('na0', '?')}")
