"""
prep_boltz_fasta.py
-------------------
Generates Boltz-2-compatible input from bar_index_snapshot.json.

Boltz-2 requires a DIRECTORY of individual FASTA files — one .fasta per
protein. A multi-sequence FASTA is NOT supported. IDs must also be clean
(no underscores), so bar_0 -> b0, bar_1 -> b1, etc.

Outputs:
  outputs/boltz_inputs/b0.fasta  b1.fasta  ...  b84.fasta
  data/boltz_id_map.json                           <- b0 -> bar_0 mapping

Colab usage:
  !boltz predict boltz_inputs/ \\
      --use_msa_server \\
      --output_format pdb \\
      --diffusion_samples 5 \\
      --out_dir boltz_outputs/

Usage:
  python pipeline/prep_boltz_fasta.py
"""

import json
from pathlib import Path

SNAPSHOT    = Path("data/bar_index_snapshot.json")
OUT_DIR     = Path("outputs/boltz_inputs")
ID_MAP      = Path("data/boltz_id_map.json")
FASTA_WIDTH = 60

with open(SNAPSHOT) as f:
    snapshot = json.load(f)

# Sort by bar_n to preserve original ordering
bars = sorted(snapshot.items(), key=lambda x: x[1]["bar_n"])

OUT_DIR.mkdir(parents=True, exist_ok=True)

id_map   = {}   # boltz_id -> bar_id
written  = 0
skipped  = 0

for i, (bar_id, data) in enumerate(bars):
    boltz_id = f"b{i}"
    id_map[boltz_id] = bar_id

    seq = data.get("fasta_seq_concordance", "")
    if not seq:
        print(f"  [SKIP] {bar_id} — no concordance sequence")
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
    print(f"Skipped {skipped} bars (no concordance sequence)")
print(f"Wrote {ID_MAP}  ({len(id_map)} mappings)")
print(f"\nUpload boltz_inputs/ directory to Colab, then run:")
print(f"  !boltz predict boltz_inputs/ \\")
print(f"      --use_msa_server \\")
print(f"      --output_format pdb \\")
print(f"      --diffusion_samples 5 \\")
print(f"      --out_dir boltz_outputs/")
print(f"\nExample mapping:  b0 -> {id_map.get('b0', '?')}")
