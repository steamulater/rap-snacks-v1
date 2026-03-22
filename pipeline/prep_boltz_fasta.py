"""
prep_boltz_fasta.py
-------------------
Generates a Boltz-2-compatible FASTA from bar_index_snapshot.json.

Boltz-2 rejects record IDs containing underscores. This script rewrites
bar_0, bar_1, ... -> b0, b1, ... and writes a mapping file so
03_parse_boltz.py can translate back to internal bar IDs.

Outputs:
  outputs/fastas/bars_v2_concordance_boltz.fasta   <- upload this to Colab
  data/boltz_id_map.json                           <- b0 -> bar_0 mapping

Usage:
  python pipeline/prep_boltz_fasta.py
"""

import json
from pathlib import Path

SNAPSHOT   = Path("data/bar_index_snapshot.json")
OUT_FASTA  = Path("outputs/fastas/bars_v2_concordance_boltz.fasta")
ID_MAP     = Path("data/boltz_id_map.json")
FASTA_WIDTH = 60

with open(SNAPSHOT) as f:
    snapshot = json.load(f)

# Sort by bar_n to preserve original ordering
bars = sorted(snapshot.items(), key=lambda x: x[1]["bar_n"])

entries = []
id_map  = {}   # boltz_id -> bar_id

for i, (bar_id, data) in enumerate(bars):
    boltz_id = f"b{i}"
    id_map[boltz_id] = bar_id

    seq = data.get("fasta_seq_concordance", "")
    if not seq:
        print(f"  [SKIP] {bar_id} — no concordance sequence")
        continue

    # Minimal header — only the boltz_id (no spaces, no pipes)
    header = f">{boltz_id}"
    wrapped = "\n".join(seq[j:j+FASTA_WIDTH] for j in range(0, len(seq), FASTA_WIDTH))
    entries.append(f"{header}\n{wrapped}")

OUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
OUT_FASTA.write_text("\n\n".join(entries) + "\n")

ID_MAP.write_text(json.dumps(id_map, indent=2))

print(f"Wrote {OUT_FASTA}  ({len(entries)} sequences)")
print(f"Wrote {ID_MAP}  ({len(id_map)} mappings)")
print(f"\nUpload to Colab:  {OUT_FASTA}")
print(f"Example mapping:  b0 -> {id_map.get('b0', '?')}")
