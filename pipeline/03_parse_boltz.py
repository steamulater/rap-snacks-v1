"""
03_parse_boltz.py
-----------------
Extract pLDDT, pTM, confidence, and pDE from Boltz-2 output files.

Expects outputs/boltz/predictions/ with the directory structure Boltz writes:
  predictions/
    bar_0/
      bar_0_model_0.cif          (or .pdb if --output_format pdb was used)
      bar_0_confidence_model_0.json

Outputs:
  outputs/boltz/boltz_confidence_scores.csv

Also verifies that bar_N IDs in Boltz output match bar_index_snapshot.json.
Flags any bars whose Boltz sequence length doesn't match the snapshot — the
data integrity failure mode that bit v1.

Usage:
  python pipeline/03_parse_boltz.py
  python pipeline/03_parse_boltz.py --boltz-dir outputs/boltz/predictions
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path

SNAPSHOT_JSON = Path("data/bar_index_snapshot.json")
BOLTZ_PREDICTIONS = Path("outputs/boltz/predictions")
OUT_CSV = Path("outputs/boltz/boltz_confidence_scores.csv")


def extract_seq_from_cif(cif_path: Path) -> str | None:
    """Pull one-letter sequence from _entity_poly.pdbx_seq_one_letter_code in CIF."""
    THREE = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    text = cif_path.read_text(errors="replace")
    match = re.search(
        r"pdbx_seq_one_letter_code\n_entity_poly.*?\n;(.*?);", text, re.DOTALL
    )
    if match:
        raw = match.group(1).strip()
        aa_codes = re.findall(r"\(([A-Z]+)\)", raw)
        if aa_codes:
            return "".join(THREE.get(a, "?") for a in aa_codes)
    # Fallback: look for a raw one-letter sequence block
    m2 = re.search(r"\n([ACDEFGHIKLMNPQRSTVWY]{10,})\n", text)
    if m2:
        return m2.group(1)
    return None


def extract_seq_from_pdb(pdb_path: Path) -> str | None:
    """Extract CA-residue sequence from PDB ATOM records."""
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }
    seq = []
    last_resnum = None
    for line in pdb_path.read_text(errors="replace").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            resnum = line[22:26].strip()
            if resnum != last_resnum:
                resname = line[17:20].strip()
                seq.append(three_to_one.get(resname, "?"))
                last_resnum = resnum
    return "".join(seq) if seq else None


def parse_confidence_json(json_path: Path) -> dict:
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)
    return {
        "plddt": data.get("plddt") or data.get("mean_plddt"),
        "ptm": data.get("ptm"),
        "confidence": data.get("confidence_score") or data.get("confidence"),
        "pde": data.get("pde") or data.get("mean_pde"),
    }


def parse_args():
    p = argparse.ArgumentParser(description="Parse Boltz-2 outputs into CSV")
    p.add_argument("--boltz-dir", default=str(BOLTZ_PREDICTIONS))
    p.add_argument("--output", default=str(OUT_CSV))
    return p.parse_args()


def main():
    args = parse_args()
    boltz_dir = Path(args.boltz_dir)
    out_csv = Path(args.output)

    if not boltz_dir.exists():
        print(f"[ERROR] Boltz predictions dir not found: {boltz_dir}", file=sys.stderr)
        print("  Run Boltz-2 on Colab first, then download outputs/ to this directory.")
        sys.exit(1)

    if not SNAPSHOT_JSON.exists():
        print(f"[ERROR] {SNAPSHOT_JSON} not found. Run 01_convert.py first.", file=sys.stderr)
        sys.exit(1)

    with open(SNAPSHOT_JSON, encoding="utf-8") as f:
        snapshot = json.load(f)

    results = []
    mismatches = []

    # Find all bar_* subdirs
    bar_dirs = sorted(boltz_dir.glob("bar_*"))
    if not bar_dirs:
        print(f"[ERROR] No bar_* directories found in {boltz_dir}")
        sys.exit(1)

    print(f"[03_parse_boltz.py] Found {len(bar_dirs)} bar directories")

    for bar_dir in bar_dirs:
        bar_id = bar_dir.name  # e.g. "bar_0"

        # Find structure file (CIF or PDB)
        struct_files = list(bar_dir.glob("*_model_0.cif")) + list(bar_dir.glob("*_model_0.pdb"))
        json_files = list(bar_dir.glob("*_confidence_model_0.json"))

        if not struct_files:
            print(f"  [WARN] No structure file for {bar_id}")
            continue
        if not json_files:
            print(f"  [WARN] No confidence JSON for {bar_id}")
            continue

        struct_path = struct_files[0]
        json_path = json_files[0]

        # Extract sequence from structure
        if struct_path.suffix == ".cif":
            boltz_seq = extract_seq_from_cif(struct_path)
        else:
            boltz_seq = extract_seq_from_pdb(struct_path)

        # Parse confidence metrics
        conf = parse_confidence_json(json_path)

        # Verify against snapshot
        snap = snapshot.get(bar_id)
        if snap is None:
            print(f"  [WARN] {bar_id} not in snapshot")
            snap_seq = None
            snap_len = None
            snap_bar = None
            snap_song = None
        else:
            snap_seq = snap.get("fasta_seq_concordance", "")
            snap_len = snap.get("lyric_cleaned_len")
            snap_bar = snap.get("canonical_bar", "")
            snap_song = snap.get("genius_song_title", "")

        boltz_len = len(boltz_seq) if boltz_seq else None

        # Flag length mismatches (the v1 integrity failure mode)
        length_match = (boltz_len == snap_len) if (boltz_len is not None and snap_len is not None) else None
        if length_match is False:
            mismatches.append({
                "bar_id": bar_id,
                "boltz_len": boltz_len,
                "snapshot_len": snap_len,
                "song": snap_song,
            })

        results.append({
            "bar_id": bar_id,
            "canonical_bar": snap_bar or "",
            "song": snap_song or "",
            "boltz_seq_len": boltz_len or "",
            "snapshot_seq_len": snap_len or "",
            "length_match": length_match if length_match is not None else "unknown",
            "plddt": conf["plddt"] or "",
            "ptm": conf["ptm"] or "",
            "confidence": conf["confidence"] or "",
            "pde": conf["pde"] or "",
        })

    # Write output
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    if results:
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
        print(f"\n  Wrote {out_csv}  ({len(results)} bars)")

    # Report mismatches
    if mismatches:
        print(f"\n[WARN] {len(mismatches)} length mismatches — lyric attributions may be wrong:")
        for m in mismatches:
            print(
                f"  {m['bar_id']:10s}  boltz_len={m['boltz_len']}  "
                f"snapshot_len={m['snapshot_len']}  song={m['song']}"
            )
        print("\n  These bars should be excluded from lyric-level analysis.")
    else:
        print(f"\n[OK] All {len(results)} bars: Boltz sequence length matches snapshot.")

    print(f"\nNext: python pipeline/05_enrich_csv.py")


if __name__ == "__main__":
    main()
