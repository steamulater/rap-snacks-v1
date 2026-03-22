"""
03_parse_boltz.py
-----------------
Extract pLDDT, pTM, confidence, and pDE from Boltz-2 output files.

Handles multiple diffusion samples per bar (--diffusion_samples 5).
For each bar, parses all model_N files and reports:
  - Per-model metrics (one row per bar per model in detailed CSV)
  - Per-bar summary: mean, SD, best model by confidence (in summary CSV)

Expects outputs/boltz/predictions/ with the directory structure Boltz writes:
  predictions/
    bar_0/
      bar_0_model_0.pdb   bar_0_model_1.pdb  ... bar_0_model_4.pdb
      bar_0_confidence_model_0.json  ...  bar_0_confidence_model_4.json

Also verifies that each bar's sequence length matches bar_index_snapshot.json.
Length mismatches are flagged — this is the v1 data integrity failure mode.

Outputs:
  outputs/boltz/boltz_models.csv       one row per bar × model
  outputs/boltz/boltz_summary.csv      one row per bar (mean/SD/best across models)

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
BOLTZ_ID_MAP  = Path("data/boltz_id_map.json")   # b0 -> bar_0
BOLTZ_PREDICTIONS = Path("outputs/boltz/predictions")
MODELS_CSV = Path("outputs/boltz/boltz_models.csv")
SUMMARY_CSV = Path("outputs/boltz/boltz_summary.csv")


# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

def extract_seq_from_pdb(pdb_path: Path) -> str | None:
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
                seq.append(three_to_one.get(line[17:20].strip(), "?"))
                last_resnum = resnum
    return "".join(seq) if seq else None


def extract_seq_from_cif(cif_path: Path) -> str | None:
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
    m2 = re.search(r"\n([ACDEFGHIKLMNPQRSTVWY]{10,})\n", text)
    return m2.group(1) if m2 else None


# ---------------------------------------------------------------------------
# Confidence JSON parsing
# ---------------------------------------------------------------------------

def parse_confidence_json(json_path: Path) -> dict:
    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)
    return {
        "plddt":      data.get("plddt") or data.get("mean_plddt") or "",
        "ptm":        data.get("ptm") or "",
        "confidence": data.get("confidence_score") or data.get("confidence") or "",
        "pde":        data.get("pde") or data.get("mean_pde") or "",
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Parse Boltz-2 multi-sample outputs")
    p.add_argument("--boltz-dir", default=str(BOLTZ_PREDICTIONS))
    return p.parse_args()


def safe_float(val) -> float | None:
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def main():
    args = parse_args()
    boltz_dir = Path(args.boltz_dir)

    if not boltz_dir.exists():
        print(f"[ERROR] Boltz predictions dir not found: {boltz_dir}", file=sys.stderr)
        print("  Run Boltz-2 on Colab first, then download outputs/boltz/ here.")
        sys.exit(1)

    if not SNAPSHOT_JSON.exists():
        print(f"[ERROR] {SNAPSHOT_JSON} not found. Run 01_convert.py first.", file=sys.stderr)
        sys.exit(1)

    with open(SNAPSHOT_JSON, encoding="utf-8") as f:
        snapshot = json.load(f)

    # Load Boltz ID map if present (b0 -> bar_0)
    id_map = {}
    if BOLTZ_ID_MAP.exists():
        with open(BOLTZ_ID_MAP, encoding="utf-8") as f:
            id_map = json.load(f)
        print(f"  Boltz ID map loaded: {len(id_map)} entries (e.g. b0 -> {id_map.get('b0','?')})")

    bar_dirs = sorted(
        list(boltz_dir.glob("bar_*")) + list(boltz_dir.glob("b[0-9]*"))
    )
    if not bar_dirs:
        print(f"[ERROR] No bar_* directories found in {boltz_dir}")
        sys.exit(1)

    print(f"[03_parse_boltz.py] Found {len(bar_dirs)} bar directories in {boltz_dir}")

    all_model_rows = []
    all_summary_rows = []
    mismatches = []
    length_verified = 0

    for bar_dir in bar_dirs:
        boltz_dir_id = bar_dir.name          # e.g. "b0" or "bar_0"
        bar_id = id_map.get(boltz_dir_id, boltz_dir_id)  # resolve to internal bar_N
        snap = snapshot.get(bar_id, {})
        snap_len = snap.get("lyric_cleaned_len")
        snap_bar = snap.get("canonical_bar", "")
        snap_song = snap.get("genius_song_title", "")
        snap_seq  = snap.get("fasta_seq_concordance", "")

        # Find all model files for this bar
        struct_files = sorted(
            list(bar_dir.glob("*_model_*.pdb")) +
            list(bar_dir.glob("*_model_*.cif"))
        )
        json_files_map = {
            f.stem.split("_confidence_")[-1]: f
            for f in bar_dir.glob("*_confidence_model_*.json")
        }

        if not struct_files:
            print(f"  [WARN] No structure files for {bar_id}")
            continue

        model_metrics = []

        for struct_path in struct_files:
            # Extract model index from filename e.g. bar_0_model_2.pdb -> model_2
            m = re.search(r"(model_\d+)", struct_path.stem)
            model_key = m.group(1) if m else "model_0"
            model_n = int(model_key.split("_")[1]) if m else 0

            # Extract sequence (for integrity check — only needed once per bar)
            boltz_seq = None
            if model_n == 0:
                if struct_path.suffix == ".cif":
                    boltz_seq = extract_seq_from_cif(struct_path)
                else:
                    boltz_seq = extract_seq_from_pdb(struct_path)

                boltz_len = len(boltz_seq) if boltz_seq else None
                length_ok = (boltz_len == snap_len) if (boltz_len and snap_len) else None

                if length_ok is False:
                    mismatches.append({
                        "bar_id": bar_id,
                        "boltz_len": boltz_len,
                        "snapshot_len": snap_len,
                        "song": snap_song,
                    })
                elif length_ok:
                    length_verified += 1

            # Parse confidence JSON
            conf_path = json_files_map.get(model_key)
            if conf_path:
                conf = parse_confidence_json(conf_path)
            else:
                conf = {"plddt": "", "ptm": "", "confidence": "", "pde": ""}

            row = {
                "bar_id": bar_id,        # always internal bar_0 format
                "boltz_id": boltz_dir_id, # b0 or bar_0 as Boltz named it
                "model": model_n,
                "canonical_bar": snap_bar,
                "song": snap_song,
                "attribution": snap.get("attribution", ""),
                "iconicity": snap.get("aggregate_iconicity", ""),
                "seq_len": snap_len or "",
                "plddt": conf["plddt"],
                "ptm": conf["ptm"],
                "confidence": conf["confidence"],
                "pde": conf["pde"],
            }
            all_model_rows.append(row)
            model_metrics.append(conf)

        # Per-bar summary: mean, SD, best model by confidence
        conf_vals   = [safe_float(m["confidence"]) for m in model_metrics if safe_float(m["confidence"]) is not None]
        plddt_vals  = [safe_float(m["plddt"])      for m in model_metrics if safe_float(m["plddt"])      is not None]
        ptm_vals    = [safe_float(m["ptm"])         for m in model_metrics if safe_float(m["ptm"])        is not None]

        def mean(vals): return sum(vals) / len(vals) if vals else None
        def sd(vals):
            if len(vals) < 2: return None
            m = mean(vals)
            return (sum((v - m) ** 2 for v in vals) / len(vals)) ** 0.5

        best_model = None
        if conf_vals:
            best_idx = conf_vals.index(max(conf_vals))
            best_model = best_idx

        n_models = len(struct_files)
        summary = {
            "bar_id": bar_id,
            "canonical_bar": snap_bar,
            "song": snap_song,
            "attribution": snap.get("attribution", ""),
            "iconicity": snap.get("aggregate_iconicity", ""),
            "seq_len": snap_len or "",
            "n_models": n_models,
            "best_model": best_model if best_model is not None else "",
            "plddt_mean":       f"{mean(plddt_vals):.6f}"  if plddt_vals  else "",
            "plddt_sd":         f"{sd(plddt_vals):.6f}"    if plddt_vals and len(plddt_vals) > 1 else "",
            "plddt_best":       f"{max(plddt_vals):.6f}"   if plddt_vals  else "",
            "ptm_mean":         f"{mean(ptm_vals):.6f}"    if ptm_vals    else "",
            "ptm_sd":           f"{sd(ptm_vals):.6f}"      if ptm_vals and len(ptm_vals) > 1 else "",
            "ptm_best":         f"{max(ptm_vals):.6f}"     if ptm_vals    else "",
            "confidence_mean":  f"{mean(conf_vals):.6f}"   if conf_vals   else "",
            "confidence_best":  f"{max(conf_vals):.6f}"    if conf_vals   else "",
        }

        # Structural class based on best-model values
        try:
            plddt_b = float(summary["plddt_best"])
            ptm_b   = float(summary["ptm_best"])
            if plddt_b >= 0.6 and ptm_b >= 0.4:
                summary["structural_class"] = "confident_protein_like"
            elif plddt_b >= 0.6 and ptm_b < 0.4:
                summary["structural_class"] = "confident_alien"
            elif plddt_b < 0.6 and ptm_b >= 0.4:
                summary["structural_class"] = "uncertain_protein_like"
            else:
                summary["structural_class"] = "disordered"
        except (ValueError, TypeError):
            summary["structural_class"] = ""

        all_summary_rows.append(summary)

    # Write outputs
    MODELS_CSV.parent.mkdir(parents=True, exist_ok=True)

    if all_model_rows:
        with open(MODELS_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_model_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_model_rows)
        print(f"\n  Wrote {MODELS_CSV}  ({len(all_model_rows)} rows)")

    if all_summary_rows:
        with open(SUMMARY_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_summary_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_summary_rows)
        print(f"  Wrote {SUMMARY_CSV}  ({len(all_summary_rows)} bars)")

    # Integrity report
    print(f"\n  Sequence length verified (model_0): {length_verified} bars")
    if mismatches:
        print(f"\n  [WARN] {len(mismatches)} length mismatches — lyric attributions may be wrong:")
        for m in mismatches:
            print(f"    {m['bar_id']:10s}  boltz={m['boltz_len']}  snapshot={m['snapshot_len']}  {m['song']}")
        print("\n  Exclude these bars from lyric-level analysis.")
    else:
        print(f"  [OK] All sequence lengths match snapshot.")

    # Quick stats
    if all_summary_rows:
        valid = [r for r in all_summary_rows if r["plddt_mean"]]
        plddt_means = [float(r["plddt_mean"]) for r in valid]
        ptm_means   = [float(r["ptm_mean"])   for r in valid if r["ptm_mean"]]
        classes = {}
        for r in all_summary_rows:
            c = r.get("structural_class", "")
            classes[c] = classes.get(c, 0) + 1

        print(f"\n  Summary stats ({len(valid)} bars):")
        print(f"    mean pLDDT (best model): {max([float(r['plddt_best']) for r in valid if r['plddt_best']]):.4f} max")
        print(f"    mean pLDDT (mean model): {sum(plddt_means)/len(plddt_means):.4f} avg")
        if ptm_means:
            print(f"    mean pTM   (mean model): {sum(ptm_means)/len(ptm_means):.4f} avg")
        print(f"\n  Structural class breakdown:")
        for cls, n in sorted(classes.items(), key=lambda x: -x[1]):
            print(f"    {cls or '(unknown)':30s}  {n}")

    print(f"\nNext: python pipeline/05_enrich_csv.py")


if __name__ == "__main__":
    main()
