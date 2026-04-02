"""
09_proteinmpnn_design.py
------------------------
Phase 2 Step 2 — ProteinMPNN soft sequence design.

Two modes:

  PREPARE (default):
    For each top-N candidate bar:
      1. Find best Boltz concordance PDB
      2. Fix 3-char chain ID bug (Boltz writes b11/b17 → shift to 'A')
      3. Write clean PDB to outputs/proteinmpnn/pdbs/{bar_id}.pdb
      4. Build fixed_positions.jsonl (lyric AA positions fixed; BJOZXU positions designable)
    Prints Colab commands to run ProteinMPNN.

  FILTER (--filter, run after Colab returns results):
    For each designed sequence in outputs/proteinmpnn/designed/:
      1. Fold with ESMFold
      2. Compute TM-score vs original Boltz backbone
      3. Keep sequences with TM-score >= 0.70
      4. Write outputs/proteinmpnn/filtered_seqs.fasta + filtered_results.csv

Inputs:
  data/phase2_candidates.csv
  outputs/boltz/boltz_summary.csv
  outputs/masks/mask_v2_concordance.json
  data/boltz_id_map.json
  outputs/boltz_outputs/boltz_results_boltz_inputs/predictions/

Outputs (PREPARE):
  outputs/proteinmpnn/pdbs/{bar_id}.pdb       <- chain-fixed PDB for MPNN
  data/phase2_fixed_positions.jsonl           <- fixed positions per bar
  outputs/proteinmpnn/colab_commands.txt      <- paste into Colab

Outputs (FILTER):
  outputs/proteinmpnn/filtered_results.csv
  outputs/proteinmpnn/filtered_seqs.fasta

Usage:
  python analysis/09_proteinmpnn_design.py                   # prepare
  python analysis/09_proteinmpnn_design.py --top-n 12        # top 12 only
  python analysis/09_proteinmpnn_design.py --filter          # self-consistency filter
  python analysis/09_proteinmpnn_design.py --filter --tm-min 0.60
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

CANDIDATES_CSV    = Path("data/phase2_candidates.csv")
BOLTZ_SUMMARY     = Path("outputs/boltz/boltz_summary.csv")
MASK_JSON         = Path("outputs/masks/mask_v2_concordance.json")
BOLTZ_ID_MAP      = Path("data/boltz_id_map.json")
BOLTZ_PREDS       = Path("outputs/boltz_outputs/boltz_results_boltz_inputs/predictions")

MPNN_OUT          = Path("outputs/proteinmpnn")
PDB_OUT           = MPNN_OUT / "pdbs"
DESIGNED_DIR      = MPNN_OUT / "designed"          # ProteinMPNN writes here on Colab
FIXED_POS_JSONL   = Path("data/phase2_fixed_positions.jsonl")
COLAB_CMD_FILE    = MPNN_OUT / "colab_commands.txt"

FILTERED_CSV      = MPNN_OUT / "filtered_results.csv"
FILTERED_FASTA    = MPNN_OUT / "filtered_seqs.fasta"

ESM_URL           = "https://api.esmatlas.com/foldSequence/v1/pdb/"
REQUEST_DELAY     = 1.5
DEFAULT_TOP_N     = 12
DEFAULT_TM_MIN    = 0.70
NUM_SEQS_PER_BAR  = 50
SAMPLING_TEMP     = 0.1    # lower = closer to most probable AA at each position


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def fix_boltz_chain(pdb_text: str) -> str:
    """Replace 3-char Boltz chain IDs (b11, b17, etc.) with single 'A'."""
    out = []
    for line in pdb_text.splitlines(keepends=True):
        if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 32:
            line = line[:21] + "A" + line[24:]
        out.append(line)
    return "".join(out)


def bar_to_boltz_id(bar_id: str, id_map: dict) -> str | None:
    """Return the Boltz dir ID (b0, b1, …) for a given bar_id."""
    rev = {v: k for k, v in id_map.items()}
    return rev.get(bar_id)


def get_best_pdb_path(bar_id: str, id_map: dict, boltz_summary: pd.DataFrame) -> Path | None:
    boltz_id = bar_to_boltz_id(bar_id, id_map)
    if boltz_id is None:
        return None
    row = boltz_summary[boltz_summary["bar_id"] == bar_id]
    best_model = int(row["best_model"].iloc[0]) if not row.empty else 0
    pdb = BOLTZ_PREDS / boltz_id / f"{boltz_id}_model_{best_model}.pdb"
    return pdb if pdb.exists() else None


# ---------------------------------------------------------------------------
# PREPARE mode
# ---------------------------------------------------------------------------

def prepare(top_n: int) -> None:
    PDB_OUT.mkdir(parents=True, exist_ok=True)
    MPNN_OUT.mkdir(parents=True, exist_ok=True)

    candidates  = pd.read_csv(CANDIDATES_CSV).head(top_n)
    boltz_sum   = pd.read_csv(BOLTZ_SUMMARY)
    id_map      = json.load(open(BOLTZ_ID_MAP))
    mask_list   = json.load(open(MASK_JSON))
    mask_by_bar = {m["bar_id"]: m for m in mask_list}

    print(f"Preparing ProteinMPNN inputs for top {top_n} candidates\n")

    fixed_pos_records = []
    prepared = []

    for _, row in candidates.iterrows():
        bar_id = row["bar_id"]
        song   = str(row.get("genius_song_title", ""))[:30]

        # --- PDB ---
        src_pdb = get_best_pdb_path(bar_id, id_map, boltz_sum)
        if src_pdb is None:
            print(f"  [SKIP] {bar_id} — PDB not found")
            continue

        pdb_text = src_pdb.read_text()
        fixed_pdb = fix_boltz_chain(pdb_text)
        out_pdb = PDB_OUT / f"{bar_id}.pdb"
        out_pdb.write_text(fixed_pdb)

        # --- Fixed positions (1-indexed for ProteinMPNN) ---
        m = mask_by_bar.get(bar_id)
        if m is None:
            print(f"  [WARN] {bar_id} — no mask entry, all positions designable")
            bjozxu_0idx = set()
            seq_len = row["fasta_seq_len"]
        else:
            bjozxu_0idx = set(m["mutation_mask"])   # 0-indexed BJOZXU positions
            seq_len = m["sequence_length"]

        # Fixed = all positions NOT in bjozxu set, converted to 1-indexed
        fixed_1idx = [i + 1 for i in range(seq_len) if i not in bjozxu_0idx]
        designable_1idx = [i + 1 for i in bjozxu_0idx]

        fixed_pos_records.append({
            "bar_id":          bar_id,
            "pdb_name":        bar_id,
            "chain":           "A",
            "fixed_positions": fixed_1idx,
            "designable_positions": designable_1idx,
            "n_fixed":         len(fixed_1idx),
            "n_designable":    len(designable_1idx),
            "seq_len":         seq_len,
        })
        prepared.append(bar_id)

        pct = 100 * len(designable_1idx) / seq_len
        print(f"  {bar_id:8s}  {song:30s}  "
              f"len={seq_len}  fixed={len(fixed_1idx)}  "
              f"designable={len(designable_1idx)} ({pct:.0f}%)")

    # Write fixed_positions.jsonl — one JSON object per line, ProteinMPNN format
    # ProteinMPNN expects: {"pdb_name": {"chain_id": [pos, ...]}}
    with open(FIXED_POS_JSONL, "w") as f:
        for rec in fixed_pos_records:
            entry = {rec["pdb_name"]: {rec["chain"]: rec["fixed_positions"]}}
            f.write(json.dumps(entry) + "\n")

    print(f"\nWrote {len(prepared)} PDB files → {PDB_OUT}/")
    print(f"Wrote fixed positions → {FIXED_POS_JSONL}")

    # --- Colab commands ---
    pdb_names = " ".join(prepared)
    colab = f"""
# ============================================================
# ProteinMPNN — Colab commands
# Run after uploading outputs/proteinmpnn/pdbs/ and
# data/phase2_fixed_positions.jsonl to Google Drive
# ============================================================

# 1. Install
!pip install protein-mpnn -q

# 2. Upload PDB folder and fixed positions to Drive, then mount
from google.colab import drive
drive.mount('/content/drive')

# 3. Run ProteinMPNN
import os
PDB_DIR       = "/content/drive/MyDrive/proteinmpnn_pdbs"
FIXED_POS     = "/content/drive/MyDrive/phase2_fixed_positions.jsonl"
OUT_DIR       = "/content/drive/MyDrive/proteinmpnn_designed"

os.makedirs(OUT_DIR, exist_ok=True)

!python /usr/local/lib/python3.*/dist-packages/protein_mpnn/protein_mpnn_run.py \\
    --pdb_path {{PDB_DIR}} \\
    --out_folder {{OUT_DIR}} \\
    --num_seq_per_target {NUM_SEQS_PER_BAR} \\
    --sampling_temp "{SAMPLING_TEMP}" \\
    --fixed_positions_jsonl {{FIXED_POS}} \\
    --batch_size 1

# 4. Zip and download results
import shutil
shutil.make_archive("/content/proteinmpnn_designed", "zip", OUT_DIR)
from google.colab import files
files.download("/content/proteinmpnn_designed.zip")

# ============================================================
# Bars submitted: {", ".join(prepared)}
# Sequences per bar: {NUM_SEQS_PER_BAR}
# Sampling temperature: {SAMPLING_TEMP}
# Fixed positions: lyric AA positions (non-BJOZXU)
# Designable: BJOZXU-derived positions only
# ============================================================
"""
    COLAB_CMD_FILE.write_text(colab)
    print(f"Wrote Colab commands → {COLAB_CMD_FILE}")
    print(f"\nNext: upload {PDB_OUT}/ and {FIXED_POS_JSONL} to Google Drive, run Colab commands.")
    print(f"Then: python analysis/09_proteinmpnn_design.py --filter")


# ---------------------------------------------------------------------------
# FILTER mode — self-consistency with ESMFold
# ---------------------------------------------------------------------------

def fold_with_esm(seq: str) -> tuple[float | None, str | None]:
    """Returns (mean_plddt, pdb_text) or (None, None) on failure."""
    try:
        resp = requests.post(ESM_URL, data=seq, timeout=60,
                             headers={"Content-Type": "application/x-www-form-urlencoded"})
        if resp.status_code != 200:
            return None, None
        pdb_text = resp.text
        b_factors = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and len(line) >= 62:
                try:
                    b_factors.append(float(line[60:66]))
                except ValueError:
                    pass
        if not b_factors:
            return None, None
        return np.mean(b_factors) / 100.0, pdb_text
    except Exception as e:
        print(f"    ESMFold error: {e}")
        return None, None


def tm_score_vs_backbone(designed_pdb: str, ref_pdb_path: Path) -> float | None:
    """
    Compute TM-score between designed ESMFold structure and original Boltz backbone.
    Uses US-align via subprocess if available, otherwise returns None with a warning.
    """
    try:
        import subprocess, tempfile, os
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as tmp:
            tmp.write(designed_pdb)
            tmp_path = tmp.name
        result = subprocess.run(
            ["USalign", tmp_path, str(ref_pdb_path), "-outfmt", "2"],
            capture_output=True, text=True, timeout=30
        )
        os.unlink(tmp_path)
        for line in result.stdout.splitlines():
            if line.startswith("TM-score=") or "TM-score" in line:
                parts = line.split()
                for p in parts:
                    try:
                        val = float(p)
                        if 0 < val <= 1:
                            return val
                    except ValueError:
                        pass
        return None
    except FileNotFoundError:
        return None   # USalign not installed — skip TM-score filter
    except Exception:
        return None


def filter_designs(tm_min: float) -> None:
    """Self-consistency filter on ProteinMPNN designed sequences."""
    if not DESIGNED_DIR.exists():
        print(f"[ERROR] {DESIGNED_DIR} not found. Download ProteinMPNN output from Colab first.")
        return

    candidates  = pd.read_csv(CANDIDATES_CSV)
    boltz_sum   = pd.read_csv(BOLTZ_SUMMARY)
    id_map      = json.load(open(BOLTZ_ID_MAP))

    # Find all designed FASTA files
    fasta_files = list(DESIGNED_DIR.rglob("*.fa")) + list(DESIGNED_DIR.rglob("*.fasta"))
    print(f"Found {len(fasta_files)} FASTA files in {DESIGNED_DIR}")

    rows = []
    passing_seqs = []

    for fa_path in sorted(fasta_files):
        bar_id = fa_path.stem  # filename = bar_id
        ref_pdb = get_best_pdb_path(bar_id, id_map, boltz_sum)
        if ref_pdb is None:
            print(f"  [SKIP] {bar_id} — no reference PDB")
            continue

        # Parse sequences from FASTA
        seqs = []
        current_header = None
        current_seq = []
        for line in fa_path.read_text().splitlines():
            if line.startswith(">"):
                if current_header:
                    seqs.append((current_header, "".join(current_seq)))
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_header:
            seqs.append((current_header, "".join(current_seq)))

        print(f"\n{bar_id}: {len(seqs)} designed sequences")

        for i, (header, seq) in enumerate(seqs):
            print(f"  seq {i+1}/{len(seqs)}", end="  ")
            time.sleep(REQUEST_DELAY)
            plddt, pdb_text = fold_with_esm(seq)

            if plddt is None:
                print("ESMFold FAILED")
                tm = None
            else:
                print(f"pLDDT={plddt:.3f}", end="  ")
                tm = tm_score_vs_backbone(pdb_text, ref_pdb) if pdb_text else None
                if tm is not None:
                    print(f"TM={tm:.3f}", end="  ")
                else:
                    print("TM=N/A (USalign not installed)", end="  ")

            passes = (tm is not None and tm >= tm_min) or (tm is None and plddt is not None and plddt >= 0.40)
            print("PASS" if passes else "FAIL")

            rows.append({
                "bar_id":      bar_id,
                "design_idx":  i,
                "header":      header,
                "sequence":    seq,
                "esm_plddt":   plddt,
                "tm_score":    tm,
                "tm_min":      tm_min,
                "passes":      passes,
            })
            if passes:
                passing_seqs.append((f">{bar_id}_design_{i:02d} | {header}", seq))

    results = pd.DataFrame(rows)
    results.to_csv(FILTERED_CSV, index=False)
    print(f"\nSaved {len(results)} rows → {FILTERED_CSV}")
    print(f"Passing designs: {results['passes'].sum()} / {len(results)}")

    # Write passing sequences to FASTA
    with open(FILTERED_FASTA, "w") as f:
        for header, seq in passing_seqs:
            f.write(f"{header}\n{seq}\n")
    print(f"Saved {len(passing_seqs)} passing sequences → {FILTERED_FASTA}")

    # Summary per bar
    print("\n--- Passing designs per bar ---")
    summary = results.groupby("bar_id")["passes"].sum()
    for bar_id, n in summary.items():
        total = (results["bar_id"] == bar_id).sum()
        print(f"  {bar_id:8s}  {int(n):2d} / {total} pass")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n",   type=int,   default=DEFAULT_TOP_N,
                        help=f"Number of top candidates to prepare (default {DEFAULT_TOP_N})")
    parser.add_argument("--filter",  action="store_true",
                        help="Run self-consistency filter on ProteinMPNN results")
    parser.add_argument("--tm-min",  type=float, default=DEFAULT_TM_MIN,
                        help=f"TM-score threshold for self-consistency filter (default {DEFAULT_TM_MIN})")
    args = parser.parse_args()

    if args.filter:
        filter_designs(args.tm_min)
    else:
        prepare(args.top_n)
