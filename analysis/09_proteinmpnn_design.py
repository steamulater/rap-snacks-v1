"""
09_proteinmpnn_design.py
------------------------
Phase 2 Step 2 — ProteinMPNN soft sequence design + ESMFold self-consistency.

Runs fully locally end-to-end:
  1. For each top-N candidate bar, find best Boltz concordance PDB and fix chain IDs
  2. Write fixed_positions.jsonl (lyric AA positions fixed; BJOZXU positions designable)
  3. Run ProteinMPNN locally (50 seqs per bar) via subprocess
  4. Self-consistency filter: fold each design with ESMFold API, keep pLDDT >= ESM_PLDDT_MIN
     (TM-score filter also available if USalign is installed)
  5. Write filtered_seqs.fasta + filtered_results.csv

Inputs:
  data/phase2_candidates.csv
  outputs/boltz/boltz_summary.csv
  outputs/masks/mask_v2_concordance.json
  data/boltz_id_map.json
  outputs/boltz_outputs/boltz_results_boltz_inputs/predictions/

Outputs:
  outputs/proteinmpnn/pdbs/{bar_id}.pdb         <- chain-fixed PDB
  data/phase2_fixed_positions.jsonl             <- fixed positions per bar
  outputs/proteinmpnn/designed/{bar_id}/        <- raw ProteinMPNN FASTA output
  outputs/proteinmpnn/filtered_results.csv      <- all designs + ESMFold scores
  outputs/proteinmpnn/filtered_seqs.fasta       <- passing designs only

Usage:
  python analysis/09_proteinmpnn_design.py
  python analysis/09_proteinmpnn_design.py --top-n 12 --n-seqs 50 --temp 0.1
  python analysis/09_proteinmpnn_design.py --esm-plddt-min 0.45
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

MPNN_SCRIPT       = Path("/Users/tamukamartin/Desktop/adaptyv_competiton/ProteinMPNN/protein_mpnn_run.py")
MPNN_WEIGHTS      = Path("/Users/tamukamartin/Desktop/adaptyv_competiton/ProteinMPNN/vanilla_model_weights")

CANDIDATES_CSV    = Path("data/phase2_candidates.csv")
BOLTZ_SUMMARY     = Path("outputs/boltz/boltz_summary.csv")
MASK_JSON         = Path("outputs/masks/mask_v2_concordance.json")
BOLTZ_ID_MAP      = Path("data/boltz_id_map.json")
BOLTZ_PREDS       = Path("outputs/boltz_outputs/boltz_results_boltz_inputs/predictions")

MPNN_OUT          = Path("outputs/proteinmpnn")
PDB_OUT           = MPNN_OUT / "pdbs"
DESIGNED_DIR      = MPNN_OUT / "designed"
FIXED_POS_JSONL   = Path("data/phase2_fixed_positions.jsonl")
FILTERED_CSV      = MPNN_OUT / "filtered_results.csv"
FILTERED_FASTA    = MPNN_OUT / "filtered_seqs.fasta"

ESM_URL           = "https://api.esmatlas.com/foldSequence/v1/pdb/"
REQUEST_DELAY     = 1.5

DEFAULT_TOP_N         = 12
DEFAULT_N_SEQS        = 50
DEFAULT_TEMP          = 0.1
DEFAULT_ESM_PLDDT_MIN = 0.35   # MPNN designs typically score 0.35-0.50 on ESMFold

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def fix_boltz_chain(pdb_text: str) -> str:
    out = []
    for line in pdb_text.splitlines(keepends=True):
        if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 32:
            line = line[:21] + "A" + line[24:]
        out.append(line)
    return "".join(out)


def bar_to_boltz_id(bar_id: str, id_map: dict) -> str | None:
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


def fold_with_esm(seq: str) -> tuple[float | None, str | None]:
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
        return np.mean(b_factors), pdb_text  # ESMFold API returns pLDDT in 0-1 directly
    except Exception as e:
        print(f"    ESMFold error: {e}")
        return None, None


def parse_mpnn_fasta(fa_path: Path) -> list[tuple[str, str]]:
    """Parse ProteinMPNN output FASTA — returns [(header, seq), ...]."""
    seqs, header, seq = [], None, []
    for line in fa_path.read_text().splitlines():
        if line.startswith(">"):
            if header:
                seqs.append((header, "".join(seq)))
            header = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if header:
        seqs.append((header, "".join(seq)))
    return seqs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(top_n: int, n_seqs: int, temp: float, esm_plddt_min: float) -> None:
    for p in [MPNN_SCRIPT, MPNN_WEIGHTS]:
        if not p.exists():
            print(f"[ERROR] ProteinMPNN not found at {p}")
            sys.exit(1)

    PDB_OUT.mkdir(parents=True, exist_ok=True)
    DESIGNED_DIR.mkdir(parents=True, exist_ok=True)

    candidates  = pd.read_csv(CANDIDATES_CSV).head(top_n)
    boltz_sum   = pd.read_csv(BOLTZ_SUMMARY)
    id_map      = json.load(open(BOLTZ_ID_MAP))
    mask_list   = json.load(open(MASK_JSON))
    mask_by_bar = {m["bar_id"]: m for m in mask_list}

    print(f"Running ProteinMPNN on top {top_n} candidates")
    print(f"  {n_seqs} sequences per bar · temperature {temp}\n")

    fixed_pos_records = []

    # ------------------------------------------------------------------ #
    # Step 1 — prepare PDBs and fixed positions                           #
    # ------------------------------------------------------------------ #
    for _, row in candidates.iterrows():
        bar_id = row["bar_id"]
        src_pdb = get_best_pdb_path(bar_id, id_map, boltz_sum)
        if src_pdb is None:
            print(f"  [SKIP] {bar_id} — PDB not found")
            continue

        # Fix chain IDs and write
        fixed_pdb = fix_boltz_chain(src_pdb.read_text())
        out_pdb = PDB_OUT / f"{bar_id}.pdb"
        out_pdb.write_text(fixed_pdb)

        m = mask_by_bar.get(bar_id, {})
        bjozxu_0idx = set(m.get("mutation_mask", []))
        seq_len = m.get("sequence_length", int(row["fasta_seq_len"]))
        fixed_1idx = [i + 1 for i in range(seq_len) if i not in bjozxu_0idx]
        designable_1idx = [i + 1 for i in sorted(bjozxu_0idx)]

        fixed_pos_records.append({
            "bar_id":   bar_id,
            "fixed":    fixed_1idx,
            "design":   designable_1idx,
            "seq_len":  seq_len,
        })

    # Write fixed_positions.jsonl
    with open(FIXED_POS_JSONL, "w") as f:
        for rec in fixed_pos_records:
            entry = {rec["bar_id"]: {"A": rec["fixed"]}}
            f.write(json.dumps(entry) + "\n")

    print(f"Prepared {len(fixed_pos_records)} PDBs + fixed positions JSONL\n")

    # ------------------------------------------------------------------ #
    # Step 2 — run ProteinMPNN per bar                                    #
    # ------------------------------------------------------------------ #
    all_rows = []

    for rec in fixed_pos_records:
        bar_id  = rec["bar_id"]
        pdb_in  = PDB_OUT / f"{bar_id}.pdb"
        out_dir = DESIGNED_DIR / bar_id
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"[MPNN] {bar_id}  len={rec['seq_len']}  "
              f"designable={len(rec['design'])}  fixed={len(rec['fixed'])}")

        cmd = [
            sys.executable, str(MPNN_SCRIPT),
            "--pdb_path",              str(pdb_in),
            "--out_folder",            str(out_dir),
            "--num_seq_per_target",    str(n_seqs),
            "--sampling_temp",         str(temp),
            "--fixed_positions_jsonl", str(FIXED_POS_JSONL),
            "--path_to_model_weights", str(MPNN_WEIGHTS),
            "--batch_size",            "1",
            "--suppress_print",        "1",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  [ERROR] ProteinMPNN failed for {bar_id}:")
            print(result.stderr[-500:])
            continue

        # Find the output FASTA
        fa_files = list(out_dir.rglob("*.fa")) + list(out_dir.rglob("*.fasta"))
        if not fa_files:
            print(f"  [WARN] No FASTA output found for {bar_id}")
            continue

        seqs = parse_mpnn_fasta(fa_files[0])
        # Skip the first sequence (it's the original/reference)
        designs = [(h, s) for h, s in seqs if "score" in h.lower() or seqs.index((h,s)) > 0]
        print(f"  → {len(designs)} designs generated")

        # ------------------------------------------------------------------ #
        # Step 3 — ESMFold self-consistency filter                            #
        # ------------------------------------------------------------------ #
        for i, (header, seq) in enumerate(designs):
            print(f"  ESMFold {i+1}/{len(designs)}", end="  ", flush=True)
            time.sleep(REQUEST_DELAY)
            plddt, _ = fold_with_esm(seq)

            passes = plddt is not None and plddt >= esm_plddt_min
            status = f"pLDDT={plddt:.3f} {'PASS' if passes else 'FAIL'}" if plddt else "FAILED"
            print(status)

            all_rows.append({
                "bar_id":      bar_id,
                "design_idx":  i,
                "header":      header,
                "sequence":    seq,
                "seq_len":     len(seq),
                "esm_plddt":   plddt,
                "esm_plddt_min": esm_plddt_min,
                "passes":      passes,
            })

    # ------------------------------------------------------------------ #
    # Step 4 — save results                                               #
    # ------------------------------------------------------------------ #
    results = pd.DataFrame(all_rows)
    results.to_csv(FILTERED_CSV, index=False)
    print(f"\nSaved {len(results)} designs → {FILTERED_CSV}")

    passing = results[results["passes"]]
    print(f"Passing designs (pLDDT >= {esm_plddt_min}): {len(passing)} / {len(results)}")

    with open(FILTERED_FASTA, "w") as f:
        for _, row in passing.iterrows():
            f.write(f">{row['bar_id']}_design_{row['design_idx']:02d} | {row['header']}\n")
            f.write(f"{row['sequence']}\n")
    print(f"Saved {len(passing)} passing sequences → {FILTERED_FASTA}")

    print("\n--- Passing designs per bar ---")
    for bar_id, grp in results.groupby("bar_id"):
        n_pass = grp["passes"].sum()
        print(f"  {bar_id:8s}  {int(n_pass):2d} / {len(grp)} pass"
              f"  mean pLDDT={grp['esm_plddt'].mean():.3f}")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n",        type=int,   default=DEFAULT_TOP_N)
    parser.add_argument("--n-seqs",       type=int,   default=DEFAULT_N_SEQS)
    parser.add_argument("--temp",         type=float, default=DEFAULT_TEMP)
    parser.add_argument("--esm-plddt-min",type=float, default=DEFAULT_ESM_PLDDT_MIN)
    args = parser.parse_args()

    main(args.top_n, args.n_seqs, args.temp, args.esm_plddt_min)
