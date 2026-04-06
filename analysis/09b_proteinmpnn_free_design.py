"""
09b_proteinmpnn_free_design.py
------------------------------
Phase 2 Step 2b — ProteinMPNN FREE sequence design + ESMFold self-consistency.

No fixed positions — MPNN optimizes the full sequence for the Boltz concordance
backbone. All positions are designable.

Purpose:
  1. Upper bound on foldability — if free_design also scores low pLDDT, the
     backbone itself is the problem (not the BJOZXU masking constraint).
  2. Wet-lab comparison — free_design vs masked_BJOZXU directly tests whether
     preserving the lyric AAs costs foldability.

Contrast with 09_proteinmpnn_design.py (masked_BJOZXU):
  - masked_BJOZXU: lyric AA positions fixed; only BJOZXU positions (~5-18) redesigned
  - free_design:   all positions designable; MPNN picks the most foldable sequence

Inputs:
  data/phase2_candidates.csv
  outputs/boltz/boltz_summary.csv
  data/boltz_id_map.json
  outputs/boltz_outputs/boltz_results_boltz_inputs/predictions/
  outputs/proteinmpnn/pdbs/{bar_id}.pdb   <- reuse chain-fixed PDBs from run 1

Outputs:
  outputs/proteinmpnn_free/designed/{bar_id}/   <- raw ProteinMPNN FASTA output
  outputs/proteinmpnn_free/filtered_results.csv <- all designs + ESMFold scores
  outputs/proteinmpnn_free/filtered_seqs.fasta  <- passing designs only

Usage:
  python analysis/09b_proteinmpnn_free_design.py
  python analysis/09b_proteinmpnn_free_design.py --top-n 12 --n-seqs 50 --temp 0.1
  python analysis/09b_proteinmpnn_free_design.py --esm-plddt-min 0.35
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

MPNN_SCRIPT   = Path("/Users/tamukamartin/Desktop/adaptyv_competiton/ProteinMPNN/protein_mpnn_run.py")
MPNN_WEIGHTS  = Path("/Users/tamukamartin/Desktop/adaptyv_competiton/ProteinMPNN/vanilla_model_weights")

CANDIDATES_CSV = Path("data/phase2_candidates.csv")
BOLTZ_SUMMARY  = Path("outputs/boltz/boltz_summary.csv")
BOLTZ_ID_MAP   = Path("data/boltz_id_map.json")
BOLTZ_PREDS    = Path("outputs/boltz_outputs/boltz_results_boltz_inputs/predictions")

# Reuse chain-fixed PDBs written by run 1 (09_proteinmpnn_design.py)
PDB_DIR        = Path("outputs/proteinmpnn/pdbs")

MPNN_OUT       = Path("outputs/proteinmpnn_free")
DESIGNED_DIR   = MPNN_OUT / "designed"
FILTERED_CSV   = MPNN_OUT / "filtered_results.csv"
FILTERED_FASTA = MPNN_OUT / "filtered_seqs.fasta"

ESM_URL        = "https://api.esmatlas.com/foldSequence/v1/pdb/"
REQUEST_DELAY  = 1.5

DEFAULT_TOP_N         = 12
DEFAULT_N_SEQS        = 50
DEFAULT_TEMP          = 0.1
DEFAULT_ESM_PLDDT_MIN = 0.35

# ---------------------------------------------------------------------------
# Helpers (identical to run 1)
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
        return np.mean(b_factors), pdb_text
    except Exception as e:
        print(f"    ESMFold error: {e}")
        return None, None


def parse_mpnn_fasta(fa_path: Path) -> list[tuple[str, str]]:
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

    MPNN_OUT.mkdir(parents=True, exist_ok=True)
    DESIGNED_DIR.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(CANDIDATES_CSV).head(top_n)
    boltz_sum  = pd.read_csv(BOLTZ_SUMMARY)
    id_map     = json.load(open(BOLTZ_ID_MAP))

    print(f"ProteinMPNN FREE DESIGN on top {top_n} candidates")
    print(f"  {n_seqs} sequences per bar · temperature {temp} · NO fixed positions\n")

    # Load any existing results for resume
    existing   = pd.read_csv(FILTERED_CSV) if FILTERED_CSV.exists() else pd.DataFrame()
    scored_bars = set()
    if not existing.empty:
        for bar_id, grp in existing.groupby("bar_id"):
            if grp["esm_plddt"].notna().all():
                scored_bars.add(bar_id)
        if scored_bars:
            print(f"Resuming — skipping {len(scored_bars)} already-scored bars: {sorted(scored_bars)}\n")

    all_rows = list(existing[existing["bar_id"].isin(scored_bars)].to_dict("records")) if not existing.empty else []

    for _, row in candidates.iterrows():
        bar_id = row["bar_id"]

        if bar_id in scored_bars:
            print(f"  [SKIP] {bar_id} — already scored")
            continue

        # Prefer pre-fixed PDB from run 1; fall back to re-fixing from Boltz
        pdb_in = PDB_DIR / f"{bar_id}.pdb"
        if not pdb_in.exists():
            src_pdb = get_best_pdb_path(bar_id, id_map, boltz_sum)
            if src_pdb is None:
                print(f"  [SKIP] {bar_id} — PDB not found")
                continue
            pdb_in.parent.mkdir(parents=True, exist_ok=True)
            pdb_in.write_text(fix_boltz_chain(src_pdb.read_text()))

        out_dir = DESIGNED_DIR / bar_id
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"[MPNN free] {bar_id}")

        # No --fixed_positions_jsonl — all positions designable
        cmd = [
            sys.executable, str(MPNN_SCRIPT),
            "--pdb_path",           str(pdb_in),
            "--out_folder",         str(out_dir),
            "--num_seq_per_target", str(n_seqs),
            "--sampling_temp",      str(temp),
            "--path_to_model_weights", str(MPNN_WEIGHTS),
            "--batch_size",         "1",
            "--suppress_print",     "1",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  [ERROR] ProteinMPNN failed for {bar_id}:")
            print(result.stderr[-500:])
            continue

        fa_files = list(out_dir.rglob("*.fa")) + list(out_dir.rglob("*.fasta"))
        if not fa_files:
            print(f"  [WARN] No FASTA output found for {bar_id}")
            continue

        seqs = parse_mpnn_fasta(fa_files[0])
        designs = [(h, s) for h, s in seqs if "score" in h.lower() or seqs.index((h, s)) > 0]
        print(f"  → {len(designs)} designs generated")

        for i, (header, seq) in enumerate(designs):
            print(f"  ESMFold {i+1}/{len(designs)}", end="  ", flush=True)
            time.sleep(REQUEST_DELAY)
            plddt, _ = fold_with_esm(seq)

            passes = plddt is not None and plddt >= esm_plddt_min
            status = f"pLDDT={plddt:.3f} {'PASS' if passes else 'FAIL'}" if plddt else "FAILED"
            print(status)

            all_rows.append({
                "bar_id":        bar_id,
                "design_idx":    i,
                "design_type":   "free_design",
                "header":        header,
                "sequence":      seq,
                "seq_len":       len(seq),
                "esm_plddt":     plddt,
                "esm_plddt_min": esm_plddt_min,
                "passes":        passes,
            })

    results = pd.DataFrame(all_rows)
    results.to_csv(FILTERED_CSV, index=False)
    print(f"\nSaved {len(results)} designs → {FILTERED_CSV}")

    passing = results[results["passes"]]
    print(f"Passing designs (pLDDT >= {esm_plddt_min}): {len(passing)} / {len(results)}")

    with open(FILTERED_FASTA, "w") as f:
        for _, r in passing.iterrows():
            f.write(f">{r['bar_id']}_free_{r['design_idx']:02d} | {r['header']}\n")
            f.write(f"{r['sequence']}\n")
    print(f"Saved {len(passing)} passing sequences → {FILTERED_FASTA}")

    print("\n--- Passing designs per bar ---")
    for bar_id, grp in results.groupby("bar_id"):
        n_pass = grp["passes"].sum()
        print(f"  {bar_id:8s}  {int(n_pass):2d} / {len(grp)} pass"
              f"  mean pLDDT={grp['esm_plddt'].mean():.3f}")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n",         type=int,   default=DEFAULT_TOP_N)
    parser.add_argument("--n-seqs",        type=int,   default=DEFAULT_N_SEQS)
    parser.add_argument("--temp",          type=float, default=DEFAULT_TEMP)
    parser.add_argument("--esm-plddt-min", type=float, default=DEFAULT_ESM_PLDDT_MIN)
    args = parser.parse_args()

    main(args.top_n, args.n_seqs, args.temp, args.esm_plddt_min)
