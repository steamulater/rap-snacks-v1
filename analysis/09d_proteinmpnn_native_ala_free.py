"""
09d_proteinmpnn_native_ala_free.py
-----------------------------------
Phase 2 Step 4 — ProteinMPNN free design on native_ala Boltz backbones.

Input backbone: best Boltz-2 model for each bar's native_ala sequence
  (BJOZXU → Ala; all standard AAs; lyric text readable in sequence)

Design strategy: no fixed positions — MPNN optimises the full sequence for the
native_ala backbone geometry.

Rationale vs free_design (Run 2):
  Run 2  used the concordance Boltz backbone  (raw lyric-derived geometry)
  Run 4  uses the native_ala Boltz backbone   (Ala-substituted, more foldable)

  native_ala backbones have mean Boltz pLDDT 0.543 vs concordance 0.441 — a
  better starting geometry. Designing on a more foldable backbone should yield
  sequences that (a) score higher in Boltz validation and (b) may occupy
  distinct structural space from Run 2.

  Additionally, bars with native_ala-specific structural novelty (bar_3, bar_8,
  bar_9: 3–6 FoldSeek hits; bar_27: 0 hits) will have MPNN exploring sequence
  space anchored to that unusual fold — the designs may be novel *and* foldable.

Inputs (local):
  data/phase2_candidates.csv
  outputs/boltz_validation/boltz_outputs/.../{bar_id}_native_ala/{bar_id}_native_ala_model_0.pdb

Outputs (local — upload to Drive before running boltz_validation_v3.ipynb):
  outputs/proteinmpnn_native_ala_free/pdbs/{bar_id}.pdb       ← chain-fixed backbone PDB
  outputs/proteinmpnn_native_ala_free/designed/{bar_id}/      ← raw MPNN FASTA output
  outputs/proteinmpnn_native_ala_free/filtered_results.csv    ← all designs + ESMFold scores
  outputs/proteinmpnn_native_ala_free/filtered_seqs.fasta     ← passing designs (pLDDT ≥ 0.35)

Usage:
  python analysis/09d_proteinmpnn_native_ala_free.py
  python analysis/09d_proteinmpnn_native_ala_free.py --n-seqs 50 --temp 0.1 --esm-plddt-min 0.35
"""

import argparse, subprocess, sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

MPNN_SCRIPT  = Path("/Users/tamukamartin/Desktop/adaptyv_competiton/ProteinMPNN/protein_mpnn_run.py")
MPNN_WEIGHTS = Path("/Users/tamukamartin/Desktop/adaptyv_competiton/ProteinMPNN/vanilla_model_weights")

ROOT           = Path(__file__).parent.parent
CANDIDATES_CSV = ROOT / "data/phase2_candidates.csv"

# native_ala PDBs live inside boltz_validation outputs
BOLTZ_VAL_ROOT = ROOT / "outputs/boltz_validation/boltz_outputs"

OUT_DIR        = ROOT / "outputs/proteinmpnn_native_ala_free"
PDB_DIR        = OUT_DIR / "pdbs"          # chain-fixed backbone PDBs
DESIGNED_DIR   = OUT_DIR / "designed"      # raw MPNN output
FILTERED_CSV   = OUT_DIR / "filtered_results.csv"
FILTERED_FASTA = OUT_DIR / "filtered_seqs.fasta"

ESM_URL        = "https://api.esmatlas.com/foldSequence/v1/pdb/"
REQUEST_DELAY  = 1.5

DEFAULT_N_SEQS        = 50
DEFAULT_TEMP          = 0.1
DEFAULT_ESM_PLDDT_MIN = 0.35

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_preds_dir(root: Path) -> Path:
    hits = list(root.rglob("predictions"))
    if not hits:
        raise RuntimeError(f"No predictions/ dir under {root}")
    return hits[0]


def find_native_ala_pdb(bar_id: str, preds_dir: Path) -> Path | None:
    """Return model_0 PDB for bar_id_native_ala."""
    seq_dir = preds_dir / f"{bar_id}_native_ala"
    pdb = seq_dir / f"{bar_id}_native_ala_model_0.pdb"
    return pdb if pdb.exists() else None


def fix_boltz_chain(pdb_text: str) -> str:
    """Boltz writes 2-char chain IDs; fix to single-char 'A' for ProteinMPNN."""
    out = []
    for line in pdb_text.splitlines(keepends=True):
        if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 32:
            line = line[:21] + "A" + line[24:]
        out.append(line)
    return "".join(out)


def fold_with_esm(seq: str) -> float | None:
    try:
        resp = requests.post(
            ESM_URL, data=seq, timeout=60,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
        )
        if resp.status_code != 200:
            return None
        b_factors = []
        for line in resp.text.splitlines():
            if line.startswith("ATOM") and len(line) >= 62:
                try:
                    b_factors.append(float(line[60:66]))
                except ValueError:
                    pass
        return float(np.mean(b_factors)) if b_factors else None
    except Exception as e:
        print(f"    ESMFold error: {e}")
        return None


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

def main(n_seqs: int, temp: float, esm_plddt_min: float) -> None:
    for p in [MPNN_SCRIPT, MPNN_WEIGHTS]:
        if not p.exists():
            print(f"[ERROR] ProteinMPNN not found at {p}")
            sys.exit(1)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    DESIGNED_DIR.mkdir(parents=True, exist_ok=True)

    preds_dir  = find_preds_dir(BOLTZ_VAL_ROOT)
    candidates = pd.read_csv(CANDIDATES_CSV)
    bar_ids    = list(candidates["bar_id"])

    print(f"ProteinMPNN FREE DESIGN on native_ala backbones")
    print(f"  {n_seqs} seqs/bar · temp={temp} · no fixed positions · pLDDT≥{esm_plddt_min}\n")

    # Resume support
    existing = pd.read_csv(FILTERED_CSV) if FILTERED_CSV.exists() else pd.DataFrame()
    scored_bars = set()
    if not existing.empty:
        for bar_id, grp in existing.groupby("bar_id"):
            if grp["esm_plddt"].notna().all():
                scored_bars.add(bar_id)
        if scored_bars:
            print(f"Resuming — skipping {len(scored_bars)} bars: {sorted(scored_bars)}\n")

    all_rows = list(existing[existing["bar_id"].isin(scored_bars)].to_dict("records")) if not existing.empty else []

    for bar_id in bar_ids:
        if bar_id in scored_bars:
            print(f"  [SKIP] {bar_id}")
            continue

        src_pdb = find_native_ala_pdb(bar_id, preds_dir)
        if src_pdb is None:
            print(f"  [SKIP] {bar_id} — native_ala model_0 PDB not found")
            continue

        # Write chain-fixed backbone PDB
        pdb_in = PDB_DIR / f"{bar_id}.pdb"
        pdb_in.write_text(fix_boltz_chain(src_pdb.read_text()))

        out_dir = DESIGNED_DIR / bar_id
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"[MPNN native_ala_free] {bar_id}  backbone={src_pdb.name}")

        cmd = [
            sys.executable, str(MPNN_SCRIPT),
            "--pdb_path",              str(pdb_in),
            "--out_folder",            str(out_dir),
            "--num_seq_per_target",    str(n_seqs),
            "--sampling_temp",         str(temp),
            "--path_to_model_weights", str(MPNN_WEIGHTS),
            "--batch_size",            "1",
            "--suppress_print",        "1",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  [ERROR] MPNN failed:\n{result.stderr[-400:]}")
            continue

        fa_files = list(out_dir.rglob("*.fa")) + list(out_dir.rglob("*.fasta"))
        if not fa_files:
            print(f"  [WARN] No FASTA output for {bar_id}")
            continue

        seqs = parse_mpnn_fasta(fa_files[0])
        # Skip first entry (the input sequence, score=0.0000)
        designs = [(h, s) for h, s in seqs if seqs.index((h, s)) > 0]
        print(f"  → {len(designs)} designs")

        for i, (header, seq) in enumerate(designs):
            print(f"  ESMFold {i+1}/{len(designs)}", end="  ", flush=True)
            time.sleep(REQUEST_DELAY)
            plddt = fold_with_esm(seq)
            passes = plddt is not None and plddt >= esm_plddt_min
            print(f"pLDDT={'N/A' if plddt is None else f'{plddt:.3f}'} {'PASS' if passes else 'FAIL'}")

            all_rows.append({
                "bar_id":      bar_id,
                "design_idx":  i,
                "bucket":      "native_ala_free",
                "header":      header,
                "sequence":    seq,
                "seq_len":     len(seq),
                "esm_plddt":   plddt,
                "passes":      passes,
            })

    results = pd.DataFrame(all_rows)
    results.to_csv(FILTERED_CSV, index=False)
    print(f"\nSaved {len(results)} designs → {FILTERED_CSV}")

    passing = results[results["passes"] == True]
    print(f"Passing (pLDDT ≥ {esm_plddt_min}): {len(passing)} / {len(results)}")

    with open(FILTERED_FASTA, "w") as f:
        for _, r in passing.iterrows():
            f.write(f">{r['bar_id']}_naf_{int(r['design_idx']):03d}\n{r['sequence']}\n")
    print(f"Saved {len(passing)} sequences → {FILTERED_FASTA}")

    print("\n--- Per-bar summary ---")
    for bar_id, grp in results.groupby("bar_id"):
        n_pass = int(grp["passes"].sum())
        mean_p = grp["esm_plddt"].mean()
        print(f"  {bar_id:8s}  {n_pass:2d}/{len(grp)} pass  mean pLDDT={mean_p:.3f}")

    print(f"\nNext step: upload {OUT_DIR.name}/ to Drive → rap_snacks/inputs/outputs/")
    print("Then run notebooks/boltz_validation_v3.ipynb on Colab.")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--n-seqs",        type=int,   default=DEFAULT_N_SEQS)
    p.add_argument("--temp",          type=float, default=DEFAULT_TEMP)
    p.add_argument("--esm-plddt-min", type=float, default=DEFAULT_ESM_PLDDT_MIN)
    args = p.parse_args()
    main(args.n_seqs, args.temp, args.esm_plddt_min)
