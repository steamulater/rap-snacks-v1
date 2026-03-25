"""
11_scrambled_control.py
-----------------------
In silico scrambled sequence control for Phase 2.

For each candidate bar, generates N shuffled variants of BOTH the concordance
and native_ala sequences (same AA composition, random order per source),
then folds them with ESMFold.

Three buckets compared per bar:
  1. concordance   — original lyric-derived sequence (pLDDT from enriched CSV)
  2. native_ala    — alanine-substituted baseline (pLDDT from enriched CSV)
  3. scrambled     — N shuffles each of concordance + native_ala (new ESMFold calls)

Original pLDDT values are read from the enriched CSV (no extra API calls needed).
Only the scrambled sequences hit the ESMFold API.

Expected result: scrambled → pLDDT < 0.3, disordered — confirming that
structure is encoded by sequence, not just amino acid composition.

Input:  data/aggregated_lines_v2_enriched.csv
        data/phase2_candidates.csv             (if exists, restricts to shortlist)
Output: outputs/scrambled/scrambled_results.csv
        outputs/scrambled/fig_scrambled_plddt.png

Usage:
    python analysis/11_scrambled_control.py
    python analysis/11_scrambled_control.py --n-shuffles 5 --all-bars
"""

import argparse
import json
import random
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ENRICHED_CSV      = Path("data/aggregated_lines_v2_enriched.csv")
CANDIDATES_CSV    = Path("data/phase2_candidates.csv")
OUT_DIR           = Path("outputs/scrambled")
RESULTS_CSV       = OUT_DIR / "scrambled_results.csv"
FIG_PATH          = OUT_DIR / "fig_scrambled_plddt.png"

ESM_URL           = "https://api.esmatlas.com/foldSequence/v1/pdb/"
SEED              = 42
DEFAULT_N_SHUFFLE = 3      # shuffled variants per bar
REQUEST_DELAY     = 1.5    # seconds between ESMFold calls (rate limit)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def shuffle_sequence(seq: str, rng: random.Random) -> str:
    """Return seq with amino acids shuffled — same composition, random order."""
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def fold_with_esm(seq: str) -> float | None:
    """
    Call ESMFold API, return mean pLDDT (0–1 scale).
    Returns None on failure.
    """
    try:
        resp = requests.post(ESM_URL, data=seq, timeout=60,
                             headers={"Content-Type": "application/x-www-form-urlencoded"})
        if resp.status_code != 200:
            print(f"    ESMFold HTTP {resp.status_code}")
            return None
        pdb_text = resp.text
        b_factors = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and len(line) >= 62:
                try:
                    b_factors.append(float(line[60:66]))
                except ValueError:
                    pass
        if not b_factors:
            return None
        # ESMFold reports pLDDT as 0–100 in B-factor column
        return np.mean(b_factors) / 100.0
    except Exception as e:
        print(f"    ESMFold error: {e}")
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(n_shuffles: int, all_bars: bool) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    df = pd.read_csv(ENRICHED_CSV)

    # Restrict to candidate shortlist if available and --all-bars not set
    if not all_bars and CANDIDATES_CSV.exists():
        candidates = pd.read_csv(CANDIDATES_CSV)
        bar_ids = set(candidates["bar_id"])
        df = df[df["bar_id"].isin(bar_ids)]
        print(f"Restricting to {len(df)} candidate bars from {CANDIDATES_CSV}")
    else:
        print(f"Running on all {len(df)} bars (--all-bars or no candidates CSV)")

    df = df.dropna(subset=["fasta_seq_concordance"]).reset_index(drop=True)
    print(f"  {len(df)} bars with concordance sequence\n")

    rows = []

    for _, row in df.iterrows():
        bar_id = row["bar_id"]
        song   = row.get("genius_song_title", "")

        conc_seq   = row.get("fasta_seq_concordance", "")
        na_seq     = row.get("fasta_seq_native_alanine", "")
        conc_plddt = row.get("esm_plddt_concordance_mean", None)
        na_plddt   = row.get("esm_plddt_native_alanine_mean", None)

        if not conc_seq and not na_seq:
            print(f"  [SKIP] {bar_id} — no sequences")
            continue

        conc_str = f"{conc_plddt:.3f}" if pd.notna(conc_plddt) else "N/A"
        na_str   = f"{na_plddt:.3f}"   if pd.notna(na_plddt)   else "N/A"
        print(f"Bar {bar_id}  concordance pLDDT={conc_str}  native_ala pLDDT={na_str}")

        # Store originals (no API call — values already in enriched CSV)
        rows.append({
            "bar_id":       bar_id,
            "song":         song,
            "source":       "concordance",
            "shuffle_idx":  -1,
            "seq":          conc_seq,
            "plddt":        conc_plddt,
            "seq_len":      len(conc_seq) if conc_seq else None,
        })
        rows.append({
            "bar_id":       bar_id,
            "song":         song,
            "source":       "native_ala",
            "shuffle_idx":  -1,
            "seq":          na_seq,
            "plddt":        na_plddt,
            "seq_len":      len(na_seq) if na_seq else None,
        })

        # Scramble concordance sequence
        for i in range(n_shuffles):
            if not conc_seq:
                continue
            shuffled = shuffle_sequence(conc_seq, rng)
            print(f"  conc shuffle {i+1}/{n_shuffles}  {shuffled[:20]}...")
            time.sleep(REQUEST_DELAY)
            plddt = fold_with_esm(shuffled)
            print(f"    → pLDDT = {plddt:.3f}" if plddt is not None else "    → FAILED")
            rows.append({
                "bar_id":       bar_id,
                "song":         song,
                "source":       "scrambled_concordance",
                "shuffle_idx":  i,
                "seq":          shuffled,
                "plddt":        plddt,
                "seq_len":      len(conc_seq),
            })

        # Scramble native_ala sequence
        for i in range(n_shuffles):
            if not na_seq:
                continue
            shuffled = shuffle_sequence(na_seq, rng)
            print(f"  na shuffle {i+1}/{n_shuffles}  {shuffled[:20]}...")
            time.sleep(REQUEST_DELAY)
            plddt = fold_with_esm(shuffled)
            print(f"    → pLDDT = {plddt:.3f}" if plddt is not None else "    → FAILED")
            rows.append({
                "bar_id":       bar_id,
                "song":         song,
                "source":       "scrambled_native_ala",
                "shuffle_idx":  i,
                "seq":          shuffled,
                "plddt":        plddt,
                "seq_len":      len(na_seq),
            })

    results = pd.DataFrame(rows)
    results.to_csv(RESULTS_CSV, index=False)
    print(f"\nSaved {len(results)} rows → {RESULTS_CSV}")

    _plot(results)


def _plot(results: pd.DataFrame) -> None:
    """
    Fig: per-bar pLDDT across 3 buckets:
      concordance   (cyan)   — original lyric-derived sequence
      native_ala    (orange) — alanine-substituted baseline
      scrambled     (grey)   — shuffled variants of concordance + native_ala combined

    Each bar gets one point per original sequence and N scatter points for scrambles.
    """
    COLORS = {
        "concordance":          "#00d4ff",   # cyan
        "native_ala":           "#ff9500",   # orange
        "scrambled_concordance": "#444455",  # dark grey-blue
        "scrambled_native_ala":  "#554444",  # dark grey-red
    }
    LABELS = {
        "concordance":           "concordance (original)",
        "native_ala":            "native_ala (original)",
        "scrambled_concordance": "scrambled concordance",
        "scrambled_native_ala":  "scrambled native_ala",
    }

    bar_ids = list(results["bar_id"].unique())
    x_pos   = {b: i for i, b in enumerate(bar_ids)}

    fig, ax = plt.subplots(figsize=(max(10, len(bar_ids) * 0.7), 5), facecolor="#0e0e0e")
    ax.set_facecolor("#0e0e0e")

    plotted_labels = set()

    for _, row in results.iterrows():
        src   = row["source"]
        x     = x_pos[row["bar_id"]]
        plddt = row["plddt"]
        if plddt is None or pd.isna(plddt):
            continue

        color  = COLORS.get(src, "#888888")
        label  = LABELS.get(src, src) if src not in plotted_labels else ""
        is_orig = row["shuffle_idx"] == -1
        size   = 70 if is_orig else 25
        zorder = 4 if is_orig else 2
        marker = "D" if is_orig else "o"

        ax.scatter(x, plddt, color=color, s=size, zorder=zorder,
                   marker=marker, alpha=0.85, label=label)
        plotted_labels.add(src)

    ax.axhline(0.3, color="#ff4444", linestyle="--", linewidth=1,
               label="pLDDT = 0.3 (disorder threshold)")

    ax.set_xticks(range(len(bar_ids)))
    ax.set_xticklabels(bar_ids, rotation=45, ha="right", fontsize=7, color="white")
    ax.set_ylabel("mean pLDDT", color="white")
    ax.set_ylim(0, 1)
    ax.set_title(
        "pLDDT: Concordance vs Native_Ala vs Scrambled\n(in silico control — shuffled = same composition, random order)",
        color="white", fontsize=11,
    )
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#333333")

    # Deduplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        if l and l not in seen:
            seen[l] = h
    ax.legend(seen.values(), seen.keys(), facecolor="#1a1a1a",
              labelcolor="white", framealpha=0.8, fontsize=8)

    plt.tight_layout()
    plt.savefig(FIG_PATH, dpi=150, facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved figure → {FIG_PATH}")

    # Summary stats
    print("\n--- pLDDT summary by bucket ---")
    for src in ["concordance", "native_ala", "scrambled_concordance", "scrambled_native_ala"]:
        vals = results.loc[results["source"] == src, "plddt"].dropna()
        if len(vals):
            print(f"  {src:25s}  n={len(vals):3d}  mean={vals.mean():.3f}  sd={vals.std():.3f}"
                  f"  <0.3: {(vals < 0.3).sum()}/{len(vals)}")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-shuffles", type=int, default=DEFAULT_N_SHUFFLE,
                        help="Number of shuffled variants per bar (default 3)")
    parser.add_argument("--all-bars", action="store_true",
                        help="Run on all 85 bars, not just phase2 candidates")
    args = parser.parse_args()

    main(args.n_shuffles, args.all_bars)
