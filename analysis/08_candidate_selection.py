"""
08_candidate_selection.py
--------------------------
Phase 2 candidate selection: filter 85 bars to a shortlist of 8-12 for
CFPS submission, and generate a per-bar concordance vs native_ala comparison.

Selection filters (applied across BOTH concordance and native_ala):
  - Boltz pTM (mean) >= 0.35  [at least one condition]
  - Boltz pLDDT (mean) >= 0.40  [at least one condition]
  - Sequence length 80-150 AA
  - Max consecutive hydrophobic residues <= 6

Ranking:
  Score = max(boltz_na_ptm, boltz_ptm) * 0.5
        + max(boltz_na_plddt, boltz_plddt) * 0.3
        + iconicity * 0.2

Outputs:
  data/phase2_candidates.csv       <- ranked shortlist
  outputs/figures/fig28_conc_vs_na_plddt.png
  outputs/figures/fig29_conc_vs_na_ptm.png
  outputs/figures/fig30_candidate_table.png

Usage:
  python analysis/08_candidate_selection.py
  python analysis/08_candidate_selection.py --ptm-min 0.30 --plddt-min 0.35
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

ENRICHED_CSV   = Path("data/aggregated_lines_v2_enriched.csv")
CANDIDATES_CSV = Path("data/phase2_candidates.csv")
FIG_DIR        = Path("outputs/figures")

HYDROPHOBIC    = set("VILMFYWCA")

DEFAULT_PTM_MIN   = 0.35
DEFAULT_PLDDT_MIN = 0.40
DEFAULT_LEN_MIN   = 80
DEFAULT_LEN_MAX   = 150
DEFAULT_HYDRO_MAX = 6

CLASS_ORDER = ["confident_protein_like", "uncertain_protein_like", "disordered"]
CLASS_COLOR = {
    "confident_protein_like":  "#00ff99",
    "uncertain_protein_like":  "#ffcc00",
    "disordered":              "#ff4444",
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def max_hydrophobic_run(seq: str) -> int:
    if not seq:
        return 0
    max_run = run = 0
    for aa in seq.upper():
        if aa in HYDROPHOBIC:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return max_run


def best_class(c_class, na_class):
    """Return whichever structural class is better."""
    for cls in CLASS_ORDER:
        if c_class == cls or na_class == cls:
            return cls
    return "disordered"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(ptm_min, plddt_min, len_min, len_max, hydro_max):
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    full = pd.read_csv(ENRICHED_CSV)
    df = full.dropna(subset=["boltz_ptm", "boltz_na_ptm"]).copy()
    print(f"Loaded {len(df)} bars with both concordance + native_ala Boltz data\n")

    # Hydrophobic run for both sequences
    df["hydro_run_conc"] = df["fasta_seq_concordance"].fillna("").apply(max_hydrophobic_run)
    df["hydro_run_na"]   = df["fasta_seq_native_alanine"].fillna("").apply(max_hydrophobic_run)
    df["hydro_run_min"]  = df[["hydro_run_conc", "hydro_run_na"]].min(axis=1)

    # Best pTM / pLDDT across both conditions
    df["best_ptm"]   = df[["boltz_ptm",   "boltz_na_ptm"]].max(axis=1)
    df["best_plddt"] = df[["boltz_plddt", "boltz_na_plddt"]].max(axis=1)
    df["best_class"] = df.apply(
        lambda r: best_class(r["boltz_structural_class"], r["boltz_na_structural_class"]),
        axis=1
    )

    # ---------- Apply filters ----------
    f_ptm    = df["best_ptm"]   >= ptm_min
    f_plddt  = df["best_plddt"] >= plddt_min
    f_len    = df["fasta_seq_len"].between(len_min, len_max)
    f_hydro  = df["hydro_run_min"] <= hydro_max

    df["pass_ptm"]   = f_ptm
    df["pass_plddt"] = f_plddt
    df["pass_len"]   = f_len
    df["pass_hydro"] = f_hydro
    df["pass_all"]   = f_ptm & f_plddt & f_len & f_hydro

    print("Filter results:")
    print(f"  pTM >= {ptm_min}:           {f_ptm.sum():3d} / {len(df)}")
    print(f"  pLDDT >= {plddt_min}:         {f_plddt.sum():3d} / {len(df)}")
    print(f"  Length {len_min}-{len_max} AA:       {f_len.sum():3d} / {len(df)}")
    print(f"  Hydrophobic run <= {hydro_max}:   {f_hydro.sum():3d} / {len(df)}")
    print(f"  ALL filters:           {df['pass_all'].sum():3d} / {len(df)}")

    candidates = df[df["pass_all"]].copy()

    # Ranking score
    candidates["rank_score"] = (
        candidates["best_ptm"]   * 0.5 +
        candidates["best_plddt"] * 0.3 +
        candidates["aggregate_iconicity"].fillna(0) * 0.2
    )
    candidates = candidates.sort_values("rank_score", ascending=False).reset_index(drop=True)
    candidates["rank"] = candidates.index + 1

    # Which condition is better per bar
    candidates["better_condition"] = candidates.apply(
        lambda r: "native_ala" if r["boltz_na_ptm"] >= r["boltz_ptm"] else "concordance",
        axis=1
    )

    # Save shortlist
    out_cols = [
        "rank", "bar_id", "genius_song_title", "aggregate_iconicity",
        "fasta_seq_len",
        "boltz_ptm", "boltz_plddt", "boltz_structural_class",
        "boltz_na_ptm", "boltz_na_plddt", "boltz_na_structural_class",
        "best_ptm", "best_plddt", "best_class",
        "better_condition", "rank_score",
        "hydro_run_conc", "hydro_run_na",
        "fasta_seq_concordance", "fasta_seq_native_alanine",
    ]
    candidates[out_cols].to_csv(CANDIDATES_CSV, index=False)
    print(f"\nSaved {len(candidates)} candidates → {CANDIDATES_CSV}")

    print("\nShortlist:")
    for _, r in candidates.iterrows():
        print(f"  #{r['rank']:2d}  {r['bar_id']:8s}  {str(r['genius_song_title'])[:30]:30s}"
              f"  pTM={r['best_ptm']:.3f}  pLDDT={r['best_plddt']:.3f}"
              f"  [{r['best_class']}]  better={r['better_condition']}")

    # ---------- Figures ----------
    _plot_plddt_comparison(df, candidates)
    _plot_ptm_comparison(df, candidates)
    _plot_candidate_table(candidates)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _scatter_base(ax):
    ax.set_facecolor("#0e0e0e")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#333333")


def _plot_plddt_comparison(df, candidates):
    """Fig 28: Concordance vs native_ala pLDDT per bar — lollipop + dots."""
    fig, ax = plt.subplots(figsize=(14, 5), facecolor="#0e0e0e")
    _scatter_base(ax)

    bars = df.sort_values("boltz_na_plddt", ascending=False).reset_index(drop=True)
    cand_ids = set(candidates["bar_id"])

    for i, row in bars.iterrows():
        conc = row["boltz_plddt"]
        na   = row["boltz_na_plddt"]
        lo, hi = min(conc, na), max(conc, na)
        color = "#444444"
        if row["bar_id"] in cand_ids:
            color = "#ff9500"

        # Connector line
        ax.plot([i, i], [lo, hi], color=color, linewidth=1.2, zorder=2)
        # Concordance (cyan)
        ax.scatter(i, conc, color="#00d4ff", s=25, zorder=3, marker="o")
        # Native_ala (orange)
        ax.scatter(i, na, color="#ff9500", s=25, zorder=3, marker="D")

    ax.axhline(0.40, color="#ffffff", linestyle="--", linewidth=0.8, alpha=0.4,
               label="pLDDT = 0.40 (filter threshold)")
    ax.axhline(0.35, color="#aaaaaa", linestyle=":", linewidth=0.8, alpha=0.4)

    ax.set_xlim(-1, len(bars))
    ax.set_ylim(0, 1)
    ax.set_xlabel("bars (sorted by native_ala pLDDT)", color="white")
    ax.set_ylabel("Boltz mean pLDDT", color="white")
    ax.set_title("Concordance vs Native_Ala — Boltz pLDDT per bar\n"
                 "orange connector = candidate passes filters",
                 color="white", fontsize=11)
    ax.set_xticks([])

    legend_els = [
        mpatches.Patch(color="#00d4ff", label="concordance"),
        mpatches.Patch(color="#ff9500", label="native_ala"),
        mpatches.Patch(color="#444444", label="non-candidate"),
        mpatches.Patch(color="#ff9500", label="candidate bar"),
    ]
    ax.legend(handles=legend_els, facecolor="#1a1a1a", labelcolor="white",
              framealpha=0.8, fontsize=8, ncol=2)

    plt.tight_layout()
    path = FIG_DIR / "fig28_conc_vs_na_plddt.png"
    plt.savefig(path, dpi=150, facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved → {path}")


def _plot_ptm_comparison(df, candidates):
    """Fig 29: 2D scatter — concordance pTM (x) vs native_ala pTM (y)."""
    fig, ax = plt.subplots(figsize=(7, 7), facecolor="#0e0e0e")
    _scatter_base(ax)

    cand_ids = set(candidates["bar_id"])

    for _, row in df.iterrows():
        is_cand = row["bar_id"] in cand_ids
        cls = row["best_class"] if "best_class" in row else "disordered"
        color = CLASS_COLOR.get(cls, "#888888")
        size  = 80 if is_cand else 30
        alpha = 1.0 if is_cand else 0.5
        marker = "D" if is_cand else "o"
        ax.scatter(row["boltz_ptm"], row["boltz_na_ptm"],
                   color=color, s=size, alpha=alpha, marker=marker, zorder=3 if is_cand else 2)
        if is_cand:
            ax.annotate(row["bar_id"], (row["boltz_ptm"], row["boltz_na_ptm"]),
                        fontsize=6, color="white", xytext=(4, 4),
                        textcoords="offset points")

    # Diagonal y=x
    lim = max(df["boltz_ptm"].max(), df["boltz_na_ptm"].max()) + 0.05
    ax.plot([0, lim], [0, lim], color="#555555", linestyle="--", linewidth=1, label="y = x")
    ax.axhline(DEFAULT_PTM_MIN, color="#aaaaaa", linestyle=":", linewidth=0.8)
    ax.axvline(DEFAULT_PTM_MIN, color="#aaaaaa", linestyle=":", linewidth=0.8)

    ax.set_xlabel("concordance pTM", color="white")
    ax.set_ylabel("native_ala pTM", color="white")
    ax.set_title("Concordance vs Native_Ala — Boltz pTM\n"
                 "above diagonal = native_ala folds better",
                 color="white", fontsize=11)

    legend_els = [
        mpatches.Patch(color=CLASS_COLOR["confident_protein_like"],  label="confident_protein_like"),
        mpatches.Patch(color=CLASS_COLOR["uncertain_protein_like"],   label="uncertain_protein_like"),
        mpatches.Patch(color=CLASS_COLOR["disordered"],               label="disordered"),
        mpatches.Patch(color="#888888", label="diamond = candidate"),
    ]
    ax.legend(handles=legend_els, facecolor="#1a1a1a", labelcolor="white",
              framealpha=0.8, fontsize=8)

    plt.tight_layout()
    path = FIG_DIR / "fig29_conc_vs_na_ptm.png"
    plt.savefig(path, dpi=150, facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved → {path}")


def _plot_candidate_table(candidates):
    """Fig 30: Table of candidates with pLDDT and pTM for both conditions."""
    cols = ["rank", "bar_id", "genius_song_title",
            "boltz_plddt", "boltz_na_plddt",
            "boltz_ptm",   "boltz_na_ptm",
            "best_class",  "better_condition", "fasta_seq_len"]
    col_labels = ["#", "bar_id", "song",
                  "conc pLDDT", "na pLDDT",
                  "conc pTM",   "na pTM",
                  "best class", "better cond.", "len"]

    tbl = candidates[cols].copy()
    tbl["genius_song_title"] = tbl["genius_song_title"].str[:22]
    for c in ["boltz_plddt", "boltz_na_plddt", "boltz_ptm", "boltz_na_ptm"]:
        tbl[c] = tbl[c].map("{:.3f}".format)

    n = len(tbl)
    fig_h = max(2.5, 0.35 * n + 1.5)
    fig, ax = plt.subplots(figsize=(14, fig_h), facecolor="#0e0e0e")
    ax.set_facecolor("#0e0e0e")
    ax.axis("off")

    cell_text = tbl.values.tolist()
    cell_colors = []
    for _, row in candidates.iterrows():
        row_c = []
        for col in cols:
            val = row[col]
            if col == "best_class":
                row_c.append(CLASS_COLOR.get(str(val), "#222222") + "55")
            elif col == "better_condition":
                row_c.append("#ff950022" if val == "native_ala" else "#00d4ff22")
            else:
                row_c.append("#1a1a1a")
        cell_colors.append(row_c)

    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellColours=cell_colors,
        cellLoc="center",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#333333")
        if row == 0:
            cell.set_facecolor("#2a2a2a")
            cell.set_text_props(color="white", fontweight="bold")
        else:
            cell.set_text_props(color="white")

    ax.set_title("Phase 2 Candidate Shortlist — Concordance vs Native_Ala",
                 color="white", fontsize=11, pad=12)

    plt.tight_layout()
    path = FIG_DIR / "fig30_candidate_table.png"
    plt.savefig(path, dpi=150, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close()
    print(f"Saved → {path}")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ptm-min",   type=float, default=DEFAULT_PTM_MIN)
    parser.add_argument("--plddt-min", type=float, default=DEFAULT_PLDDT_MIN)
    parser.add_argument("--len-min",   type=int,   default=DEFAULT_LEN_MIN)
    parser.add_argument("--len-max",   type=int,   default=DEFAULT_LEN_MAX)
    parser.add_argument("--hydro-max", type=int,   default=DEFAULT_HYDRO_MAX)
    args = parser.parse_args()

    main(args.ptm_min, args.plddt_min, args.len_min, args.len_max, args.hydro_max)
