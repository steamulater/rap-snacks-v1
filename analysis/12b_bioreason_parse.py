"""
12b_bioreason_parse.py
-----------------------
Parse BioReason-Pro output TSV and cross-reference with project metadata.

Answers:
  1. Does BioReason predict confident GO function for lyric-derived sequences?
  2. Does predicted GO confidence correlate with Boltz pLDDT or ESMFold pLDDT?
  3. Does iconicity correlate with predicted functional confidence?
  4. Do MPNN-designed sequences (free_design / native_ala_free) get more
     confident / different GO predictions than concordance sequences?
  5. Which GO terms / protein families does BioReason hallucinate for rap lyrics?

Usage:
  python analysis/12b_bioreason_parse.py \
    --results outputs/bioreason/bioreason_results.tsv

BioReason-Pro output columns (expected):
  protein_id, GO_term, GO_description, confidence, reasoning_trace
  (adjust RESULT_COLS below if different)
"""

import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT       = Path(__file__).parent.parent
ENRICHED   = ROOT / "data/aggregated_lines_v2_enriched.csv"
BOLTZ_CSV  = ROOT / "outputs/boltz_validation/boltz_validation_results.csv"
NAF_CSV    = ROOT / "outputs/proteinmpnn_native_ala_free/filtered_results.csv"
FREE_CSV   = ROOT / "outputs/proteinmpnn_free/filtered_results.csv"
OUT_DIR    = ROOT / "outputs/bioreason"
FIG_DIR    = ROOT / "outputs/figures"

DARK_BG = "#0e0e0e"
COLORS  = {
    "concordance":     "#00d4ff",
    "native_ala":      "#ff9500",
    "free_design":     "#a855f7",
    "native_ala_free": "#39d353",
    "breadth":         "#888888",
}

# Adjust if BioReason output columns differ
RESULT_COLS = {
    "id":         "protein_id",
    "go":         "GO_term",
    "go_desc":    "GO_description",
    "confidence": "confidence",
    "trace":      "reasoning_trace",
}


def parse_protein_id(pid: str) -> dict:
    """Reverse the make_id() encoding from 12_bioreason_prep.py."""
    parts = pid.split("__")
    return {
        "bar_id":   parts[0] if len(parts) > 0 else pid,
        "bucket":   parts[1] if len(parts) > 1 else "unknown",
        "song":     parts[2].replace("_", " ") if len(parts) > 2 else "",
        "iconicity": float(parts[3].replace("ic", "").split("_")[0]) if len(parts) > 3 else float("nan"),
    }


def main(results_path: Path) -> None:
    if not results_path.exists():
        print(f"[ERROR] Results not found: {results_path}")
        print("Run BioReason-Pro first:")
        print("  python predict.py --input outputs/bioreason/bioreason_all.tsv \\")
        print("                    --output outputs/bioreason/bioreason_results.tsv")
        return

    results = pd.read_csv(results_path, sep="\t")
    print(f"Loaded {len(results)} BioReason-Pro predictions from {results_path}")
    print(f"Columns: {list(results.columns)}")

    # Parse protein_id metadata
    meta = results[RESULT_COLS["id"]].apply(parse_protein_id).apply(pd.Series)
    results = pd.concat([results, meta], axis=1)

    # Load Boltz pLDDT for cross-reference
    boltz = pd.read_csv(BOLTZ_CSV) if BOLTZ_CSV.exists() else pd.DataFrame()
    enriched = pd.read_csv(ENRICHED)

    # Merge Boltz pLDDT by bar_id + bucket
    if not boltz.empty:
        boltz_mean = boltz.groupby(["bar_id", "bucket"])["boltz_plddt"].mean().reset_index()
        boltz_mean.rename(columns={"boltz_plddt": "boltz_plddt_mean"}, inplace=True)
        results = results.merge(boltz_mean, on=["bar_id", "bucket"], how="left")

    conf_col = RESULT_COLS["confidence"]
    if conf_col not in results.columns:
        print(f"[WARN] '{conf_col}' column not found. Available: {list(results.columns)}")
        conf_col = results.select_dtypes(include="number").columns[0]
        print(f"  Using '{conf_col}' as confidence proxy.")

    results["_conf"] = pd.to_numeric(results[conf_col], errors="coerce")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. Per-bucket confidence distribution ─────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 5), facecolor=DARK_BG)
    ax.set_facecolor(DARK_BG)
    buckets = [b for b in ["concordance", "native_ala", "free_design", "native_ala_free"]
               if b in results["bucket"].unique()]
    for xi, bucket in enumerate(buckets):
        sub = results[results["bucket"] == bucket]["_conf"].dropna()
        if sub.empty:
            continue
        parts = ax.violinplot([sub.values], positions=[xi], widths=0.7, showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor(COLORS.get(bucket, "#888")); pc.set_alpha(0.75)
        parts["cmedians"].set_colors(COLORS.get(bucket, "#888"))
        for k in ["cbars", "cmins", "cmaxes"]:
            parts[k].set_colors("#555")
    ax.set_xticks(range(len(buckets)))
    ax.set_xticklabels(buckets, rotation=15, ha="right", color="white", fontsize=10)
    ax.set_ylabel("BioReason-Pro confidence", color="white")
    ax.set_title("BioReason GO prediction confidence — rap lyrics as protein", color="white", fontsize=12)
    ax.tick_params(colors="white")
    for s in ax.spines.values(): s.set_edgecolor("#333")
    plt.tight_layout()
    fig_path = FIG_DIR / "fig_bioreason_confidence_violin.png"
    plt.savefig(fig_path, dpi=150, facecolor=DARK_BG, bbox_inches="tight")
    print(f"Saved {fig_path.name}")

    # ── 2. Confidence vs Boltz pLDDT scatter ─────────────────────────────
    if "boltz_plddt_mean" in results.columns:
        fig, ax = plt.subplots(figsize=(8, 6), facecolor=DARK_BG)
        ax.set_facecolor(DARK_BG)
        for bucket in buckets:
            sub = results[results["bucket"] == bucket].dropna(subset=["_conf", "boltz_plddt_mean"])
            ax.scatter(sub["boltz_plddt_mean"], sub["_conf"],
                       color=COLORS.get(bucket, "#888"), alpha=0.6, s=30, label=bucket)
        ax.set_xlabel("Boltz-2 pLDDT", color="white")
        ax.set_ylabel("BioReason confidence", color="white")
        ax.set_title("BioReason confidence vs Boltz pLDDT", color="white", fontsize=11)
        ax.tick_params(colors="white")
        for s in ax.spines.values(): s.set_edgecolor("#333")
        ax.legend(facecolor="#1a1a1a", labelcolor="white", fontsize=9)
        # Correlation
        sub_all = results.dropna(subset=["_conf", "boltz_plddt_mean"])
        if len(sub_all) > 2:
            r = np.corrcoef(sub_all["boltz_plddt_mean"], sub_all["_conf"])[0, 1]
            ax.text(0.05, 0.95, f"r = {r:.3f}", transform=ax.transAxes,
                    color="white", fontsize=10, va="top")
        plt.tight_layout()
        plt.savefig(FIG_DIR / "fig_bioreason_vs_plddt.png", dpi=150, facecolor=DARK_BG, bbox_inches="tight")
        print("Saved fig_bioreason_vs_plddt.png")

    # ── 3. Confidence vs iconicity ────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5), facecolor=DARK_BG)
    ax.set_facecolor(DARK_BG)
    conc = results[results["bucket"] == "concordance"].dropna(subset=["_conf", "iconicity"])
    if not conc.empty:
        ax.scatter(conc["iconicity"], conc["_conf"], color=COLORS["concordance"], alpha=0.7, s=40)
        if len(conc) > 2:
            r = np.corrcoef(conc["iconicity"], conc["_conf"])[0, 1]
            ax.text(0.05, 0.95, f"r = {r:.3f} (concordance only)", transform=ax.transAxes,
                    color="white", fontsize=10, va="top")
    ax.set_xlabel("Cultural iconicity score", color="white")
    ax.set_ylabel("BioReason confidence", color="white")
    ax.set_title("BioReason confidence vs iconicity — concordance sequences", color="white", fontsize=11)
    ax.tick_params(colors="white")
    for s in ax.spines.values(): s.set_edgecolor("#333")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "fig_bioreason_vs_iconicity.png", dpi=150, facecolor=DARK_BG, bbox_inches="tight")
    print("Saved fig_bioreason_vs_iconicity.png")

    # ── 4. Top GO terms predicted ─────────────────────────────────────────
    go_col = RESULT_COLS["go_desc"]
    if go_col in results.columns:
        top_go = (results.groupby(go_col)["_conf"]
                  .agg(["count", "mean"])
                  .sort_values("mean", ascending=False)
                  .head(20))
        print("\n--- Top 20 GO terms predicted for rap lyric sequences ---")
        print(top_go.to_string())
        top_go.to_csv(OUT_DIR / "bioreason_top_go_terms.csv")

    # ── 5. Per-bar summary table ──────────────────────────────────────────
    summary = (results.groupby(["bar_id", "bucket"])["_conf"]
               .agg(n="count", mean_conf="mean", max_conf="max")
               .reset_index()
               .sort_values(["bar_id", "bucket"]))
    summary.to_csv(OUT_DIR / "bioreason_per_bar_summary.csv", index=False)
    print("\n--- Per-bar BioReason confidence summary ---")
    print(summary.to_string(index=False))

    # ── 6. Reasoning trace sample ─────────────────────────────────────────
    trace_col = RESULT_COLS["trace"]
    if trace_col in results.columns:
        print("\n--- Sample reasoning traces (top 3 by confidence) ---")
        top3 = results.nlargest(3, "_conf")
        for _, row in top3.iterrows():
            print(f"\n{row[RESULT_COLS['id']]} (conf={row['_conf']:.3f})")
            print(f"  GO: {row.get(RESULT_COLS['go_desc'], 'N/A')}")
            trace = str(row.get(trace_col, ""))[:500]
            print(f"  Trace: {trace}...")

    print(f"\nAll outputs → outputs/bioreason/ and outputs/figures/")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--results", type=Path,
                   default=ROOT / "outputs/bioreason/bioreason_results.tsv")
    args = p.parse_args()
    main(args.results)
