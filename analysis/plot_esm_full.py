"""
plot_esm_full.py
----------------
Full ESMFold analysis figures — all 85 bars, all 5 conditions.

Reads outputs/esm/plddt_scores.csv (complete run).

Figures produced (outputs/figures/):
  fig11_full_plddt_by_condition.png      Mean pLDDT ± SE all 85 bars
  fig12_full_violin.png                  Violin all conditions
  fig13_native_vs_concordance.png        Per-bar paired scatter native vs concordance
  fig14_plddt_vs_length_all_conds.png    pLDDT vs length, faceted by condition
  fig15_iconicity_vs_plddt_full.png      Iconicity vs concordance pLDDT (n=85)
  fig16_concordance_vs_alanine_full.png  Per-bar delta conc - alanine (n=85)
  fig17_plddt_sd_full.png                SD per bar across seeds (concordance, n=85)
  fig18_condition_correlation.png        Pairwise pLDDT correlation matrix
  fig19_top_bottom_bars.png              Top/bottom 10 bars by concordance pLDDT
  fig20_length_bucket_heatmap.png        pLDDT heatmap by length bucket × condition

Run:
  python analysis/plot_esm_full.py
"""

import json
import math
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PLDDT_CSV  = Path("outputs/esm/plddt_scores.csv")
SNAPSHOT   = Path("data/bar_index_snapshot.json")
FIG_DIR    = Path("outputs/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)

COND_ORDER  = ["concordance", "alanine", "random", "native", "native_alanine"]
COND_COLORS = {
    "concordance":    "#2196F3",
    "alanine":        "#4CAF50",
    "random":         "#9E9E9E",
    "native":         "#FF5722",
    "native_alanine": "#FF9800",
}
COND_LABELS = {
    "concordance":    "Concordance",
    "alanine":        "Alanine",
    "random":         "Random",
    "native":         "Native",
    "native_alanine": "Native-Ala",
}

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv(PLDDT_CSV)
df = df[df["status"] == "ok"].copy()
df["plddt"] = df["plddt"].astype(float)
df["seq_len"] = df["sequence_len"].astype(int)
df["iconicity"] = df["iconicity"].astype(float)

with open(SNAPSHOT) as f:
    snapshot = json.load(f)

n_bars = df["bar_id"].nunique()
print(f"Loaded {len(df)} rows — {n_bars} bars, conditions: {df['condition'].unique().tolist()}")

# Per-bar per-condition mean
bar_cond_mean = df.groupby(["bar_id", "condition"])["plddt"].mean().reset_index()
bar_cond_mean.columns = ["bar_id", "condition", "plddt_mean"]

# Bar metadata
meta = {}
for bar_id, data in snapshot.items():
    meta[bar_id] = {
        "seq_len":   data.get("lyric_cleaned_len", 0),
        "iconicity": data.get("aggregate_iconicity", 0),
        "song":      data.get("genius_song_title", ""),
        "bar":       data.get("canonical_bar", ""),
    }

# ---------------------------------------------------------------------------
# Fig 11 — Mean pLDDT ± SE all 85 bars per condition
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
stats = df.groupby("condition")["plddt"].agg(["mean", "sem"]).reindex(COND_ORDER).dropna()
bars = ax.bar(
    [COND_LABELS[c] for c in stats.index],
    stats["mean"],
    yerr=stats["sem"],
    color=[COND_COLORS[c] for c in stats.index],
    capsize=5, edgecolor="white", linewidth=0.5,
)
ax.set_ylabel("Mean pLDDT")
ax.set_title(f"Mean pLDDT by condition (n={n_bars} bars)")
ax.set_ylim(0, 0.55)
ax.axhline(0.3, color="grey", lw=0.8, ls="--", alpha=0.5)
for bar, (_, row) in zip(bars, stats.iterrows()):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + row["sem"] + 0.005,
            f"{row['mean']:.3f}", ha="center", va="bottom", fontsize=8)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig11_full_plddt_by_condition.png", dpi=150)
plt.close()
print("fig11 done")

# ---------------------------------------------------------------------------
# Fig 12 — Violin all conditions
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))
present = [c for c in COND_ORDER if c in df["condition"].unique()]
data_violin = [df[df["condition"] == c]["plddt"].values for c in present]
parts = ax.violinplot(data_violin, positions=range(len(present)), showmedians=True, showextrema=False)
for i, (pc, c) in enumerate(zip(parts["bodies"], present)):
    pc.set_facecolor(COND_COLORS[c])
    pc.set_alpha(0.7)
parts["cmedians"].set_color("black")
parts["cmedians"].set_linewidth(1.5)
ax.set_xticks(range(len(present)))
ax.set_xticklabels([COND_LABELS[c] for c in present])
ax.set_ylabel("pLDDT")
ax.set_title(f"pLDDT distribution by condition (n={n_bars} bars)")
ax.axhline(0.3, color="grey", lw=0.8, ls="--", alpha=0.5, label="pLDDT=0.3")
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig12_full_violin.png", dpi=150)
plt.close()
print("fig12 done")

# ---------------------------------------------------------------------------
# Fig 13 — Native vs Concordance per-bar paired scatter
# ---------------------------------------------------------------------------
native_df = bar_cond_mean[bar_cond_mean["condition"] == "native"].set_index("bar_id")["plddt_mean"]
conc_df   = bar_cond_mean[bar_cond_mean["condition"] == "concordance"].set_index("bar_id")["plddt_mean"]
both = native_df.index.intersection(conc_df.index)

fig, ax = plt.subplots(figsize=(6, 6))
x = conc_df[both].values
y = native_df[both].values
ax.scatter(x, y, alpha=0.7, color="#2196F3", edgecolors="white", s=50)
lims = [min(x.min(), y.min()) - 0.01, max(x.max(), y.max()) + 0.01]
ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
ax.set_xlabel("Concordance pLDDT")
ax.set_ylabel("Native pLDDT")
ax.set_title(f"Native vs Concordance pLDDT (n={len(both)} bars)")
r = np.corrcoef(x, y)[0, 1]
ax.text(0.05, 0.93, f"r = {r:.3f}", transform=ax.transAxes, fontsize=10)
n_above = (y > x).sum()
ax.text(0.05, 0.87, f"native > conc: {n_above}/{len(both)}", transform=ax.transAxes, fontsize=9, color="grey")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig13_native_vs_concordance.png", dpi=150)
plt.close()
print("fig13 done")

# ---------------------------------------------------------------------------
# Fig 14 — pLDDT vs length, faceted by condition
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, len(present), figsize=(4 * len(present), 4), sharey=True)
for ax, cond in zip(axes, present):
    sub = df[df["condition"] == cond].groupby("bar_id").agg(
        plddt=("plddt", "mean"), seq_len=("seq_len", "first")
    ).reset_index()
    ax.scatter(sub["seq_len"], sub["plddt"], alpha=0.6, color=COND_COLORS[cond], s=30)
    r = np.corrcoef(sub["seq_len"], sub["plddt"])[0, 1]
    ax.set_title(f"{COND_LABELS[cond]}\nr={r:.3f}", fontsize=9)
    ax.set_xlabel("Seq length (AA)")
    if ax == axes[0]:
        ax.set_ylabel("Mean pLDDT")
    # trend line
    z = np.polyfit(sub["seq_len"], sub["plddt"], 1)
    xr = np.linspace(sub["seq_len"].min(), sub["seq_len"].max(), 50)
    ax.plot(xr, np.poly1d(z)(xr), "k--", lw=0.8, alpha=0.6)
fig.suptitle("pLDDT vs sequence length by condition", y=1.02)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig14_plddt_vs_length_all_conds.png", dpi=150, bbox_inches="tight")
plt.close()
print("fig14 done")

# ---------------------------------------------------------------------------
# Fig 15 — Iconicity vs concordance pLDDT (n=85)
# ---------------------------------------------------------------------------
conc_by_bar = df[df["condition"] == "concordance"].groupby("bar_id").agg(
    plddt=("plddt", "mean"),
    iconicity=("iconicity", "first"),
    seq_len=("seq_len", "first"),
).reset_index()

fig, ax = plt.subplots(figsize=(7, 5))
sc = ax.scatter(conc_by_bar["iconicity"], conc_by_bar["plddt"],
                c=conc_by_bar["seq_len"], cmap="viridis", alpha=0.8, s=50, edgecolors="white")
plt.colorbar(sc, ax=ax, label="Seq length (AA)")
r = np.corrcoef(conc_by_bar["iconicity"], conc_by_bar["plddt"])[0, 1]
ax.set_xlabel("Aggregate iconicity")
ax.set_ylabel("Mean concordance pLDDT")
ax.set_title(f"Iconicity vs pLDDT — concordance (n={len(conc_by_bar)} bars)\nr = {r:.3f}")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig15_iconicity_vs_plddt_full.png", dpi=150)
plt.close()
print("fig15 done")

# ---------------------------------------------------------------------------
# Fig 16 — Concordance vs Alanine delta (n=85)
# ---------------------------------------------------------------------------
ala_df = bar_cond_mean[bar_cond_mean["condition"] == "alanine"].set_index("bar_id")["plddt_mean"]
both16 = conc_df.index.intersection(ala_df.index)
delta = (conc_df[both16] - ala_df[both16]).sort_values()

colors16 = ["#FF5722" if d > 0 else "#2196F3" for d in delta]
fig, ax = plt.subplots(figsize=(14, 4))
ax.bar(range(len(delta)), delta.values, color=colors16, edgecolor="white", linewidth=0.3)
ax.axhline(0, color="black", lw=0.8)
ax.set_xticks([])
ax.set_ylabel("Δ pLDDT (concordance − alanine)")
ax.set_title(f"Concordance vs Alanine per bar (n={len(delta)}, orange=conc higher, blue=ala higher)")
n_pos = (delta > 0).sum()
ax.text(0.01, 0.95, f"conc > ala: {n_pos}/{len(delta)}", transform=ax.transAxes, fontsize=9)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig16_concordance_vs_alanine_full.png", dpi=150)
plt.close()
print("fig16 done")

# ---------------------------------------------------------------------------
# Fig 17 — pLDDT SD per bar (concordance, structural sensitivity)
# ---------------------------------------------------------------------------
conc_sd = df[df["condition"] == "concordance"].groupby("bar_id")["plddt"].std().sort_values(ascending=False)
fig, ax = plt.subplots(figsize=(14, 4))
ax.bar(range(len(conc_sd)), conc_sd.values, color="#2196F3", edgecolor="white", linewidth=0.3)
ax.set_xticks([])
ax.set_ylabel("SD of pLDDT (concordance, 15 seeds)")
ax.set_title(f"Structural sensitivity per bar (n={len(conc_sd)})")
ax.axhline(conc_sd.mean(), color="red", lw=0.8, ls="--", label=f"mean SD={conc_sd.mean():.4f}")
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(FIG_DIR / "fig17_plddt_sd_full.png", dpi=150)
plt.close()
print("fig17 done")

# ---------------------------------------------------------------------------
# Fig 18 — Pairwise pLDDT correlation matrix across conditions
# ---------------------------------------------------------------------------
pivot = bar_cond_mean.pivot(index="bar_id", columns="condition", values="plddt_mean")
present_cols = [c for c in COND_ORDER if c in pivot.columns]
corr = pivot[present_cols].corr()

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(corr.values, vmin=0, vmax=1, cmap="RdYlGn")
plt.colorbar(im, ax=ax)
labels = [COND_LABELS[c] for c in present_cols]
ax.set_xticks(range(len(labels)))
ax.set_yticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=9)
ax.set_yticklabels(labels, fontsize=9)
for i in range(len(labels)):
    for j in range(len(labels)):
        ax.text(j, i, f"{corr.values[i,j]:.2f}", ha="center", va="center", fontsize=9,
                color="black" if corr.values[i,j] > 0.5 else "white")
ax.set_title("Pairwise pLDDT correlation (per-bar means)")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig18_condition_correlation.png", dpi=150)
plt.close()
print("fig18 done")

# ---------------------------------------------------------------------------
# Fig 19 — Top / bottom 10 bars by concordance pLDDT
# ---------------------------------------------------------------------------
conc_ranked = conc_by_bar.sort_values("plddt", ascending=False)
top10    = conc_ranked.head(10)
bottom10 = conc_ranked.tail(10)

def short_label(bar_id, song):
    return f"{bar_id}\n{song[:22]}"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
for ax, subset, title, color in [
    (ax1, top10,    "Top 10",    "#2196F3"),
    (ax2, bottom10, "Bottom 10", "#FF5722"),
]:
    ax.barh(
        [short_label(r["bar_id"], meta.get(r["bar_id"], {}).get("song", "")) for _, r in subset.iterrows()],
        subset["plddt"].values,
        color=color, edgecolor="white",
    )
    ax.set_xlabel("Mean concordance pLDDT")
    ax.set_title(title)
    ax.invert_yaxis()
fig.suptitle("Top / Bottom 10 bars by concordance pLDDT")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig19_top_bottom_bars.png", dpi=150)
plt.close()
print("fig19 done")

# ---------------------------------------------------------------------------
# Fig 20 — pLDDT heatmap by length bucket × condition
# ---------------------------------------------------------------------------
df["len_bucket"] = pd.cut(df["seq_len"], bins=[79, 100, 130, 160, 200, 300],
                           labels=["80-100", "101-130", "131-160", "161-200", "201-300"])
heatmap_df = df.groupby(["len_bucket", "condition"])["plddt"].mean().unstack()
heatmap_df = heatmap_df[[c for c in COND_ORDER if c in heatmap_df.columns]]

fig, ax = plt.subplots(figsize=(8, 4))
im = ax.imshow(heatmap_df.values, aspect="auto", cmap="YlOrRd", vmin=0.25, vmax=0.50)
plt.colorbar(im, ax=ax, label="Mean pLDDT")
ax.set_xticks(range(len(heatmap_df.columns)))
ax.set_xticklabels([COND_LABELS[c] for c in heatmap_df.columns], rotation=20, ha="right", fontsize=9)
ax.set_yticks(range(len(heatmap_df.index)))
ax.set_yticklabels(heatmap_df.index, fontsize=9)
for i in range(len(heatmap_df.index)):
    for j in range(len(heatmap_df.columns)):
        v = heatmap_df.values[i, j]
        if not np.isnan(v):
            ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=8)
ax.set_title("Mean pLDDT: length bucket × condition")
plt.tight_layout()
plt.savefig(FIG_DIR / "fig20_length_bucket_heatmap.png", dpi=150)
plt.close()
print("fig20 done")

print(f"\nAll figures written to {FIG_DIR}/")
print(f"Key stats:")
print(f"  Bars:       {n_bars}")
print(f"  Conditions: {len(present)}")
cond_stats = df.groupby("condition")["plddt"].agg(["mean","std"]).reindex(COND_ORDER).dropna()
for cond, row in cond_stats.iterrows():
    print(f"  {COND_LABELS[cond]:<15}: mean={row['mean']:.4f}  SD={row['std']:.4f}")
