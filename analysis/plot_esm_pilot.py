"""
plot_esm_pilot.py
-----------------
Generates all summary figures from the ESMFold pilot data.

Figures produced (outputs/figures/):
  fig01_length_distribution.png       Length histogram for 85 bars
  fig02_iconicity_distribution.png    Iconicity score distribution by badge
  fig03_mean_plddt_by_condition.png   Bar chart: mean pLDDT ± SE per condition
  fig04_plddt_violin.png              Violin: pLDDT distribution per condition
  fig05_plddt_vs_length.png           Scatter: pLDDT vs sequence length (concordance)
  fig06_plddt_by_condition_bucket.png Grouped bar: pLDDT by condition × length bucket
  fig07_concordance_vs_alanine.png    Per-bar delta: concordance minus alanine
  fig08_iconicity_vs_plddt.png        Scatter: iconicity vs mean concordance pLDDT
  fig09_plddt_sd_per_bar.png          Structural sensitivity: SD across 15 seeds (concordance)
  fig10_factorial_effects.png         2×2 factorial effect sizes

Run:
  python analysis/plot_esm_pilot.py
"""

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PLDDT_CSV   = Path("outputs/esm/plddt_scores.csv")
SNAPSHOT    = Path("data/bar_index_snapshot.json")
FIGURES_DIR = Path("outputs/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

CONDITION_ORDER = ["concordance", "alanine", "random", "native", "native_alanine"]
CONDITION_COLORS = {
    "concordance":    "#4E79A7",
    "alanine":        "#59A14F",
    "random":         "#F28E2B",
    "native":         "#B07AA1",
    "native_alanine": "#76B7B2",
}
BADGE_COLORS = {"Viral": "#E15759", "Aligned": "#4E79A7", "Deep Cut": "#59A14F"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "figure.dpi": 150,
})

def save(fig, name):
    path = FIGURES_DIR / name
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Saved {path}")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

print("Loading data...")
esm = pd.read_csv(PLDDT_CSV)
esm = esm[esm["status"] == "ok"].copy()

with open(SNAPSHOT) as f:
    snap_raw = json.load(f)

snap = pd.DataFrame([{"bar_id": k, **v} for k, v in snap_raw.items()])
snap["aggregate_iconicity"] = snap["aggregate_iconicity"].astype(float)
snap["lyric_cleaned_len"]   = snap["lyric_cleaned_len"].astype(int)

# Per-bar summary (concordance only, 15 seeds)
conc = esm[esm["condition"] == "concordance"]
bar_conc = conc.groupby("bar_id").agg(
    plddt_mean=("plddt", "mean"),
    plddt_sd=("plddt", "std"),
    plddt_min=("plddt", "min"),
    plddt_max=("plddt", "max"),
    n=("plddt", "count"),
    sequence_len=("sequence_len", "first"),
    iconicity=("iconicity", "first"),
    song=("song", "first"),
).reset_index()
bar_conc = bar_conc.merge(snap[["bar_id","divergence_badge","lyric_cleaned_len"]], on="bar_id", how="left")

# Per-bar per-condition means
bar_cond = esm.groupby(["bar_id", "condition"]).agg(
    plddt_mean=("plddt", "mean"),
    sequence_len=("sequence_len", "first"),
    iconicity=("iconicity", "first"),
).reset_index()
bar_cond = bar_cond.merge(snap[["bar_id","lyric_cleaned_len"]], on="bar_id", how="left")

def length_bucket(length):
    if length <= 100:
        return "80–100"
    elif length <= 120:
        return "101–120"
    else:
        return "121–155"

bar_cond["bucket"] = bar_cond["sequence_len"].apply(length_bucket)
bar_conc["bucket"] = bar_conc["sequence_len"].apply(length_bucket)
BUCKET_ORDER = ["80–100", "101–120", "121–155"]

print(f"  ESMFold rows (ok): {len(esm)}")
print(f"  Bars in snapshot:  {len(snap)}")
print(f"  Top-25 bars (concordance): {len(bar_conc)}")

# ---------------------------------------------------------------------------
# Fig 01 — Length distribution of 85 bars
# ---------------------------------------------------------------------------

print("\nFig 01 — Length distribution")
fig, ax = plt.subplots(figsize=(8, 4))

lengths = snap["lyric_cleaned_len"].values
bins = range(80, 165, 5)
ax.hist(lengths, bins=bins, color="#4E79A7", edgecolor="white", linewidth=0.5)

ax.axvline(100, color="#E15759", linestyle="--", linewidth=1.2, label="100 AA")
ax.axvline(120, color="#F28E2B", linestyle="--", linewidth=1.2, label="120 AA")

ax.set_xlabel("Sequence length (AA)")
ax.set_ylabel("Number of bars")
ax.set_title("Fig 01 — Sequence length distribution (85 bars, 80–300 AA filter)")
ax.legend(frameon=False, fontsize=9)

stats_txt = f"n=85  min={lengths.min()}  median={int(np.median(lengths))}  max={lengths.max()}"
ax.text(0.98, 0.95, stats_txt, transform=ax.transAxes,
        ha="right", va="top", fontsize=9, color="#555")
save(fig, "fig01_length_distribution.png")

# ---------------------------------------------------------------------------
# Fig 02 — Iconicity distribution by badge
# ---------------------------------------------------------------------------

print("Fig 02 — Iconicity distribution")
fig, ax = plt.subplots(figsize=(8, 4))

# Only bars in the pilot (top 25 by iconicity, in snap)
top25_ids = bar_conc["bar_id"].tolist()
snap_top = snap[snap["bar_id"].isin(top25_ids)].copy()
snap_all = snap.copy()

for badge, grp in snap_all.groupby("divergence_badge"):
    ax.hist(grp["aggregate_iconicity"], bins=15,
            color=BADGE_COLORS.get(badge, "#aaa"), alpha=0.7,
            label=badge, edgecolor="white", linewidth=0.4)

# Mark top-25 cutoff
cutoff = sorted(snap["aggregate_iconicity"], reverse=True)[24]
ax.axvline(cutoff, color="black", linestyle=":", linewidth=1.2,
           label=f"Top-25 cutoff ({cutoff:.3f})")

ax.set_xlabel("Aggregate iconicity score")
ax.set_ylabel("Number of bars")
ax.set_title("Fig 02 — Iconicity score distribution (85 bars, coloured by divergence badge)")
ax.legend(frameon=False, fontsize=9)
save(fig, "fig02_iconicity_distribution.png")

# ---------------------------------------------------------------------------
# Fig 03 — Mean pLDDT by condition (bar chart ± SE)
# ---------------------------------------------------------------------------

print("Fig 03 — Mean pLDDT by condition")
fig, ax = plt.subplots(figsize=(8, 5))

cond_stats = (
    esm.groupby(["bar_id", "condition"])["plddt"].mean()
       .reset_index()
       .groupby("condition")["plddt"]
       .agg(["mean", "sem", "count"])
       .reindex(CONDITION_ORDER)
)

bars = ax.bar(
    range(len(CONDITION_ORDER)),
    cond_stats["mean"],
    yerr=cond_stats["sem"],
    color=[CONDITION_COLORS[c] for c in CONDITION_ORDER],
    capsize=4, error_kw={"linewidth": 1.2},
    edgecolor="white", linewidth=0.5,
)

for i, (cond, row) in enumerate(cond_stats.iterrows()):
    ax.text(i, row["mean"] + row["sem"] + 0.004,
            f"{row['mean']:.4f}", ha="center", va="bottom", fontsize=9)
    ax.text(i, 0.01, f"n={int(row['count'])}", ha="center", va="bottom",
            fontsize=8, color="#666")

ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels(CONDITION_ORDER, rotation=15, ha="right")
ax.set_ylabel("Mean pLDDT")
ax.set_ylim(0, max(cond_stats["mean"]) * 1.25)
ax.set_title("Fig 03 — Mean pLDDT by condition ± SE\n(top 25 bars, bar-level means then averaged)")
save(fig, "fig03_mean_plddt_by_condition.png")

# ---------------------------------------------------------------------------
# Fig 04 — pLDDT violin by condition
# ---------------------------------------------------------------------------

print("Fig 04 — pLDDT violin by condition")
fig, ax = plt.subplots(figsize=(9, 5))

# Bar-level means per condition (not per seed)
bar_means = (
    esm.groupby(["bar_id", "condition"])["plddt"].mean()
       .reset_index()
)

data_by_cond = [bar_means[bar_means["condition"] == c]["plddt"].values
                for c in CONDITION_ORDER]

vp = ax.violinplot(data_by_cond, positions=range(len(CONDITION_ORDER)),
                   showmedians=True, showextrema=True)

for i, (pc, cond) in enumerate(zip(vp["bodies"], CONDITION_ORDER)):
    pc.set_facecolor(CONDITION_COLORS[cond])
    pc.set_alpha(0.75)

for part in ("cmedians", "cmins", "cmaxes", "cbars"):
    vp[part].set_color("#333")
    vp[part].set_linewidth(1.2)

# Overlay individual bar points
for i, cond in enumerate(CONDITION_ORDER):
    vals = bar_means[bar_means["condition"] == cond]["plddt"].values
    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals,
               color=CONDITION_COLORS[cond], s=20, alpha=0.6, zorder=3)

ax.set_xticks(range(len(CONDITION_ORDER)))
ax.set_xticklabels(CONDITION_ORDER, rotation=15, ha="right")
ax.set_ylabel("pLDDT (bar-level mean)")
ax.set_title("Fig 04 — pLDDT distribution by condition\n(bar-level means, n=25 per condition, dots = individual bars)")
save(fig, "fig04_plddt_violin.png")

# ---------------------------------------------------------------------------
# Fig 05 — pLDDT vs sequence length (concordance)
# ---------------------------------------------------------------------------

print("Fig 05 — pLDDT vs sequence length")
fig, ax = plt.subplots(figsize=(8, 5))

sc = ax.scatter(
    bar_conc["sequence_len"], bar_conc["plddt_mean"],
    c=bar_conc["iconicity"], cmap="YlOrRd", s=60,
    edgecolors="#333", linewidths=0.4, zorder=3,
)
cb = plt.colorbar(sc, ax=ax, shrink=0.8)
cb.set_label("Iconicity score", fontsize=9)

# Linear trend
m, b = np.polyfit(bar_conc["sequence_len"], bar_conc["plddt_mean"], 1)
xs = np.linspace(bar_conc["sequence_len"].min(), bar_conc["sequence_len"].max(), 100)
ax.plot(xs, m * xs + b, "--", color="#555", linewidth=1.2, label=f"trend (slope={m:.5f})")

r = np.corrcoef(bar_conc["sequence_len"], bar_conc["plddt_mean"])[0, 1]
ax.text(0.97, 0.95, f"r = {r:.3f}", transform=ax.transAxes,
        ha="right", va="top", fontsize=10, color="#333")

# Label top 3 and bottom 3
for _, row in bar_conc.nlargest(3, "plddt_mean").iterrows():
    ax.annotate(row["bar_id"], (row["sequence_len"], row["plddt_mean"]),
                xytext=(4, 4), textcoords="offset points", fontsize=7, color="#333")
for _, row in bar_conc.nsmallest(2, "plddt_mean").iterrows():
    ax.annotate(row["bar_id"], (row["sequence_len"], row["plddt_mean"]),
                xytext=(4, -8), textcoords="offset points", fontsize=7, color="#888")

ax.set_xlabel("Sequence length (AA)")
ax.set_ylabel("Mean pLDDT (concordance, 15 seeds)")
ax.set_title("Fig 05 — pLDDT vs sequence length (concordance condition)\nColour = iconicity score")
ax.legend(frameon=False, fontsize=9)
save(fig, "fig05_plddt_vs_length.png")

# ---------------------------------------------------------------------------
# Fig 06 — Mean pLDDT by condition × length bucket
# ---------------------------------------------------------------------------

print("Fig 06 — pLDDT by condition × length bucket")
fig, ax = plt.subplots(figsize=(10, 5))

bucket_cond = (
    bar_cond.groupby(["condition", "bucket"])["plddt_mean"]
            .agg(["mean", "sem", "count"])
            .reset_index()
)

n_conds = len(CONDITION_ORDER)
n_buckets = len(BUCKET_ORDER)
width = 0.15
x = np.arange(n_buckets)

for i, cond in enumerate(CONDITION_ORDER):
    sub = bucket_cond[bucket_cond["condition"] == cond].set_index("bucket")
    means = [sub.loc[b, "mean"] if b in sub.index else np.nan for b in BUCKET_ORDER]
    sems  = [sub.loc[b, "sem"]  if b in sub.index else 0      for b in BUCKET_ORDER]
    ns    = [int(sub.loc[b, "count"]) if b in sub.index else 0 for b in BUCKET_ORDER]
    offset = (i - n_conds / 2 + 0.5) * width
    bars = ax.bar(x + offset, means, width, yerr=sems,
                  color=CONDITION_COLORS[cond], label=cond,
                  capsize=3, error_kw={"linewidth": 1}, alpha=0.85)
    for j, (bar, n) in enumerate(zip(bars, ns)):
        if n > 0:
            ax.text(bar.get_x() + bar.get_width()/2, 0.005,
                    f"n={n}", ha="center", va="bottom", fontsize=6.5, color="#555",
                    rotation=90)

ax.set_xticks(x)
ax.set_xticklabels(BUCKET_ORDER)
ax.set_xlabel("Sequence length bucket (AA)")
ax.set_ylabel("Mean pLDDT")
ax.set_ylim(0, 0.65)
ax.set_title("Fig 06 — Mean pLDDT by condition × length bucket\n(top 25 bars, bar-level means)")
ax.legend(frameon=False, fontsize=9, ncol=2)
save(fig, "fig06_plddt_by_condition_bucket.png")

# ---------------------------------------------------------------------------
# Fig 07 — Per-bar delta: concordance minus alanine
# ---------------------------------------------------------------------------

print("Fig 07 — Concordance vs alanine delta")
pivot = bar_cond.pivot(index="bar_id", columns="condition", values="plddt_mean")
pivot = pivot.join(bar_cond.drop_duplicates("bar_id").set_index("bar_id")[["sequence_len", "bucket"]])
pivot["delta"] = pivot["concordance"] - pivot["alanine"]
pivot = pivot.dropna(subset=["delta"]).sort_values("delta")

fig, ax = plt.subplots(figsize=(10, 5))
colors = [CONDITION_COLORS["concordance"] if d >= 0 else CONDITION_COLORS["alanine"]
          for d in pivot["delta"]]
ax.bar(range(len(pivot)), pivot["delta"], color=colors, edgecolor="white", linewidth=0.4)
ax.axhline(0, color="#333", linewidth=1)

ax.set_xticks(range(len(pivot)))
ax.set_xticklabels(pivot.index, rotation=90, fontsize=7)
ax.set_xlabel("Bar ID (sorted by delta)")
ax.set_ylabel("pLDDT delta (concordance − alanine)")
ax.set_title("Fig 07 — Per-bar pLDDT delta: concordance minus alanine\nBlue = concordance higher  |  Green = alanine higher")

conc_patch = mpatches.Patch(color=CONDITION_COLORS["concordance"], label="concordance higher")
alanine_patch = mpatches.Patch(color=CONDITION_COLORS["alanine"], label="alanine higher")
ax.legend(handles=[conc_patch, alanine_patch], frameon=False, fontsize=9)

mean_delta = pivot["delta"].mean()
ax.text(0.98, 0.95, f"mean Δ = {mean_delta:.4f}", transform=ax.transAxes,
        ha="right", va="top", fontsize=10)
save(fig, "fig07_concordance_vs_alanine.png")

# ---------------------------------------------------------------------------
# Fig 08 — Iconicity vs mean pLDDT
# ---------------------------------------------------------------------------

print("Fig 08 — Iconicity vs mean pLDDT")
fig, ax = plt.subplots(figsize=(8, 5))

for badge, grp in bar_conc.groupby("divergence_badge"):
    ax.scatter(grp["iconicity"], grp["plddt_mean"],
               label=badge, color=BADGE_COLORS.get(badge, "#aaa"),
               s=55, alpha=0.8, edgecolors="#333", linewidths=0.4, zorder=3)

m, b = np.polyfit(bar_conc["iconicity"], bar_conc["plddt_mean"], 1)
xs = np.linspace(bar_conc["iconicity"].min(), bar_conc["iconicity"].max(), 100)
ax.plot(xs, m * xs + b, "--", color="#555", linewidth=1.2)

r = np.corrcoef(bar_conc["iconicity"], bar_conc["plddt_mean"])[0, 1]
ax.text(0.97, 0.95, f"r = {r:.3f}", transform=ax.transAxes,
        ha="right", va="top", fontsize=11, color="#333",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#ccc"))

# Label top and bottom pLDDT bars
for _, row in bar_conc.nlargest(2, "plddt_mean").iterrows():
    ax.annotate(f"{row['bar_id']}\n{row['song'][:20]}",
                (row["iconicity"], row["plddt_mean"]),
                xytext=(6, 4), textcoords="offset points", fontsize=7)

ax.set_xlabel("Aggregate iconicity score")
ax.set_ylabel("Mean pLDDT (concordance, 15 seeds)")
ax.set_title("Fig 08 — Iconicity vs structural confidence\n(top 25 bars, concordance condition)")
ax.legend(frameon=False, fontsize=9)
save(fig, "fig08_iconicity_vs_plddt.png")

# ---------------------------------------------------------------------------
# Fig 09 — pLDDT SD per bar (structural sensitivity)
# ---------------------------------------------------------------------------

print("Fig 09 — pLDDT SD per bar")
fig, ax = plt.subplots(figsize=(10, 5))

bar_conc_sorted = bar_conc.sort_values("plddt_sd", ascending=False).reset_index(drop=True)
bucket_palette = {"80–100": "#4E79A7", "101–120": "#F28E2B", "121–155": "#76B7B2"}
colors = [bucket_palette.get(b, "#aaa") for b in bar_conc_sorted["bucket"]]

ax.bar(range(len(bar_conc_sorted)), bar_conc_sorted["plddt_sd"],
       color=colors, edgecolor="white", linewidth=0.4)

ax.set_xticks(range(len(bar_conc_sorted)))
ax.set_xticklabels(bar_conc_sorted["bar_id"], rotation=90, fontsize=7)
ax.set_xlabel("Bar ID (sorted by SD, highest first)")
ax.set_ylabel("pLDDT standard deviation (15 seeds)")
ax.set_title("Fig 09 — Structural sensitivity: pLDDT SD across seeds (concordance condition)\nHigh SD = BOJUXZ substitutions strongly affect fold confidence")

patches = [mpatches.Patch(color=c, label=b) for b, c in bucket_palette.items()]
ax.legend(handles=patches, frameon=False, fontsize=9, title="Length bucket")

mean_sd = bar_conc_sorted["plddt_sd"].mean()
ax.axhline(mean_sd, color="#333", linestyle=":", linewidth=1.2)
ax.text(len(bar_conc_sorted) - 0.5, mean_sd + 0.001,
        f"mean SD = {mean_sd:.4f}", ha="right", fontsize=9, color="#333")
save(fig, "fig09_plddt_sd_per_bar.png")

# ---------------------------------------------------------------------------
# Fig 10 — 2×2 factorial effect sizes
# ---------------------------------------------------------------------------

print("Fig 10 — Factorial effect sizes")

# Compute effects using bar-level means
bp = bar_cond.pivot(index="bar_id", columns="condition", values="plddt_mean").dropna()

effects = {
    "conc − native\n(freq remapping, softmax BOJUXZ)":
        (bp["concordance"] - bp["native"]).mean(),
    "alanine − native_alanine\n(freq remapping, alanine BOJUXZ)":
        (bp["alanine"] - bp["native_alanine"]).mean(),
    "conc − alanine\n(softmax vs alanine, freq-mapped)":
        (bp["concordance"] - bp["alanine"]).mean(),
    "native − native_alanine\n(softmax vs alanine, native)":
        (bp["native"] - bp["native_alanine"]).mean(),
    "conc − random\n(softmax vs random)":
        (bp["concordance"] - bp["random"]).mean(),
    "alanine − random\n(alanine vs random)":
        (bp["alanine"] - bp["random"]).mean(),
}

fig, ax = plt.subplots(figsize=(10, 5))
labels = list(effects.keys())
values = list(effects.values())
colors = ["#59A14F" if v >= 0 else "#E15759" for v in values]

bars = ax.barh(range(len(labels)), values, color=colors, edgecolor="white", linewidth=0.5)
ax.axvline(0, color="#333", linewidth=1)

for i, (bar, val) in enumerate(zip(bars, values)):
    ax.text(val + (0.0005 if val >= 0 else -0.0005), i,
            f"{val:+.4f}", va="center",
            ha="left" if val >= 0 else "right", fontsize=9)

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel("Mean pLDDT delta (row condition A − condition B)")
ax.set_title("Fig 10 — 2×2 factorial effect sizes\nGreen = first condition higher  |  Red = second condition higher")

pos_patch = mpatches.Patch(color="#59A14F", label="first condition higher")
neg_patch = mpatches.Patch(color="#E15759", label="second condition higher")
ax.legend(handles=[pos_patch, neg_patch], frameon=False, fontsize=9, loc="lower right")
save(fig, "fig10_factorial_effects.png")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

print(f"\n--- Done. {len(list(FIGURES_DIR.glob('fig*.png')))} figures written to {FIGURES_DIR} ---")
for f in sorted(FIGURES_DIR.glob("fig*.png")):
    print(f"  {f.name}")
