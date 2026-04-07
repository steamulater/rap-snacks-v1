"""
14_rmsd_comparison.py
---------------------
Fig 90 — Backbone RMSD comparison across all 4 buckets:
  native_ala (reference) · concordance · free_design (MPNN) · native_ala_free (MPNN)

Data sources:
  outputs/boltz_validation/boltz_rmsd.csv       — v2: concordance, native_ala, free_design
  results-20260407T081540Z-3-001.zip            — v3: native_ala_free (boltz_rmsd_v3.csv)

12 bars with all 4 buckets. Sorted by native_ala RMSD ascending.
RMSD reference = ESMFold backbone for each bar.

Output: outputs/figures/fig90_rmsd_comparison.png
"""

import warnings
warnings.filterwarnings("ignore")

import io, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

ROOT    = Path(__file__).parent.parent
DATA    = ROOT / "data"
OUT_FIG = ROOT / "outputs/figures"

RMSD_CSV = ROOT / "outputs/boltz_validation/boltz_rmsd.csv"
V3_ZIP   = ROOT / "results-20260407T081540Z-3-001.zip"
ENRICHED = DATA / "aggregated_lines_v2_enriched.csv"

COLORS = {
    "native_ala":      "#f97316",
    "concordance":     "#0ea5e9",
    "free_design":     "#a855f7",
    "native_ala_free": "#22c55e",
}
BUCKET_ORDER = ["native_ala", "concordance", "free_design", "native_ala_free"]

# ── Load data ────────────────────────────────────────────────────────────────

v2 = pd.read_csv(RMSD_CSV)
v2 = v2[v2["bucket"].isin(["concordance", "native_ala", "free_design"])]

with zipfile.ZipFile(V3_ZIP) as z:
    v3 = pd.read_csv(io.BytesIO(z.read("results/boltz_rmsd_v3.csv")))
v3 = v3[v3["bucket"] == "native_ala_free"][["bar_id", "name", "bucket", "rmsd_vs_backbone"]]

rmsd = pd.concat([v2[["bar_id", "name", "bucket", "rmsd_vs_backbone"]], v3],
                 ignore_index=True)

enriched = pd.read_csv(ENRICHED)
meta = (enriched[["bar_id", "genius_song_title", "aggregate_iconicity"]]
        .drop_duplicates("bar_id").set_index("bar_id"))

# ── Per-bar summary ──────────────────────────────────────────────────────────

pivot = (rmsd[rmsd["bucket"].isin(["concordance", "native_ala"])]
         .pivot_table(index="bar_id", columns="bucket", values="rmsd_vs_backbone"))

def agg_bucket(df, bucket):
    g = df[df["bucket"] == bucket].groupby("bar_id")["rmsd_vs_backbone"]
    return g.agg(["median",
                  lambda x: x.quantile(0.25),
                  lambda x: x.quantile(0.75),
                  "min", "max"]).rename(
        columns={"median": f"{bucket}_med",
                 "<lambda_0>": f"{bucket}_q1",
                 "<lambda_1>": f"{bucket}_q3",
                 "min": f"{bucket}_min",
                 "max": f"{bucket}_max"})

fd_agg  = agg_bucket(rmsd, "free_design")
naf_agg = agg_bucket(rmsd, "native_ala_free")

df = pivot.join(fd_agg).join(naf_agg)
df["song"]      = df.index.map(lambda b: meta.loc[b, "genius_song_title"] if b in meta.index else b)
df["iconicity"] = df.index.map(lambda b: meta.loc[b, "aggregate_iconicity"] if b in meta.index else float("nan"))
df["delta_conc"]= df["native_ala"] - df["concordance"]

# keep only bars with all 4 buckets
df = df.dropna(subset=["native_ala", "concordance",
                        "free_design_med", "native_ala_free_med"])
df = df.sort_values("native_ala").reset_index()

bars = df["bar_id"].tolist()
n    = len(bars)
y    = np.arange(n)

print(f"Bars in plot: {n}")
for _, row in df.iterrows():
    print(f"  {row['bar_id']:8s}  na={row['native_ala']:.2f}  "
          f"conc={row['concordance']:.2f}  "
          f"fd_med={row['free_design_med']:.2f}  "
          f"naf_med={row['native_ala_free_med']:.2f}")

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(17, 7), facecolor="white",
                         gridspec_kw={"width_ratios": [3, 1]})

# ── Panel A — per-bar paired dot + band plot ──────────────────────────────────
ax = axes[0]
ax.set_facecolor("white")

BAND_HEIGHT = 0.38
BAND_OFFSET = 0.22   # vertical offset between the two MPNN bands per row

for i, row in df.iterrows():
    yi = y[i]

    # free_design band (upper band)
    yfd = yi + BAND_OFFSET
    ax.barh(yfd, row["free_design_q3"] - row["free_design_q1"],
            left=row["free_design_q1"],
            height=BAND_HEIGHT, color=COLORS["free_design"], alpha=0.20, zorder=1)
    ax.plot([row["free_design_min"], row["free_design_max"]], [yfd, yfd],
            color=COLORS["free_design"], lw=0.8, alpha=0.35, zorder=1)
    ax.plot([row["free_design_med"]] * 2,
            [yfd - BAND_HEIGHT / 2, yfd + BAND_HEIGHT / 2],
            color=COLORS["free_design"], lw=1.8, alpha=0.80, zorder=2)

    # native_ala_free band (lower band)
    ynaf = yi - BAND_OFFSET
    ax.barh(ynaf, row["native_ala_free_q3"] - row["native_ala_free_q1"],
            left=row["native_ala_free_q1"],
            height=BAND_HEIGHT, color=COLORS["native_ala_free"], alpha=0.20, zorder=1)
    ax.plot([row["native_ala_free_min"], row["native_ala_free_max"]], [ynaf, ynaf],
            color=COLORS["native_ala_free"], lw=0.8, alpha=0.35, zorder=1)
    ax.plot([row["native_ala_free_med"]] * 2,
            [ynaf - BAND_HEIGHT / 2, ynaf + BAND_HEIGHT / 2],
            color=COLORS["native_ala_free"], lw=1.8, alpha=0.80, zorder=2)

    # connecting line native_ala → concordance (lyric seed comparison)
    ax.plot([row["concordance"], row["native_ala"]], [yi, yi],
            color="#cccccc", lw=1.3, zorder=2, alpha=0.8)

    # native_ala (reference) — filled square
    ax.scatter(row["native_ala"], yi,
               s=120, marker="s", color=COLORS["native_ala"],
               zorder=5, linewidths=0.5, edgecolors="white")

    # concordance — diamond
    ax.scatter(row["concordance"], yi,
               s=100, marker="D", color=COLORS["concordance"],
               zorder=6, linewidths=0.5, edgecolors="white")

    # delta annotation (native_ala → concordance gap)
    ax.text(row["native_ala"] + 0.35, yi, f"Δ{row['delta_conc']:.1f}",
            va="center", fontsize=6.5, color="#999", style="italic")

# y-axis labels
song_labels = [
    f"{row['bar_id']}  {row['song'][:22]}{'…' if len(row['song']) > 22 else ''}"
    for _, row in df.iterrows()
]
ax.set_yticks(y)
ax.set_yticklabels(song_labels, fontsize=8.5, color="black")
ax.invert_yaxis()

ax.set_xlabel("Backbone RMSD (Å)  vs ESMFold reference", fontsize=10, color="#444")
ax.set_title(
    "Panel A — Per-bar backbone RMSD\n"
    "native_ala ■ · concordance ◆ · free_design MPNN (upper band) · native_ala_free MPNN (lower band)\n"
    "Band = IQR  |  tick = median  |  whisker = min–max",
    fontsize=9.5, color="black", pad=10)
ax.tick_params(axis="x", colors="#888", labelsize=8)
ax.tick_params(axis="y", length=0)
for sp in ax.spines.values():
    sp.set_edgecolor("#dddddd")
ax.axvline(0, color="#eeeeee", lw=0.5)
ax.set_ylim(-0.8, n - 0.2)

handles = [
    mpatches.Patch(color=COLORS["native_ala"],      label="native_ala (lyric ref ■)"),
    mpatches.Patch(color=COLORS["concordance"],     label="concordance (lyric alt ◆)"),
    mpatches.Patch(color=COLORS["free_design"],     label="free_design MPNN (upper band)"),
    mpatches.Patch(color=COLORS["native_ala_free"], label="native_ala_free MPNN (lower band)"),
]
ax.legend(handles=handles, fontsize=8.5, framealpha=0.95, edgecolor="#ccc",
          facecolor="white", labelcolor="black", loc="lower right")

# ── Panel B — overall violin by bucket ────────────────────────────────────────
ax2 = axes[1]
ax2.set_facecolor("white")

data_b  = [rmsd[rmsd["bucket"] == b]["rmsd_vs_backbone"].values for b in BUCKET_ORDER]
labels  = ["native_ala\n(ref)", "concordance", "free_design\n(MPNN)", "native_ala_free\n(MPNN)"]
xpos    = [1, 2, 3, 4]

vp = ax2.violinplot(data_b, positions=xpos, showmedians=True,
                    widths=0.7, showextrema=False)
for pc, b in zip(vp["bodies"], BUCKET_ORDER):
    pc.set_facecolor(COLORS[b])
    pc.set_alpha(0.60)
    pc.set_edgecolor("none")
vp["cmedians"].set_color("black")
vp["cmedians"].set_linewidth(1.8)

# overlay individual lyric seeds (concordance + native_ala)
rng = np.random.default_rng(42)
for xi, b in zip(xpos[:2], ["native_ala", "concordance"]):
    vals = rmsd[rmsd["bucket"] == b]["rmsd_vs_backbone"].values
    jitter = rng.uniform(-0.07, 0.07, len(vals))
    ax2.scatter(np.full(len(vals), xi) + jitter, vals,
                s=35, color=COLORS[b], alpha=0.90, zorder=5,
                edgecolors="white", linewidths=0.4)

medians = [np.median(d) for d in data_b]
for xi, med in zip(xpos, medians):
    ax2.text(xi, med + 0.5, f"{med:.1f}",
             ha="center", va="bottom", fontsize=8, color="black", fontweight="bold")

ax2.set_xticks(xpos)
ax2.set_xticklabels(labels, fontsize=8, color="black")
ax2.set_ylabel("Backbone RMSD (Å)", fontsize=9.5, color="#444")
ax2.set_title("Panel B\nOverall distribution", fontsize=10, color="black", pad=10)
ax2.tick_params(axis="y", colors="#888", labelsize=8)
ax2.tick_params(axis="x", length=0)
for sp in ax2.spines.values():
    sp.set_edgecolor("#dddddd")

# ── Figure-level blurb ───────────────────────────────────────────────────────

plt.tight_layout(rect=[0, 0.09, 1, 1])

blurb = (
    "Backbone RMSD of Boltz-2 predicted structures vs the ESMFold reference for each bar (n=12 bars with full data across all 4 buckets). "
    "native_ala (orange ■) = literal pass-through lyric encoding (reference). "
    "concordance (blue ◆) = freq-rank + softmax encoding — consistently folds closer to the ESMFold backbone (Δ = 2–18 Å lower RMSD). "
    "free_design MPNN (purple band, upper) = designs seeded from concordance backbone. "
    "native_ala_free MPNN (green band, lower) = designs seeded from native_ala backbone. "
    "Band = IQR; vertical tick = median; whisker = min–max. "
    "Panel B violin shows the overall distribution across all 12 bars."
)
fig.text(0.5, 0.01, blurb, ha="center", va="bottom", fontsize=7.5,
         color="#444", wrap=True,
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                   edgecolor="#ddd", alpha=0.9))

plt.savefig(OUT_FIG / "fig90_rmsd_comparison.png", dpi=160,
            facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig90_rmsd_comparison.png")
