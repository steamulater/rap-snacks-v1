"""
14_rmsd_comparison.py
---------------------
Fig 90 — Pairwise Boltz-structure RMSD across buckets, per bar.

Computes Cα RMSD (Kabsch alignment) between Boltz-2 predicted structures
for these pairs, across 12 bars:

  Pair 1: native_ala  vs  concordance         (lyric encoding comparison)
  Pair 2: native_ala  vs  native_ala_free     (MPNN rescue of na seed)
  Pair 3: concordance vs  free_design         (MPNN rescue of conc seed)
  Pair 4: free_design vs  native_ala_free     (MPNN strategy comparison)

All structures are Boltz-2 model_0 predictions.
For MPNN buckets (free_design, native_ala_free) RMSD is computed for every
design vs the corresponding lyric seed; median + IQR shown per bar.

Outputs:
  outputs/figures/fig90_rmsd_pairs.png   — 2×2 per-pair subplots
  outputs/figures/fig91_rmsd_combined.png — all 4 pairs in one panel
"""

import warnings
warnings.filterwarnings("ignore")

import io, zipfile, os
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

V3_ZIP   = ROOT / "results-20260407T081540Z-3-001.zip"
ENRICHED = DATA / "aggregated_lines_v2_enriched.csv"

PRED_BASE = ROOT / "outputs/boltz_validation/boltz_outputs/boltz_results_boltz_inputs/predictions"
NA_BASE   = ROOT / "boltz_results_boltz_inputs_native_ala/predictions"

COLORS = {
    "native_ala":      "#f97316",
    "concordance":     "#0ea5e9",
    "free_design":     "#a855f7",
    "native_ala_free": "#22c55e",
}

NA_MAP = {
    "bar_0":  "na0",  "bar_3":  "na3",  "bar_6":  "na6",  "bar_8":  "na8",
    "bar_9":  "na9",  "bar_11": "na11", "bar_13": "na13", "bar_17": "na17",
    "bar_27": "na27", "bar_32": "na32", "bar_46": "na46", "bar_77": "na77",
}

# ── Kabsch RMSD ───────────────────────────────────────────────────────────────

def parse_ca(pdb_source):
    """Return (N,3) Cα coords from a PDB file path or text string.
    Handles non-standard multi-char chain IDs by splitting on whitespace."""
    import re
    if isinstance(pdb_source, (str, Path)):
        lines = Path(pdb_source).read_text().splitlines()
    else:
        lines = pdb_source.decode().splitlines()
    coords = []
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        # atom name is at cols 12-16 — check for CA
        # also handle lines where chain ID shifts columns
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            # try splitting — some Boltz PDBs use multi-char chain IDs
            parts = line.split()
            if len(parts) >= 9 and parts[2] == "CA":
                try:
                    # x,y,z are the 3 floats after resnum/icode
                    # find them: last 6 tokens are x y z occ bfac [element]
                    floats = [float(p) for p in parts if re.match(r'^-?\d+\.\d+$', p)]
                    if len(floats) >= 3:
                        coords.append(floats[:3])
                except (ValueError, IndexError):
                    pass
            continue
        try:
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            # fallback for shifted columns (multi-char chain ID)
            parts = line.split()
            floats = [float(p) for p in parts if re.match(r'^-?\d+\.\d+$', p)]
            if len(floats) >= 3:
                coords.append(floats[:3])
    return np.array(coords)

def kabsch_rmsd(P, Q):
    """Kabsch-aligned Cα RMSD between two (N,3) arrays."""
    assert P.shape == Q.shape, f"Shape mismatch: {P.shape} vs {Q.shape}"
    p = P - P.mean(0)
    q = Q - Q.mean(0)
    H  = p.T @ q
    U, _, Vt = np.linalg.svd(H)
    d  = np.linalg.det(Vt.T @ U.T)
    D  = np.diag([1, 1, d])
    R  = Vt.T @ D @ U.T
    p_rot = p @ R.T
    return float(np.sqrt(((p_rot - q) ** 2).sum(1).mean()))

# ── Load structures ───────────────────────────────────────────────────────────

enriched = pd.read_csv(ENRICHED)
meta = (enriched[["bar_id", "genius_song_title", "aggregate_iconicity"]]
        .drop_duplicates("bar_id").set_index("bar_id"))

# Open v3 zip once
z = zipfile.ZipFile(V3_ZIP)
naf_zip_map = {}  # bar_id → list of zip paths (model_0 PDBs)
for name in z.namelist():
    if "boltz_outputs_naf/" in name and name.endswith("_model_0.pdb"):
        bar = "_".join(name.split("/")[-2].split("_")[:2])  # bar_XX
        naf_zip_map.setdefault(bar, []).append(name)

fd_paths_map = {}  # bar_id → list of PDB paths (model_0)
for item in PRED_BASE.iterdir():
    if "_free_" in item.name:
        bar = "_".join(item.name.split("_")[:2])
        pdb = item / f"{item.name}_model_0.pdb"
        if pdb.exists():
            fd_paths_map.setdefault(bar, []).append(pdb)

bars = sorted(NA_MAP.keys(), key=lambda b: int(b.split("_")[1]))

# ── Compute pairwise RMSDs ────────────────────────────────────────────────────

records = []

for bar in bars:
    na_id   = NA_MAP[bar]
    na_pdb  = NA_BASE / na_id / f"{na_id}_model_0.pdb"
    conc_pdb = PRED_BASE / f"{bar}_concordance" / f"{bar}_concordance_model_0.pdb"

    if not na_pdb.exists() or not conc_pdb.exists():
        print(f"  {bar}: missing na or conc PDB — skip")
        continue

    ca_na   = parse_ca(na_pdb)
    ca_conc = parse_ca(conc_pdb)
    assert len(ca_na) == len(ca_conc), f"{bar}: length mismatch na={len(ca_na)} conc={len(ca_conc)}"

    # Pair 1: na vs conc
    p1 = kabsch_rmsd(ca_na, ca_conc)

    # Pair 2: na vs each native_ala_free design
    naf_rmsds = []
    for zpath in naf_zip_map.get(bar, []):
        ca_naf = parse_ca(z.read(zpath))
        if len(ca_naf) == len(ca_na):
            naf_rmsds.append(kabsch_rmsd(ca_na, ca_naf))

    # Pair 3: conc vs each free_design
    fd_rmsds = []
    for pdb_path in fd_paths_map.get(bar, []):
        ca_fd = parse_ca(pdb_path)
        if len(ca_fd) == len(ca_conc):
            fd_rmsds.append(kabsch_rmsd(ca_conc, ca_fd))

    # Pair 4: free_design vs native_ala_free (conc-seed vs na-seed MPNN designs)
    #   compute median fd structure vs median naf structure as representative
    #   approximate: use pairwise between all fd and all naf, take median
    cross_rmsds = []
    for pdb_fd in fd_paths_map.get(bar, []):
        ca_fd = parse_ca(pdb_fd)
        if len(ca_fd) != len(ca_na):
            continue
        for zpath in naf_zip_map.get(bar, []):
            ca_naf = parse_ca(z.read(zpath))
            if len(ca_naf) == len(ca_fd):
                cross_rmsds.append(kabsch_rmsd(ca_fd, ca_naf))
        # break after a few fd to keep compute tractable
        if len(cross_rmsds) >= 50:
            break

    song = meta.loc[bar, "genius_song_title"] if bar in meta.index else bar
    records.append({
        "bar_id": bar,
        "song":   song,
        # pair 1
        "p1":     p1,
        # pair 2
        "p2_med": np.median(naf_rmsds) if naf_rmsds else np.nan,
        "p2_q1":  np.quantile(naf_rmsds, 0.25) if naf_rmsds else np.nan,
        "p2_q3":  np.quantile(naf_rmsds, 0.75) if naf_rmsds else np.nan,
        "p2_all": naf_rmsds,
        # pair 3
        "p3_med": np.median(fd_rmsds) if fd_rmsds else np.nan,
        "p3_q1":  np.quantile(fd_rmsds, 0.25) if fd_rmsds else np.nan,
        "p3_q3":  np.quantile(fd_rmsds, 0.75) if fd_rmsds else np.nan,
        "p3_all": fd_rmsds,
        # pair 4
        "p4_med": np.median(cross_rmsds) if cross_rmsds else np.nan,
        "p4_q1":  np.quantile(cross_rmsds, 0.25) if cross_rmsds else np.nan,
        "p4_q3":  np.quantile(cross_rmsds, 0.75) if cross_rmsds else np.nan,
    })
    print(f"  {bar:8s}  p1={p1:.2f}  "
          f"p2_med={records[-1]['p2_med']:.2f}  "
          f"p3_med={records[-1]['p3_med']:.2f}  "
          f"p4_med={records[-1]['p4_med']:.2f}")

z.close()
df = pd.DataFrame(records).sort_values("p1").reset_index(drop=True)
n  = len(df)
y  = np.arange(n)

# ── Plot helpers ──────────────────────────────────────────────────────────────

def style(ax, title):
    ax.set_facecolor("white")
    ax.set_title(title, fontsize=9.5, color="black", pad=8)
    ax.set_xlabel("Cα RMSD (Å)  [Kabsch-aligned, Boltz model_0]", fontsize=8.5, color="#555")
    ax.tick_params(axis="x", colors="#888", labelsize=7.5)
    ax.tick_params(axis="y", length=0)
    for sp in ax.spines.values():
        sp.set_edgecolor("#dddddd")

def song_label(row):
    s = row["song"]
    return f"{row['bar_id']}  {s[:20]}{'…' if len(s)>20 else ''}"

def draw_pair_single(ax, vals_a, vals_b, color_a, color_b, label_a, label_b,
                     is_dist_b=False):
    """vals_a = scalar per bar; vals_b = scalar or (med,q1,q3) per bar."""
    for i, row in df.iterrows():
        yi = y[i]
        va = vals_a[i]
        if is_dist_b:
            vb_med, vb_q1, vb_q3 = vals_b[i]
        else:
            vb_med = vals_b[i]; vb_q1 = vb_q3 = vb_med

        # connecting line
        ax.plot([va, vb_med], [yi, yi], color="#cccccc", lw=1.1, zorder=1, alpha=0.75)
        # IQR band for distribution
        if is_dist_b and not np.isnan(vb_q1):
            ax.barh(yi, vb_q3 - vb_q1, left=vb_q1, height=0.45,
                    color=color_b, alpha=0.18, zorder=2)
        # markers
        ax.scatter(va, yi, s=90, marker="s", color=color_a, zorder=5,
                   linewidths=0.5, edgecolors="white")
        ax.scatter(vb_med, yi, s=70, marker="D", color=color_b, zorder=6,
                   linewidths=0.5, edgecolors="white")

    labels = [song_label(r) for _, r in df.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8, color="black")
    ax.invert_yaxis()
    handles = [
        mpatches.Patch(color=color_a, label=f"{label_a} (■)"),
        mpatches.Patch(color=color_b, label=f"{label_b} (◆, {'median + IQR' if is_dist_b else 'single'})"),
    ]
    ax.legend(handles=handles, fontsize=7.5, framealpha=0.9,
              edgecolor="#ccc", facecolor="white", loc="lower right")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 90 — 2×2 per-pair subplots
# ══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(2, 2, figsize=(18, 12), facecolor="white")
axes = axes.flatten()

# Pair 1: na vs conc (single vs single)
ax = axes[0]
style(ax, "Pair 1 — native_ala  vs  concordance\n(same lyric, different encoding)")
draw_pair_single(ax,
    vals_a=df["p1"].values, vals_b=df["p1"].values,   # placeholder — draw manually
    color_a=COLORS["native_ala"], color_b=COLORS["concordance"],
    label_a="native_ala", label_b="concordance")

# redo pair 1 manually (both are single values, connect from na side to conc side)
ax.cla()
style(ax, "Pair 1 — native_ala  vs  concordance\n(same lyric, different encoding strategy)")
na_v2  = pd.read_csv(ROOT/"outputs/boltz_validation/boltz_rmsd.csv")
na_v2  = na_v2[na_v2["bucket"].isin(["concordance","native_ala"])]
pivot12 = na_v2.pivot_table(index="bar_id", columns="bucket", values="rmsd_vs_backbone")
# actually pair 1 IS the na vs conc Boltz-structure RMSD = df["p1"]
# let me show it as a single-value comparison: both points on the bar
for i, row in df.iterrows():
    yi = y[i]
    # p1 is the RMSD between the two structures — show as a single horizontal bar from 0
    ax.barh(yi, row["p1"], left=0, height=0.55,
            color="#888", alpha=0.25, zorder=1)
    ax.scatter(row["p1"], yi, s=90, marker="o",
               color="#555", zorder=5, linewidths=0.5, edgecolors="white")

labels = [song_label(r) for _, r in df.iterrows()]
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8, color="black")
ax.invert_yaxis()
ax.set_xlabel("Cα RMSD (Å)  [na structure vs conc structure, Boltz model_0]",
              fontsize=8.5, color="#555")
ax.legend(handles=[mpatches.Patch(color="#555", label="RMSD(native_ala, concordance)")],
          fontsize=8, facecolor="white", edgecolor="#ccc", loc="lower right")

# Pair 2: na vs native_ala_free
ax = axes[1]
style(ax, "Pair 2 — native_ala  vs  native_ala_free MPNN\n(how much does MPNN redesign the seed?)")
for i, row in df.iterrows():
    yi = y[i]
    ax.barh(yi, row["p2_q3"]-row["p2_q1"], left=row["p2_q1"],
            height=0.45, color=COLORS["native_ala_free"], alpha=0.20, zorder=2)
    ax.scatter(row["p2_med"], yi, s=75, marker="D",
               color=COLORS["native_ala_free"], zorder=5, linewidths=0.5, edgecolors="white")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8, color="black")
ax.invert_yaxis()
ax.legend(handles=[mpatches.Patch(color=COLORS["native_ala_free"],
                                   label="RMSD(native_ala, naf_mpnn) — median + IQR")],
          fontsize=8, facecolor="white", edgecolor="#ccc", loc="lower right")

# Pair 3: conc vs free_design
ax = axes[2]
style(ax, "Pair 3 — concordance  vs  free_design MPNN\n(how much does MPNN redesign the conc seed?)")
for i, row in df.iterrows():
    yi = y[i]
    ax.barh(yi, row["p3_q3"]-row["p3_q1"], left=row["p3_q1"],
            height=0.45, color=COLORS["free_design"], alpha=0.20, zorder=2)
    ax.scatter(row["p3_med"], yi, s=75, marker="D",
               color=COLORS["free_design"], zorder=5, linewidths=0.5, edgecolors="white")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8, color="black")
ax.invert_yaxis()
ax.legend(handles=[mpatches.Patch(color=COLORS["free_design"],
                                   label="RMSD(concordance, fd_mpnn) — median + IQR")],
          fontsize=8, facecolor="white", edgecolor="#ccc", loc="lower right")

# Pair 4: free_design vs native_ala_free
ax = axes[3]
style(ax, "Pair 4 — free_design MPNN  vs  native_ala_free MPNN\n(concordance-seeded vs native_ala-seeded designs)")
for i, row in df.iterrows():
    yi = y[i]
    ax.barh(yi, row["p4_q3"]-row["p4_q1"], left=row["p4_q1"],
            height=0.45, color="#94a3b8", alpha=0.25, zorder=2)
    ax.scatter(row["p4_med"], yi, s=75, marker="D",
               color="#475569", zorder=5, linewidths=0.5, edgecolors="white")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8, color="black")
ax.invert_yaxis()
ax.legend(handles=[mpatches.Patch(color="#475569",
                                   label="RMSD(fd_mpnn, naf_mpnn) — median + IQR")],
          fontsize=8, facecolor="white", edgecolor="#ccc", loc="lower right")

plt.suptitle("Pairwise Boltz-structure Cα RMSD — 4 pairs × 12 bars",
             fontsize=13, color="black", y=1.01)
plt.tight_layout(rect=[0, 0.05, 1, 1])

blurb90 = (
    "Cα RMSD (Kabsch-aligned) between Boltz-2 predicted structures. All values are "
    "structure-vs-structure (not vs ESMFold backbone). model_0 used for lyric seeds; "
    "median + IQR shown for MPNN design distributions. "
    "Pair 1: how different are the two lyric-seed structures from each other. "
    "Pairs 2–3: how much MPNN redesigns each seed. "
    "Pair 4: how structurally similar are the two MPNN strategy outputs."
)
fig.text(0.5, 0.01, blurb90, ha="center", va="bottom", fontsize=7.5, color="#444",
         wrap=True, bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                               edgecolor="#ddd", alpha=0.9))
fig.patch.set_facecolor("white")
plt.savefig(OUT_FIG / "fig90_rmsd_pairs.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig90_rmsd_pairs.png")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 91 — All 4 pairs in one panel (sorted by pair 1 RMSD)
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(13, 8), facecolor="white")
ax.set_facecolor("white")

PAIR_COLORS = ["#555555", COLORS["native_ala_free"], COLORS["free_design"], "#475569"]
PAIR_LABELS = [
    "Pair 1: na ↔ conc",
    "Pair 2: na ↔ naf MPNN",
    "Pair 3: conc ↔ fd MPNN",
    "Pair 4: fd MPNN ↔ naf MPNN",
]
MARKERS     = ["o", "D", "s", "^"]
OFFSETS     = [-0.28, -0.09, 0.09, 0.28]  # vertical jitter per pair

for i, row in df.iterrows():
    yi = y[i]
    pairs_med = [row["p1"], row["p2_med"], row["p3_med"], row["p4_med"]]
    pairs_q1  = [row["p1"], row["p2_q1"],  row["p3_q1"],  row["p4_q1"]]
    pairs_q3  = [row["p1"], row["p2_q3"],  row["p3_q3"],  row["p4_q3"]]

    for j, (med, q1, q3, col, mk, off) in enumerate(
            zip(pairs_med, pairs_q1, pairs_q3, PAIR_COLORS, MARKERS, OFFSETS)):
        if np.isnan(med):
            continue
        yj = yi + off
        # IQR band
        if j > 0:
            ax.barh(yj, q3 - q1, left=q1, height=0.15,
                    color=col, alpha=0.22, zorder=1)
        ax.scatter(med, yj, s=55, marker=mk, color=col, zorder=5,
                   linewidths=0.4, edgecolors="white", alpha=0.92)

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=8.5, color="black")
ax.invert_yaxis()
ax.set_xlabel("Cα RMSD (Å)  [Boltz-2 structure vs structure, Kabsch-aligned]",
              fontsize=10, color="#444")
ax.set_title("All 4 pairwise RMSD comparisons — Boltz-2 structures × 12 bars",
             fontsize=11, color="black", pad=10)
ax.tick_params(axis="x", colors="#888", labelsize=8)
ax.tick_params(axis="y", length=0)
for sp in ax.spines.values():
    sp.set_edgecolor("#dddddd")

handles = [mpatches.Patch(color=c, label=l)
           for c, l in zip(PAIR_COLORS, PAIR_LABELS)]
ax.legend(handles=handles, fontsize=8.5, framealpha=0.95,
          edgecolor="#ccc", facecolor="white", labelcolor="black",
          loc="lower right")

plt.tight_layout(rect=[0, 0.08, 1, 1])

blurb91 = (
    "All 4 pairwise Boltz Cα RMSD comparisons overlaid per bar (sorted by Pair 1 RMSD). "
    "Circle = Pair 1 (na vs conc): structural divergence from encoding strategy alone. "
    "Diamond = Pair 2 (na vs naf MPNN): how far MPNN redesigns from na seed. "
    "Square = Pair 3 (conc vs fd MPNN): how far MPNN redesigns from conc seed. "
    "Triangle = Pair 4 (fd vs naf MPNN): structural divergence between the two MPNN strategy outputs. "
    "IQR band shown for distribution pairs."
)
fig.text(0.5, 0.01, blurb91, ha="center", va="bottom", fontsize=7.5, color="#444",
         wrap=True, bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9",
                               edgecolor="#ddd", alpha=0.9))
fig.patch.set_facecolor("white")
plt.savefig(OUT_FIG / "fig91_rmsd_combined.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig91_rmsd_combined.png")
