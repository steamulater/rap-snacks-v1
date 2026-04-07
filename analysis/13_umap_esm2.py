"""
13_umap_esm2.py
---------------
ESM-2 embeddings → UMAP variants (white background).

Outputs:
  fig83_umap_all_buckets.png       — all 4 buckets, colour by bucket
  fig84_umap_native_ala.png        — native_ala + native_ala_free only
  fig85_umap_concordance.png       — concordance + free_design only
  fig86_umap_barbie_dangerous.png  — bar_32 only, all 4 buckets
  fig87_umap_lyric_seeds.png       — concordance + native_ala only (no MPNN)
  fig88_umap_plddt.png             — all seqs coloured by Boltz pLDDT
  fig89_umap_bar_XX.png × 37      — per-bar, all 4 buckets

  outputs/embeddings/esm2_embeddings.csv  (cached)
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial import ConvexHull
from pathlib import Path
import umap as umap_lib

ROOT    = Path(__file__).parent.parent
DATA    = ROOT / "data"
OUT_FIG = ROOT / "outputs/figures"
OUT_EMB = ROOT / "outputs/embeddings"
OUT_FIG.mkdir(exist_ok=True)
OUT_EMB.mkdir(exist_ok=True)

ENRICHED   = DATA / "aggregated_lines_v2_enriched.csv"
CANDIDATES = DATA / "phase2_candidates.csv"
FREE_CSV   = ROOT / "outputs/proteinmpnn_free/filtered_results.csv"
NAF_CSV    = ROOT / "outputs/proteinmpnn_native_ala_free/filtered_results.csv"
V2_CSV     = ROOT / "outputs/boltz_validation/boltz_validation_results.csv"
FS2_CSV    = ROOT / "outputs/foldseek_phase2/foldseek_phase2_summary.csv"
EMB_CACHE  = OUT_EMB / "esm2_embeddings.csv"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

COLORS = {
    "concordance":     "#0ea5e9",   # sky blue
    "native_ala":      "#f97316",   # orange
    "free_design":     "#a855f7",   # purple
    "native_ala_free": "#22c55e",   # green
}
BUCKET_ORDER = ["concordance", "native_ala", "free_design", "native_ala_free"]

# ── Load metadata ──────────────────────────────────────────────────────────

enriched   = pd.read_csv(ENRICHED)
candidates = pd.read_csv(CANDIDATES)
bar_ids    = list(candidates["bar_id"])

meta = (enriched[enriched["bar_id"].isin(bar_ids)]
        [["bar_id", "genius_song_title", "aggregate_iconicity"]]
        .drop_duplicates("bar_id")
        .set_index("bar_id"))

fs_hits = {}
if FS2_CSV.exists():
    fs2 = pd.read_csv(FS2_CSV)
    for _, row in fs2.iterrows():
        fs_hits[(row["bar_id"], row["bucket"])] = int(row["n_hits"])

boltz_plddt = {}
if V2_CSV.exists():
    v2 = pd.read_csv(V2_CSV)
    for _, row in v2.iterrows():
        boltz_plddt[(row["bar_id"], row["bucket"])] = row["boltz_plddt"]

# ── Collect sequences ──────────────────────────────────────────────────────

rows = []

def add(bar_id, bucket, seq, name):
    seq = str(seq).upper().strip()
    if not seq or not set(seq) <= VALID_AA:
        return
    song  = meta.loc[bar_id, "genius_song_title"] if bar_id in meta.index else bar_id
    icon  = meta.loc[bar_id, "aggregate_iconicity"] if bar_id in meta.index else float("nan")
    rows.append({
        "name":      name,
        "bar_id":    bar_id,
        "bucket":    bucket,
        "song":      song,
        "iconicity": icon,
        "sequence":  seq,
        "n_hits":    fs_hits.get((bar_id, bucket), 0),
        "boltz_plddt": boltz_plddt.get((bar_id, bucket), float("nan")),
    })

for bar_id in bar_ids:
    r = enriched[enriched["bar_id"] == bar_id].iloc[0]
    add(bar_id, "concordance", r.get("fasta_seq_concordance"),    f"{bar_id}__concordance")
    add(bar_id, "native_ala",  r.get("fasta_seq_native_alanine"), f"{bar_id}__native_ala")

if FREE_CSV.exists():
    fd = pd.read_csv(FREE_CSV)
    for i, row in fd.iterrows():
        if row.get("bar_id") not in bar_ids: continue
        add(row["bar_id"], "free_design", row["sequence"],
            f"{row['bar_id']}__fd_{int(row.get('design_idx', i)):03d}")

if NAF_CSV.exists():
    naf = pd.read_csv(NAF_CSV)
    for i, row in naf.iterrows():
        if row.get("bar_id") not in bar_ids: continue
        add(row["bar_id"], "native_ala_free", row["sequence"],
            f"{row['bar_id']}__naf_{int(row.get('design_idx', i)):03d}")

df = pd.DataFrame(rows).drop_duplicates("name").reset_index(drop=True)
print(f"Total sequences: {len(df)}")
for bucket, g in df.groupby("bucket"):
    print(f"  {bucket:22s}  n={len(g)}")

# ── ESM-2 embeddings ────────────────────────────────────────────────────────

if EMB_CACHE.exists():
    print(f"\nLoading cached embeddings...")
    emb_df  = pd.read_csv(EMB_CACHE, index_col=0)
    missing = set(df["name"]) - set(emb_df.index)
    if missing:
        print(f"  {len(missing)} sequences missing from cache — recomputing all")
        EMB_CACHE.unlink()

if not EMB_CACHE.exists():
    print("\nComputing ESM-2 embeddings (esm2_t6_8M_UR50D)...")
    import esm as esmlib, torch
    model, alphabet = esmlib.pretrained.esm2_t6_8M_UR50D()
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    BATCH = 32
    all_embs, names = [], list(df["name"])
    seqs = list(df["sequence"])
    for i in range(0, len(seqs), BATCH):
        data = [(n, s[:1022]) for n, s in zip(names[i:i+BATCH], seqs[i:i+BATCH])]
        _, _, tokens = batch_converter(data)
        with torch.no_grad():
            reps = model(tokens, repr_layers=[6])["representations"][6]
        for j, s in enumerate(seqs[i:i+BATCH]):
            L = min(len(s), 1022)
            all_embs.append(reps[j, 1:L+1].mean(0).cpu().numpy())
        if (i // BATCH) % 5 == 0:
            print(f"  {min(i+BATCH, len(seqs))}/{len(seqs)}")
    emb_matrix = np.stack(all_embs)
    emb_df = pd.DataFrame(emb_matrix, index=names)
    emb_df.to_csv(EMB_CACHE)
    print(f"Saved {EMB_CACHE.name}  shape={emb_matrix.shape}")

emb_matrix = emb_df.loc[df["name"]].values

# ── UMAP (full) ─────────────────────────────────────────────────────────────

print("\nRunning UMAP (full, 1286 seqs)...")
reducer = umap_lib.UMAP(n_neighbors=15, min_dist=0.1, metric="cosine",
                        random_state=42, verbose=False)
coords = reducer.fit_transform(emb_matrix)
df["umap_x"] = coords[:, 0]
df["umap_y"] = coords[:, 1]
print("Done.")

# ── Helpers ──────────────────────────────────────────────────────────────────

def draw_hull(ax, pts, color, alpha=0.06, lw=0.8):
    if len(pts) < 3:
        return
    try:
        hull = ConvexHull(pts)
        verts = pts[hull.vertices]
        verts = np.vstack([verts, verts[0]])
        ax.fill(verts[:, 0], verts[:, 1], color=color, alpha=alpha, zorder=1)
        ax.plot(verts[:, 0], verts[:, 1], color=color, alpha=0.4, lw=lw, zorder=2)
    except Exception:
        pass

def label_bars(ax, sub_df, bar_ids_to_label, fontsize=6.5, color="black"):
    for bar_id in bar_ids_to_label:
        g = sub_df[sub_df["bar_id"] == bar_id]
        if g.empty: continue
        cx, cy = g["umap_x"].mean(), g["umap_y"].mean()
        song = meta.loc[bar_id, "genius_song_title"] if bar_id in meta.index else bar_id
        short = song[:20] + "…" if len(song) > 20 else song
        ax.text(cx, cy, f"{bar_id}\n{short}", ha="center", va="center",
                fontsize=fontsize, color=color, alpha=0.9,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="#cccccc", alpha=0.85))

def style_ax(ax, title):
    ax.set_facecolor("white")
    ax.set_title(title, fontsize=11, color="black", pad=8)
    ax.set_xlabel("UMAP 1", fontsize=9, color="#555")
    ax.set_ylabel("UMAP 2", fontsize=9, color="#555")
    ax.tick_params(colors="#999", labelsize=7)
    for s in ax.spines.values():
        s.set_edgecolor("#dddddd")

def legend(ax, buckets):
    handles = [mpatches.Patch(color=COLORS[b], label=b.replace("_", " ")) for b in buckets]
    ax.legend(handles=handles, fontsize=8, framealpha=0.9,
              edgecolor="#cccccc", facecolor="white", labelcolor="black")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 83 — All 4 buckets
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(14, 10), facecolor="white")
style_ax(ax, "ESM-2 UMAP — all 4 buckets  (concordance · native_ala · free_design · native_ala_free)")

for bar_id in bar_ids:
    g = df[df["bar_id"] == bar_id]
    pts = g[["umap_x", "umap_y"]].values
    draw_hull(ax, pts, color="#aaaaaa", alpha=0.05)

for bucket in BUCKET_ORDER:
    g = df[df["bucket"] == bucket]
    is_lyric = bucket in ("concordance", "native_ala")
    ax.scatter(g["umap_x"], g["umap_y"],
               c=COLORS[bucket],
               s=55 if is_lyric else 10,
               alpha=0.9 if is_lyric else 0.35,
               zorder=4 if is_lyric else 2,
               linewidths=0.4 if is_lyric else 0,
               edgecolors="white" if is_lyric else "none",
               label=bucket)

label_bars(ax, df, bar_ids)
legend(ax, BUCKET_ORDER)
plt.tight_layout()
plt.savefig(OUT_FIG / "fig83_umap_all_buckets.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig83_umap_all_buckets.png")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 84 — native_ala + native_ala_free only
# ══════════════════════════════════════════════════════════════════════════════

sub = df[df["bucket"].isin(["native_ala", "native_ala_free"])].copy()

# Re-run UMAP on subset for cleaner layout
print("Running UMAP (native_ala subset)...")
emb_sub = emb_df.loc[sub["name"]].values
coords_sub = umap_lib.UMAP(n_neighbors=12, min_dist=0.08, metric="cosine",
                           random_state=42, verbose=False).fit_transform(emb_sub)
sub["umap_x"] = coords_sub[:, 0]
sub["umap_y"] = coords_sub[:, 1]

fig, ax = plt.subplots(figsize=(14, 10), facecolor="white")
style_ax(ax, "native_ala sequences + native_ala_free MPNN designs — ESM-2 UMAP")

for bar_id in bar_ids:
    g = sub[sub["bar_id"] == bar_id]
    pts = g[["umap_x", "umap_y"]].values
    draw_hull(ax, pts, color=COLORS["native_ala"], alpha=0.07)

ax.scatter(sub[sub["bucket"] == "native_ala_free"]["umap_x"],
           sub[sub["bucket"] == "native_ala_free"]["umap_y"],
           c=COLORS["native_ala_free"], s=12, alpha=0.4, zorder=2,
           linewidths=0, label="native_ala_free (MPNN)")

na = sub[sub["bucket"] == "native_ala"]
ax.scatter(na["umap_x"], na["umap_y"],
           c=COLORS["native_ala"], s=90, alpha=0.95, zorder=5,
           linewidths=0.6, edgecolors="white", label="native_ala (lyric seed)")

label_bars(ax, sub, bar_ids)
legend(ax, ["native_ala", "native_ala_free"])
plt.tight_layout()
plt.savefig(OUT_FIG / "fig84_umap_native_ala.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig84_umap_native_ala.png")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 85 — concordance + free_design only
# ══════════════════════════════════════════════════════════════════════════════

sub = df[df["bucket"].isin(["concordance", "free_design"])].copy()

print("Running UMAP (concordance subset)...")
emb_sub = emb_df.loc[sub["name"]].values
coords_sub = umap_lib.UMAP(n_neighbors=12, min_dist=0.08, metric="cosine",
                           random_state=42, verbose=False).fit_transform(emb_sub)
sub["umap_x"] = coords_sub[:, 0]
sub["umap_y"] = coords_sub[:, 1]

fig, ax = plt.subplots(figsize=(14, 10), facecolor="white")
style_ax(ax, "concordance sequences + free_design MPNN designs — ESM-2 UMAP")

for bar_id in bar_ids:
    g = sub[sub["bar_id"] == bar_id]
    pts = g[["umap_x", "umap_y"]].values
    draw_hull(ax, pts, color=COLORS["concordance"], alpha=0.07)

ax.scatter(sub[sub["bucket"] == "free_design"]["umap_x"],
           sub[sub["bucket"] == "free_design"]["umap_y"],
           c=COLORS["free_design"], s=12, alpha=0.4, zorder=2,
           linewidths=0, label="free_design (MPNN)")

conc = sub[sub["bucket"] == "concordance"]
ax.scatter(conc["umap_x"], conc["umap_y"],
           c=COLORS["concordance"], s=90, alpha=0.95, zorder=5,
           linewidths=0.6, edgecolors="white", label="concordance (lyric seed)")

label_bars(ax, sub, bar_ids)
legend(ax, ["concordance", "free_design"])
plt.tight_layout()
plt.savefig(OUT_FIG / "fig85_umap_concordance.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig85_umap_concordance.png")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 86 — bar_32 (Barbie Dangerous) only
# ══════════════════════════════════════════════════════════════════════════════

sub = df[df["bar_id"] == "bar_32"].copy()
print(f"\nbar_32 sequences: {len(sub)} ({sub['bucket'].value_counts().to_dict()})")

print("Running UMAP (bar_32 only)...")
emb_sub = emb_df.loc[sub["name"]].values
# fewer neighbors for small set
coords_sub = umap_lib.UMAP(n_neighbors=8, min_dist=0.05, metric="cosine",
                           random_state=42, verbose=False).fit_transform(emb_sub)
sub["umap_x"] = coords_sub[:, 0]
sub["umap_y"] = coords_sub[:, 1]

fig, ax = plt.subplots(figsize=(12, 9), facecolor="white")
style_ax(ax, "bar_32 — Barbie Dangerous — all 4 buckets  (ESM-2 UMAP)")

for bucket in BUCKET_ORDER:
    g = sub[sub["bucket"] == bucket]
    if g.empty: continue
    is_lyric = bucket in ("concordance", "native_ala")
    ax.scatter(g["umap_x"], g["umap_y"],
               c=COLORS[bucket],
               s=120 if is_lyric else 25,
               alpha=0.95 if is_lyric else 0.55,
               zorder=5 if is_lyric else 3,
               linewidths=0.8 if is_lyric else 0,
               edgecolors="white" if is_lyric else "none",
               label=bucket.replace("_", " "))
    # label the single lyric sequences
    if is_lyric:
        for _, row in g.iterrows():
            ax.annotate(bucket.replace("_", " "),
                        (row["umap_x"], row["umap_y"]),
                        xytext=(8, 6), textcoords="offset points",
                        fontsize=8, color=COLORS[bucket],
                        fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color=COLORS[bucket],
                                        lw=0.8, alpha=0.6))

song = meta.loc["bar_32", "genius_song_title"] if "bar_32" in meta.index else "Barbie Dangerous"
ax.set_title(f"bar_32 — {song} — all 4 buckets  (ESM-2 UMAP)\n"
             f"concordance (1) · native_ala (1) · free_design ({len(sub[sub['bucket']=='free_design'])}) "
             f"· native_ala_free ({len(sub[sub['bucket']=='native_ala_free'])})",
             fontsize=10, color="black", pad=8)
legend(ax, [b for b in BUCKET_ORDER if b in sub["bucket"].values])
plt.tight_layout()
plt.savefig(OUT_FIG / "fig86_umap_barbie_dangerous.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig86_umap_barbie_dangerous.png")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 87 — lyric seeds only (concordance + native_ala, no MPNN)
# ══════════════════════════════════════════════════════════════════════════════

sub87 = df[df["bucket"].isin(["concordance", "native_ala"])].copy()

print(f"\nRunning UMAP (lyric seeds only, n={len(sub87)})...")
emb_sub87 = emb_df.loc[sub87["name"]].values
coords87 = umap_lib.UMAP(n_neighbors=12, min_dist=0.12, metric="cosine",
                         random_state=42, verbose=False).fit_transform(emb_sub87)
sub87["umap_x"] = coords87[:, 0]
sub87["umap_y"] = coords87[:, 1]

fig, ax = plt.subplots(figsize=(14, 10), facecolor="white")
style_ax(ax, "ESM-2 UMAP — lyric seeds only  (concordance vs native_ala, n=74)")

for bucket in ["concordance", "native_ala"]:
    g = sub87[sub87["bucket"] == bucket]
    ax.scatter(g["umap_x"], g["umap_y"],
               c=COLORS[bucket], s=110, alpha=0.92, zorder=4,
               linewidths=0.6, edgecolors="white",
               label=bucket.replace("_", " "))

label_bars(ax, sub87, bar_ids, fontsize=6)
legend(ax, ["concordance", "native_ala"])
plt.tight_layout(rect=[0, 0.06, 1, 1])

blurb87 = (
    "Each point is a single bar's lyric-derived protein seed — concordance (blue, freq-rank + softmax draw) "
    "vs native_ala (orange, literal pass-through + Ala for non-standard residues). "
    "No MPNN designs included. Clusters reflect shared amino-acid composition driven by rap phonology."
)
fig.text(0.5, 0.01, blurb87, ha="center", va="bottom", fontsize=7.5,
         color="#444", wrap=True,
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9", edgecolor="#ddd", alpha=0.9))

plt.savefig(OUT_FIG / "fig87_umap_lyric_seeds.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig87_umap_lyric_seeds.png")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 88 — all sequences coloured by Boltz pLDDT
# ══════════════════════════════════════════════════════════════════════════════

# Build per-sequence pLDDT lookup from V2_CSV (concordance / native_ala / free_design)
# V2_CSV names are "bar_6_concordance" — join on (bar_id, bucket) instead to avoid
# name-format mismatches with df names ("bar_6__concordance")
seq_plddt = {}
if V2_CSV.exists():
    v2_seq = pd.read_csv(V2_CSV)
    for _, row in v2_seq.iterrows():
        seq_plddt[(row["bar_id"], row["bucket"])] = row["boltz_plddt"]

# Map onto df — native_ala_free will stay NaN (no Boltz run yet)
df["plddt_val"] = df.apply(
    lambda r: seq_plddt.get((r["bar_id"], r["bucket"]), float("nan")), axis=1)

import matplotlib.cm as cm
from matplotlib.colors import Normalize

vmin, vmax = 0.30, 0.95
norm = Normalize(vmin=vmin, vmax=vmax)
cmap = cm.RdYlGn

fig, ax = plt.subplots(figsize=(14, 10), facecolor="white")
style_ax(ax, "ESM-2 UMAP — all sequences coloured by Boltz pLDDT  (n=1286)")

# Grey background for sequences with no pLDDT (native_ala_free)
no_plddt = df[df["plddt_val"].isna()]
ax.scatter(no_plddt["umap_x"], no_plddt["umap_y"],
           c="#cccccc", s=8, alpha=0.25, zorder=1,
           linewidths=0, label="native_ala_free (pLDDT not yet run)")

# Colour all other sequences by pLDDT — draw each bucket separately to control size/alpha
SIZE_MAP  = {"concordance": 90, "native_ala": 90, "free_design": 14}
ALPHA_MAP = {"concordance": 0.95, "native_ala": 0.95, "free_design": 0.55}
LW_MAP    = {"concordance": 0.6, "native_ala": 0.6, "free_design": 0}

has_plddt = df[df["plddt_val"].notna()]
sc = None
for bucket in ["free_design", "concordance", "native_ala"]:
    g = has_plddt[has_plddt["bucket"] == bucket]
    if g.empty:
        continue
    sc = ax.scatter(g["umap_x"], g["umap_y"],
                    c=g["plddt_val"], cmap=cmap, norm=norm,
                    s=SIZE_MAP[bucket], alpha=ALPHA_MAP[bucket],
                    linewidths=LW_MAP[bucket], edgecolors="white",
                    zorder=3 if bucket == "free_design" else 5)

if sc is None:
    import matplotlib.cm as cm2
    sc = ax.scatter([], [], c=[], cmap=cmap, norm=norm)  # dummy for colorbar
cbar = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label("Boltz pLDDT", fontsize=9, color="#555")
cbar.ax.tick_params(labelsize=7, colors="#555")

ax.legend(handles=[mpatches.Patch(color="#cccccc", label="native_ala_free (pLDDT unavailable)")],
          fontsize=8, framealpha=0.9, edgecolor="#ccc", facecolor="white")

plt.tight_layout(rect=[0, 0.07, 1, 1])

blurb88 = (
    "All 1286 sequences overlaid on the shared UMAP embedding, coloured by Boltz-2 mean pLDDT "
    "(green = high confidence, red = low). Lyric seeds (large dots) and free MPNN designs (small dots) "
    "share the same embedding space. Grey = native_ala_free MPNN designs (Boltz run pending). "
    "Higher pLDDT concentrates in the upper clusters, away from the disordered lyric-seed cloud."
)
fig.text(0.5, 0.01, blurb88, ha="center", va="bottom", fontsize=7.5,
         color="#444", wrap=True,
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9", edgecolor="#ddd", alpha=0.9))

plt.savefig(OUT_FIG / "fig88_umap_plddt.png", dpi=160, facecolor="white", bbox_inches="tight")
plt.close()
print("Saved fig88_umap_plddt.png")

# ══════════════════════════════════════════════════════════════════════════════
# Fig 89 — per-bar UMAPs  (one plot per bar, all 4 buckets)
# ══════════════════════════════════════════════════════════════════════════════

print("\nGenerating per-bar UMAPs (fig89_bar_XX)...")

for bar_id in bar_ids:
    sub_b = df[df["bar_id"] == bar_id].copy()
    if len(sub_b) < 4:
        print(f"  {bar_id}: skipping (only {len(sub_b)} seqs)")
        continue

    nn = min(8, len(sub_b) - 1)
    emb_b = emb_df.loc[sub_b["name"]].values
    coords_b = umap_lib.UMAP(n_neighbors=nn, min_dist=0.05, metric="cosine",
                             random_state=42, verbose=False).fit_transform(emb_b)
    sub_b["umap_x"] = coords_b[:, 0]
    sub_b["umap_y"] = coords_b[:, 1]

    song  = meta.loc[bar_id, "genius_song_title"] if bar_id in meta.index else bar_id
    icon  = meta.loc[bar_id, "aggregate_iconicity"] if bar_id in meta.index else float("nan")
    bucks = sub_b["bucket"].value_counts().to_dict()

    fig, ax = plt.subplots(figsize=(10, 8), facecolor="white")
    ax.set_facecolor("white")

    for bucket in BUCKET_ORDER:
        g = sub_b[sub_b["bucket"] == bucket]
        if g.empty: continue
        is_lyric = bucket in ("concordance", "native_ala")
        ax.scatter(g["umap_x"], g["umap_y"],
                   c=COLORS[bucket],
                   s=130 if is_lyric else 20,
                   alpha=0.95 if is_lyric else 0.55,
                   zorder=5 if is_lyric else 3,
                   linewidths=0.8 if is_lyric else 0,
                   edgecolors="white" if is_lyric else "none",
                   label=f"{bucket.replace('_',' ')} (n={len(g)})")
        if is_lyric:
            for _, row in g.iterrows():
                ax.annotate(bucket.replace("_", " "),
                            (row["umap_x"], row["umap_y"]),
                            xytext=(8, 5), textcoords="offset points",
                            fontsize=8, color=COLORS[bucket], fontweight="bold",
                            arrowprops=dict(arrowstyle="-", color=COLORS[bucket],
                                            lw=0.7, alpha=0.6))

    title_str = (f"{bar_id} — {song}\n"
                 f"iconicity={icon:.3f}  |  "
                 + "  ·  ".join(f"{b.replace('_',' ')} {bucks.get(b,0)}" for b in BUCKET_ORDER
                                if bucks.get(b, 0) > 0))
    ax.set_title(title_str, fontsize=9.5, color="black", pad=8)
    ax.set_xlabel("UMAP 1", fontsize=9, color="#555")
    ax.set_ylabel("UMAP 2", fontsize=9, color="#555")
    ax.tick_params(colors="#999", labelsize=7)
    for s in ax.spines.values():
        s.set_edgecolor("#dddddd")

    present = [b for b in BUCKET_ORDER if b in sub_b["bucket"].values]
    handles = [mpatches.Patch(color=COLORS[b], label=f"{b.replace('_',' ')} (n={bucks.get(b,0)})")
               for b in present]
    ax.legend(handles=handles, fontsize=8, framealpha=0.9,
              edgecolor="#cccccc", facecolor="white", labelcolor="black")

    plt.tight_layout(rect=[0, 0.09, 1, 1])

    blurb89 = (
        f"Per-bar UMAP for {bar_id} ({song}). "
        "Lyric seeds (concordance/native_ala, large dots) anchor the space; "
        "MPNN designs (small dots) scatter around them. "
        "Tight clusters = MPNN designs converging on a shared fold family. "
        "Spread = diverse sequence space exploration."
    )
    fig.text(0.5, 0.01, blurb89, ha="center", va="bottom", fontsize=7,
             color="#444", wrap=True,
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#f9f9f9", edgecolor="#ddd", alpha=0.9))

    bar_num = bar_id.replace("bar_", "")
    fname = f"fig89_umap_bar_{bar_num.zfill(2)}.png"
    plt.savefig(OUT_FIG / fname, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close()
    print(f"  Saved {fname}")

print("\nAll done.")
print("  fig83_umap_all_buckets.png")
print("  fig84_umap_native_ala.png")
print("  fig85_umap_concordance.png")
print("  fig86_umap_barbie_dangerous.png")
print("  fig87_umap_lyric_seeds.png")
print("  fig88_umap_plddt.png")
print("  fig89_umap_bar_XX.png × 37")
