"""
13_umap_esm2.py
---------------
ESM-2 embeddings → UMAP for the 4 core buckets:
  concordance, native_ala, free_design, native_ala_free

One dot per sequence. Layout:
  - Colour   : bucket
  - Hull     : bar_id (convex hull around each bar's designs)
  - Dot size : FoldSeek hit count (concordance/native_ala only; designs = 40)
  - Labels   : bar_id + song title on hull centroid

Panel A : coloured by bucket
Panel B : same layout, coloured by Boltz pLDDT

Output:
  outputs/figures/fig83_umap_esm2.png
  outputs/embeddings/esm2_embeddings.csv   (cached — skip recompute on re-run)
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from scipy.spatial import ConvexHull
from pathlib import Path

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
FOLDSEEK   = ROOT / "outputs/foldseek/foldseek_hits.csv"
FS2_CSV    = ROOT / "outputs/foldseek_phase2/foldseek_phase2_summary.csv"
EMB_CACHE  = OUT_EMB / "esm2_embeddings.csv"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

DARK_BG = "#0e0e0e"
COLORS  = {
    "concordance":     "#00d4ff",
    "native_ala":      "#ff9500",
    "free_design":     "#a855f7",
    "native_ala_free": "#39d353",
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

# FoldSeek hit counts (phase 2 summary if available, else phase 1)
fs_hits = {}
if FS2_CSV.exists():
    fs2 = pd.read_csv(FS2_CSV)
    for _, row in fs2.iterrows():
        fs_hits[(row["bar_id"], row["bucket"])] = int(row["n_hits"])
elif FOLDSEEK.exists():
    fs1 = pd.read_csv(FOLDSEEK)
    for bar_id, g in fs1.groupby("bar_id"):
        n = len(g[g["foldseek_result"] != "no_hits_novel"])
        fs_hits[(bar_id, "concordance")] = n

# Boltz pLDDT lookup
boltz_plddt = {}
if V2_CSV.exists():
    v2 = pd.read_csv(V2_CSV)
    for _, row in v2.iterrows():
        boltz_plddt[(row["name"], row["bucket"])] = row["boltz_plddt"]

# ── Collect sequences ──────────────────────────────────────────────────────

rows = []

def add(bar_id, bucket, seq, name=None, plddt=None):
    seq = str(seq).upper().strip()
    if not seq or not set(seq) <= VALID_AA:
        return
    song  = meta.loc[bar_id, "genius_song_title"] if bar_id in meta.index else ""
    icon  = meta.loc[bar_id, "aggregate_iconicity"] if bar_id in meta.index else float("nan")
    n_hits = fs_hits.get((bar_id, bucket), 0)
    rows.append({
        "name":     name or f"{bar_id}__{bucket}",
        "bar_id":   bar_id,
        "bucket":   bucket,
        "song":     song,
        "iconicity": icon,
        "sequence": seq,
        "n_hits":   n_hits,
        "boltz_plddt": plddt,
    })

# concordance + native_ala (one per bar from enriched)
for bar_id in bar_ids:
    r = enriched[enriched["bar_id"] == bar_id].iloc[0]
    add(bar_id, "concordance", r.get("fasta_seq_concordance"),
        name=f"{bar_id}__concordance",
        plddt=boltz_plddt.get((f"{bar_id}__concordance", "concordance")))
    add(bar_id, "native_ala", r.get("fasta_seq_native_alanine"),
        name=f"{bar_id}__native_ala",
        plddt=boltz_plddt.get((f"{bar_id}__native_ala", "native_ala")))

# free_design — all passing designs
if FREE_CSV.exists():
    fd = pd.read_csv(FREE_CSV)
    for _, row in fd.iterrows():
        if row.get("bar_id") not in bar_ids:
            continue
        plddt = boltz_plddt.get((row["name"], "free_design")) if "name" in row else None
        add(row["bar_id"], "free_design", row["sequence"],
            name=row.get("name", f"{row['bar_id']}__free_design"),
            plddt=plddt)

# native_ala_free — all passing designs
if NAF_CSV.exists():
    naf = pd.read_csv(NAF_CSV)
    for _, row in naf.iterrows():
        if row.get("bar_id") not in bar_ids:
            continue
        add(row["bar_id"], "native_ala_free", row["sequence"],
            name=row.get("name", f"{row['bar_id']}__native_ala_free"),
            plddt=None)

df = pd.DataFrame(rows).drop_duplicates("name").reset_index(drop=True)
print(f"Total sequences: {len(df)}")
for bucket, g in df.groupby("bucket"):
    print(f"  {bucket:22s}  n={len(g)}")

# ── ESM-2 embeddings ────────────────────────────────────────────────────────

if EMB_CACHE.exists():
    print(f"\nLoading cached embeddings from {EMB_CACHE.name}")
    emb_df = pd.read_csv(EMB_CACHE, index_col=0)
    # align to df order
    missing = set(df["name"]) - set(emb_df.index)
    if missing:
        print(f"  {len(missing)} new sequences not in cache — recomputing all")
        EMB_CACHE.unlink()
    else:
        emb_matrix = emb_df.loc[df["name"]].values
        print(f"  {emb_matrix.shape[1]}-dim embeddings loaded")

if not EMB_CACHE.exists():
    print("\nComputing ESM-2 embeddings (esm2_t6_8M_UR50D)...")
    import esm as esmlib
    model, alphabet = esmlib.pretrained.esm2_t6_8M_UR50D()
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    import torch
    BATCH = 32
    all_embs = []
    seqs = list(df["sequence"])
    names = list(df["name"])

    for i in range(0, len(seqs), BATCH):
        batch_seqs = seqs[i:i+BATCH]
        batch_names = names[i:i+BATCH]
        # truncate to 1022 (ESM-2 limit)
        data = [(n, s[:1022]) for n, s in zip(batch_names, batch_seqs)]
        _, _, tokens = batch_converter(data)
        with torch.no_grad():
            results = model(tokens, repr_layers=[6], return_contacts=False)
        reps = results["representations"][6]  # (B, L+2, D)
        for j, (n, s) in enumerate(zip(batch_names, batch_seqs)):
            L = min(len(s), 1022)
            emb = reps[j, 1:L+1].mean(0).cpu().numpy()  # mean pool over residues
            all_embs.append(emb)
        if (i // BATCH) % 5 == 0:
            print(f"  {i+len(batch_seqs)}/{len(seqs)} done")

    emb_matrix = np.stack(all_embs)
    emb_df = pd.DataFrame(emb_matrix, index=names)
    emb_df.to_csv(EMB_CACHE)
    print(f"Embeddings saved → {EMB_CACHE.name}  shape={emb_matrix.shape}")

# ── UMAP ────────────────────────────────────────────────────────────────────

print("\nRunning UMAP...")
import umap as umap_lib
reducer = umap_lib.UMAP(n_neighbors=15, min_dist=0.1, metric="cosine",
                        random_state=42, verbose=False)
coords = reducer.fit_transform(emb_matrix)
df["umap_x"] = coords[:, 0]
df["umap_y"] = coords[:, 1]
print("UMAP done.")

# ── Plot ─────────────────────────────────────────────────────────────────────

def draw_hull(ax, pts, color, alpha=0.08):
    """Draw convex hull polygon around a set of 2D points."""
    if len(pts) < 3:
        return
    try:
        hull = ConvexHull(pts)
        verts = pts[hull.vertices]
        verts = np.vstack([verts, verts[0]])
        ax.fill(verts[:, 0], verts[:, 1], color=color, alpha=alpha, zorder=1)
        ax.plot(verts[:, 0], verts[:, 1], color=color, alpha=0.3, lw=0.8, zorder=2)
    except Exception:
        pass

fig, axes = plt.subplots(1, 2, figsize=(20, 9), facecolor=DARK_BG)
fig.suptitle("ESM-2 → UMAP: Rap Snacks sequence space\n"
             "concordance · native_ala · free_design · native_ala_free",
             color="white", fontsize=13, y=1.01)

for ax_i, ax in enumerate(axes):
    ax.set_facecolor(DARK_BG)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(colors="#555")
    ax.set_xlabel("UMAP 1", color="#666", fontsize=9)
    ax.set_ylabel("UMAP 2", color="#666", fontsize=9)

# ── Panel A: colour by bucket ─────────────────────────────────────────────
ax = axes[0]
ax.set_title("Coloured by bucket", color="white", fontsize=11)

# Draw hulls per bar (using all 4 buckets combined)
for bar_id in bar_ids:
    g = df[df["bar_id"] == bar_id]
    pts = g[["umap_x", "umap_y"]].values
    # hull colour = average of bucket colours (just use white/dim)
    draw_hull(ax, pts, color="#ffffff", alpha=0.04)

# Plot points
for bucket in BUCKET_ORDER:
    g = df[df["bucket"] == bucket]
    is_lyric = bucket in ("concordance", "native_ala")
    sizes = np.clip(g["n_hits"].fillna(0).values * 0.15 + 20, 20, 120) if is_lyric else np.full(len(g), 12)
    ax.scatter(g["umap_x"], g["umap_y"],
               c=COLORS[bucket], s=sizes, alpha=0.75 if is_lyric else 0.35,
               zorder=3 if is_lyric else 2, label=bucket, linewidths=0)

# Label bar centroids
for bar_id in bar_ids:
    g = df[df["bar_id"] == bar_id]
    cx, cy = g["umap_x"].mean(), g["umap_y"].mean()
    song = meta.loc[bar_id, "genius_song_title"] if bar_id in meta.index else bar_id
    short = song[:18] + "…" if len(song) > 18 else song
    ax.text(cx, cy, f"{bar_id}\n{short}", ha="center", va="center",
            fontsize=5.5, color="white", alpha=0.85,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#111", edgecolor="none", alpha=0.6))

legend_handles = [mpatches.Patch(color=COLORS[b], label=b) for b in BUCKET_ORDER]
ax.legend(handles=legend_handles, facecolor="#1a1a1a", labelcolor="white",
          fontsize=8, loc="upper left", framealpha=0.8)

# ── Panel B: colour by Boltz pLDDT ───────────────────────────────────────
ax = axes[1]
ax.set_title("Coloured by Boltz pLDDT (concordance/native_ala/free_design only)",
             color="white", fontsize=10)

# Draw hulls per bar
for bar_id in bar_ids:
    g = df[df["bar_id"] == bar_id]
    pts = g[["umap_x", "umap_y"]].values
    draw_hull(ax, pts, color="#ffffff", alpha=0.04)

# Plot MPNN designs (no pLDDT) as dim grey background
mpnn = df[df["bucket"].isin(["free_design", "native_ala_free"])]
ax.scatter(mpnn["umap_x"], mpnn["umap_y"], c="#333333", s=8, alpha=0.3, zorder=1, linewidths=0)

# Plot lyric sequences coloured by pLDDT
lyric = df[df["bucket"].isin(["concordance", "native_ala"])].dropna(subset=["boltz_plddt"])
sc = ax.scatter(lyric["umap_x"], lyric["umap_y"],
                c=lyric["boltz_plddt"], cmap="plasma",
                vmin=0.3, vmax=0.9, s=80, alpha=0.9, zorder=4,
                marker="o", linewidths=0.5, edgecolors="white")
plt.colorbar(sc, ax=ax, label="Boltz pLDDT", fraction=0.03, pad=0.02)

# Label bar centroids
for bar_id in bar_ids:
    g = df[df["bar_id"] == bar_id]
    cx, cy = g["umap_x"].mean(), g["umap_y"].mean()
    song = meta.loc[bar_id, "genius_song_title"] if bar_id in meta.index else bar_id
    short = song[:18] + "…" if len(song) > 18 else song
    ax.text(cx, cy, f"{bar_id}\n{short}", ha="center", va="center",
            fontsize=5.5, color="white", alpha=0.85,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#111", edgecolor="none", alpha=0.6))

for s in ax.spines.values():
    s.set_visible(False)

plt.tight_layout()
out_path = OUT_FIG / "fig83_umap_esm2.png"
plt.savefig(out_path, dpi=160, facecolor=DARK_BG, bbox_inches="tight")
plt.close()
print(f"\nSaved → {out_path}")
print(f"Embeddings cached → {EMB_CACHE}")
