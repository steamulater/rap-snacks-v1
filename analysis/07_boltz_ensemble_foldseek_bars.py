"""
analysis/07_boltz_ensemble_foldseek_bars.py
--------------------------------------------
Deep dive into the 6 FoldSeek-searched bars (pTM >= 0.4).
Compares all 5 Boltz-2 diffusion models per bar across three figures.

Figure 1 — Per-residue pLDDT profiles
  2×3 grid. Each subplot: one bar, 5 model lines + mean ± SD ribbon.
  Shows WHERE in the sequence confidence is high/low and whether
  diffusion samples agree (tight ribbon = well-determined fold region).

Figure 2 — Within-bar model agreement: pairwise RMSD matrix (5×5)
  6 small heatmaps side by side (one per bar). Each cell = Kabsch RMSD
  between two diffusion models. Answers: are the 5 structures consistent?

Figure 3 — UMAP of all 30 structures (5 models × 6 bars)
  30×30 pairwise RMSD → UMAP embedding. Color by bar, marker by model index.
  Shows inter-bar separation vs intra-bar structural spread.

Outputs:
  outputs/pairwise/fig4_plddt_profiles.png
  outputs/pairwise/fig5_intrabar_rmsd.png
  outputs/pairwise/fig6_umap_30structs.png

Usage:
  python analysis/07_boltz_ensemble_foldseek_bars.py
"""

import csv
import json
import re
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import umap

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------
ROOT          = Path(__file__).resolve().parent.parent
ID_MAP_JSON   = ROOT / "data/boltz_id_map.json"
SNAPSHOT_JSON = ROOT / "data/bar_index_snapshot.json"
BOLTZ_SUMMARY = ROOT / "outputs/boltz/boltz_summary.csv"
BOLTZ_MODELS  = ROOT / "outputs/boltz/boltz_models.csv"
PRED_DIR      = ROOT / "outputs/boltz_outputs/boltz_results_boltz_inputs/predictions"
OUT_DIR       = ROOT / "outputs/pairwise"

FOLDSEEK_BARS = ["bar_11", "bar_17", "bar_27", "bar_38", "bar_49", "bar_53"]
N_MODELS      = 5

BAR_COLORS = {
    "bar_11": "#e74c3c",  # red
    "bar_17": "#3498db",  # blue
    "bar_27": "#2ecc71",  # green  (confident_protein_like)
    "bar_38": "#f39c12",  # orange
    "bar_49": "#9b59b6",  # purple
    "bar_53": "#1abc9c",  # teal
}

THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

# ---------------------------------------------------------------------------
# PDB parsing — Boltz 2-char chain IDs via regex
# ---------------------------------------------------------------------------
_CA_RE = re.compile(
    r"^ATOM\s+\d+\s+CA\s+(\w{3})\s+\S+\s+(-?\d+)\s+"
    r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+[\d.]+\s+([\d.]+)"
)

def parse_ca(pdb_path: Path) -> tuple[str, np.ndarray, np.ndarray]:
    """Returns (sequence, coords [N,3], per_residue_plddt [N])."""
    seq, coords, plddt = [], [], []
    last_resnum = None
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        m = _CA_RE.match(line)
        if m:
            resname, resnum, x, y, z, bf = m.groups()
            if resnum != last_resnum:
                seq.append(THREE_TO_ONE.get(resname, "X"))
                coords.append([float(x), float(y), float(z)])
                plddt.append(float(bf) / 100.0)   # B-factor is pLDDT × 100
                last_resnum = resnum
    return "".join(seq), np.array(coords, dtype=np.float64), np.array(plddt)


# ---------------------------------------------------------------------------
# Kabsch RMSD (no sequence alignment needed — same bar, same sequence)
# ---------------------------------------------------------------------------
def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    n = min(len(P), len(Q))
    P, Q = P[:n] - P[:n].mean(0), Q[:n] - Q[:n].mean(0)
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    diff = (P @ R.T) - Q
    return float(np.sqrt((diff ** 2).sum(1).mean()))


# ---------------------------------------------------------------------------
# Load all data
# ---------------------------------------------------------------------------
def load_all():
    id_map   = json.load(open(ID_MAP_JSON))
    inv_map  = {v: k for k, v in id_map.items()}
    snap     = json.load(open(SNAPSHOT_JSON))
    summary  = {r["bar_id"]: r for r in csv.DictReader(open(BOLTZ_SUMMARY))}

    data = {}   # bar_id -> {model_i -> {seq, coords, plddt_per_res, plddt_mean, ptm}}
    for bar_id in FOLDSEEK_BARS:
        boltz_id = inv_map[bar_id]
        models = {}
        for mi in range(N_MODELS):
            pdb_path = PRED_DIR / boltz_id / f"{boltz_id}_model_{mi}.pdb"
            seq, coords, plddt_res = parse_ca(pdb_path)
            # Per-model scalar metrics from boltz_models.csv (loaded below)
            models[mi] = {
                "seq":       seq,
                "coords":    coords,
                "plddt_res": plddt_res,
            }
        song = snap.get(bar_id, {}).get("genius_song_title", bar_id)
        data[bar_id] = {"models": models, "song": song}

    # Attach scalar pLDDT/pTM per model from boltz_models.csv
    for row in csv.DictReader(open(BOLTZ_MODELS)):
        bar_id = row["bar_id"]
        if bar_id not in data:
            continue
        mi = int(row["model"])
        if mi in data[bar_id]["models"]:
            data[bar_id]["models"][mi]["plddt_scalar"] = float(row.get("plddt") or 0)
            data[bar_id]["models"][mi]["ptm_scalar"]   = float(row.get("ptm") or 0)

    return data


# ---------------------------------------------------------------------------
# Figure 1 — Per-residue pLDDT profiles
# ---------------------------------------------------------------------------
def fig_plddt_profiles(data: dict, out_path: Path):
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
    axes = axes.flatten()

    plddt_bands = [
        (0.0, 0.5,  "#ff7f7f", "Low (<0.5)"),
        (0.5, 0.7,  "#ffe680", "Medium (0.5–0.7)"),
        (0.7, 1.01, "#c8e6c9", "High (>0.7)"),
    ]

    model_palette = cm.get_cmap("tab10")

    for ax_idx, bar_id in enumerate(FOLDSEEK_BARS):
        ax = axes[ax_idx]
        bar_data = data[bar_id]
        song = bar_data["song"]

        all_profiles = np.array([bar_data["models"][mi]["plddt_res"]
                                  for mi in range(N_MODELS)])
        n_res = all_profiles.shape[1]
        x = np.arange(1, n_res + 1)
        mean_profile = all_profiles.mean(0)
        sd_profile   = all_profiles.std(0)

        # pLDDT confidence band shading
        for lo, hi, color, _ in plddt_bands:
            ax.axhspan(lo, hi, alpha=0.12, color=color, zorder=0)

        # Individual model lines (thin, semi-transparent)
        for mi in range(N_MODELS):
            ax.plot(x, all_profiles[mi], color=model_palette(mi),
                    lw=0.8, alpha=0.55, zorder=2)

        # Mean ± SD ribbon
        ax.fill_between(x, mean_profile - sd_profile, mean_profile + sd_profile,
                        alpha=0.25, color=BAR_COLORS[bar_id], zorder=3)
        ax.plot(x, mean_profile, color=BAR_COLORS[bar_id], lw=2.0, zorder=4,
                label="mean")

        # Threshold lines
        ax.axhline(0.7, color="#888888", lw=0.7, ls="--", zorder=1)
        ax.axhline(0.5, color="#bbbbbb", lw=0.7, ls=":", zorder=1)

        ax.set_title(f"{bar_id}  ·  {song[:28]}", fontsize=9, fontweight="bold",
                     color=BAR_COLORS[bar_id])
        ax.set_xlabel("Residue position", fontsize=8)
        if ax_idx % 3 == 0:
            ax.set_ylabel("pLDDT", fontsize=8)
        ax.set_ylim(0, 1)
        ax.set_xlim(1, n_res)
        ax.tick_params(labelsize=7)

        # Annotation: mean pLDDT ± SD
        mean_val = mean_profile.mean()
        sd_val   = sd_profile.mean()
        ax.text(0.97, 0.06, f"mean={mean_val:.3f} ± {sd_val:.3f}",
                transform=ax.transAxes, ha="right", fontsize=7, color="#333333")

    # Shared legend for model lines
    model_lines = [plt.Line2D([0], [0], color=model_palette(mi), lw=1.2, alpha=0.7,
                               label=f"Model {mi}") for mi in range(N_MODELS)]
    mean_patch = plt.Line2D([0], [0], color="black", lw=2, label="Mean")
    sd_patch   = mpatches.Patch(color="grey", alpha=0.3, label="± SD ribbon")
    fig.legend(handles=model_lines + [mean_patch, sd_patch],
               loc="lower center", ncol=7, fontsize=8, frameon=True,
               bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "Per-residue pLDDT Profiles — 5 Diffusion Models per Bar\n"
        "Shaded ribbon = SD across models · Dashed lines at 0.5 and 0.7",
        fontsize=12, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Figure 2 — Within-bar pairwise model RMSD heatmaps (5×5 × 6 bars)
# ---------------------------------------------------------------------------
def fig_intrabar_rmsd(data: dict, out_path: Path):
    fig, axes = plt.subplots(1, 6, figsize=(17, 3.2))

    model_labels = [f"M{i}" for i in range(N_MODELS)]
    vmax_global = 0.0

    # Pre-compute all matrices to get consistent color scale
    matrices = {}
    for bar_id in FOLDSEEK_BARS:
        mat = np.zeros((N_MODELS, N_MODELS))
        for i in range(N_MODELS):
            for j in range(i + 1, N_MODELS):
                r = kabsch_rmsd(data[bar_id]["models"][i]["coords"],
                                data[bar_id]["models"][j]["coords"])
                mat[i, j] = mat[j, i] = r
        matrices[bar_id] = mat
        vmax_global = max(vmax_global, mat.max())

    import matplotlib.colors as mcolors
    cmap = plt.get_cmap("YlOrRd")

    for ax_idx, bar_id in enumerate(FOLDSEEK_BARS):
        ax = axes[ax_idx]
        mat = matrices[bar_id]
        song = data[bar_id]["song"]

        im = ax.imshow(mat, cmap=cmap, vmin=0, vmax=vmax_global,
                       aspect="auto", interpolation="nearest")

        # Annotate cells
        for i in range(N_MODELS):
            for j in range(N_MODELS):
                val = mat[i, j]
                text_color = "white" if val > vmax_global * 0.6 else "#333333"
                ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                        fontsize=7.5, color=text_color, fontweight="bold")

        ax.set_xticks(range(N_MODELS))
        ax.set_yticks(range(N_MODELS))
        ax.set_xticklabels(model_labels, fontsize=7)
        ax.set_yticklabels(model_labels, fontsize=7)
        ax.set_title(
            f"{bar_id}\n{song[:20]}",
            fontsize=8, fontweight="bold", color=BAR_COLORS[bar_id],
        )

        # Mean off-diagonal RMSD
        off_diag = mat[np.triu_indices(N_MODELS, k=1)]
        ax.set_xlabel(f"mean={off_diag.mean():.2f}Å  max={off_diag.max():.2f}Å",
                      fontsize=7)

    # Shared colorbar
    cbar = fig.colorbar(
        im, ax=axes, orientation="vertical",
        fraction=0.02, pad=0.02,
        label="Cα RMSD between models (Å)",
    )
    cbar.ax.tick_params(labelsize=8)

    fig.suptitle(
        "Within-bar Pairwise Model RMSD (Å) — 5 Diffusion Samples × 6 FoldSeek Bars\n"
        "Low values = diffusion samples converge to same fold · "
        "High values = structural plasticity",
        fontsize=11, y=1.04,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Figure 3 — UMAP of all 30 structures
# ---------------------------------------------------------------------------
def fig_umap(data: dict, out_path: Path):
    # Build 30×30 pairwise RMSD matrix
    bar_model_pairs = [(b, mi) for b in FOLDSEEK_BARS for mi in range(N_MODELS)]
    n = len(bar_model_pairs)

    rmsd_mat = np.zeros((n, n))
    for i, (b1, m1) in enumerate(bar_model_pairs):
        for j, (b2, m2) in enumerate(bar_model_pairs):
            if i >= j:
                continue
            # Cross-bar: need sequence-guided alignment
            s1 = data[b1]["models"][m1]["seq"]
            s2 = data[b2]["models"][m2]["seq"]
            c1 = data[b1]["models"][m1]["coords"]
            c2 = data[b2]["models"][m2]["coords"]
            if b1 == b2:
                # Same bar: same sequence, direct Kabsch
                r = kabsch_rmsd(c1, c2)
            else:
                # Different bar: sequence-guided alignment
                from Bio.Align import PairwiseAligner
                aligner = PairwiseAligner()
                aligner.mode = "global"
                aligner.match_score, aligner.mismatch_score = 1, 0
                aligner.open_gap_score, aligner.extend_gap_score = -1, -0.1
                aln = aligner.align(s1, s2)[0]
                lines = str(aln).splitlines()
                idx1, idx2 = [], []
                ii = jj = 0
                for a, b_ch in zip(lines[0], lines[2]):
                    if a != "-" and b_ch != "-":
                        idx1.append(ii)
                        idx2.append(jj)
                    if a != "-": ii += 1
                    if b_ch != "-": jj += 1
                if len(idx1) >= 4:
                    r = kabsch_rmsd(c1[idx1], c2[idx2])
                else:
                    r = np.nan
            if np.isnan(r):
                r = 20.0   # fallback for near-zero alignment
            rmsd_mat[i, j] = rmsd_mat[j, i] = r

    # UMAP on precomputed distance
    reducer = umap.UMAP(
        n_components=2,
        metric="precomputed",
        n_neighbors=8,
        min_dist=0.3,
        random_state=42,
        n_epochs=500,
    )
    embedding = reducer.fit_transform(rmsd_mat)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 8))

    model_markers = ["o", "s", "^", "D", "P"]   # model 0–4
    model_sizes   = [110, 90, 90, 90, 90]

    for i, (bar_id, mi) in enumerate(bar_model_pairs):
        ax.scatter(
            embedding[i, 0], embedding[i, 1],
            c=BAR_COLORS[bar_id],
            marker=model_markers[mi],
            s=model_sizes[mi],
            alpha=0.85,
            edgecolors="white",
            linewidths=0.6,
            zorder=3,
        )

    # Draw convex hull / ellipse per bar to show within-bar spread
    from matplotlib.patches import Ellipse
    from scipy.spatial import ConvexHull

    for bar_id in FOLDSEEK_BARS:
        idxs = [i for i, (b, _) in enumerate(bar_model_pairs) if b == bar_id]
        pts  = embedding[idxs]
        if len(pts) >= 3:
            try:
                hull = ConvexHull(pts)
                hull_pts = np.append(pts[hull.vertices], [pts[hull.vertices[0]]], axis=0)
                ax.fill(hull_pts[:, 0], hull_pts[:, 1],
                        alpha=0.08, color=BAR_COLORS[bar_id], zorder=1)
                ax.plot(hull_pts[:, 0], hull_pts[:, 1],
                        color=BAR_COLORS[bar_id], lw=1.0, alpha=0.4, zorder=2)
            except Exception:
                pass

    # Label centroid of each bar's cluster
    for bar_id in FOLDSEEK_BARS:
        idxs = [i for i, (b, _) in enumerate(bar_model_pairs) if b == bar_id]
        cx, cy = embedding[idxs].mean(0)
        song = data[bar_id]["song"]
        ax.annotate(
            f"{bar_id.replace('bar_', '')} · {song[:20]}",
            (cx, cy),
            xytext=(cx + 0.15, cy + 0.15),
            fontsize=8.5,
            fontweight="bold",
            color=BAR_COLORS[bar_id],
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                      edgecolor=BAR_COLORS[bar_id], alpha=0.85, lw=1),
        )

    # Legend: bars (color) + models (shape)
    bar_patches = [mpatches.Patch(color=BAR_COLORS[b], label=b) for b in FOLDSEEK_BARS]
    model_handles = [
        plt.scatter([], [], marker=model_markers[mi], c="#666666",
                    s=model_sizes[mi], label=f"Model {mi}")
        for mi in range(N_MODELS)
    ]
    l1 = ax.legend(handles=bar_patches, title="Bar", loc="upper left",
                   fontsize=8, title_fontsize=8)
    ax.legend(handles=model_handles, title="Diffusion model", loc="upper right",
              fontsize=8, title_fontsize=8)
    ax.add_artist(l1)

    ax.set_title(
        "UMAP of 30 Structures (5 Diffusion Models × 6 FoldSeek Bars)\n"
        "Distance = pairwise Cα RMSD · Shaded hull = within-bar spread",
        fontsize=12,
    )
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.grid(True, alpha=0.15)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("[07_boltz_ensemble_foldseek_bars.py] Loading 30 PDB files (5 models × 6 bars)...")
    data = load_all()

    print("  Figure 1: per-residue pLDDT profiles...")
    fig_plddt_profiles(data, OUT_DIR / "fig4_plddt_profiles.png")

    print("  Figure 2: within-bar model RMSD heatmaps...")
    fig_intrabar_rmsd(data, OUT_DIR / "fig5_intrabar_rmsd.png")

    print("  Figure 3: UMAP of 30 structures...")
    fig_umap(data, OUT_DIR / "fig6_umap_30structs.png")

    # Print within-bar RMSD summary
    print("\n  Within-bar model consistency (mean Cα RMSD between diffusion samples):")
    for bar_id in FOLDSEEK_BARS:
        coords_list = [data[bar_id]["models"][mi]["coords"] for mi in range(N_MODELS)]
        rmsds = []
        for i in range(N_MODELS):
            for j in range(i + 1, N_MODELS):
                rmsds.append(kabsch_rmsd(coords_list[i], coords_list[j]))
        song = data[bar_id]["song"]
        print(f"    {bar_id:10s} ({song[:25]:25s})  "
              f"mean={np.mean(rmsds):.2f}Å  max={np.max(rmsds):.2f}Å  "
              f"min={np.min(rmsds):.2f}Å")

    print(f"\nDone. Figures saved to {OUT_DIR}/")


if __name__ == "__main__":
    main()
