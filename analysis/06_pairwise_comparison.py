"""
analysis/06_pairwise_comparison.py
------------------------------------
Pairwise sequence identity and structural RMSD for all 85 bars.

For each bar we use the best-confidence Boltz-2 model (by confidence_mean).
Sequence identity: global alignment (Needleman-Wunsch via Bio.Align.PairwiseAligner),
  % identity = matches / alignment_length.
Structural RMSD: sequence-guided Cα superposition (Kabsch algorithm in numpy).
  Only aligned (non-gap) residue pairs are superimposed — handles variable lengths.

Outputs:
  outputs/pairwise/seq_identity.csv          85×85 sequence identity matrix
  outputs/pairwise/struct_rmsd.csv           85×85 structural RMSD matrix (Angstroms)
  outputs/pairwise/fig1_seq_heatmap.png
  outputs/pairwise/fig2_struct_heatmap.png
  outputs/pairwise/fig3_novelty_scatter.png

Usage:
  python analysis/06_pairwise_comparison.py
  python analysis/06_pairwise_comparison.py --skip-rmsd   # seq only (fast)
"""

import argparse
import csv
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from Bio.Align import PairwiseAligner

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
SNAPSHOT_JSON   = ROOT / "data/bar_index_snapshot.json"
ID_MAP_JSON     = ROOT / "data/boltz_id_map.json"
BOLTZ_SUMMARY   = ROOT / "outputs/boltz/boltz_summary.csv"
PRED_DIR        = ROOT / "outputs/boltz_outputs/boltz_results_boltz_inputs/predictions"
OUT_DIR         = ROOT / "outputs/pairwise"

THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

# ---------------------------------------------------------------------------
# PDB parsing — handles Boltz 2-char chain IDs via regex
# ---------------------------------------------------------------------------
_CA_RE = re.compile(
    r"^ATOM\s+\d+\s+CA\s+(\w{3})\s+\S+\s+(-?\d+)\s+"
    r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)"
)

def parse_ca_atoms(pdb_path: Path) -> tuple[str, np.ndarray]:
    """Return (sequence_str, coords_array shape [N,3]) for Cα atoms."""
    seq = []
    coords = []
    last_resnum = None
    for line in pdb_path.read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        m = _CA_RE.match(line)
        if m:
            resname, resnum, x, y, z = m.groups()
            if resnum != last_resnum:
                seq.append(THREE_TO_ONE.get(resname, "X"))
                coords.append([float(x), float(y), float(z)])
                last_resnum = resnum
    return "".join(seq), np.array(coords, dtype=np.float64)


# ---------------------------------------------------------------------------
# Sequence identity
# ---------------------------------------------------------------------------
_aligner = PairwiseAligner()
_aligner.mode = "global"
_aligner.match_score = 1
_aligner.mismatch_score = 0
_aligner.open_gap_score = -1
_aligner.extend_gap_score = -0.1

def seq_identity(s1: str, s2: str) -> float:
    """Global alignment % identity = matches / alignment_length."""
    aln = _aligner.align(s1, s2)[0]
    aligned = str(aln).splitlines()
    # Lines: seq1, match_row, seq2
    if len(aligned) < 3:
        return 0.0
    match_row = aligned[1]
    matches = match_row.count("|")
    aln_len = len(match_row)
    return matches / aln_len if aln_len else 0.0


def aligned_indices(s1: str, s2: str) -> tuple[list[int], list[int]]:
    """Return index pairs (i1, i2) for aligned (non-gap) positions."""
    aln = _aligner.align(s1, s2)[0]
    i1_list, i2_list = [], []
    i1 = i2 = 0
    for a1_char, a2_char in zip(*[str(aln).splitlines()[k] for k in [0, 2]]):
        in_gap1 = a1_char == "-"
        in_gap2 = a2_char == "-"
        if not in_gap1 and not in_gap2:
            i1_list.append(i1)
            i2_list.append(i2)
        if not in_gap1:
            i1 += 1
        if not in_gap2:
            i2 += 1
    return i1_list, i2_list


# ---------------------------------------------------------------------------
# Kabsch RMSD
# ---------------------------------------------------------------------------
def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD after optimal superposition of equal-length coord sets (Kabsch)."""
    if len(P) == 0:
        return np.nan
    P = P - P.mean(axis=0)
    Q = Q - Q.mean(axis=0)
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    P_rot = P @ R.T
    diff = P_rot - Q
    return float(np.sqrt((diff ** 2).sum(axis=1).mean()))


def struct_rmsd(coords1: np.ndarray, seq1: str,
                coords2: np.ndarray, seq2: str) -> float:
    """Sequence-guided Cα RMSD: superimpose only aligned residue pairs."""
    idx1, idx2 = aligned_indices(seq1, seq2)
    if len(idx1) < 4:
        return np.nan
    P = coords1[idx1]
    Q = coords2[idx2]
    return kabsch_rmsd(P, Q)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--skip-rmsd", action="store_true",
                   help="Compute sequence identity only (skip RMSD)")
    p.add_argument("--load-cached", action="store_true",
                   help="Load previously saved matrices from outputs/pairwise/")
    return p.parse_args()


def main():
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- Load metadata ---
    snap = json.load(open(SNAPSHOT_JSON))
    id_map = json.load(open(ID_MAP_JSON))          # boltz_id -> bar_id
    inv_map = {v: k for k, v in id_map.items()}    # bar_id -> boltz_id

    summary_rows = {r["bar_id"]: r
                    for r in csv.DictReader(open(BOLTZ_SUMMARY))}

    # Sort bars by bar number for consistent ordering
    bar_ids = sorted(snap.keys(), key=lambda b: int(b.split("_")[1]))
    N = len(bar_ids)
    print(f"[06_pairwise_comparison.py] {N} bars")

    # --- Build bar metadata lookup ---
    meta = {}
    for bar_id in bar_ids:
        s = snap.get(bar_id, {})
        r = summary_rows.get(bar_id, {})
        meta[bar_id] = {
            "song":     s.get("genius_song_title", ""),
            "bar":      s.get("canonical_bar", ""),
            "iconicity": float(s.get("aggregate_iconicity") or 0),
            "plddt":    float(r.get("plddt_mean") or 0),
            "ptm":      float(r.get("ptm_mean") or 0),
            "struct_class": r.get("structural_class", "disordered"),
        }

    # FoldSeek-searched bars (pTM >= 0.4)
    foldseek_bars = {b for b, m in meta.items() if m["ptm"] >= 0.4}

    # --- Load sequences and PDB coords ---
    seqs = {}
    coords_map = {}

    print("  Loading PDB files...")
    for bar_id in bar_ids:
        boltz_id = inv_map.get(bar_id, bar_id)
        best_model = summary_rows.get(bar_id, {}).get("best_model", "0")
        pdb_path = PRED_DIR / boltz_id / f"{boltz_id}_model_{best_model}.pdb"
        if not pdb_path.exists():
            # fallback to model_0
            pdb_path = PRED_DIR / boltz_id / f"{boltz_id}_model_0.pdb"
        seq, coords = parse_ca_atoms(pdb_path)
        seqs[bar_id] = seq
        coords_map[bar_id] = coords

    # --- Compute or load matrices ---
    seq_csv = OUT_DIR / "seq_identity.csv"
    rmsd_csv = OUT_DIR / "struct_rmsd.csv"

    if args.load_cached and seq_csv.exists():
        print("  Loading cached sequence identity matrix...")
        seq_mat = np.loadtxt(seq_csv, delimiter=",")
    else:
        print(f"  Computing {N}×{N} sequence identity matrix ({N*(N-1)//2} pairs)...")
        seq_mat = np.eye(N, dtype=np.float64)
        for i in range(N):
            for j in range(i + 1, N):
                sid = seq_identity(seqs[bar_ids[i]], seqs[bar_ids[j]])
                seq_mat[i, j] = sid
                seq_mat[j, i] = sid
            if (i + 1) % 10 == 0:
                print(f"    seq: {i+1}/{N}")
        np.savetxt(seq_csv, seq_mat, delimiter=",", fmt="%.6f")
        print(f"  Saved {seq_csv}")

    if args.skip_rmsd:
        rmsd_mat = None
    elif args.load_cached and rmsd_csv.exists():
        print("  Loading cached structural RMSD matrix...")
        rmsd_mat = np.loadtxt(rmsd_csv, delimiter=",")
    else:
        print(f"  Computing {N}×{N} structural RMSD matrix ({N*(N-1)//2} pairs)...")
        rmsd_mat = np.zeros((N, N), dtype=np.float64)
        for i in range(N):
            for j in range(i + 1, N):
                rmsd = struct_rmsd(
                    coords_map[bar_ids[i]], seqs[bar_ids[i]],
                    coords_map[bar_ids[j]], seqs[bar_ids[j]],
                )
                rmsd_mat[i, j] = rmsd if not np.isnan(rmsd) else 0
                rmsd_mat[j, i] = rmsd_mat[i, j]
            if (i + 1) % 10 == 0:
                print(f"    rmsd: {i+1}/{N}")
        np.savetxt(rmsd_csv, rmsd_mat, delimiter=",", fmt="%.6f")
        print(f"  Saved {rmsd_csv}")

    # -----------------------------------------------------------------------
    # Figure helpers
    # -----------------------------------------------------------------------
    CLASS_COLORS = {
        "confident_protein_like":  "#e74c3c",   # red
        "uncertain_protein_like":  "#f39c12",   # orange
        "disordered":              "#95a5a6",   # grey
        "confident_alien":         "#9b59b6",   # purple
    }

    labels = [bid.replace("bar_", "") for bid in bar_ids]

    # -----------------------------------------------------------------------
    # Figure 1 — Sequence Identity Heatmap
    # -----------------------------------------------------------------------
    print("  Plotting Figure 1: sequence identity heatmap...")
    fig1_path = OUT_DIR / "fig1_seq_heatmap.png"

    row_colors = [CLASS_COLORS.get(meta[b]["struct_class"], "#cccccc") for b in bar_ids]

    cg = sns.clustermap(
        seq_mat,
        cmap="YlOrRd",
        vmin=0, vmax=1,
        xticklabels=labels,
        yticklabels=labels,
        figsize=(16, 14),
        row_colors=row_colors,
        col_colors=row_colors,
        dendrogram_ratio=0.12,
        cbar_pos=(0.02, 0.85, 0.03, 0.12),
    )
    cg.ax_heatmap.set_xlabel("Bar ID", fontsize=10)
    cg.ax_heatmap.set_ylabel("Bar ID", fontsize=10)
    cg.ax_heatmap.tick_params(labelsize=5)
    cg.figure.suptitle(
        "Pairwise Sequence Identity — All 85 Bars (concordance mapping)\n"
        "Clustered by similarity · Color bar = structural class",
        y=1.01, fontsize=11,
    )

    legend_patches = [mpatches.Patch(color=c, label=cls.replace("_", " "))
                      for cls, c in CLASS_COLORS.items() if cls != "confident_alien"]
    cg.ax_heatmap.legend(handles=legend_patches, title="Structural class",
                         loc="upper right", fontsize=7, title_fontsize=7,
                         bbox_to_anchor=(1.25, 1.0))
    cg.figure.savefig(fig1_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved {fig1_path}")

    # -----------------------------------------------------------------------
    # Figure 2 — Structural RMSD Heatmap
    # -----------------------------------------------------------------------
    if rmsd_mat is not None:
        print("  Plotting Figure 2: structural RMSD heatmap...")
        fig2_path = OUT_DIR / "fig2_struct_heatmap.png"

        cg2 = sns.clustermap(
            rmsd_mat,
            cmap="viridis_r",
            xticklabels=labels,
            yticklabels=labels,
            figsize=(16, 14),
            row_colors=row_colors,
            col_colors=row_colors,
            dendrogram_ratio=0.12,
            cbar_pos=(0.02, 0.85, 0.03, 0.12),
        )
        cg2.ax_heatmap.set_xlabel("Bar ID", fontsize=10)
        cg2.ax_heatmap.set_ylabel("Bar ID", fontsize=10)
        cg2.ax_heatmap.tick_params(labelsize=5)
        cg2.figure.suptitle(
            "Pairwise Structural RMSD (Å) — All 85 Bars · Sequence-guided Cα superposition\n"
            "Clustered by similarity · Color bar = structural class",
            y=1.01, fontsize=11,
        )
        cg2.ax_heatmap.legend(handles=legend_patches, title="Structural class",
                              loc="upper right", fontsize=7, title_fontsize=7,
                              bbox_to_anchor=(1.25, 1.0))
        cg2.figure.savefig(fig2_path, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"  Saved {fig2_path}")

    # -----------------------------------------------------------------------
    # Figure 3 — Novelty Scatter (seq novelty × structural novelty)
    # -----------------------------------------------------------------------
    if rmsd_mat is not None:
        print("  Plotting Figure 3: novelty scatter...")
        fig3_path = OUT_DIR / "fig3_novelty_scatter.png"

        # Per-bar novelty scores
        # Sequence novelty: 1 - max seq identity to any OTHER bar (higher = more novel)
        # Structural novelty: min RMSD to any OTHER bar (higher = more structurally distinct)
        seq_novelty = []
        struct_novelty = []
        for i in range(N):
            other_seq = np.concatenate([seq_mat[i, :i], seq_mat[i, i+1:]])
            other_rmsd = np.concatenate([rmsd_mat[i, :i], rmsd_mat[i, i+1:]])
            seq_novelty.append(1.0 - float(np.max(other_seq)))
            struct_novelty.append(float(np.percentile(other_rmsd[other_rmsd > 0], 10))
                                  if np.any(other_rmsd > 0) else 0.0)

        seq_novelty = np.array(seq_novelty)
        struct_novelty = np.array(struct_novelty)

        fig3, ax = plt.subplots(figsize=(12, 9))

        colors = [CLASS_COLORS.get(meta[b]["struct_class"], "#cccccc") for b in bar_ids]
        sizes  = [40 + 120 * meta[b]["iconicity"] for b in bar_ids]

        sc = ax.scatter(seq_novelty, struct_novelty,
                        c=colors, s=sizes, alpha=0.75, linewidths=0.5,
                        edgecolors="white", zorder=3)

        # --- Novelty quadrant threshold lines ---
        seq_thresh   = float(np.percentile(seq_novelty, 50))
        rmsd_thresh  = float(np.percentile(struct_novelty, 50))
        ax.axvline(seq_thresh,  color="#bdc3c7", lw=1, ls="--", zorder=1)
        ax.axhline(rmsd_thresh, color="#bdc3c7", lw=1, ls="--", zorder=1)

        # Shade "novel in both" quadrant (top-right)
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_xlim(xlim); ax.set_ylim(ylim)
        ax.fill_between(
            [seq_thresh, max(seq_novelty) * 1.05],
            rmsd_thresh, max(struct_novelty) * 1.08,
            alpha=0.07, color="#27ae60", zorder=0,
        )
        ax.text(max(seq_novelty) * 0.99, max(struct_novelty) * 1.04,
                "Novel in both", ha="right", va="top",
                fontsize=8, color="#27ae60", style="italic")

        # --- Annotate notable bars ---
        # Always label: foldseek bars, bar_27, top/bottom iconicity extremes
        top_icon   = sorted(bar_ids, key=lambda b: meta[b]["iconicity"], reverse=True)[:3]
        bottom_icon = sorted(bar_ids, key=lambda b: meta[b]["iconicity"])[:2]
        top_seq    = sorted(range(N), key=lambda i: seq_novelty[i], reverse=True)[:3]
        top_rmsd   = sorted(range(N), key=lambda i: struct_novelty[i], reverse=True)[:3]

        annotate_set = (
            foldseek_bars
            | {"bar_27"}
            | set(top_icon)
            | set(bottom_icon)
            | {bar_ids[i] for i in top_seq}
            | {bar_ids[i] for i in top_rmsd}
        )

        for i, bar_id in enumerate(bar_ids):
            if bar_id not in annotate_set:
                continue
            song = meta[bar_id]["song"]
            label = f"{bar_id.replace('bar_','')}·{song[:18]}" if song else bar_id
            xoff = 0.005
            yoff = 0.3
            ax.annotate(
                label,
                (seq_novelty[i], struct_novelty[i]),
                xytext=(seq_novelty[i] + xoff, struct_novelty[i] + yoff),
                fontsize=6.5,
                color="#2c3e50",
                arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.6),
            )

        ax.set_xlabel("Sequence Novelty  (1 − max pairwise identity to any other bar)", fontsize=11)
        ax.set_ylabel("Structural Novelty  (10th-pct Cα RMSD to nearest neighbors, Å)", fontsize=11)
        ax.set_title(
            "Sequence × Structural Novelty — All 85 Bars\n"
            "Size ∝ iconicity · Color = structural class · Green zone = novel in both axes",
            fontsize=12,
        )

        legend_patches2 = [
            mpatches.Patch(color=c, label=cls.replace("_", " "))
            for cls, c in CLASS_COLORS.items() if cls != "confident_alien"
        ]
        size_legend = [
            plt.scatter([], [], s=40 + 120 * v, c="#888888", alpha=0.6, label=f"iconicity={v:.1f}")
            for v in [0.0, 0.5, 1.0]
        ]
        l1 = ax.legend(handles=legend_patches2, title="Structural class",
                       loc="lower left", fontsize=8, title_fontsize=8)
        ax.legend(handles=size_legend, title="Iconicity (size)",
                  loc="lower right", fontsize=8, title_fontsize=8)
        ax.add_artist(l1)

        ax.grid(True, alpha=0.2)
        fig3.tight_layout()
        fig3.savefig(fig3_path, dpi=150, bbox_inches="tight")
        plt.close("all")
        print(f"  Saved {fig3_path}")

        # --- Print summary stats ---
        novel_mask = (seq_novelty > seq_thresh) & (struct_novelty > rmsd_thresh)
        print(f"\n  Novel in both (top-right quadrant): {novel_mask.sum()} bars")
        for i in np.where(novel_mask)[0]:
            b = bar_ids[i]
            print(f"    {b:10s}  seq_nov={seq_novelty[i]:.3f}  "
                  f"struct_nov={struct_novelty[i]:.2f}Å  "
                  f"song={meta[b]['song'][:40]}  "
                  f"ptm={meta[b]['ptm']:.3f}  plddt={meta[b]['plddt']:.3f}")

    print(f"\nDone. Figures in {OUT_DIR}/")
    print("Seq identity stats:")
    mask = ~np.eye(N, dtype=bool)
    print(f"  mean={seq_mat[mask].mean():.3f}  "
          f"min={seq_mat[mask].min():.3f}  "
          f"max={seq_mat[mask].max():.3f}")
    if rmsd_mat is not None:
        rmsd_vals = rmsd_mat[mask & (rmsd_mat > 0)]
        print("Structural RMSD stats (Å):")
        print(f"  mean={rmsd_vals.mean():.2f}  "
              f"min={rmsd_vals.min():.2f}  "
              f"max={rmsd_vals.max():.2f}")


if __name__ == "__main__":
    main()
