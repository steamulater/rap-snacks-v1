"""
Visualize how the 3 scrambles of each native_ala sequence differ.

For each bar, plots:
  Row 0: native_ala  (original)
  Row 1: scramble_0
  Row 2: scramble_1
  Row 3: scramble_2

Each cell = one amino acid, coloured by physicochemical class.
Right panel = pairwise identity matrix (4×4).
Bottom panel = AA composition (must be identical across all 4 rows — sanity check).
"""

import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
from pathlib import Path

ROOT      = Path(__file__).parent.parent
ENRICHED  = ROOT / "data/aggregated_lines_v2_enriched.csv"
CANDS_CSV = ROOT / "data/phase2_candidates.csv"
FIG_DIR   = ROOT / "outputs/figures"
FIG_DIR.mkdir(exist_ok=True)

DARK_BG   = "#0e0e0e"
N_SCRAMBLES = 3
SEED        = 42

# Physicochemical AA colour classes
AA_COLORS = {
    # Hydrophobic — amber
    'A':'#f59e0b','V':'#f59e0b','L':'#f59e0b','I':'#f59e0b',
    'M':'#f59e0b','F':'#f59e0b','W':'#f59e0b','P':'#f59e0b',
    # Polar uncharged — teal
    'S':'#14b8a6','T':'#14b8a6','C':'#14b8a6','Y':'#14b8a6',
    'N':'#14b8a6','Q':'#14b8a6','G':'#14b8a6',
    # Positive — cyan
    'R':'#00d4ff','K':'#00d4ff','H':'#00d4ff',
    # Negative — orange-red
    'D':'#ef4444','E':'#ef4444',
}

def shuffle_seq(seq, seed):
    rng = random.Random(seed)
    s = list(seq); rng.shuffle(s)
    return "".join(s)

def hamming_identity(a, b):
    matches = sum(x == y for x, y in zip(a, b))
    return matches / max(len(a), len(b))

enriched   = pd.read_csv(ENRICHED)
candidates = pd.read_csv(CANDS_CSV)
bar_ids    = list(candidates["bar_id"])
sub        = enriched[enriched["bar_id"].isin(bar_ids)].set_index("bar_id")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# ── Per-bar 4-panel figure ─────────────────────────────────────────────────

for bar_id in bar_ids:
    na_seq = str(sub.loc[bar_id, "fasta_seq_native_alanine"]).upper().strip()
    if not set(na_seq) <= VALID_AA:
        continue

    variants = [("native_ala", na_seq)] + \
               [(f"scramble_{i}", shuffle_seq(na_seq, SEED + i)) for i in range(N_SCRAMBLES)]
    n_vars   = len(variants)
    seq_len  = len(na_seq)

    fig = plt.figure(figsize=(max(14, seq_len * 0.18 + 4), 7), facecolor=DARK_BG)
    gs  = fig.add_gridspec(3, 2, width_ratios=[4, 1],
                           height_ratios=[n_vars, 0.6, 0.6],
                           hspace=0.4, wspace=0.25)

    # ── Left: sequence grid ────────────────────────────────────────────────
    ax_seq = fig.add_subplot(gs[0, 0])
    ax_seq.set_facecolor(DARK_BG)
    ax_seq.set_xlim(0, seq_len)
    ax_seq.set_ylim(-0.5, n_vars - 0.5)

    for ri, (label, seq) in enumerate(variants):
        y = n_vars - 1 - ri
        for ci, aa in enumerate(seq):
            color = AA_COLORS.get(aa, "#888888")
            rect  = mpatches.FancyBboxPatch(
                (ci + 0.05, y - 0.4), 0.9, 0.8,
                boxstyle="round,pad=0.05",
                facecolor=color, edgecolor="none", alpha=0.85,
            )
            ax_seq.add_patch(rect)
            if seq_len <= 60:
                ax_seq.text(ci + 0.5, y, aa, ha="center", va="center",
                            fontsize=5.5, color="black", fontweight="bold")

    ax_seq.set_yticks(range(n_vars))
    ax_seq.set_yticklabels([v[0] for v in reversed(variants)],
                           color="white", fontsize=9)
    ax_seq.set_xlabel("Position", color="white", fontsize=9)
    ax_seq.tick_params(colors="white")
    ax_seq.set_title(f"{bar_id} — native_ala sequence vs 3 scrambles\n"
                     f"(same AA composition, different order — len={seq_len})",
                     color="white", fontsize=10, pad=6)
    for s in ax_seq.spines.values(): s.set_visible(False)

    # ── Right: pairwise identity matrix ───────────────────────────────────
    ax_mat = fig.add_subplot(gs[0, 1])
    ax_mat.set_facecolor(DARK_BG)
    labels = [v[0][:9] for v in variants]
    seqs   = [v[1] for v in variants]
    mat    = np.array([[hamming_identity(seqs[i], seqs[j])
                        for j in range(n_vars)]
                       for i in range(n_vars)])
    im = ax_mat.imshow(mat, cmap="Blues", vmin=0, vmax=1)
    for i in range(n_vars):
        for j in range(n_vars):
            ax_mat.text(j, i, f"{mat[i,j]:.2f}",
                        ha="center", va="center", fontsize=8,
                        color="white" if mat[i,j] < 0.5 else "black")
    ax_mat.set_xticks(range(n_vars))
    ax_mat.set_yticks(range(n_vars))
    ax_mat.set_xticklabels(labels, rotation=45, ha="right", fontsize=7, color="white")
    ax_mat.set_yticklabels(labels, fontsize=7, color="white")
    ax_mat.set_title("Pairwise\nidentity", color="white", fontsize=8)
    ax_mat.tick_params(colors="white")
    plt.colorbar(im, ax=ax_mat, fraction=0.046, pad=0.04)

    # ── Bottom left: AA composition (should be identical) ─────────────────
    ax_comp = fig.add_subplot(gs[1:, 0])
    ax_comp.set_facecolor(DARK_BG)
    aa_order = sorted(VALID_AA)
    x        = np.arange(len(aa_order))
    width    = 0.2
    colors_v = ["white", "#00d4ff", "#ff9500", "#a855f7"]
    for vi, (label, seq) in enumerate(variants):
        cnt = Counter(seq)
        vals = [cnt.get(aa, 0) / len(seq) for aa in aa_order]
        ax_comp.bar(x + (vi - 1.5) * width, vals, width=width * 0.9,
                    color=colors_v[vi], alpha=0.75, label=label)
    ax_comp.set_xticks(x)
    ax_comp.set_xticklabels(aa_order, fontsize=7, color="white")
    ax_comp.set_ylabel("Frequency", color="white", fontsize=8)
    ax_comp.set_title("AA composition — identical across all variants (sanity check)",
                      color="white", fontsize=8)
    ax_comp.tick_params(colors="white")
    ax_comp.legend(facecolor="#1a1a1a", labelcolor="white", fontsize=7,
                   ncol=4, loc="upper right")
    for s in ax_comp.spines.values(): s.set_edgecolor("#333")

    # ── Legend for AA classes ──────────────────────────────────────────────
    class_patches = [
        mpatches.Patch(color="#f59e0b", label="Hydrophobic (AVLIMFWP)"),
        mpatches.Patch(color="#14b8a6", label="Polar (STCYNGQ)"),
        mpatches.Patch(color="#00d4ff", label="Positive (RKH)"),
        mpatches.Patch(color="#ef4444", label="Negative (DE)"),
    ]
    ax_seq.legend(handles=class_patches, loc="lower right",
                  facecolor="#1a1a1a", labelcolor="white", fontsize=7,
                  framealpha=0.8, ncol=2)

    plt.savefig(FIG_DIR / f"fig_scramble_{bar_id}.png",
                dpi=140, facecolor=DARK_BG, bbox_inches="tight")
    plt.close()
    print(f"  {bar_id}  identity(na, sc0)={mat[0,1]:.3f}  "
          f"identity(sc0,sc1)={mat[1,2]:.3f}")

# ── Summary: average pairwise identity across bars ────────────────────────
print("\n--- Average pairwise identity across bars ---")
print("  native_ala vs scramble: ~0.05 expected (random shuffle)")
print("  scramble vs scramble:   ~0.05 expected (independent shuffles)")
print(f"\nFigures saved → outputs/figures/fig_scramble_bar_XX.png")
print("Each figure shows:")
print("  - Coloured sequence grid (position × variant)")
print("  - 4×4 pairwise identity matrix")
print("  - AA composition bar chart (identical across all 4 — confirms only order changes)")
