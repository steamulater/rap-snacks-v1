"""
Stage 15 — Sequence composition analysis of the 24 selected proteins.
Flags: long AA runs, low complexity, high charge, cysteine count, biased composition.
Outputs: fig92_seq_composition.png, fig93_aa_heatmap.png
"""
import csv, re
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

OUT_DIR = "outputs/figures"

# ── load sequences ────────────────────────────────────────────────
seqs = []
with open('outputs/selected_proteins.csv') as f:
    for row in csv.DictReader(f):
        seqs.append({'name': row['name'], 'seq': row['sequence'],
                     'group': row['group'], 'bucket': row['bucket']})

# ── helpers ────────────────────────────────────────────────────────
def find_runs(seq, min_len=4):
    runs = []
    for m in re.finditer(r'(.)\1{%d,}' % (min_len - 1), seq):
        runs.append((m.group()[0], len(m.group()), m.start()))
    return runs

def low_complexity_frac(seq, win=12, thresh=4):
    bad = sum(1 for i in range(len(seq) - win + 1) if len(set(seq[i:i+win])) <= thresh)
    return bad / max(1, len(seq) - win + 1)

def net_charge(seq):
    return seq.count('R') + seq.count('K') - seq.count('D') - seq.count('E')

AAS = list('ACDEFGHIKLMNPQRSTVWY')

results = []
for s in seqs:
    seq = s['seq']
    runs  = find_runs(seq, min_len=4)
    lc    = low_complexity_frac(seq)
    nc    = net_charge(seq)
    cys   = seq.count('C')
    pct   = {aa: seq.count(aa) / len(seq) * 100 for aa in AAS}
    results.append({**s, 'runs': runs, 'lc': lc, 'nc': nc,
                    'cys': cys, 'pct': pct, 'len': len(seq)})

# ── figure 92 — per-sequence flag matrix ─────────────────────────
GROUP_COLOR = {'A': '#2166ac', 'B': '#4dac26', 'C': '#d01c8b', 'D': '#b35806'}

names  = [r['name'] for r in results]
n      = len(names)

# metric columns
run_lens    = [max((x[1] for x in r['runs']), default=0) for r in results]
lc_fracs    = [r['lc'] for r in results]
cys_counts  = [r['cys'] for r in results]
abs_charges = [abs(r['nc']) for r in results]
ala_pcts    = [r['pct']['A'] for r in results]
lys_pcts    = [r['pct']['K'] for r in results]
glu_pcts    = [r['pct']['E'] for r in results]

fig, axes = plt.subplots(1, 7, figsize=(18, 10),
                         gridspec_kw={'wspace': 0.05})

metrics = [
    ('Max repeat\nrun (aa)',   run_lens,    4,    12,  '#d62728', 'Homopolymer runs ≥4'),
    ('Low-complexity\nfraction', lc_fracs, 0.05, 0.5,  '#ff7f0e', 'Low diversity windows'),
    ('Cysteines',              cys_counts,  2,    8,   '#9467bd', 'Cys → disulfide risk'),
    ('|Net charge|',           abs_charges, 12,  25,  '#e377c2', 'Charge imbalance'),
    ('Ala %',                  ala_pcts,   30,   60,  '#8c564b', 'Ala-biased'),
    ('Lys %',                  lys_pcts,   20,   40,  '#1f77b4', 'Lys-biased'),
    ('Glu %',                  glu_pcts,   20,   40,  '#17becf', 'Glu-biased'),
]

short_names = []
for r in results:
    n2 = r['name']
    if n2 == '6E5C_positive_control':
        n2 = '6E5C (+ctrl)'
    elif n2.startswith('bar_'):
        parts = n2.split('_')
        bid = parts[0]+'_'+parts[1]
        n2 = bid + '_' + '_'.join(parts[2:]) if len(parts) > 2 else bid
        if len(n2) > 22:
            n2 = n2[:22]
    short_names.append(n2)

y = np.arange(n)

for col, (label, vals, warn_thr, crit_thr, color, _desc) in enumerate(metrics):
    ax = axes[col]
    norm = np.array(vals, dtype=float)
    vmax = max(crit_thr * 1.1, max(norm) * 1.05)
    
    colors = []
    for v in norm:
        if v >= crit_thr:
            colors.append('#d62728')
        elif v >= warn_thr:
            colors.append('#ff7f0e')
        else:
            colors.append('#2ca02c')
    
    ax.barh(y, norm, color=colors, alpha=0.85, height=0.7)
    ax.axvline(warn_thr, color='orange', lw=1, ls='--', alpha=0.7)
    ax.axvline(crit_thr, color='red', lw=1, ls='--', alpha=0.7)
    ax.set_xlim(0, vmax)
    ax.set_xlabel(label, fontsize=8)
    ax.set_yticks(y)
    
    # Draw group colour stripe on left axis only
    if col == 0:
        ax.set_yticklabels(short_names, fontsize=7.5)
        for tick, r in zip(ax.get_yticklabels(), results):
            tick.set_color(GROUP_COLOR[r['group']])
    else:
        ax.set_yticklabels([])
    
    ax.tick_params(axis='x', labelsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Annotate value
    for i, v in enumerate(norm):
        if v > 0:
            ax.text(min(v + vmax*0.02, vmax*0.95), i, f'{v:.0f}' if v >= 1 else f'{v:.2f}',
                    va='center', fontsize=6.5)

# Horizontal group separators (A=12, B=5, C=2, D=5)
boundaries = [12, 17, 19]  # indices where group changes
for ax in axes:
    for b in boundaries:
        ax.axhline(b - 0.5, color='grey', lw=0.8, ls=':')

# Legend
legend_handles = [
    mpatches.Patch(color='#2ca02c', label='OK'),
    mpatches.Patch(color='#ff7f0e', label='Warning'),
    mpatches.Patch(color='#d62728', label='Critical'),
]
group_handles = [mpatches.Patch(color=c, label=f'Group {g}')
                 for g, c in GROUP_COLOR.items()]
fig.legend(handles=legend_handles + group_handles,
           loc='lower center', ncol=7, fontsize=8,
           bbox_to_anchor=(0.5, -0.01))

fig.suptitle('Fig 92 — Sequence composition flags: 24 selected proteins',
             fontsize=12, fontweight='bold', y=1.005)
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/fig92_seq_composition.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved fig92')

# ── figure 93 — AA composition heatmap ────────────────────────────
mat = np.array([[r['pct'][aa] for aa in AAS] for r in results])

fig, ax = plt.subplots(figsize=(16, 9))
im = ax.imshow(mat, aspect='auto', cmap='YlOrRd', vmin=0, vmax=30)

ax.set_xticks(range(len(AAS)))
ax.set_xticklabels(AAS, fontsize=9)
ax.set_yticks(range(n))
ax.set_yticklabels(short_names, fontsize=8)
for tick, r in zip(ax.get_yticklabels(), results):
    tick.set_color(GROUP_COLOR[r['group']])

# Annotate cells ≥ 15%
for i in range(n):
    for j, aa in enumerate(AAS):
        v = mat[i, j]
        if v >= 15:
            ax.text(j, i, f'{v:.0f}', ha='center', va='center',
                    fontsize=6.5, color='black' if v < 22 else 'white', fontweight='bold')

# Group lines
for b in boundaries:
    ax.axhline(b - 0.5, color='white', lw=1.5)

plt.colorbar(im, ax=ax, shrink=0.6, label='Residue % composition')

# Group labels on right
group_info = [('Group A\n(naf top)', 0, 11), ('Group B\n(fd)', 12, 16),
              ('Group C\n(replicate)', 17, 18), ('Group D\n(seeds)', 19, 23)]
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks([])
for label, lo, hi in group_info:
    mid = (lo + hi) / 2
    ax2.text(1.01, mid / (n-1), label, transform=ax2.transAxes,
             va='center', fontsize=8, color=GROUP_COLOR[label[6]])

ax.set_title('Fig 93 — Amino acid composition heatmap: 24 selected proteins\n'
             '(bold values ≥ 15%; colour = group)', fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/fig93_aa_heatmap.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved fig93')

# ── print clean flag report ────────────────────────────────────────
print()
print("="*90)
print("DETAILED FLAG REPORT")
print("="*90)
for r in results:
    flags = []
    for aa, length, pos in r['runs']:
        flags.append(f"REPEAT {aa}×{length} at pos {pos}")
    if r['lc'] > 0.05:
        flags.append(f"low-complexity {r['lc']:.0%} of windows")
    if r['cys'] >= 2:
        flags.append(f"{r['cys']} cysteines (disulfide/folding risk)")
    if abs(r['nc']) >= 12:
        flags.append(f"net charge {r['nc']:+d} (high)")
    if r['pct']['A'] > 30:
        flags.append(f"Ala {r['pct']['A']:.0f}% (low-complexity coil risk)")
    if r['pct']['K'] > 20:
        flags.append(f"Lys {r['pct']['K']:.0f}% (aggregation/charge risk)")
    if r['pct']['E'] > 20:
        flags.append(f"Glu {r['pct']['E']:.0f}% (charge risk)")
    
    severity = "CRITICAL" if any("REPEAT" in f or "cystein" in f for f in flags) or r['lc'] > 0.3 else \
               "WARN" if flags else "OK"
    
    print(f"[{r['group']}] {r['name']:<32s}  {severity}")
    for f in flags:
        print(f"       • {f}")
    if not flags:
        print(f"       • no issues detected")
