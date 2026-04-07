#!/usr/bin/env python3
"""
Stage 16 — Sequence composition audit (same style as adaptyv_competition project).
9-panel audit figure + 6-panel cysteine figure.
Outputs: fig94_seq_audit.png, fig95_cysteine_analysis.png
"""
import re, math, csv
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

OUT_DIR  = "outputs/figures"
CSV_PATH = "outputs/selected_proteins.csv"

# SwissProt average frequencies
SWISSPROT = {'A':8.25,'R':5.53,'N':4.06,'D':5.45,'C':1.37,'Q':3.93,'E':6.75,
             'G':7.07,'H':2.27,'I':5.96,'L':9.66,'K':5.84,'M':2.42,'F':3.86,
             'P':4.70,'S':6.56,'T':5.34,'W':1.08,'Y':2.92,'V':6.87}

# ── load data ─────────────────────────────────────────────────────
records = []
with open(CSV_PATH) as f:
    for row in csv.DictReader(f):
        records.append(row)

def safe_float(v, default=0.0):
    try:
        return float(v)
    except (ValueError, TypeError):
        return default

# ── flag functions (copied from competition script) ───────────────

def find_homopolymers(seq, min_run=4):
    hits = []
    for m in re.finditer(r'(.)\1{%d,}' % (min_run - 1), seq):
        hits.append({'aa': m.group(1), 'start': m.start(), 'end': m.end(),
                     'len': len(m.group()), 'text': m.group()})
    return hits

def find_dipeptide_repeats(seq, min_reps=3):
    hits = []
    for m in re.finditer(r'(..)(\1{%d,})' % (min_reps - 1), seq):
        full = m.group(1) + m.group(2)
        hits.append({'motif': m.group(1), 'start': m.start(), 'end': m.end(),
                     'reps': len(full) // 2, 'text': full})
    return hits

def polyelectrostatic_stretch(seq, window=10, min_net=7):
    hits = []
    for i in range(len(seq) - window + 1):
        w = seq[i:i + window]
        pos = sum(c in 'KR' for c in w)
        neg = sum(c in 'DE' for c in w)
        if pos >= min_net:
            hits.append({'type': 'poly-K/R', 'start': i, 'end': i + window, 'text': w, 'count': pos})
        if neg >= min_net:
            hits.append({'type': 'poly-D/E', 'start': i, 'end': i + window, 'text': w, 'count': neg})
    merged = {}
    for h in hits:
        key = (h['type'], h['start'] // 5)
        if key not in merged or h['count'] > merged[key]['count']:
            merged[key] = h
    return list(merged.values())

def low_complexity_score(seq, window=20):
    min_entropy = 100
    worst_window = ''
    for i in range(max(1, len(seq) - window + 1)):
        w = seq[i:i + window]
        cnt = Counter(w)
        h = -sum((v / window) * math.log2(v / window) for v in cnt.values())
        if h < min_entropy:
            min_entropy = h
            worst_window = w
    return min_entropy, worst_window

def net_charge_at_pH7(seq):
    return seq.count('K') + seq.count('R') + seq.count('H') * 0.1 - seq.count('D') - seq.count('E')

# ── run audit ─────────────────────────────────────────────────────
results = []
for i, row in enumerate(records):
    seq   = row['sequence']
    name  = row['name']
    group = row['group']
    rank  = int(row['rank'])
    plddt = safe_float(row['boltz_plddt'])
    flags = []

    hp = find_homopolymers(seq, min_run=4)
    for h in hp:
        level = 'CRITICAL' if h['len'] >= 6 else 'WARNING'
        flags.append({'type': f'homopolymer_{h["aa"]}x{h["len"]}',
                      'level': level,
                      'detail': f'{h["aa"]}×{h["len"]} at pos {h["start"] + 1}'})

    dp = find_dipeptide_repeats(seq, min_reps=3)
    for d in dp:
        flags.append({'type': f'dipeptide_{d["motif"]}', 'level': 'WARNING',
                      'detail': f'({d["motif"]})×{d["reps"]} at pos {d["start"] + 1}'})

    pe = polyelectrostatic_stretch(seq)
    for p in pe:
        flags.append({'type': f'poly_{p["type"]}', 'level': 'WARNING',
                      'detail': f'{p["type"]} at pos {p["start"] + 1}: {p["text"]}'})

    ent, ww = low_complexity_score(seq, window=min(20, len(seq)))
    if ent < 2.5:
        level = 'CRITICAL' if ent < 1.8 else 'WARNING'
        flags.append({'type': 'low_complexity', 'level': level,
                      'detail': f'entropy={ent:.2f} (worst: {ww})'})

    n        = len(seq)
    pct_A    = seq.count('A') / n * 100
    pct_E    = seq.count('E') / n * 100
    pct_K    = seq.count('K') / n * 100
    pct_EK   = pct_E + pct_K
    pct_hydr = sum(seq.count(aa) for aa in 'VILMFYW') / n * 100
    cys      = seq.count('C')

    if pct_A > 25:
        flags.append({'type': 'extreme_A_bias', 'level': 'CRITICAL',
                      'detail': f'Ala={pct_A:.1f}% (expected ~8.3%)'})
    if pct_EK > 40:
        flags.append({'type': 'extreme_EK_bias', 'level': 'WARNING',
                      'detail': f'E+K={pct_EK:.1f}% (expected ~12.6%)'})
    if pct_hydr < 10:
        flags.append({'type': 'no_hydrophobic_core', 'level': 'WARNING',
                      'detail': f'Hydrophobic={pct_hydr:.1f}% (expected ~38%)'})
    if cys >= 4:
        flags.append({'type': 'high_cysteine', 'level': 'CRITICAL',
                      'detail': f'{cys} cysteines (disulfide/misfolding risk)'})
    elif cys >= 2:
        flags.append({'type': 'cysteine_pair', 'level': 'WARNING',
                      'detail': f'{cys} cysteines'})

    nc = net_charge_at_pH7(seq)
    if abs(nc) > 15:
        flags.append({'type': 'extreme_charge', 'level': 'WARNING',
                      'detail': f'net charge {nc:+.0f}'})

    results.append({
        'name': name, 'group': group, 'rank': rank, 'seq': seq,
        'length': n, 'plddt': plddt,
        'flags': flags,
        'n_critical': sum(1 for f in flags if f['level'] == 'CRITICAL'),
        'n_warning':  sum(1 for f in flags if f['level'] == 'WARNING'),
        'n_info':     0,
        'pct_A': pct_A, 'pct_EK': pct_EK, 'pct_hydr': pct_hydr,
        'entropy_min': ent, 'net_charge': nc, 'cys': cys,
        'pct_cys': cys / n * 100,
    })

# sorted by rank
results = sorted(results, key=lambda x: x['rank'])
n_seqs  = len(results)

# ── colours ───────────────────────────────────────────────────────
GROUP_COLOR = {'A': '#2166ac', 'B': '#4dac26', 'C': '#d01c8b', 'D': '#b35806'}

def rank_color(r):
    if r <= 8:  return '#E53935'
    if r <= 16: return '#FF9800'
    return '#9E9E9E'

# ══════════════════════════════════════════════════════════════════
# FIGURE 94 — 9-panel composition audit
# ══════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(18, 14))
gs  = GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.38)
fig.suptitle('Sequence Composition Audit — 24 Selected Rap-Protein Sequences',
             fontsize=13, fontweight='bold')

ranks      = [r['rank'] for r in results]
colors_r   = [rank_color(r['rank']) for r in results]
group_cols = [GROUP_COLOR[r['group']] for r in results]

# ── A: %Ala by rank ──
ax = fig.add_subplot(gs[0, 0])
ax.scatter(ranks, [r['pct_A'] for r in results],
           c=colors_r, s=55, alpha=0.85, edgecolors='white', lw=0.5, zorder=3)
ax.axhline(25, color='red', ls='--', lw=1.3, label='25% CRITICAL')
ax.axhline(8.25, color='grey', ls=':', lw=1, label='SwissProt avg 8.3%')
ax.set_xlabel('Selection rank'); ax.set_ylabel('% Alanine')
ax.set_title('Alanine content by rank', fontsize=10)
ax.legend(fontsize=7)
for r in results:
    if r['pct_A'] > 20:
        short = r['name'].replace('bar_', 'b').replace('_naf_', '_n').replace('_free_', '_f')[:12]
        ax.annotate(short, (r['rank'], r['pct_A']), fontsize=5.5, ha='left', va='bottom')

# ── B: %E+K by rank ──
ax2 = fig.add_subplot(gs[0, 1])
ax2.scatter(ranks, [r['pct_EK'] for r in results],
            c=colors_r, s=55, alpha=0.85, edgecolors='white', lw=0.5, zorder=3)
ax2.axhline(40, color='red', ls='--', lw=1.3, label='40% threshold')
ax2.axhline(12.6, color='grey', ls=':', lw=1, label='SwissProt avg 12.6%')
ax2.set_xlabel('Selection rank'); ax2.set_ylabel('% E + K')
ax2.set_title('Glutamate + Lysine content by rank', fontsize=10)
ax2.legend(fontsize=7)
for r in results:
    if r['pct_EK'] > 40:
        short = r['name'].replace('bar_', 'b')[:12]
        ax2.annotate(short, (r['rank'], r['pct_EK']), fontsize=5.5, ha='left', va='bottom')

# ── C: min Shannon entropy ──
ax3 = fig.add_subplot(gs[0, 2])
ax3.scatter(ranks, [r['entropy_min'] for r in results],
            c=colors_r, s=55, alpha=0.85, edgecolors='white', lw=0.5, zorder=3)
ax3.axhline(1.8, color='red', ls='--', lw=1.3, label='CRITICAL <1.8 bits')
ax3.axhline(2.5, color='orange', ls='--', lw=1.3, label='WARNING <2.5 bits')
ax3.set_xlabel('Selection rank'); ax3.set_ylabel('Min window entropy (bits)')
ax3.set_title('Sequence complexity (Shannon entropy)', fontsize=10)
ax3.set_ylim(0, 5)
ax3.legend(fontsize=7)
for r in results:
    if r['entropy_min'] < 2.5:
        short = r['name'].replace('bar_', 'b')[:12]
        ax3.annotate(short, (r['rank'], r['entropy_min']), fontsize=5.5, ha='left', va='top')

# ── D: % hydrophobic ──
ax4 = fig.add_subplot(gs[1, 0])
ax4.scatter(ranks, [r['pct_hydr'] for r in results],
            c=colors_r, s=55, alpha=0.85, edgecolors='white', lw=0.5, zorder=3)
ax4.axhline(10, color='red', ls='--', lw=1.3, label='<10% WARNING')
ax4.axhline(38, color='grey', ls=':', lw=1, label='SwissProt avg 38%')
ax4.set_xlabel('Selection rank'); ax4.set_ylabel('% hydrophobic (V,I,L,M,F,Y,W)')
ax4.set_title('Hydrophobic residue content', fontsize=10)
ax4.legend(fontsize=7)

# ── E: net charge ──
ax5 = fig.add_subplot(gs[1, 1])
ax5.scatter(ranks, [r['net_charge'] for r in results],
            c=colors_r, s=55, alpha=0.85, edgecolors='white', lw=0.5, zorder=3)
ax5.axhline(0, color='black', lw=0.8)
ax5.axhline(15, color='red', ls=':', lw=0.9, label='|charge|>15 risky')
ax5.axhline(-15, color='red', ls=':', lw=0.9)
ax5.set_xlabel('Selection rank'); ax5.set_ylabel('Net charge at pH 7')
ax5.set_title('Net charge (pos=K/R/H, neg=D/E)', fontsize=10)
ax5.legend(fontsize=7)
for r in results:
    if abs(r['net_charge']) > 12:
        short = r['name'].replace('bar_', 'b')[:12]
        ax5.annotate(short, (r['rank'], r['net_charge']), fontsize=5.5, ha='left')

# ── F: flag counts stacked bar (all 24) ──
ax6 = fig.add_subplot(gs[1, 2])
xs     = ranks
crits  = [r['n_critical'] for r in results]
warns  = [r['n_warning']  for r in results]
infos  = [r['n_info']     for r in results]
ax6.bar(xs, crits, color='#E53935', label='Critical', width=0.7)
ax6.bar(xs, warns, bottom=crits, color='#FF9800', label='Warning', width=0.7)
ax6.bar(xs, infos, bottom=[c + w for c, w in zip(crits, warns)],
        color='#90A4AE', label='Info', width=0.7)
ax6.set_xlabel('Selection rank'); ax6.set_ylabel('Flag count')
ax6.set_title('Issue flags — all 24 designs', fontsize=10)
ax6.legend(fontsize=7)
# annotate group
for r in results:
    if r['n_critical'] + r['n_warning'] > 0:
        ax6.text(r['rank'], r['n_critical'] + r['n_warning'] + 0.05,
                 r['group'], ha='center', va='bottom', fontsize=6,
                 color=GROUP_COLOR[r['group']], fontweight='bold')

# ── G: max homopolymer run bar ──
ax7 = fig.add_subplot(gs[2, 0])
max_runs = []
run_aas  = []
for r in results:
    hp = find_homopolymers(r['seq'], min_run=1)
    if hp:
        best = max(hp, key=lambda h: h['len'])
        max_runs.append(best['len'])
        run_aas.append(best['aa'])
    else:
        max_runs.append(1)
        run_aas.append('-')
colors_hp = ['#E53935' if l >= 6 else '#FF9800' if l >= 4 else '#4CAF50' for l in max_runs]
ax7.bar(ranks, max_runs, color=colors_hp, width=0.7)
ax7.axhline(6, color='red', ls='--', lw=1.3, label='≥6: CRITICAL')
ax7.axhline(4, color='orange', ls='--', lw=1.3, label='≥4: WARNING')
for x, l, aa in zip(ranks, max_runs, run_aas):
    if l >= 4:
        ax7.text(x, l + 0.15, aa, ha='center', va='bottom', fontsize=7, fontweight='bold')
ax7.set_xlabel('Selection rank'); ax7.set_ylabel('Max homopolymer run (AA)')
ax7.set_title('Longest homopolymer run', fontsize=10)
ax7.legend(fontsize=7)

# ── H: AA composition heatmap (all 24) ──
ax8 = fig.add_subplot(gs[2, 1:])
AAS = list('ACDEFGHIKLMNPQRSTVWY')
mat = np.zeros((n_seqs, len(AAS)))
for i, r in enumerate(results):
    n = len(r['seq'])
    for j, aa in enumerate(AAS):
        mat[i, j] = r['seq'].count(aa) / n * 100

short_labels = []
for r in results:
    nm = r['name']
    if nm == '6E5C_positive_control':
        short_labels.append(f"#{r['rank']} 6E5C(+ctrl)")
    else:
        parts = nm.split('_')
        short = f"#{r['rank']} {parts[0]}_{parts[1]}_{('_'.join(parts[2:]))[:8]}"
        short_labels.append(short)

im = ax8.imshow(mat, aspect='auto', cmap='YlOrRd', vmin=0, vmax=35)
ax8.set_xticks(range(len(AAS))); ax8.set_xticklabels(list(AAS), fontsize=8)
ax8.set_yticks(range(n_seqs)); ax8.set_yticklabels(short_labels, fontsize=7)
for tick, r in zip(ax8.get_yticklabels(), results):
    tick.set_color(GROUP_COLOR[r['group']])
# Annotate high values
for i in range(n_seqs):
    for j in range(len(AAS)):
        if mat[i, j] >= 20:
            ax8.text(j, i, f'{mat[i,j]:.0f}', ha='center', va='center',
                     fontsize=5.5, color='black' if mat[i,j] < 28 else 'white', fontweight='bold')
ax8.set_title('AA composition heatmap — all 24 sequences (% per residue)', fontsize=10)
plt.colorbar(im, ax=ax8, fraction=0.015, pad=0.02, label='%')

# Group boundary lines
for b in [12, 17, 19]:  # after rank 12, 17, 19
    ax8.axhline(b - 0.5, color='white', lw=1.5)

# Rank legend
patches = [mpatches.Patch(color='#E53935', label='Rank 1–8'),
           mpatches.Patch(color='#FF9800', label='Rank 9–16'),
           mpatches.Patch(color='#9E9E9E', label='Rank 17–24')]
group_patches = [mpatches.Patch(color=c, label=f'Group {g}')
                 for g, c in GROUP_COLOR.items()]
fig.legend(handles=patches + group_patches, loc='lower left', fontsize=8, ncol=7,
           bbox_to_anchor=(0.01, -0.01))

plt.savefig(f'{OUT_DIR}/fig94_seq_audit.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved fig94')

# ══════════════════════════════════════════════════════════════════
# FIGURE 95 — 6-panel cysteine analysis
# ══════════════════════════════════════════════════════════════════
fig2, axes = plt.subplots(2, 3, figsize=(15, 9))
fig2.suptitle('Cysteine Composition in 24 Selected Rap-Protein Sequences',
              fontsize=13, fontweight='bold')

groups = ['A', 'B', 'C', 'D']
gdata  = {g: [r for r in results if r['group'] == g] for g in groups}
gcols  = [GROUP_COLOR[g] for g in groups]
jitter = np.random.default_rng(42)

# ── panel 1: Cys% per design (strip plot by group) ──
ax = axes[0, 0]
for gi, g in enumerate(groups):
    ys = [r['pct_cys'] for r in gdata[g]]
    xs = jitter.uniform(gi - 0.25, gi + 0.25, len(ys))
    ax.scatter(xs, ys, color=GROUP_COLOR[g], s=45, alpha=0.8, zorder=3)
    ax.hlines(np.mean(ys), gi - 0.3, gi + 0.3, colors=GROUP_COLOR[g], lw=2, zorder=4)
    ax.text(gi, np.mean(ys) + 0.05, f'avg {np.mean(ys):.2f}%', ha='center',
            fontsize=7, color=GROUP_COLOR[g])
ax.axhline(1.37, color='grey', ls=':', lw=1, label='SwissProt avg 1.37%')
ax.set_xticks(range(4)); ax.set_xticklabels([f'Group {g}' for g in groups])
ax.set_ylabel('% Cysteine'); ax.set_title('Cys% per design (bar=mean)', fontsize=10)
ax.legend(fontsize=7)

# ── panel 2: Cysteine count distribution ──
ax2 = axes[0, 1]
from collections import Counter as C2
max_cys = max(r['cys'] for r in results)
xs_all  = range(max_cys + 1)
bottom = np.zeros(max_cys + 1)
for g in groups:
    counts = C2(r['cys'] for r in gdata[g])
    heights = np.array([counts.get(x, 0) for x in xs_all], dtype=float)
    n_g = len(gdata[g])
    ax2.bar(xs_all, heights / n_g, bottom=bottom[:len(xs_all)], color=GROUP_COLOR[g],
            alpha=0.85, label=f'Group {g}', width=0.6)
    bottom[:len(xs_all)] += heights / n_g
ax2.set_xlabel('Number of cysteines per design'); ax2.set_ylabel('Fraction')
ax2.set_title('Cysteine count distribution', fontsize=10)
ax2.legend(fontsize=7)

# ── panel 3: Cys% vs pLDDT ──
ax3 = axes[0, 2]
for g in groups:
    xs3 = [r['pct_cys'] for r in gdata[g]]
    ys3 = [r['plddt'] for r in gdata[g] if r['plddt'] > 0]
    xs3 = [r['pct_cys'] for r in gdata[g] if r['plddt'] > 0]
    ax3.scatter(xs3, ys3, color=GROUP_COLOR[g], s=45, alpha=0.8, label=f'Group {g}', zorder=3)
ax3.set_xlabel('% Cysteine'); ax3.set_ylabel('Boltz-2 pLDDT')
ax3.set_title('Cys% vs pLDDT', fontsize=10)
ax3.legend(fontsize=7)

# ── panel 4: Cys count vs pLDDT ──
ax4 = axes[1, 0]
for g in groups:
    xs4 = [r['cys'] + jitter.uniform(-0.15, 0.15) for r in gdata[g] if r['plddt'] > 0]
    ys4 = [r['plddt'] for r in gdata[g] if r['plddt'] > 0]
    ax4.scatter(xs4, ys4, color=GROUP_COLOR[g], s=45, alpha=0.8, label=f'Group {g}', zorder=3)
ax4.set_xlabel('Cysteine count'); ax4.set_ylabel('Boltz-2 pLDDT')
ax4.set_title('Cys count vs pLDDT', fontsize=10)
ax4.legend(fontsize=7)

# ── panel 5: Mean Cys% ± SD by group ──
ax5 = axes[1, 1]
means = [np.mean([r['pct_cys'] for r in gdata[g]]) for g in groups]
sds   = [np.std([r['pct_cys'] for r in gdata[g]]) for g in groups]
bars5 = ax5.bar(groups, means, color=gcols, alpha=0.8,
                yerr=sds, capsize=5, error_kw={'lw': 1.5})
ax5.axhline(1.37, color='grey', ls=':', lw=1, label='SwissProt avg')
for i, (m, g) in enumerate(zip(means, groups)):
    ax5.text(i, m + sds[i] + 0.05, f'{m:.2f}%', ha='center', fontsize=8)
ax5.set_ylabel('Mean Cys%'); ax5.set_title('Mean Cys% ±SD by group', fontsize=10)
ax5.set_xticklabels([f'Group {g}' for g in groups])
ax5.legend(fontsize=7)

# ── panel 6: Fraction with ≥1 cysteine ──
ax6 = axes[1, 2]
fracs = [sum(1 for r in gdata[g] if r['cys'] >= 1) / len(gdata[g]) * 100 for g in groups]
bars6 = ax6.bar(groups, fracs, color=gcols, alpha=0.8)
for i, (f, g) in enumerate(zip(fracs, groups)):
    ax6.text(i, f + 1, f'{f:.0f}%', ha='center', fontsize=9, fontweight='bold')
ax6.set_ylim(0, 115)
ax6.set_ylabel('% designs with ≥1 Cys')
ax6.set_title('Fraction with any cysteine', fontsize=10)
ax6.set_xticklabels([f'Group {g}' for g in groups])

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/fig95_cysteine_analysis.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved fig95')

# ── print summary ─────────────────────────────────────────────────
print()
critical = [r for r in results if r['n_critical'] > 0]
warned   = [r for r in results if r['n_warning'] > 0 and r['n_critical'] == 0]
clean    = [r for r in results if r['n_critical'] == 0 and r['n_warning'] == 0]
print(f"{'='*70}")
print(f"AUDIT SUMMARY — 24 sequences")
print(f"{'='*70}")
print(f"  CRITICAL: {len(critical)}  |  WARNING: {len(warned)}  |  CLEAN: {len(clean)}")
print()
for r in sorted(results, key=lambda x: -(x['n_critical']*10 + x['n_warning'])):
    if not r['flags']: continue
    print(f"[{r['group']}] #{r['rank']:>2} {r['name']:<32} — {r['n_critical']}C {r['n_warning']}W")
    for f in r['flags']:
        marker = '🚨' if f['level'] == 'CRITICAL' else '⚠️ '
        print(f"       {marker} {f['type']}: {f['detail']}")
