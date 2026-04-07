"""
10_foldseek_phase2.py
---------------------
FoldSeek structural homolog search for Phase 2 validation structures.

Submits one PDB per bucket per bar to search.foldseek.com:
  - concordance  : best Boltz model from validation run
  - native_ala   : best Boltz model from validation run
  - free_design  : highest boltz_plddt design per bar

Searches: pdb100, afdb-swissprot, mgnify_esm30
Compares hit profiles across the three buckets.

This script is designed to run as a Colab cell (PDBs are on Drive).
Paste the content between the dashed lines into a new Colab cell.

Outputs (saved to Drive):
  results/foldseek_phase2/foldseek_phase2_hits.csv
  results/foldseek_phase2/foldseek_phase2_summary.csv
  results/figures/fig_foldseek_phase2_hits.png
  results/figures/fig_foldseek_phase2_db.png
"""

# ── COLAB CELL — paste everything below ───────────────────────────────────

import json, time, urllib.request, urllib.error, shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DARK_BG = '#0e0e0e'

# Paths (Drive)
DRIVE_ROOT   = Path('/content/drive/MyDrive/rap_snacks')
DRIVE_RES    = DRIVE_ROOT / 'results'
DRIVE_FIGS   = DRIVE_RES  / 'figures'
FS_DIR       = DRIVE_RES  / 'foldseek_phase2'
FS_DIR.mkdir(parents=True, exist_ok=True)

BOLTZ_OUT    = Path('/content/scratch/boltz_outputs')
RESULTS_CSV  = DRIVE_RES  / 'boltz_validation_results.csv'

# Restore Boltz outputs from Drive if scratch empty
if not BOLTZ_OUT.exists() or not any(BOLTZ_OUT.rglob('*.pdb')):
    print('Restoring Boltz outputs from Drive...')
    shutil.copytree(DRIVE_RES / 'boltz_outputs', BOLTZ_OUT, dirs_exist_ok=True)
    print(f'  {sum(1 for _ in BOLTZ_OUT.rglob("*.pdb"))} PDBs restored')

meta = pd.read_csv(RESULTS_CSV)
preds_dirs = list(BOLTZ_OUT.rglob('predictions'))
if not preds_dirs:
    raise RuntimeError('No predictions dir found. Restore from Drive first.')
PREDS_DIR = preds_dirs[0]
print(f'Predictions dir: {PREDS_DIR}')

# ── Pick one PDB per bucket per bar ───────────────────────────────────────
DATABASES = ['pdb100', 'afdb-swissprot', 'mgnify_esm30']
BUCKETS   = ['concordance', 'native_ala', 'free_design']
COLORS    = {'concordance': '#00d4ff', 'native_ala': '#ff9500', 'free_design': '#a855f7'}

def best_pdb(bar_id, bucket, preds_dir, meta_df):
    """Return path to best PDB for this bar/bucket."""
    grp = meta_df[(meta_df['bar_id'] == bar_id) & (meta_df['bucket'] == bucket)]
    if grp.empty:
        return None
    # Pick highest boltz_plddt; fall back to first row
    grp = grp.sort_values('boltz_plddt', ascending=False)
    for _, row in grp.iterrows():
        seq_dir  = preds_dir / row['name']
        pdb_list = sorted(seq_dir.glob('*.pdb')) if seq_dir.exists() else []
        if pdb_list:
            return pdb_list[0]   # model_0 is highest confidence
    return None

submissions = []  # (bar_id, bucket, pdb_path)
bar_ids = sorted(meta['bar_id'].unique())

for bar_id in bar_ids:
    for bucket in BUCKETS:
        pdb = best_pdb(bar_id, bucket, PREDS_DIR, meta)
        if pdb:
            submissions.append((bar_id, bucket, pdb))
        else:
            print(f'  [SKIP] {bar_id} {bucket} — no PDB found')

print(f'\n{len(submissions)} structures to submit to FoldSeek')
print(f'  ({len(bar_ids)} bars × {len(BUCKETS)} buckets)')

# ── FoldSeek REST API helpers ──────────────────────────────────────────────
SUBMIT_URL = 'https://search.foldseek.com/api/ticket'
RESULT_URL = 'https://search.foldseek.com/api/result/{ticket}/0'

def submit_foldseek(pdb_path: Path, databases: list) -> str | None:
    boundary = '----FoldSeekBoundary'
    body = b''
    file_data = pdb_path.read_bytes()

    body += f'--{boundary}\r\n'.encode()
    body += (
        f'Content-Disposition: form-data; name="q"; filename="{pdb_path.name}"\r\n'
        f'Content-Type: application/octet-stream\r\n\r\n'
    ).encode()
    body += file_data + b'\r\n'

    for db in databases:
        body += f'--{boundary}\r\n'.encode()
        body += f'Content-Disposition: form-data; name="database[]"\r\n\r\n{db}\r\n'.encode()

    body += f'--{boundary}\r\n'.encode()
    body += f'Content-Disposition: form-data; name="mode"\r\n\r\n3diaa\r\n'.encode()
    body += f'--{boundary}--\r\n'.encode()

    req = urllib.request.Request(
        SUBMIT_URL, data=body, method='POST',
        headers={'Content-Type': f'multipart/form-data; boundary={boundary}'},
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.loads(r.read())['id']
    except Exception as e:
        print(f'    Submit error: {e}')
        return None


def poll_foldseek(ticket: str, timeout=300, interval=5) -> dict | None:
    url = RESULT_URL.format(ticket=ticket)
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = json.loads(r.read())
                if data.get('status') == 'COMPLETE':
                    return data
                if data.get('status') == 'ERROR':
                    print(f'    FoldSeek error for ticket {ticket}')
                    return None
        except Exception:
            pass
        time.sleep(interval)
    print(f'    Timeout for ticket {ticket}')
    return None


def parse_hits(result: dict, bar_id: str, bucket: str) -> list[dict]:
    rows = []
    if not result or 'results' not in result:
        return rows
    for db_result in result['results']:
        db = db_result.get('db', '')
        for alignment in db_result.get('alignments', [[]])[0]:
            target = alignment.get('target', '')
            parts  = target.split(' ', 1)
            acc    = parts[0]
            desc   = parts[1] if len(parts) > 1 else ''
            rows.append({
                'bar_id':   bar_id,
                'bucket':   bucket,
                'db':       db,
                'accession': acc,
                'description': desc,
                'prob':     float(alignment.get('prob', 0)),
                'evalue':   float(alignment.get('evalue', 999)),
                'score':    float(alignment.get('score', 0)),
                'seqId':    float(alignment.get('seqId', 0)),
                'alnLen':   int(alignment.get('alnLen', 0)),
            })
    return rows


# ── Submit all, with rate-limit backoff ───────────────────────────────────
RAW_DIR = FS_DIR / 'raw'
RAW_DIR.mkdir(exist_ok=True)

all_hits = []
DELAY    = 8.0   # seconds between submissions (avoid 429)

for i, (bar_id, bucket, pdb_path) in enumerate(submissions):
    label     = f'{bar_id}_{bucket}'
    cache_path = RAW_DIR / f'{label}.json'

    if cache_path.exists():
        print(f'  [CACHED] {label}')
        result = json.loads(cache_path.read_text())
    else:
        print(f'  [{i+1}/{len(submissions)}] Submitting {label} ...', end=' ', flush=True)
        ticket = None
        for attempt in range(4):
            ticket = submit_foldseek(pdb_path, DATABASES)
            if ticket:
                break
            wait = DELAY * (2 ** attempt)
            print(f'retry in {wait:.0f}s...', end=' ', flush=True)
            time.sleep(wait)

        if not ticket:
            print('FAILED — skipping')
            continue

        result = poll_foldseek(ticket)
        if result:
            cache_path.write_text(json.dumps(result))
            print('done')
        else:
            print('no result')
            continue

        time.sleep(DELAY)

    hits = parse_hits(result, bar_id, bucket)
    all_hits.extend(hits)
    print(f'    → {len(hits)} hits')

# ── Save results ───────────────────────────────────────────────────────────
hits_df = pd.DataFrame(all_hits)
hits_csv = FS_DIR / 'foldseek_phase2_hits.csv'
hits_df.to_csv(hits_csv, index=False)
print(f'\nSaved {len(hits_df)} total hits → {hits_csv}')

# Summary: hits per bar per bucket, mean prob of top hit
summary_rows = []
for (bar_id, bucket), grp in hits_df.groupby(['bar_id', 'bucket']):
    pdb_hits = grp[grp['db'].str.contains('pdb', case=False)]
    top_prob = grp['prob'].max() if not grp.empty else 0
    top_desc = grp.sort_values('prob', ascending=False).iloc[0]['description'][:60] if not grp.empty else ''
    summary_rows.append({
        'bar_id': bar_id, 'bucket': bucket,
        'n_hits': len(grp), 'n_pdb_hits': len(pdb_hits),
        'top_prob': round(top_prob, 3), 'top_hit': top_desc,
    })

summary_df = pd.DataFrame(summary_rows)
summary_csv = FS_DIR / 'foldseek_phase2_summary.csv'
summary_df.to_csv(summary_csv, index=False)
print(f'Saved summary → {summary_csv}')

print('\n--- Top hits per bar per bucket ---')
for _, r in summary_df.sort_values(['bar_id', 'bucket']).iterrows():
    print(f'  {r["bar_id"]:8s}  {r["bucket"]:12s}  '
          f'n={r["n_hits"]:3d}  top_prob={r["top_prob"]:.3f}  {r["top_hit"][:50]}')

# ── Fig F: top hit probability per bar — concordance vs native_ala vs free_design ──
fig, axes = plt.subplots(1, 2, figsize=(15, 5), facecolor=DARK_BG)
bar_ids_s = sorted(summary_df['bar_id'].unique())
x = np.arange(len(bar_ids_s))
width = 0.25

# Fig F1: top prob (pdb100 only)
ax = axes[0]
ax.set_facecolor(DARK_BG)
pdb_hits = hits_df[hits_df['db'].str.contains('pdb100', case=False)]

for oi, bucket in enumerate(BUCKETS):
    vals = []
    for bar_id in bar_ids_s:
        grp = pdb_hits[(pdb_hits['bar_id'] == bar_id) & (pdb_hits['bucket'] == bucket)]
        vals.append(grp['prob'].max() if not grp.empty else 0)
    ax.bar(x + (oi - 1) * width, vals, width=width * 0.9,
           color=COLORS[bucket], alpha=0.85, label=bucket)

ax.set_xticks(x)
ax.set_xticklabels(bar_ids_s, rotation=45, ha='right', fontsize=8, color='white')
ax.set_ylabel('FoldSeek probability (pdb100)', color='white')
ax.set_ylim(0, 1)
ax.set_title('Top FoldSeek hit probability — pdb100', color='white', fontsize=11)
ax.tick_params(colors='white')
ax.axhline(0.5, color='white', linestyle='--', linewidth=0.8, alpha=0.4)
for s in ax.spines.values(): s.set_edgecolor('#333')
ax.legend(facecolor='#1a1a1a', labelcolor='white', fontsize=9)

# Fig F2: n_hits per bar per bucket (all databases)
ax2 = axes[1]
ax2.set_facecolor(DARK_BG)
for oi, bucket in enumerate(BUCKETS):
    vals = []
    for bar_id in bar_ids_s:
        grp = hits_df[(hits_df['bar_id'] == bar_id) & (hits_df['bucket'] == bucket)]
        vals.append(len(grp))
    ax2.bar(x + (oi - 1) * width, vals, width=width * 0.9,
            color=COLORS[bucket], alpha=0.85, label=bucket)

ax2.set_xticks(x)
ax2.set_xticklabels(bar_ids_s, rotation=45, ha='right', fontsize=8, color='white')
ax2.set_ylabel('Number of FoldSeek hits (all DBs)', color='white')
ax2.set_title('FoldSeek hit count — all databases', color='white', fontsize=11)
ax2.tick_params(colors='white')
for s in ax2.spines.values(): s.set_edgecolor('#333')
ax2.legend(facecolor='#1a1a1a', labelcolor='white', fontsize=9)

plt.tight_layout()
fig_f = Path('/content/scratch') / 'fig_foldseek_phase2.png'
plt.savefig(fig_f, dpi=150, facecolor=DARK_BG, bbox_inches='tight')
shutil.copy2(fig_f, DRIVE_FIGS / 'fig_foldseek_phase2.png')
plt.show()
print('Fig F saved → Drive: results/figures/fig_foldseek_phase2.png')

# ── Fig G: per-bar hit breakdown (which protein families?) ────────────────
# Show top 3 hit descriptions per bar per bucket as an annotation table
print('\n--- Top PDB hits per bar ---')
for bar_id in bar_ids_s:
    print(f'\n{bar_id}:')
    for bucket in BUCKETS:
        grp = pdb_hits[(pdb_hits['bar_id'] == bar_id) & (pdb_hits['bucket'] == bucket)]
        grp = grp.sort_values('prob', ascending=False).head(3)
        if grp.empty:
            print(f'  {bucket:12s} — no hits')
        else:
            for _, r in grp.iterrows():
                print(f'  {bucket:12s}  prob={r["prob"]:.3f}  {r["accession"]}  {r["description"][:55]}')
