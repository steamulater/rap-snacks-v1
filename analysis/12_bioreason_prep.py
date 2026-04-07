"""
12_bioreason_prep.py
--------------------
Prepare BioReason-Pro input TSV from all four sequence buckets across the 12
Phase 2 candidate bars, plus a breadth panel sampling the full 85-bar set at
high / mid / low iconicity.

Columns (BioReason-Pro format):
  protein_id   — unique identifier encoding bar, bucket, song, iconicity
  organism     — "Homo sapiens" (consistent; triggers human annotation)
  sequence     — standard-AA protein sequence

Buckets included:
  concordance        — raw lyric-derived (probabilistic BJOZXU substitution)
  native_ala         — BJOZXU → Ala (lyric text literally readable)
  free_design        — MPNN, concordance backbone, no fixed positions
  native_ala_free    — MPNN, native_ala backbone, no fixed positions (Run 4)

Outputs:
  outputs/bioreason/bioreason_candidates.tsv   ← 12-bar × 4-bucket panel (48 seqs)
  outputs/bioreason/bioreason_breadth.tsv      ← full-range iconicity sample (20 seqs)
  outputs/bioreason/bioreason_all.tsv          ← both combined + controls

Usage:
  python analysis/12_bioreason_prep.py

Then run BioReason-Pro:
  python predict.py --input outputs/bioreason/bioreason_all.tsv \\
                    --output outputs/bioreason/bioreason_results.tsv
"""

import pandas as pd
import numpy as np
from pathlib import Path

ROOT      = Path(__file__).parent.parent
OUT_DIR   = ROOT / "outputs/bioreason"
OUT_DIR.mkdir(parents=True, exist_ok=True)

ENRICHED    = ROOT / "data/aggregated_lines_v2_enriched.csv"
CANDIDATES  = ROOT / "data/phase2_candidates.csv"
FREE_CSV    = ROOT / "outputs/proteinmpnn_free/filtered_results.csv"
NAF_CSV     = ROOT / "outputs/proteinmpnn_native_ala_free/filtered_results.csv"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# ── Load data ─────────────────────────────────────────────────────────────

enriched    = pd.read_csv(ENRICHED)
candidates  = pd.read_csv(CANDIDATES)
free_df     = pd.read_csv(FREE_CSV)     if FREE_CSV.exists()  else pd.DataFrame()
naf_df      = pd.read_csv(NAF_CSV)      if NAF_CSV.exists()   else pd.DataFrame()

cand_ids = list(candidates["bar_id"])
cand_set = set(cand_ids)

def clean_seq(s):
    """Return uppercase sequence if all-standard-AA, else None."""
    if pd.isna(s):
        return None
    s = str(s).upper().strip()
    return s if set(s) <= VALID_AA else None

def make_id(bar_id, bucket, song, iconicity, idx=None):
    """Build a BioReason-compatible protein_id string."""
    song_slug = song.replace(" ", "_").replace("/", "-")[:30] if pd.notna(song) else "unknown"
    icon_str  = f"ic{iconicity:.3f}" if pd.notna(iconicity) else "ic?"
    suffix    = f"_{idx:03d}" if idx is not None else ""
    return f"{bar_id}__{bucket}__{song_slug}__{icon_str}{suffix}"

rows = []    # final TSV rows: protein_id, organism, sequence, _meta (dropped before save)

# ── 1. Candidate bars — all four buckets ──────────────────────────────────

for bar_id in cand_ids:
    sub = enriched[enriched["bar_id"] == bar_id]
    if sub.empty:
        continue
    r         = sub.iloc[0]
    song      = r.get("genius_song_title", "")
    iconicity = r.get("aggregate_iconicity", float("nan"))
    url       = r.get("genius_url", "")

    # concordance
    seq = clean_seq(r.get("fasta_seq_concordance"))
    if seq:
        rows.append({
            "protein_id": make_id(bar_id, "concordance", song, iconicity),
            "organism":   "Homo sapiens",
            "sequence":   seq,
            "_bar_id":    bar_id,
            "_bucket":    "concordance",
            "_iconicity": iconicity,
            "_url":       url,
            "_panel":     "candidates",
        })

    # native_ala
    seq = clean_seq(r.get("fasta_seq_native_alanine"))
    if seq:
        rows.append({
            "protein_id": make_id(bar_id, "native_ala", song, iconicity),
            "organism":   "Homo sapiens",
            "sequence":   seq,
            "_bar_id":    bar_id,
            "_bucket":    "native_ala",
            "_iconicity": iconicity,
            "_url":       url,
            "_panel":     "candidates",
        })

    # free_design — best ESMFold pLDDT design per bar
    if not free_df.empty:
        grp = free_df[free_df["bar_id"] == bar_id].copy()
        if not grp.empty:
            best = grp.sort_values("esm_plddt", ascending=False).iloc[0]
            seq  = clean_seq(best.get("sequence"))
            if seq:
                rows.append({
                    "protein_id": make_id(bar_id, "free_design", song, iconicity),
                    "organism":   "Homo sapiens",
                    "sequence":   seq,
                    "_bar_id":    bar_id,
                    "_bucket":    "free_design",
                    "_iconicity": iconicity,
                    "_url":       url,
                    "_panel":     "candidates",
                })

    # native_ala_free — best ESMFold pLDDT design per bar
    if not naf_df.empty:
        grp = naf_df[naf_df["bar_id"] == bar_id].copy()
        if not grp.empty:
            best = grp.sort_values("esm_plddt", ascending=False).iloc[0]
            seq  = clean_seq(best.get("sequence"))
            if seq:
                rows.append({
                    "protein_id": make_id(bar_id, "native_ala_free", song, iconicity),
                    "organism":   "Homo sapiens",
                    "sequence":   seq,
                    "_bar_id":    bar_id,
                    "_bucket":    "native_ala_free",
                    "_iconicity": iconicity,
                    "_url":       url,
                    "_panel":     "candidates",
                })

candidates_df = pd.DataFrame(rows)
tsv_cols = ["protein_id", "organism", "sequence"]
candidates_df[tsv_cols].to_csv(OUT_DIR / "bioreason_candidates.tsv", sep="\t", index=False)
print(f"Candidates panel: {len(candidates_df)} sequences → bioreason_candidates.tsv")

# ── 2. Breadth panel — full iconicity range, concordance only ─────────────
# Sample ~20 bars spanning high→low iconicity from the full 85-bar set,
# excluding bars already in the candidate panel.

# Filter to rows that actually have a concordance sequence
all_bars  = enriched[enriched["fasta_seq_concordance"].notna()].copy()
non_cands = all_bars[~all_bars["bar_id"].isin(cand_set)].copy()
non_cands = non_cands.dropna(subset=["aggregate_iconicity"]).sort_values(
    "aggregate_iconicity", ascending=False
)

# Evenly spaced across iconicity range
n_breadth = 20
indices = np.linspace(0, len(non_cands) - 1, n_breadth, dtype=int)
breadth_rows = []
for idx in indices:
    r         = non_cands.iloc[idx]
    bar_id    = r["bar_id"]
    song      = r.get("genius_song_title", "")
    iconicity = r.get("aggregate_iconicity", float("nan"))
    url       = r.get("genius_url", "")

    seq = clean_seq(r.get("fasta_seq_concordance"))
    if seq:
        breadth_rows.append({
            "protein_id": make_id(bar_id, "concordance", song, iconicity),
            "organism":   "Homo sapiens",
            "sequence":   seq,
            "_bar_id":    bar_id,
            "_bucket":    "concordance",
            "_iconicity": iconicity,
            "_url":       url,
            "_panel":     "breadth",
        })

breadth_df = pd.DataFrame(breadth_rows)
breadth_df[tsv_cols].to_csv(OUT_DIR / "bioreason_breadth.tsv", sep="\t", index=False)
print(f"Breadth panel:    {len(breadth_df)} sequences → bioreason_breadth.tsv")

# ── 3. Combined output ────────────────────────────────────────────────────

all_df = pd.concat([candidates_df, breadth_df], ignore_index=True)

# Drop any duplicate sequences (concordance appears in both panels for some bars)
all_df = all_df.drop_duplicates(subset=["protein_id"])

all_df[tsv_cols].to_csv(OUT_DIR / "bioreason_all.tsv", sep="\t", index=False)
print(f"Combined:         {len(all_df)} sequences → bioreason_all.tsv")

# ── 4. Summary print ──────────────────────────────────────────────────────

print("\n--- Candidate panel breakdown ---")
for bucket, grp in candidates_df.groupby("_bucket"):
    print(f"  {bucket:20s}  n={len(grp)}")

print("\n--- Breadth panel iconicity range ---")
print(f"  max:  {breadth_df['_iconicity'].max():.3f}  ({breadth_df.sort_values('_iconicity', ascending=False).iloc[0]['_bar_id']})")
print(f"  min:  {breadth_df['_iconicity'].min():.3f}  ({breadth_df.sort_values('_iconicity').iloc[0]['_bar_id']})")
print(f"  mean: {breadth_df['_iconicity'].mean():.3f}")

print("\n--- Sequence length stats ---")
all_df["_len"] = all_df["sequence"].str.len()
print(f"  min={all_df['_len'].min()}  max={all_df['_len'].max()}  mean={all_df['_len'].mean():.0f}")

print(f"\nReady for BioReason-Pro:")
print(f"  python predict.py \\")
print(f"    --input  outputs/bioreason/bioreason_all.tsv \\")
print(f"    --output outputs/bioreason/bioreason_results.tsv")
