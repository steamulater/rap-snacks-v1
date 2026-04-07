# BioReason-Pro Setup

## Install

```bash
# Python 3.11+ required
git clone https://github.com/bowang-lab/BioReason-Pro
cd BioReason-Pro
pip install -r requirements.txt
```

GPU required. If running on Colab, use A100 (same runtime as Boltz).

## Run on rap-snacks inputs

From the repo root:

```bash
python predict.py \
  --input  /path/to/rap-snacks-v1/outputs/bioreason/bioreason_all.tsv \
  --output /path/to/rap-snacks-v1/outputs/bioreason/bioreason_results.tsv
```

## Input file

`bioreason_all.tsv` — 118 sequences, 3 columns:

| protein_id | organism | sequence |
|---|---|---|
| bar_9__concordance__the_light_is_coming__ic0.707 | Homo sapiens | GEKLWLIC... |
| bar_9__native_ala__the_light_is_coming__ic0.707 | Homo sapiens | ANDEVERYTHING... |
| ... | ... | ... |

**Panels:**
- `bioreason_candidates.tsv` — 12 Phase 2 bars × 4 buckets (98 seqs)
- `bioreason_breadth.tsv`    — 20 bars spanning full iconicity range, concordance only
- `bioreason_all.tsv`        — both combined (118 seqs) ← use this

## Parse results

```bash
python analysis/12b_bioreason_parse.py \
  --results outputs/bioreason/bioreason_results.tsv
```

Outputs:
- `outputs/figures/fig_bioreason_confidence_violin.png` — confidence by bucket
- `outputs/figures/fig_bioreason_vs_plddt.png`         — confidence vs Boltz pLDDT
- `outputs/figures/fig_bioreason_vs_iconicity.png`      — confidence vs iconicity
- `outputs/bioreason/bioreason_top_go_terms.csv`        — top hallucinated GO terms
- `outputs/bioreason/bioreason_per_bar_summary.csv`     — per-bar summary

## What to look for

**If BioReason assigns HIGH confidence to lyric sequences:**
- The model is not grounded — it predicts function for any sequence regardless of foldability
- Compare confidence distribution to Boltz pLDDT: if r ≈ 0, predictions are sequence-agnostic

**If BioReason assigns LOW confidence or flags as disordered:**
- Concordance sequences structurally resemble nothing real (consistent with FoldSeek: 13 bars with 0 hits in Phase 1)
- Native_ala and free_design sequences may score higher — use as validation

**Key comparisons:**
- concordance vs free_design: does MPNN design increase BioReason confidence?
- concordance vs native_ala: does Ala substitution change predicted function?
- iconicity vs confidence: does cultural salience predict biological plausibility?
- bar_27 native_ala (0 FoldSeek hits, novel fold): does BioReason also flag as unknown?
- bar_6 concordance (TIM barrel, 2596 hits): does BioReason predict TIM barrel function?
- bar_46 (backbone failure): does BioReason predict low confidence for the dead backbone?
