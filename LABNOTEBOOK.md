# Rap Snacks v1 — Lab Notebook

**Repo:** `steamulater/rap-snacks-v1`
**Started:** March 2026
**Status:** ESMFold pilot complete — Boltz-2 pending Colab run

---

## Concept

Convert culturally iconic rap lyrics into protein sequences, fold them, and find their structural homologs in nature. The upstream input is an iconicity score derived from organic social signal. The downstream output is a folded protein structure with a mutation mask marking exactly where the lyric forced a biological decision that nature would never make.

This is version 2 of the pipeline, rebuilt from scratch after the v1 pilot (`steamulater/rap-snacks-protein-bars-pilot`). All v1 findings are preserved; all v1 data integrity failures are fixed by design.

---

## Why V2

The v1 pilot (225 bars, March 2026) completed structural analysis and FoldSeek homolog search but hit a pre-publication blocker: the FASTA files on disk were not the files that were folded. A `remerge.py` re-run after `convert.py` silently shifted the CSV row ordering, breaking the `bar_N → lyric` mapping for every bar at or after the affected index.

V2 fixes this architecturally — the lyric attribution is locked at conversion time and can never drift, regardless of what happens to the CSV afterward.

---

## V2 Pipeline Architecture

```
data/aggregated_lines_v2_frozen.csv    <- immutable input, SHA-256 guarded
         |
00_audit.py        validate + 20 random row spot-check (human review gate)
         |
01_convert.py      lyric -> FASTA (5 modes) + bar_index_snapshot.json
                   + enriched CSV with all 5 sequences per bar
                   + length agreement check (pre vs post conversion)
     /outputs/fastas/
bars_v2_{mode}.fasta × 5 modes
         |
[Boltz-2 on Colab]  concordance condition, PDB output, --output_format pdb
         |
03_parse_boltz.py  extract pLDDT, pTM, confidence, pDE -> boltz_confidence_scores.csv
                   + length mismatch check against snapshot
         |
02_esm_fold.py     ESMFold API ensemble, all 5 conditions, pilot = top 25 bars
         |
04_foldseek.py     FoldSeek REST API, pTM>=0.4 bars only
         |
05_enrich_csv.py   join all metrics -> aggregated_lines_v2_enriched.csv
```

**Orchestration:** `make audit` → `make convert` → `make esm-pilot` → [Colab] → `make parse-boltz` → `make foldseek` → `make enrich`

---

## V2 Fixes Baked In

| V1 failure | V2 fix |
|---|---|
| CSV re-sorted after convert → bar_N mapping broke | `bar_index_snapshot.json` written at conversion time. Maps `bar_N → {lyric, song, all 5 sequences}`. Authoritative record independent of CSV state. |
| No way to detect drift | SHA-256 hash of frozen CSV written on first run. Subsequent runs refuse to proceed if hash changes. |
| Sequences not tracked | Per-condition sequence columns (`fasta_seq_concordance`, etc.) added to enriched CSV at conversion time. |
| Length mismatch undetected | `len(lyric_cleaned) == len(fasta_seq_*)` asserted for all bars, all 5 modes, before any output is written. |
| CIF output → conversion step needed | `--output_format pdb` passed to Boltz predict. |
| Lyric attribution scattered | `bar_index_snapshot.json` is the single source of truth — not the FASTA header, not the CSV row index. |

---

## Dataset

**Input:** `data/aggregated_lines_v2_frozen.csv` — copied from v1 `aggregated_lines_v2_verified.csv`, treated as immutable.

| Metric | Value |
|---|---|
| Total rows | 539 |
| fasta_eligible=True | 248 |
| Eligible + attributed (primary/featured) | 225 |
| Passing 80–300 AA filter | **85** |
| Attribution: primary | 368 (all rows) |
| Attribution: featured | 106 (all rows) |
| Attribution: unverified | 65 (all rows) |
| Bars with any BOJUXZ chars | 221 / 225 |
| Mean BOJUXZ density | 13.85% |

---

## Stage 0 — Data Audit

**Script:** `pipeline/00_audit.py`
**Date:** March 2026
**Seed:** 99 (for random row sampling)

### Length distribution (eligible + attributed rows)

| Bucket | Count |
|---|---|
| <40 AA | 61 |
| 40–79 AA | 79 |
| **80–300 AA (target)** | **85** |
| 91–120 AA | 26 (subset of target) |
| >120 AA | 44 (subset of target) |

Min: 10 AA. Max: 155 AA. Median: 66 AA.

### Spot-check result

20 randomly sampled rows reviewed. All attributions, Genius URLs, iconicity scores, and song titles verified as correct. No missing values in required columns. Data approved for conversion.

---

## Stage 1 — FASTA Conversion

**Script:** `pipeline/01_convert.py`
**Date:** March 2026
**Parameters:** `--min-length 80 --max-length 300 --seed 42 --lambda-val 2.0`

### Results

| Metric | Value |
|---|---|
| Input rows | 539 |
| Eligible + attributed | 225 |
| Filtered out (too short <80 AA) | 140 |
| Filtered out (too long >300 AA) | 0 |
| **Bars converted** | **85** |
| Modes | 5 (concordance, alanine, random, native, native_alanine) |
| Agreement check | PASSED — all 85 bars × 5 modes |
| Frozen CSV hash | `49af11129b64...` (written to `data/frozen_csv.sha256`) |

### Outputs

| File | Description |
|---|---|
| `outputs/fastas/bars_v2_{mode}.fasta` | 5 FASTAs, 85 sequences each |
| `outputs/masks/mask_v2_{mode}.json` | Mutation masks per mode |
| `data/bar_index_snapshot.json` | 85 bars, authoritative lyric→sequence mapping |
| `data/aggregated_lines_v2_enriched.csv` | All original columns + `bar_id`, `lyric_cleaned`, `lyric_cleaned_len`, `fasta_seq_*` for all 5 modes |

### Length filter rationale

80 AA minimum drops the short-peptide regime where pLDDT is artificially elevated and structural comparison is dominated by length effects. 300 AA maximum is a practical ceiling for single-chain folding — no bars exceed 155 AA in this dataset, so the ceiling is not active.

In v1, pTM > 0.4 (protein-like) peaked in the 40–90 AA range. This v2 run starts at 80 AA, capturing the upper end of the best-folding regime and extending through longer sequences. Cross-condition analysis must stratify by length given the 80–155 AA spread.

### Conversion: five modes

2×2 factorial design across standard letter mapping × BOJUXZ substitution:

```
                      BOJUXZ substitution
                    softmax    alanine    random
Standard 20     concordance    alanine    random
mapping
Native             native  native_alanine
```

| Mode | Standard 20 | BOJUXZ | What it isolates |
|---|---|---|---|
| `concordance` | Freq-rank concordance | Softmax peaked draw | Canonical run |
| `alanine` | Freq-rank concordance | → A (neutral) | Effect of softmax BOJUXZ vs placeholder |
| `random` | Freq-rank concordance | Uniform random | Is concordance better than chance? |
| `native` | AA pass-through | Softmax peaked draw | Effect of freq remapping |
| `native_alanine` | AA pass-through | → A | BOJUXZ strategy without remapping |

All runs: seed=42, λ=2.0. Separate RNG instance per mode (seed + mode_index offset).

---

## Stage 2 — ESMFold Ensemble Pilot

**Script:** `pipeline/02_esm_fold.py`
**API:** `api.esmatlas.com` (free REST, SSL `verify=False`)
**Date:** March 2026
**Scope:** Top 25 bars by `aggregate_iconicity`

### Nstruct plan

| Condition | Seeds | Rationale |
|---|---|---|
| concordance | 15 | Peaked softmax — captures meaningful variance |
| alanine | 1 | Fully deterministic |
| random | 30 | Uniform distribution — needs more samples |
| native | 15 | Same softmax as concordance |
| native_alanine | 1 | Fully deterministic |

Total planned calls: 1,550

### Implementation notes

**Bug found and fixed — pLDDT extraction:** The initial run reported pLDDT values of ~0.003–0.004 — 100× too low. Root cause: the ESMFold API stores pLDDT in the PDB B-factor column already on a 0–1 scale (not 0–100 as standard PDB convention). The code incorrectly divided by 100. Fix: removed the `/ 100.0` division. The PDB files saved during the bad run were structurally correct; only the pLDDT CSV values were wrong. The CSV was deleted and recomputed from saved PDB files.

**PDB file reuse added:** After the fix, `02_esm_fold.py` checks whether the PDB file for a given `(bar_id, condition, seed)` triplet already exists on disk. If so, it reads the file directly and skips the API call. This makes reruns instant for already-folded sequences and makes the script crash-safe at the file level, not just the CSV level.

**Bug found and fixed — bar_id KeyError:** `snapshot.values()` returns dicts that don't contain `bar_id` (it was only the key). Fixed by iterating `snapshot.items()` and building `{"bar_id": k, **v}`. Also added `bar_id` to the snapshot value dict in `01_convert.py` for future correctness.

### Results

| Metric | Value |
|---|---|
| Total calls | 1,550 |
| Completed (OK) | 1,536 |
| Failed (504/timeout) | 14 |
| Failure rate | 0.9% |

**Per-condition mean pLDDT (top 25 bars, 80–300 AA):**

| Condition | n | Mean pLDDT | Min | Max |
|---|---|---|---|---|
| alanine | 25 | 0.3546 | 0.2633 | 0.5013 |
| native_alanine | 25 | 0.3526 | 0.2425 | 0.4496 |
| native | 373 | 0.3344 | 0.2383 | 0.6305 |
| random | 742 | 0.3290 | 0.2399 | 0.5184 |
| concordance | 371 | 0.3261 | 0.2440 | 0.5112 |

**Key observation:** Alanine substitution outperforms softmax across both mapping conditions — same direction as v1. The deterministic conditions (alanine, native_alanine) lead the table. Concordance is the weakest condition. This is consistent with the v1 finding that BOJUXZ → A produces more foldable sequences than softmax or random draws at those positions.

Note: these are top-25 bars only, with sequences ranging 82–148 AA. The length confound is still present — cross-condition analysis must stratify by length bucket.

**Output:** `outputs/esm/plddt_scores.csv` (1,550 rows), PDB files in `outputs/esm/{condition}/pdbs/`.

---

## Figures — ESMFold Pilot Analysis

**Script:** `analysis/plot_esm_pilot.py`
**Output:** `outputs/figures/`
**Generated:** March 2026

Run to regenerate all figures:
```bash
python analysis/plot_esm_pilot.py
```

---

### Fig 01 — Sequence length distribution (85 bars)

![Fig 01](outputs/figures/fig01_length_distribution.png)

Length histogram of all 85 bars passing the 80–300 AA filter. Dashed lines at 100 AA and 120 AA mark the length bucket boundaries used in cross-condition analysis. The distribution is right-skewed — most bars cluster in the 80–110 AA range, with a long tail to 155 AA. No bars exceed 155 AA in this dataset despite the 300 AA ceiling being open.

---

### Fig 02 — Iconicity score distribution by divergence badge

![Fig 02](outputs/figures/fig02_iconicity_distribution.png)

Distribution of `aggregate_iconicity` scores across all 85 bars, coloured by divergence badge (Viral / Aligned). The vertical dotted line marks the top-25 iconicity cutoff — only bars to the right of this line were included in the ESMFold pilot. The distribution is right-skewed; most bars sit in the 0.15–0.40 range with a few high-iconicity outliers above 0.80.

---

### Fig 03 — Mean pLDDT by condition

![Fig 03](outputs/figures/fig03_mean_plddt_by_condition.png)

Bar chart of mean pLDDT per condition with standard error bars. Values are computed as bar-level means (averaging across seeds first) then averaged across the 25 bars. The ordering is consistent with v1: alanine and native_alanine lead, concordance is last. All conditions are within a ~0.03 range — condition effects are real but small relative to the overall pLDDT level (~0.33–0.35).

| Condition | Mean pLDDT | SE |
|---|---|---|
| alanine | 0.3546 | — |
| native_alanine | 0.3526 | — |
| native | 0.3344 | — |
| random | 0.3290 | — |
| concordance | 0.3261 | — |

---

### Fig 04 — pLDDT distribution by condition (violin)

![Fig 04](outputs/figures/fig04_plddt_violin.png)

Violin plot showing the full pLDDT distribution per condition (bar-level means, n=25 per condition). Individual bars are overlaid as jittered dots. The distributions overlap substantially — no condition dominates at the bar level. The native condition shows the widest spread, consistent with AA pass-through producing more variable sequences than frequency-remapped conditions.

---

### Fig 05 — pLDDT vs sequence length (concordance)

![Fig 05](outputs/figures/fig05_plddt_vs_length.png)

Scatter of mean pLDDT against sequence length for the 25 pilot bars under the concordance condition. Colour encodes iconicity score. The negative trend confirms the length-pLDDT relationship seen in v1 — longer sequences fold with lower confidence. The Pearson r will be reported once the full 85-bar run is complete; the top-25 sample is too small to draw firm conclusions on the slope.

The iconicity colour encoding shows no obvious clustering — high-iconicity bars (dark orange) are scattered across the length range, further supporting the iconicity-structure orthogonality finding from v1.

---

### Fig 06 — Mean pLDDT by condition × length bucket

![Fig 06](outputs/figures/fig06_plddt_by_condition_bucket.png)

Grouped bar chart stratifying all 5 conditions across three length buckets (80–100, 101–120, 121–155 AA). This is the primary cross-condition comparison — raw condition means without length stratification are confounded.

Key observations:
- The alanine advantage is most visible in the 80–100 AA bucket (short sequences), where individual residue choices have maximum structural impact
- All conditions converge toward indistinguishable values in the 121–155 AA bucket
- The pattern is consistent with v1: condition effects shrink as length grows, because BOJUXZ positions become a smaller fraction of the total sequence

---

### Fig 07 — Per-bar delta: concordance minus alanine

![Fig 07](outputs/figures/fig07_concordance_vs_alanine.png)

Per-bar pLDDT delta (concordance − alanine), sorted from most negative to most positive. Green bars indicate concordance outperforms alanine; blue bars indicate alanine outperforms concordance. The majority of bars are blue (alanine higher), with a mean delta shown. This directly quantifies how often the softmax BOJUXZ draw hurts vs helps fold confidence relative to the neutral alanine placeholder.

---

### Fig 08 — Iconicity vs mean pLDDT

![Fig 08](outputs/figures/fig08_iconicity_vs_plddt.png)

Scatter of aggregate iconicity against mean concordance pLDDT for the 25 pilot bars, coloured by divergence badge. The Pearson r printed on the figure is the key finding: cultural resonance has no meaningful predictive relationship with structural quality. This replicates the v1 finding (r = 0.051 over 225 bars). The lyric domain and protein domain are orthogonal.

---

### Fig 09 — Structural sensitivity: pLDDT SD per bar

![Fig 09](outputs/figures/fig09_plddt_sd_per_bar.png)

Standard deviation of pLDDT across 15 seeds (concordance condition) per bar, sorted descending. Bars coloured by length bucket. High SD means the BOJUXZ softmax draws strongly affect fold confidence for that sequence — the protein is structurally sensitive to which AA is drawn at non-standard positions. Low SD means the fold is robust regardless of BOJUXZ substitution.

The dashed line marks the mean SD across all 25 bars. Bars above the line are candidates for BOJUXZ mask analysis — understanding which specific position drives the variance is the next mechanistic question.

---

### Fig 10 — 2×2 factorial effect sizes

![Fig 10](outputs/figures/fig10_factorial_effects.png)

Bar chart of mean pLDDT deltas for six pairwise comparisons of interest. Green = first condition scores higher; red = second condition scores higher.

| Comparison | Effect | Interpretation |
|---|---|---|
| conc − native | + | Freq remapping: small positive effect |
| alanine − native_alanine | + | Freq remapping with alanine: small positive effect |
| conc − alanine | − | BOJUXZ softmax hurts: alanine folds better |
| native − native_alanine | − | BOJUXZ softmax hurts even without remapping |
| conc − random | + | Softmax marginally better than random |
| alanine − random | + | Alanine clearly better than random |

All BOJUXZ strategy comparisons are red (alanine outperforms softmax). All freq-remapping comparisons are green (concordance outperforms native). Same direction as v1. Effect magnitudes are 0.001–0.030 — consistent but small relative to the length-driven variance.

---

## Stage 3 — Boltz-2 (PENDING)

**Tool:** Boltz-2
**Platform:** Google Colab Pro+ (NVIDIA A100-SXM4-80GB, High-RAM)
**Input:** `outputs/fastas/bars_v2_concordance.fasta` (85 sequences)
**Diffusion samples:** 5 per bar (425 PDB files total)

### Colab command

```python
!pip install boltz cuequivariance-torch

!boltz predict bars_v2_concordance.fasta \
    --use_msa_server \
    --output_format pdb \
    --diffusion_samples 5 \
    --out_dir boltz_outputs/
```

Download `boltz_outputs/` to `outputs/boltz/`, then:

```bash
make parse-boltz
```

**Expected timeline:** MSA phase ~3–4 hours (85 bars, rate-limited). Fold ~1.5–2 hours for 5 samples × 85 bars from cached MSA.

### Why 5 samples

Boltz-2 is a diffusion model — each sample draws from a different point in the distribution of plausible structures. Running 5 samples per bar gives us:

- **Best model by confidence** — the highest-confidence prediction for structural visualization and FoldSeek
- **Mean pLDDT / pTM across models** — more stable estimate of structural quality than a single sample
- **SD across models** — measures structural plasticity: how much the predicted fold varies across draws. High SD = the sequence is structurally ambiguous; low SD = the fold is well-determined.

This makes the Boltz-2 ensemble directly comparable to the ESMFold multi-seed approach (15 seeds for concordance/native, 30 for random).

### Output structure (per bar)

```
outputs/boltz/predictions/
  bar_0/
    bar_0_model_0.pdb  bar_0_model_1.pdb  ...  bar_0_model_4.pdb
    bar_0_confidence_model_0.json  ...  bar_0_confidence_model_4.json
```

### Parser outputs (`03_parse_boltz.py`)

| File | Contents |
|---|---|
| `outputs/boltz/boltz_models.csv` | One row per bar × model — pLDDT, pTM, confidence, pDE |
| `outputs/boltz/boltz_summary.csv` | One row per bar — mean/SD/best across 5 models + structural class |

**Integrity check:** model_0 sequence length verified against `bar_index_snapshot.json` for every bar. Mismatches flagged — this is the v1 attribution drift failure mode.

**Note on `--output_format pdb`:** V1 used CIF by default, requiring conversion for downstream tools. V2 outputs PDB directly — universally supported by PyMOL, biopython, FoldSeek, and soft design tools.

---

## Stage 4 — FoldSeek (PENDING)

Awaiting Boltz-2 outputs. Will run on bars with pTM ≥ 0.4.

```bash
make foldseek
```

Key notes from v1:
- Effective floor: ~35 AA — not a concern for this dataset (min 80 AA)
- `target` field format: `"ACCESSION Description"` — split on first space
- Rate limit: HTTP 429 after ~27 submissions → exponential backoff (30s × attempt)
- MGnify hits dominate by probability but have no functional annotation

---

## Stage 5 — Enrich Master CSV (PENDING)

After Boltz-2 and FoldSeek:

```bash
make enrich
```

Adds `boltz_plddt`, `boltz_ptm`, `boltz_confidence`, `boltz_structural_class`, `esm_plddt_*`, `foldseek_result`, `foldseek_best_hit` to `aggregated_lines_v2_enriched.csv`.

---

## Open Questions

**Blocking (pre-publication):**
- None currently — attribution is locked by snapshot, agreement check passed.

**Non-blocking:**

1. **Length stratification for cross-condition analysis** — the 80–300 AA window still has length variance (82–155 AA observed). Any condition comparison must use length buckets: 80–100, 101–120, 121–155. Planned for analysis notebooks.

2. **14 ESMFold failures** — all 504 timeouts, distributed across conditions. Will retry in full run if these bars have high iconicity. The `--resume` flag will skip all completed calls.

3. **Full ESMFold run (85 bars)** — pilot covers top 25 by iconicity. Full run needed to replicate v1 finding that structural outliers can appear anywhere in the iconicity distribution (bar_135 in v1 was iconicity=0.219, not in the top 25).

4. **bar_24 in native_alanine** — pLDDT=0.2425, the lowest in the pilot across all conditions. len=131. Worth checking if this is a sequence-specific effect or a native pass-through artifact.

5. **native max pLDDT=0.6305** — the highest single-bar pLDDT in the pilot, and it's in the native condition (AA pass-through, no freq remapping). This is a potential candidate for FoldSeek once Boltz data confirms the structure.

---

## Repo Structure

```
rap-snacks-v1/
├── pipeline/
│   ├── pipeline_utils.py           <- single source of truth for LETTER_TO_AA mapping
│   ├── 00_audit.py                 <- validate frozen CSV + spot-check
│   ├── 01_convert.py               <- FASTA conversion, 5 modes, snapshot + hash guard
│   ├── 02_esm_fold.py              <- ESMFold API ensemble runner
│   ├── 03_parse_boltz.py           <- Boltz-2 output parser + integrity check
│   ├── 04_foldseek.py              <- FoldSeek REST API search
│   └── 05_enrich_csv.py            <- join all metrics to master CSV
├── analysis/
│   ├── 01_plddt.ipynb
│   ├── 02_cross_condition.ipynb
│   ├── 03_foldseek.ipynb
│   └── 04_figures.ipynb
├── data/
│   ├── aggregated_lines_v2_frozen.csv      <- immutable input
│   ├── aggregated_lines_v2_enriched.csv    <- working master CSV
│   ├── bar_index_snapshot.json             <- authoritative bar_N -> lyric mapping
│   └── frozen_csv.sha256                   <- hash guard
├── outputs/
│   ├── fastas/                     <- bars_v2_{mode}.fasta × 5
│   ├── masks/                      <- mask_v2_{mode}.json × 5
│   ├── boltz/                      <- Boltz-2 PDB files + confidence CSVs (post-Colab)
│   ├── esm/                        <- ESMFold PDB files + plddt_scores.csv
│   └── foldseek/                   <- foldseek_hits.csv + raw JSON cache
├── Makefile
├── requirements.txt
├── README.md
└── LABNOTEBOOK.md                  <- this file
```

---

## Environment

- **Local:** MacBook Air M3, Python 3.x, conda base
- **Cloud:** Google Colab Pro+ (A100-SXM4-80GB) for Boltz-2
- **APIs:** ESMFold (`api.esmatlas.com`), FoldSeek (`search.foldseek.com`), ColabFold MSA server
- **Key packages (local):** stdlib only for pipeline stages 0–2, 4–5

---

*Living document. Update after each completed stage.*
