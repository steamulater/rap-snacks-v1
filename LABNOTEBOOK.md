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

## Stage 3 — Boltz-2 (PENDING)

**Tool:** Boltz-2
**Platform:** Google Colab Pro+ (NVIDIA A100-SXM4-80GB, High-RAM)
**Input:** `outputs/fastas/bars_v2_concordance.fasta` (85 sequences)
**Command:**

```bash
!pip install boltz cuequivariance-torch
!boltz predict bars_v2_concordance.fasta \
    --use_msa_server \
    --output_format pdb \
    --diffusion_samples 1 \
    --out_dir boltz_outputs/
```

Download `boltz_outputs/` to `outputs/boltz/`, then:

```bash
make parse-boltz
```

**Expected timeline:** MSA phase ~3–4 hours (85 bars). Fold ~30–45 min from cached MSA.

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
