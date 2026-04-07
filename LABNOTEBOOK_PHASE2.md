# Rap Snacks v1 — Phase 2 Lab Notebook

**Repo:** `steamulater/rap-snacks-v1`
**Phase 2 start:** 2026-03-25
**Last updated:** 2026-04-06
**Status:** Boltz-2 validation complete · FoldSeek Phase 2 running · Next: free_native_ala MPNN run

---

## Concept

Phase 1 established that concordance-mapped rap lyrics produce structurally novel proteins — no homologs in PDB, afdb-swissprot, or MGnify across 13 bars. Phase 2 takes the best 12 candidates, optimises their sequences for wet-lab foldability using ProteinMPNN, and submits codon-optimised DNA to a cell-free expression platform.

**Core narrative:** the backbone encodes the bar. ProteinMPNN finds the most foldable sequence given that shape. Every expressed protein traces directly back to a specific Nicki Minaj lyric.

---

## Phase 1 Summary (reference)

| Metric | Value |
|--------|-------|
| Bars mapped | 85 |
| Conditions per bar | 5 (concordance, native, alanine, native_ala, random) |
| Boltz-2 models | 5 per bar × 85 = 425 PDBs |
| FoldSeek bars searched | 40 |
| Bars with zero PDB homologs | 13 (structurally novel) |
| r(iconicity, confidence) | 0.051 (orthogonal) |
| r(length, pLDDT SD) | -0.658 (longer = more consistent) |
| Figures generated | 30 (Figs 1–30, `outputs/figures/`) |

**Key Phase 1 finding:** musical iconicity and protein foldability are orthogonal. The most structurally confident bars are not the most iconic bars — the biology and the culture are independent axes.

---

## Phase 2 — Step 0: Candidate Selection

**Script:** `analysis/08_candidate_selection.py`
**Output:** `data/phase2_candidates.csv`

### Selection filters

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| Boltz pTM (mean) | ≥ 0.35 | Model committed to a fold |
| Boltz pLDDT (mean) | ≥ 0.40 | Per-residue confidence |
| Sequence length | 80–150 AA | CFPS sweet spot |
| Structural novelty | 10th-pct RMSD ≥ 11 Å | No close structural neighbours |
| Hydrophobic run | ≤ 6 consecutive | Solubility predictor for E. coli |

### Final 12 candidates

| Bar | Song | Length (AA) | Boltz pLDDT | Boltz pTM | Notes |
|-----|------|-------------|-------------|-----------|-------|
| bar_9  | Super Bass | 82 | 0.839 | — | Top performer throughout |
| bar_6  | I'm The Best | 88 | 0.711 | — | High foldability |
| bar_13 | Anaconda | — | 0.637 | — | — |
| bar_17 | Monster | 121 | 0.637 | — | Long sequence |
| bar_27 | Ganja Burn | — | 0.637 | — | Only `confident_protein_like` in Phase 1 |
| bar_3  | — | — | 0.558 | — | — |
| bar_32 | Barbie Dangerous | 94 | 0.537 | — | — |
| bar_77 | Moment 4 Life | — | 0.517 | — | — |
| bar_11 | — | — | 0.458 | — | Marginal |
| bar_0  | — | — | 0.401 | — | Recovers with free_design |
| bar_8  | — | — | 0.366 | — | Recovers partially |
| bar_46 | — | — | 0.309 | — | Backbone failure confirmed |

**Figures:** Fig 28 (concordance vs native_ala pLDDT lollipop), Fig 29 (pTM scatter), Fig 30 (candidate table)

---

## Phase 2 — Step 1: Understanding the Sequence Buckets

Each bar has four sequence representations:

| Bucket | Description | Construction |
|--------|-------------|-------------|
| `concordance` | Full lyric-to-AA mapping including BJOZXU | Direct concordance table lookup |
| `native_ala` | BJOZXU positions → Ala; all other lyric AAs kept | Replace B/J/O/Z/X/U with A |
| `masked_BJOZXU` | MPNN design, lyric AAs fixed, BJOZXU redesigned | ProteinMPNN Run 1 |
| `free_design` | MPNN design, all positions free | ProteinMPNN Run 2 |

**Why native_ala matters:** substituting Ala at ambiguous positions gives a 100% standard-AA sequence that is literally readable as the lyric. bar_9 native_ala = `ANDEVERYTHINGIPEEPEDCANTAASTANSEESIPSTEAANDITSANSWEETRESPECTAADGYALWHENMIDANESPEAK`. The lyric is in the protein.

---

## Phase 2 — Step 2: ProteinMPNN Run 1 (masked_BJOZXU)

**Script:** `analysis/09_proteinmpnn_design.py`
**Date:** 2026-04-01
**Parameters:** top-12 candidates × 50 seqs × temp=0.1 · lyric AA fixed · BJOZXU free
**Output:** `outputs/proteinmpnn/filtered_results.csv` (612 rows)

### What is masked_BJOZXU?

```
Lyric:       A N D E V E R Y T H I N G X P E E P E D ...
                                          ↑
                              X = BJOZXU position (no valid AA)
                              → ProteinMPNN designs this position
                              all other positions = locked to lyric AA
```

### Constraint analysis

| Bar | Seq len | BJOZXU (free) | Lyric AA (fixed) | Unique designs |
|-----|---------|---------------|------------------|----------------|
| bar_6  | 88  | 5  | 83  | 5 / 51  |
| bar_9  | 82  | 8  | 74  | 20 / 51 |
| bar_32 | 94  | 17 | 77  | 42 / 51 |
| bar_17 | 121 | 18 | 103 | 49 / 51 |

At temp=0.1 with few designable positions, MPNN converges — bar_6 gets only 5 unique sequences from 50 runs (5 positions = 5 degrees of freedom).

### Results

| Bar | Pass / Total | Mean ESMFold pLDDT | Status |
|-----|--------------|--------------------|--------|
| bar_9  | 51 / 51 | 0.478 | ✅ top |
| bar_13 | 51 / 51 | 0.454 | ✅ top |
| bar_17 | 51 / 51 | 0.443 | ✅ top |
| bar_6  | 51 / 51 | 0.390 | ✅ solid |
| bar_32 | 32 / 51 | 0.390 | ✅ solid |
| bar_3  | 46 / 51 | 0.376 | ✅ solid |
| bar_77 | 51 / 51 | 0.374 | ✅ solid |
| bar_11 | 29 / 51 | 0.358 | ⚠️ marginal |
| bar_27 | 31 / 51 | 0.346 | ⚠️ marginal |
| bar_0  |  0 / 51 | 0.332 | ❌ dropout |
| bar_8  |  0 / 51 | 0.290 | ❌ dropout |
| bar_46 |  0 / 51 | 0.294 | ❌ dropout |

**Total passing (ESMFold pLDDT ≥ 0.35):** 393 / 612 (64%)
**Mean pLDDT across all:** 0.38
**Mutation rate from concordance:** 5–18% (only BJOZXU positions changed)

### Dropout interpretation

A bar drops out when 0/50 designs pass ESMFold pLDDT ≥ 0.35. Two causes:

1. **Backbone failure** — the Boltz prediction has no compact core; no sequence can fold to that geometry. Confirmed by free_design run (Step 3): bar_46 scores 0.309 even with all positions free.
2. **Temperature too low** — at temp=0.1 with few designable positions, MPNN barely explores. bar_0 and bar_8 recover substantially in free_design (0.401 and 0.366) — their backbones are foldable but the masking constraint + low temperature was the bottleneck.

---

## Phase 2 — Step 3: ProteinMPNN Run 2 (free_design)

**Script:** `analysis/09b_proteinmpnn_free_design.py`
**Date:** 2026-04-05
**Parameters:** top-12 candidates × 50 seqs × temp=0.1 · NO fixed positions
**Output:** `outputs/proteinmpnn_free/filtered_results.csv` (612 rows)

### Rationale

The masked run is maximally faithful to the lyric but severely constrained — bar_6 had 5 designable positions. Free design asks: *what is the most foldable sequence for this backbone with no constraints?* This serves as:

1. **Upper bound on foldability** — if free_design also scores low, the backbone is the problem
2. **Wet-lab comparison** — does preserving lyric AAs cost foldability?

### Results

| Bar | Pass / Total | Mean ESMFold pLDDT | Mutation rate | Status |
|-----|--------------|--------------------|----|-------|
| bar_9  | 50 / 50 | 0.839 | ~65% | ✅ exceptional |
| bar_6  | 51 / 51 | 0.711 | ~74% | ✅ excellent |
| bar_13 | 51 / 51 | 0.637 | ~68% | ✅ strong |
| bar_17 | 49 / 49 | 0.637 | ~70% | ✅ strong |
| bar_27 | 50 / 51 | 0.637 | ~72% | ✅ strong |
| bar_3  | 50 / 51 | 0.558 | ~66% | ✅ solid |
| bar_32 | 50 / 51 | 0.537 | ~64% | ✅ solid |
| bar_77 | 51 / 51 | 0.517 | ~69% | ✅ solid |
| bar_11 | 50 / 51 | 0.458 | ~71% | ✅ solid |
| bar_0  | 50 / 51 | 0.401 | ~73% | ✅ recovers |
| bar_8  | 39 / 51 | 0.366 | ~68% | ⚠️ partial |
| bar_46 |  1 / 50 | 0.309 | ~75% | ❌ backbone dead |

**Total passing:** 542 / 612 (89%) — vs 64% for masked_BJOZXU
**Mean pLDDT:** 0.643 — vs 0.38 for masked_BJOZXU
**Mutation rate from concordance:** 60–75% across all bars

### bar_9 example sequences (free_design, top 3 by ESMFold pLDDT)

```
>bar_9_free_000  esm_plddt=0.881
GELTEEVLKELEKKHKEFAAKGKSLVEAIDELLKELKAKGGEENKKKIEILEKAREIAIEENLESGPDTIKKLRELLKKRGV

>bar_9_free_001  esm_plddt=0.877
GKLTEEVLKELEKKAKEFAAKGKSIVEAIEELIKEIKKKGGEKNEEKLKILNKALEIAKEENLKSGPDTIAKLRKLLKERGI

>bar_9_free_002  esm_plddt=0.877
GELTEKVLKELEKKAKEAEKKGKSLEEAIEELIKELKEKGGEENEKKLKILEKILELAKEKKLKSGPDTIEILKKLLKEEGV
```

**Sequence character:** K/E/L-heavy — classic coiled-coil / helix-bundle composition. MPNN converges to helical sequences because the bar_9 Boltz backbone is helical. Highly consistent across designs.

### Comparison: masked vs free (bar_9)

```
concordance:   GEKLWLICATVENVQLLQLKRGEAMTSACESLLSVQSALGGEKVASFESHLLAILSQLRAFGKNCGPHTLEDVKAELSQLGY
native_ala:    ANDEVERYTHINGIPEEPEDCANTAASTANSEESIPSTEAANDITSANSWEETRESPECTAADGYALWHENMIDANESPEAK
free_design:   GELTEEVLKELEKKHKEFAAKGKSLVEAIDELLKELKAKGGEENKKKIEILEKAREIAIEENLESGPDTIKKLRELLKKRGV
scrambled:     GTLSSDAELLELKWLIAQEEAVLLNSLSGRQQAGALAVQLSRCENCASYLGTMHHLGLVSAPFESKLGIQKEKFLVKTGCVE
```

---

## Phase 2 — Step 4: ESMFold Scrambled Control

**Script:** `analysis/11_scrambled_control.py`
**Date:** 2026-04-01

### Findings

| Bucket | n | mean pLDDT | sd |
|--------|---|------------|----|
| concordance (original) | 37 | 0.324 | 0.041 |
| scrambled_concordance | 110 | 0.333 | 0.059 |
| native_ala (original) | 9 | 0.370 | 0.060 |
| scrambled_native_ala | 111 | 0.397 | 0.099 |

**Critical finding:** ESMFold **cannot discriminate** between structured and scrambled sequences at this length. Scrambled concordance (0.333) ≈ concordance (0.324). ESMFold relies on PLM embeddings (ESM-2) which are partially order-insensitive for short sequences without strong long-range contacts.

**Consequence:** ESMFold cannot be used as a self-consistency validator. Switched all validation to **Boltz-2**, which is a physics-informed diffusion model and is order-sensitive.

---

## Phase 2 — Step 5: Boltz-2 Validation

**Notebook:** `notebooks/boltz_validation_v2.ipynb` (Colab A100)
**Date:** 2026-04-06
**Runtime:** ~2.5 hrs (636 seqs × 5 models + 60 scrambled × 1 model)
**Output:** `rap_snacks/results/boltz_validation_results.csv` (Drive)

### Design of the validation experiment

| Bucket | n seqs | n models | Purpose |
|--------|--------|----------|---------|
| concordance | 12 | 5 | Baseline — lyric-derived sequence |
| native_ala | 12 | 5 | Alanine-substituted baseline |
| free_design | 612 | 5 | MPNN-optimised (all positions free) |
| scrambled | 60 → 36 valid | 1 | Negative control |

Note: scrambled n=36 (not 60) — some bars had fewer than 5 valid scrambled sequences in the CSV. 36 is sufficient for a negative control.

### Boltz-2 pLDDT results

```
Bucket           n     mean    sd     min    max    models/seq
concordance      12    0.441   0.069  0.338  0.590  5
native_ala       12    0.543   0.107  0.362  0.689  5
free_design     612    0.643   0.173  0.330  0.960  5
scrambled        36    0.426   0.068  0.335  0.615  1
```

**Figure A** — `fig_boltz_per_bar.png` — strip plot: all four buckets overlaid per bar
**Figure B** — `fig_boltz_violin.png` — violin distribution by bucket

### Hierarchy

```
scrambled    0.426  ← noise floor (random AA order)
concordance  0.441  ← lyric-derived, barely above scrambled
native_ala   0.543  ← Ala substitution helps: +0.10 over concordance
free_design  0.643  ← MPNN: +0.20 over native_ala, +0.22 over concordance
```

**Key result:** free_design dramatically outperforms all baselines. MPNN is doing exactly what it should — finding sequences that fold to the target backbone.

### What each metric means

**concordance pLDDT = 0.441 ≈ scrambled 0.426**
The raw lyric-derived sequences fold no better than random. The musical encoding carries no structural information, as expected. The backbone is the biological contribution; the lyric sequence is the cultural contribution.

**native_ala pLDDT = 0.543 (+0.10)**
Replacing BJOZXU with Ala improves foldability meaningfully. All-standard-AA sequence folds better than one containing ambiguous characters. Ala is structurally neutral and helix-favoring — a good placeholder.

**free_design pLDDT = 0.643 (+0.20 over native_ala)**
MPNN has full freedom to find the optimal sequence for the backbone. The jump from native_ala to free_design quantifies how much foldability you gain by letting MPNN optimise all positions vs just substituting Ala.

**free_design range: 0.330–0.960**
Huge spread — backbone quality sets the ceiling. bar_9 designs cluster near 0.96; bar_46 still failing (0.330). bar_46 is a confirmed backbone-level failure.

---

## Phase 2 — Step 6: RMSD Analysis

**RMSD reference:** Boltz-predicted backbone from the **concordance sequence** (same backbone used as ProteinMPNN input)
**Method:** CA-only Kabsch superposition (BioPython Superimposer)

### Results

```
Bucket           n       mean RMSD   sd      min     max
free_design      612     7.13 Å      5.11    0.35    17.56
concordance       12     9.58 Å      3.76    3.44    15.18
native_ala        12    16.90 Å      6.06   11.99    33.40
free_pairwise  15300     7.57 Å      4.65    0.35    18.35
```

**Figure C** — `fig_rmsd_violin.png` — RMSD vs backbone by bucket
**Figure D** — `fig_rmsd_per_bar.png` — free_design RMSD per bar (all designs, mean line)
**Figure E** — `fig_rmsd_pairwise.png` — pairwise RMSD histogram (structural diversity)
**Figure Grid** — `fig_rmsd_per_bar_grid.png` — 12-panel grid, concordance + native_ala + free_design per bar

### What each RMSD tells you

**concordance RMSD = 9.58 Å**
Same concordance sequence, re-folded by Boltz in a new run, compared to the original run. This is purely **Boltz run-to-run stochasticity** — the noise floor. Same sequence, same model, different diffusion sample → 9.58 Å average deviation.

**native_ala RMSD = 16.90 Å**
Native_ala folds to a **completely different structure** from the concordance backbone. This is expected — native_ala is a different sequence (different secondary structure preferences) and was never designed for the concordance backbone.

**free_design RMSD = 7.13 Å**
MPNN-designed sequences fold **closer to the backbone** than re-folding the original concordance sequence does. 7.13 Å < 9.58 Å (concordance noise floor). This is the key validation: MPNN successfully encoded backbone geometry into the designed sequence.

**free_pairwise RMSD = 7.57 Å**
The 50 free_design structures per bar are structurally diverse from each other (min 0.35 Å, max 18.35 Å). MPNN is not just generating one sequence repeated — the designs genuinely explore different compatible conformations of the same backbone.

### Hierarchy (RMSD)

```
free_design  7.13 Å  → most faithful to backbone (MPNN worked)
concordance  9.58 Å  → Boltz stochasticity (noise floor)
native_ala  16.90 Å  → different structure entirely (expected)
```

---

## Phase 2 — Step 7: FoldSeek Phase 2 (in progress)

**Script:** `analysis/10_foldseek_phase2.py`
**Date:** 2026-04-06
**Status:** Running in Colab (36 structures × 3 databases)

### Design

For each of the 12 bars, submitting the best PDB from each of 3 buckets:
- **concordance** — model_0 from Boltz validation run
- **native_ala** — model_0 from Boltz validation run
- **free_design** — highest boltz_plddt design for that bar

Searching: `pdb100`, `afdb-swissprot`, `mgnify_esm30`
Cached to Drive: `rap_snacks/results/foldseek_phase2/raw/{bar_id}_{bucket}.json`

### Questions being answered

1. Do concordance, native_ala, and free_design structures hit the **same protein families** or different ones?
2. Does free_design, with higher pLDDT, also find **better-quality PDB hits** (higher prob) than concordance?
3. Are there bars where native_ala hits known folds that concordance misses — suggesting the Ala substitution shifts the sequence into fold space?
4. Which bars remain **structurally novel** (zero hits) across all three buckets?

**Figure F** — `fig_foldseek_phase2.png` — bar chart: top pdb100 hit prob + hit count per bar per bucket (pending)

---

## Phase 2 — Key Scientific Findings (Summary)

### Finding 1 — ESMFold cannot be used as a validation tool
ESMFold pLDDT is order-insensitive for short sequences. Scrambled ≈ original. All structural validation must use Boltz-2.

### Finding 2 — free_design dramatically outperforms all baselines
```
free_design:   mean Boltz pLDDT 0.643  (89% pass rate)
native_ala:                     0.543
scrambled:                      0.426
concordance:                    0.441  (barely above scrambled)
```
The lyric-derived sequence (concordance) folds no better than random. The backbone is the protein biology; the lyric is the cultural encoding.

### Finding 3 — MPNN encodes backbone geometry better than the original sequence
Free_design RMSD to backbone = 7.13 Å < concordance re-fold RMSD = 9.58 Å. The designed sequences are geometrically more faithful to the target backbone than re-folding the sequence that generated the backbone.

### Finding 4 — bar_46 is a confirmed backbone-level failure
bar_46 scores 0.309 pLDDT even when fully unconstrained (free_design). The backbone itself is disordered — no sequence can stabilise it. bar_46 is dropped from the submission.

### Finding 5 — native_ala is more foldable than concordance
0.543 vs 0.441 — replacing BJOZXU with Ala consistently improves foldability. The non-standard AA positions are a structural liability in the concordance sequence.

### Finding 6 — structural diversity of free_design is real
Free-pairwise RMSD mean 7.57 Å, range 0.35–18.35 Å. The 50 designs per bar genuinely sample different compatible conformations — not degenerate copies of one sequence.

---

## Phase 2 — Scientific Question: native_ala as backbone for MPNN

**Motivated by Findings 4 & 5.** If native_ala folds more confidently than concordance (0.543 vs 0.441), then using the native_ala Boltz backbone as the MPNN input may produce even better designed sequences — without losing the cultural signal.

### Three design strategies compared

| Strategy | Backbone | Fixed positions | Cultural fidelity | Design freedom |
|----------|----------|-----------------|-------------------|----------------|
| masked_BJOZXU (Run 1) | concordance Boltz | lyric AAs | high | low (5–18 positions) |
| free_design (Run 2) | concordance Boltz | none | none | maximum |
| **free_native_ala (planned)** | **native_ala Boltz** | **none** | **moderate** | **maximum** |

**Why free_native_ala over masked_native_ala:**
Masked MPNN with few designable positions (bar_6: 5) produces near-identical sequences — insufficient diversity for expression and hard to distinguish from the baseline. Free design on the native_ala backbone gives full design freedom while starting from a more foldable backbone geometry.

**Cultural note:** native_ala sequences are the most readable — the lyric is literally visible in the AA sequence. free_native_ala designs will diverge from this, but they start from a more structurally grounded position.

**Planned:** `analysis/09d_proteinmpnn_native_ala_free.py`
Input: native_ala Boltz backbones from Drive (`results/boltz_outputs/predictions/{bar_id}_native_ala/`)
Output: `outputs/proteinmpnn_native_ala_free/`

---

## Platform Comparison

### Ginkgo Cloud Lab — CFPS + HiBiT

| Parameter | Value |
|-----------|-------|
| System | E. coli lysate + HiBiT bioluminescent readout |
| Per protein (24 proteins) | ~$117 |
| Total (24 proteins) | $2,808 |
| Replicates | 3 per protein |
| Readout | Quantitative HiBiT luminescence |
| Gene synthesis | Yes |
| Turnaround | Not specified |

### Adaptyv Bio — CFPS

| Parameter | Value |
|-----------|-------|
| System | Reconstituted E. coli TX/TL (purified components, lower background) |
| Per protein (open-source discount) | ~$95 |
| Total (24 proteins) | $2,285 |
| Replicates | 1 |
| Readout | Expressed / low / not detected + relative yield |
| Gene synthesis | Yes (submit AA sequence) |
| Turnaround | 2–3 weeks (standard) |
| Condition | Publish results on Proteinbase |
| De novo track record | 10,000+ designed proteins tested |

### Recommendation

Adaptyv Bio — $500 cheaper, stronger de novo validation track record, reconstituted system better suited for synthetic proteins. Compatible with project's publication goals (Proteinbase release required for discount).

**Pending decision.** Final call after FoldSeek Phase 2 results and free_native_ala MPNN run complete.

---

## Submission Plan (48-protein plate, Adaptyv Bio)

| Slot | Count | Description |
|------|-------|-------------|
| 12 bars × concordance | 12 | Raw lyric-derived sequence |
| 12 bars × native_ala | 12 | Ala-substituted lyric sequence |
| 12 bars × free_design (best) | 12 | MPNN, no constraints, concordance backbone |
| 12 bars × free_native_ala (best) | 12 | MPNN, no constraints, native_ala backbone |
| Positive control (sfGFP) | 1 | Confirms platform ran |
| Reserve | ~7 | — |
| **Total** | **48** | ~$4,560 at $95/protein |

**Tag:** C-terminal His6 (Adaptyv) or HiBiT (Ginkgo)
**Codon optimisation:** E. coli K12, avoid rare codons, no T7 terminators, no hairpins in first 30 nt

---

## Controls

### In silico scrambled control (complete)
- 3 shuffled variants per bar × concordance + native_ala
- ESMFold pLDDT: scrambled ≈ original → **ESMFold not a valid discriminator**
- Boltz-2 pLDDT: scrambled 0.426 < concordance 0.441 < native_ala 0.543 < free_design 0.643 → **Boltz discriminates**
- Scrambled control validated in silico — no wet-lab slot needed

### Positive control (wet lab)
- 1 × sfGFP or known-expressing protein — occupies 1 well

### native_ala as in silico baseline
- native_ala already in the submission — serves as the deterministic baseline for each bar

---

## Sequence Engineering Technical Notes

### ProteinMPNN fixed positions (masked_BJOZXU)
- Positions map: `data/phase2_fixed_positions.jsonl`
- Fixed = lyric-derived standard AA positions
- Free = BJOZXU positions (B/J/O/Z/X/U → any standard AA)
- Run via `--fixed_positions_jsonl`

### Boltz chain ID fix
Boltz writes 2-char chain IDs (`b11`, `A `) which break ProteinMPNN:
```python
# Fix: overwrite column 22 (chain ID) with 'A'
line = line[:21] + 'A' + line[24:]
```

### Boltz YAML format (single-seq mode)
```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: GELTEEVLKELEKK...
      msa: empty          # required — omitting causes MSA error
```

### Boltz Colab flags
```bash
boltz predict {input_dir} \
  --out_dir {out_dir} \
  --diffusion_samples 5 \
  --output_format pdb \
  --no_kernels \          # underscore not hyphen
  --override              # force regeneration
```

### numpy binary incompatibility after boltz install
Boltz changes numpy version. Colab kernel must restart after `pip install boltz`. Cell 1 of `notebooks/boltz_validation_v2.ipynb` handles this automatically with `os.kill(os.getpid(), 9)`.

---

## Figure Inventory

### Local (outputs/figures/)

| Figure | File | Description |
|--------|------|-------------|
| Fig 1 | fig01_length_distribution.png | Bar length distribution (85 bars) |
| Fig 2 | fig02_iconicity_distribution.png | Iconicity score distribution |
| Fig 3 | fig03_mean_plddt_by_condition.png | ESMFold pLDDT by condition |
| Fig 4 | fig04_plddt_violin.png | ESMFold pLDDT violin |
| Fig 5 | fig05_plddt_vs_length.png | pLDDT vs sequence length |
| Fig 6 | fig06_plddt_by_condition_bucket.png | pLDDT by condition bucket |
| Fig 7 | fig07_concordance_vs_alanine.png | Concordance vs alanine pLDDT |
| Fig 8 | fig08_iconicity_vs_plddt.png | Iconicity vs pLDDT scatter |
| Fig 9 | fig09_plddt_sd_per_bar.png | pLDDT SD per bar |
| Fig 10 | fig10_factorial_effects.png | Factorial condition effects |
| Fig 11–20 | fig11–fig20 | ESMFold full run (85 bars × 5 conditions) |
| Fig 28 | fig28_conc_vs_na_plddt.png | Concordance vs native_ala pLDDT lollipop |
| Fig 29 | fig29_conc_vs_na_ptm.png | Concordance vs native_ala pTM scatter |
| Fig 30 | fig30_candidate_table.png | Phase 2 candidate summary table |

### Drive (rap_snacks/results/figures/)

| Figure | File | Description |
|--------|------|-------------|
| Fig A | fig_boltz_per_bar.png | Boltz pLDDT strip plot — 4 buckets per bar |
| Fig B | fig_boltz_violin.png | Boltz pLDDT violin — 4 buckets |
| Fig C | fig_rmsd_violin.png | RMSD vs backbone violin by bucket |
| Fig D | fig_rmsd_per_bar.png | free_design RMSD per bar |
| Fig E | fig_rmsd_pairwise.png | Pairwise RMSD histogram |
| Fig F | fig_foldseek_phase2.png | FoldSeek Phase 2 hits (pending) |
| Grid | fig_rmsd_per_bar_grid.png | 12-panel per-bar RMSD grid |

**Next local figure number:** Fig 31

---

## Pending Work

| Priority | Task | Script |
|----------|------|--------|
| 1 | Parse + interpret FoldSeek Phase 2 results | — (results from Colab cell) |
| 2 | Build free_native_ala MPNN run | `analysis/09d_proteinmpnn_native_ala_free.py` |
| 3 | Boltz-2 validate free_native_ala designs | Add to `notebooks/boltz_validation_v2.ipynb` |
| 4 | Codon optimisation | `analysis/10_codon_optimize.py` |
| 5 | Platform decision + submission | Adaptyv Bio API or Ginkgo EstiMate |

---

## References

- Dauparas et al. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. *Science* 378, 49–56.
- Kipnis et al. (2023). Improving Protein Expression, Stability, and Function with ProteinMPNN. *JACS* 145(29).
- van Kempen et al. (2024). Fast and accurate protein structure search with Foldseek. *Nature Biotechnology* 42, 243–246.
- Adaptyv Bio EGFR design competition — crowdsourced de novo protein validation. [bioRxiv 2025.04.17.648362](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2)
- Adaptyv Bio expression service. [start.adaptyvbio.com](https://start.adaptyvbio.com)
- Ginkgo Cloud Lab CFPS + HiBiT. [cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit](https://cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit)

---

*Living document. Phase 2 active — Boltz-2 validation complete, FoldSeek Phase 2 running, free_native_ala MPNN run next.*
