# Rap Snacks v1 — Phase 2 Lab Notebook

**Repo:** `steamulater/rap-snacks-v1`
**Phase 2 start:** 2026-03-25
**Last updated:** 2026-04-07
**Status:** native_ala_free MPNN + Boltz v3 complete · BioReason-Pro pipeline ready · Next: codon optimisation + platform decision

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

---

### Figure 28 — Concordance vs Native-Alanine pLDDT (Boltz-2)

![Figure 28](outputs/figures/fig28_conc_vs_na_plddt.png)

**Figure 28 | Boltz-2 pLDDT comparison between concordance and native-alanine conditions across Phase 2 candidate bars.** Each bar represents one of the 12 candidate sequences scored under two encoding conditions: *concordance* (frequency-rank remapping of standard letters + softmax-peaked probabilistic draw for BJOZXU characters) and *native_ala* (lyric characters passed through literally as amino acids, BJOZXU → Ala — the more raw encoding). Bars are sorted by the concordance pLDDT score. Native-alanine outperforms concordance in the majority of cases, with a mean improvement of +0.10 pLDDT units. The concordance condition is nonetheless retained in the submission panel as it represents the canonical probabilistic mapping.

---

### Figure 29 — Concordance vs Native-Alanine pTM Scatter

![Figure 29](outputs/figures/fig29_conc_vs_na_ptm.png)

**Figure 29 | Predicted TM-score (pTM) scatter comparing concordance and native-alanine conditions.** Each point represents one bar, plotted with concordance pTM on the x-axis and native-alanine pTM on the y-axis. Points above the diagonal indicate bars where alanine substitution improves global fold confidence. The majority of bars lie above the diagonal, consistent with the pLDDT analysis in Figure 28. bar_27 (Ganja Burn) is a notable exception — one of the rare cases where concordance equals or exceeds native-alanine in pTM, consistent with its Phase 1 classification as the only `confident_protein_like` bar in the concordance condition. This scatter establishes the per-bar baseline before any sequence design intervention.

---

### Figure 30 — Phase 2 Candidate Summary Table

![Figure 30](outputs/figures/fig30_candidate_table.png)

**Figure 30 | Summary table of the 12 Phase 2 candidate bars ranked by composite score.** Columns include bar identifier, source song, sequence length, Boltz-2 pLDDT and pTM under concordance and native-alanine conditions, structural novelty (FoldSeek hit status), and hydrophobic run length (solubility predictor). Candidates were selected by applying the five filters described above; the 12 shown are the highest-ranked bars satisfying all thresholds. bar_46 is included because it passed initial filters but was subsequently confirmed as a backbone-level failure in the free_design run (see Step 3). bar_9, bar_6, and bar_27 represent the top tier based on combined structural confidence and cultural salience.

---

## Phase 2 — Step 1: Understanding the Sequence Buckets

Each bar has four sequence representations, forming a hierarchy from most culturally faithful to most computationally optimised:

| Bucket | Description | Construction |
|--------|-------------|-------------|
| `concordance` | Full lyric-to-AA mapping including BJOZXU | Direct concordance table lookup |
| `native_ala` | BJOZXU → Ala; all other lyric AAs kept | Replace B/J/O/Z/X/U with A |
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
**Mean pLDDT:** 0.38 · **Mutation rate from concordance:** 5–18%

### Dropout interpretation

A bar drops out when 0/50 designs pass ESMFold pLDDT ≥ 0.35. Two causes:

1. **Backbone failure** — the Boltz prediction has no compact core; no sequence can fold to that geometry. Confirmed by free_design run (Step 3): bar_46 scores 0.309 even with all positions free.
2. **Temperature too low** — at temp=0.1 with few designable positions, MPNN barely explores. bar_0 and bar_8 recover substantially in free_design (0.401 and 0.366).

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
**Mean pLDDT:** 0.643 — vs 0.38 for masked_BJOZXU · **Mutation rate:** 60–75%

### bar_9 sequence comparison across all buckets

```
concordance:   GEKLWLICATVENVQLLQLKRGEAMTSACESLLSVQSALGGEKVASFESHLLAILSQLRAFGKNCGPHTLEDVKAELSQLGY
native_ala:    ANDEVERYTHINGIPEEPEDCANTAASTANSEESIPSTEAANDITSANSWEETRESPECTAADGYALWHENMIDANESPEAK
free_design:   GELTEEVLKELEKKHKEFAAKGKSLVEAIDELLKELKAKGGEENKKKIEILEKAREIAIEENLESGPDTIKKLRELLKKRGV
scrambled:     GTLSSDAELLELKWLIAQEEAVLLNSLSGRQQAGALAVQLSRCENCASYLGTMHHLGLVSAPFESKLGIQKEKFLVKTGCVE
```

The free_design sequence is K/E/L-heavy — classic coiled-coil / helix-bundle composition. MPNN converges on helical sequences because the bar_9 Boltz backbone is helical. The lyric in native_ala (`ANDEVERYTHINGIPEEPED...`) is literally readable English embedded in a valid protein sequence.

---

## Phase 2 — Step 4: ESMFold Scrambled Control

**Script:** `analysis/11_scrambled_control.py`
**Date:** 2026-04-01

| Bucket | n | mean pLDDT | sd |
|--------|---|------------|----|
| concordance (original) | 37 | 0.324 | 0.041 |
| scrambled_concordance | 110 | 0.333 | 0.059 |
| native_ala (original) | 9 | 0.370 | 0.060 |
| scrambled_native_ala | 111 | 0.397 | 0.099 |

**Critical finding:** ESMFold **cannot discriminate** between structured and scrambled sequences at this length. Scrambled concordance (0.333) ≈ concordance (0.324). ESMFold relies on PLM embeddings (ESM-2) which are partially order-insensitive for short sequences without strong long-range contacts. All downstream validation switched to **Boltz-2**.

---

## Phase 2 — Step 5: Boltz-2 Validation

**Notebook:** `notebooks/boltz_validation_v2.ipynb` (Colab A100)
**Date:** 2026-04-06
**Runtime:** ~2.5 hrs (636 seqs × 5 models + 60 scrambled × 1 model)
**Data:** `outputs/boltz_validation/boltz_validation_results.csv`

### Experimental design

| Bucket | n seqs | n models | Purpose |
|--------|--------|----------|---------|
| concordance | 12 | 5 | Baseline — lyric-derived sequence |
| native_ala | 12 | 5 | Alanine-substituted baseline |
| free_design | 612 | 5 | MPNN-optimised (all positions free) |
| scrambled | 36 | 1 | Negative control |

### Boltz-2 pLDDT results

```
Bucket           n     mean    sd     min    max
concordance      12    0.441   0.069  0.338  0.590
native_ala       12    0.543   0.107  0.362  0.689
free_design     612    0.643   0.173  0.330  0.960
scrambled        36    0.426   0.068  0.335  0.615
```

---

### Figure 31 — Boltz-2 pLDDT Strip Plot Per Bar

![Figure 31](outputs/figures/fig31_boltz_plddt_per_bar.png)

**Figure 31 | Boltz-2 per-residue confidence (pLDDT) for all four sequence buckets across the 12 candidate bars.** Each vertical cluster of points represents one bar. Within each cluster, the four buckets — scrambled (dark grey), concordance (cyan), native_ala (orange), and free_design (purple) — are offset horizontally for clarity. Individual points represent single sequence predictions; horizontal tick marks denote the bucket mean. The dashed line at pLDDT = 0.5 marks the conventional threshold separating low-confidence (disordered) from moderate-confidence (structured) predictions. The progressive improvement from scrambled → concordance → native_ala → free_design is consistent across bars, confirming that the observed hierarchy is not driven by one or two outlier sequences. bar_9 and bar_6 are exceptional — their free_design clouds extend well above 0.8, approaching the confidence levels typically seen for computationally designed proteins. bar_46 is the only bar where free_design fails to rise above scrambled, confirming backbone-level structural failure.

---

### Figure 32 — Boltz-2 pLDDT Violin Plot by Bucket

![Figure 32](outputs/figures/fig32_boltz_plddt_violin.png)

**Figure 32 | Distribution of Boltz-2 pLDDT scores across all sequences in each bucket, shown as kernel density estimates (violin plots).** White horizontal lines mark the median; violin width encodes local density. The scrambled and concordance distributions are nearly identical — both centred around 0.42–0.44 — confirming that the raw lyric-derived sequences fold no better than randomly shuffled amino acid sequences of the same composition. The native_ala distribution shifts upward by approximately 0.10 pLDDT units, demonstrating the structural benefit of replacing non-standard BJOZXU characters with alanine. The free_design distribution is both higher (mean 0.643) and substantially wider (sd 0.173), reflecting the heterogeneity of the 612 MPNN-designed sequences across 12 bars with different backbone geometries. The rightward tail of the free_design violin, extending to pLDDT ≈ 0.96, represents the top-performing bar_9 designs and indicates that MPNN can recover near-native confidence levels when given an appropriate backbone.

---

## Phase 2 — Step 6: RMSD Analysis

**Method:** CA-only Kabsch superposition (BioPython Superimposer)
**Reference:** Boltz-predicted concordance backbone (same structure used as ProteinMPNN input)
**Data:** `outputs/boltz_validation/boltz_rmsd.csv`

### Results

```
Bucket           n       mean RMSD   sd      min     max
free_design      612     7.13 Å      5.11    0.35    17.56
concordance       12     9.58 Å      3.76    3.44    15.18
native_ala        12    16.90 Å      6.06   11.99    33.40
free_pairwise  15300     7.57 Å      4.65    0.35    18.35
```

### What each RMSD tells you

**concordance RMSD = 9.58 Å — the Boltz noise floor.**
The concordance sequence re-folded in a new Boltz run deviates 9.58 Å from the original run of the same sequence. This is not sequence-dependent signal — it is purely stochastic variation between diffusion samples of the same input. All subsequent RMSD comparisons must be interpreted relative to this 9.58 Å baseline.

**native_ala RMSD = 16.90 Å — a structurally unrelated fold.**
Native_ala is a chemically distinct sequence from concordance (BJOZXU→Ala at 5–18 positions per bar). Its Boltz-predicted structure lies 16.90 Å from the concordance backbone on average, confirming that native_ala does not fold to the same shape. This is expected — native_ala was not designed for the concordance backbone. The high RMSD does not indicate failure; it simply means native_ala and concordance occupy different regions of fold space, which is itself an interesting finding.

**free_design RMSD = 7.13 Å — backbone recovery below the noise floor.**
MPNN-designed sequences, when re-folded by Boltz, deviate only 7.13 Å from the concordance backbone — less than the 9.58 Å deviation seen when re-folding the original concordance sequence itself. This is the key structural validation result: MPNN has encoded the backbone geometry into the designed sequences more faithfully than the original sequence that generated the backbone. The minimum free_design RMSD of 0.35 Å indicates near-perfect backbone recovery for at least one design.

**free_pairwise RMSD = 7.57 Å — genuine structural diversity.**
The 15,300 pairwise RMSD values between free_design structures within the same bar average 7.57 Å, with a range of 0.35–18.35 Å. This confirms that the 50 MPNN designs per bar are not degenerate copies of a single solution — they genuinely sample a diverse set of backbone-compatible conformations.

---

### Figure 33 — RMSD vs Backbone Violin by Bucket

![Figure 33](outputs/figures/fig33_rmsd_violin.png)

**Figure 33 | Distribution of Cα RMSD values (Å) relative to the concordance backbone for each sequence bucket.** RMSD was computed between the Boltz-2 model_0 prediction for each sequence and the original concordance Boltz backbone used as the ProteinMPNN input. Violin width encodes density; white lines mark medians. The ordering along the RMSD axis — free_design (7.13 Å) < concordance (9.58 Å) < native_ala (16.90 Å) — reveals a counterintuitive result: the MPNN-designed sequences are geometrically closer to the target backbone than re-folding the original concordance sequence. The concordance RMSD (9.58 Å) represents the irreducible stochasticity of the Boltz-2 diffusion process for these sequences; free_design beats this baseline, confirming successful backbone encoding by MPNN. The native_ala bucket has the widest distribution (sd 6.06 Å) because different bars have structurally very different native_ala folds depending on the composition of their lyric amino acids.

---

### Figure 34 — free_design RMSD Per Bar

![Figure 34](outputs/figures/fig34_rmsd_per_bar.png)

**Figure 34 | Cα RMSD of each free_design Boltz-2 prediction relative to the concordance backbone, shown per bar.** Each point represents one MPNN-designed sequence folded by Boltz-2; the white horizontal bar marks the per-bar mean. Substantial variation exists across bars, reflecting differences in backbone quality and fold type. Bars with compact helical backbones (e.g., bar_9, bar_6) show lower RMSD values and tighter point clouds, indicating that MPNN designs for these bars consistently recover the target geometry. Bars with less defined backbones show broader distributions and higher means, consistent with Boltz-2's stochastic sampling behaviour in regions of low structural confidence. The per-bar spread within each cluster (vertical extent of points) reflects the structural diversity of the 50 MPNN designs — diverse sequences can fold to structurally distinct but backbone-compatible conformations.

---

### Figure 35 — Pairwise RMSD Histogram (free_design)

![Figure 35](outputs/figures/fig35_rmsd_pairwise.png)

**Figure 35 | Distribution of pairwise Cα RMSD values between all free_design Boltz-2 structures within the same bar (n = 15,300 pairs).** For each bar, all 50 MPNN-designed sequences were folded by Boltz-2 and pairwise RMSD was computed between every pair of structures within that bar. The distribution is broad (range 0.35–18.35 Å, mean 7.57 Å, dashed white line), demonstrating that MPNN explores genuinely diverse structural solutions for the same backbone template. The left tail (0.35–2 Å) represents near-identical structures — pairs of sequences that converge on the same low-energy conformation. The right tail (>14 Å) represents structurally distinct solutions that satisfy the backbone constraint through different conformational routes. This breadth is favourable for a wet-lab submission: the 50 designs per bar are not redundant, and selecting diverse representatives for expression increases the probability of identifying at least one that expresses well.

---

### Figure 36 — Per-Bar RMSD Grid (12 bars × 3 buckets)

![Figure 36](outputs/figures/fig36_rmsd_per_bar_grid.png)

**Figure 36 | Per-bar structural comparison grid showing Cα RMSD vs concordance backbone for three sequence buckets across all 12 candidate bars.** Each panel represents one bar. Within each panel, the x-axis shows the three buckets (concordance, native_ala, free_design) and the y-axis shows RMSD in Å. Points are jittered for free_design (multiple sequences per bar) and placed precisely for concordance and native_ala (single sequences). Coloured mean lines and numerical labels facilitate direct comparison. This panel view reveals bar-level heterogeneity not visible in the aggregate plots. Notably, the concordance–free_design RMSD gap varies substantially across bars — bars with more helical backbones (bar_9, bar_6) show tighter free_design clusters and a larger separation from native_ala, while bars with mixed or uncertain backbone geometry show more diffuse patterns. The consistently high native_ala RMSD across all 12 bars confirms that alanine substitution at BJOZXU positions reliably redirects folding away from the concordance backbone — native_ala and concordance occupy orthogonal structural subspaces.

---

## Phase 2 — Key Scientific Findings (Summary)

### Finding 1 — ESMFold cannot be used as a validation tool
ESMFold pLDDT is order-insensitive for short sequences. Scrambled ≈ original. All structural validation must use Boltz-2 or equivalent physics-informed models.

### Finding 2 — concordance folds no better than scrambled
Boltz-2 pLDDT: concordance 0.441 ≈ scrambled 0.426. The raw lyric-derived sequence carries no inherent structural information — the biological signal is in the backbone geometry extracted by Boltz-2, not in the lyric amino acid sequence itself.

### Finding 3 — native_ala is meaningfully more foldable (+0.10)
Replacing BJOZXU with alanine improves Boltz-2 pLDDT from 0.441 to 0.543 across 12 bars. BJOZXU characters are a structural liability in the concordance encoding.

### Finding 4 — free_design dramatically outperforms all baselines
free_design mean Boltz-2 pLDDT = 0.643, pass rate 89% (542/612). MPNN, given full design freedom, finds sequences that fold confidently to the target backbone. The 0.20 gap over native_ala quantifies the cost of imposing lyric constraints on sequence design.

### Finding 5 — MPNN encodes backbone geometry better than the original sequence
free_design RMSD to backbone = 7.13 Å < concordance re-fold RMSD = 9.58 Å. Designed sequences are geometrically more faithful to the target backbone than re-folding the sequence that generated it. This is the core structural validation of the MPNN design strategy.

### Finding 6 — bar_46 is a confirmed backbone-level failure
bar_46 scores 0.309 pLDDT even when fully unconstrained. The backbone itself is disordered — no sequence design strategy can salvage it. bar_46 is dropped from the wet-lab submission.

### Finding 7 — free_design structural diversity is real and favourable
Pairwise RMSD between free_design structures: mean 7.57 Å, range 0.35–18.35 Å across 15,300 pairs. The designs are not degenerate. This diversity increases the probability of finding an expressible sequence in the wet-lab screen.

---

## Phase 2 — Step 7: FoldSeek Phase 2

**Script:** `analysis/10_foldseek_phase2.py`
**Date:** 2026-04-06
**Output:** `outputs/foldseek_phase2/foldseek_phase2_hits.csv`, `foldseek_phase2_summary.csv`

### What FoldSeek is doing

FoldSeek is a structural search tool — you give it a PDB file and it searches protein structure databases to find proteins with similar 3D shapes, regardless of sequence similarity. For each of the 12 bars we submit 3 structures (concordance, native_ala, best free_design) and search against three databases:

- **pdb100** — the entire Protein Data Bank (experimentally solved structures)
- **afdb-swissprot** — AlphaFold predictions for all SwissProt proteins (well-characterised proteins)
- **mgnify_esm30** — metagenomic proteins (environmental sequences, no annotation)

### The question we are asking

Does each structure resemble anything biology has already made? And critically — do the three versions of the same bar hit the same protein families, or does each version occupy a different structural neighbourhood?

### Why each bucket matters

**concordance** — if it hits known proteins, the lyric backbone accidentally encodes a real fold. If 0 hits, it is genuinely structurally novel. Phase 1 found 13 bars with 0 hits across all databases — this is the headline novelty claim of the project.

**native_ala** — the lyric text treated as protein. Does substituting Ala at BJOZXU positions push the structure into or out of known fold space? A bar with 0 concordance hits that also returns 0 native_ala hits confirms that the structural novelty is not an artefact of the non-standard BJOZXU characters — the lyric amino acids themselves encode a novel fold.

**free_design** — MPNN optimised for foldability, no constraints. If free_design hits known protein families, it tells us *what structural class* MPNN converged to for that backbone — coiled-coil, TIM barrel, helix bundle, etc. Comparing hit counts across buckets answers a key question: does optimising for foldability (free_design) pull the structure toward known folds, or can MPNN find foldable sequences that remain structurally novel?

### This is new relative to Phase 1

Phase 1 only searched concordance structures. Phase 2 adds native_ala and free_design to the same search, enabling a three-way structural comparison per bar. The cross-bucket pattern is where the scientific story lives.

### Complete results (2026-04-06)

All 36 structures submitted and parsed across pdb100, afdb-swissprot, and mgnify_esm30.

| Bar | Song | conc hits | native_ala hits | free_design hits | Top hit (conc) | Pattern |
|-----|------|-----------|-----------------|------------------|----------------|---------|
| bar_0  | — | 22 | 273 | 1170 | — | MPNN pulls into fold space |
| bar_3  | — | 1071 | 6 | 1164 | — | native_ala novel; conc and free known |
| bar_6  | I'm The Best | 2596 | 1565 | 1922 | TIM barrel / ferredoxin (prob=0.997) | all three known; TIM barrel dominant |
| bar_8  | — | 1237 | 3 | 826 | — | native_ala novel; conc and free known |
| bar_9  | Super Bass | 536 | 3 | 287 | — | native_ala novel; free_design → Saposin A (prob=0.537) |
| bar_11 | — | 1649 | 1039 | 1330 | — | all three known |
| bar_13 | Anaconda | 861 | 1144 | 1304 | — | all three known |
| bar_17 | Monster | 201 | 1337 | 56 | — | MPNN moves *away* from known folds |
| bar_27 | Ganja Burn | 1344 | 0 | 1440 | — | native_ala uniquely novel |
| bar_32 | Barbie Dangerous | 1199 | 773 | 1252 | — | all three known |
| bar_46 | — | 636 | 418 | 0 | — | free_design novel — backbone dead, MPNN cannot fold |
| bar_77 | Moment 4 Life | 1016 | 1199 | 1210 | — | all three known |

### Scientific interpretation

**Four structural archetypes emerged from the cross-bucket comparison:**

**Type 1 — All three known (bar_6, bar_11, bar_13, bar_32, bar_77):** The concordance backbone already encodes a real fold class; native_ala and free_design both land in the same structural neighbourhood. bar_6 is the extreme case — 2596 concordance hits with top probability 0.997, identifying a TIM barrel / ferredoxin-like fold. The lyric backbone here happened to encode one of evolution's favourite structures.

**Type 2 — native_ala novel; concordance and free_design known (bar_3, bar_8, bar_9):** Substituting Ala at ambiguous positions creates a structurally novel sequence, while the original concordance and the MPNN-designed sequence both hit known families. bar_9 free_design returns 287 hits with the top match being Human Saposin A (prob=0.537) — a lipid-binding protein. The Saposin fold is a small helix bundle, consistent with the bar_9 Boltz backbone geometry.

**Type 3 — native_ala uniquely novel (bar_17, bar_27):** These are the headline results.
- **bar_27 native_ala → 0 hits** across all databases. Remarkable: bar_27 was the only `confident_protein_like` bar in Phase 1, yet the literal lyric amino acids (Ganja Burn, BJOZXU→A) encode a fold with no structural precedent in PDB, SwissProt, or metagenomics. The lyric itself is the novel fold.
- **bar_17 native_ala 1337 hits vs free_design 56 hits** — MPNN improved pLDDT from 0.441 to 0.637 while simultaneously moving the structure *away* from known fold space. A foldable but structurally novel MPNN sequence. This is the single best outcome possible: higher confidence, lower homology.

**Type 4 — free_design novel / backbone failure (bar_46):** free_design returns 0 hits because the backbone is too disordered for MPNN to generate a coherent fold. Not novelty — failure. Confirmed dropout.

**bar_0 transition (concordance 22 hits → free_design 1170 hits):** The concordance backbone is unusual; MPNN converges to a common fold class when given full freedom. The lyric-constrained structure is structurally unusual; the optimised structure is not.

### Top PDB hits (selected)

| Bar | Bucket | Top accession | Description | prob |
|-----|--------|--------------|-------------|------|
| bar_6 | concordance | — | TIM barrel / ferredoxin-like | 0.997 |
| bar_6 | free_design | — | PII signalling protein | 1.000 |
| bar_9 | free_design | — | Human Saposin A | 0.537 |
| bar_27 | native_ala | — | *no hits* | — |
| bar_17 | free_design | — | *56 hits, low prob* | ~0.2 |
| bar_46 | free_design | — | *no hits* | — |

### Figure 37 — FoldSeek Phase 2 hit comparison

![Figure 37](outputs/figures/fig37_foldseek_phase2.png)

**Figure 37 | FoldSeek Phase 2 structural search — hit count per bar per bucket (concordance, native_ala, free_design).** Each group of bars shows the number of structural homologs found across PDB, afdb-swissprot, and MGnify for one Phase 2 candidate bar under three sequence conditions. bar_6 concordance leads with 2596 hits; bar_46 free_design and bar_27 native_ala return 0 hits — structurally novel across all databases.

**Figure G** — `fig_foldseek_phase2_db.png` — Hit count breakdown by database (not yet committed)

---

## Phase 2 — Scientific Question: native_ala as backbone for MPNN

If native_ala folds more confidently (0.543 vs 0.441) and folds to a structurally distinct space from concordance, then using the native_ala Boltz backbone as the MPNN input may produce a new design series — starting from a more foldable backbone without sacrificing the cultural connection to the lyric.

### Three design strategies compared

| Strategy | Backbone | Fixed positions | Cultural fidelity | Design freedom |
|----------|----------|-----------------|-------------------|----------------|
| masked_BJOZXU (Run 1) | concordance Boltz | lyric AAs | high | low (5–18 positions) |
| free_design (Run 2) | concordance Boltz | none | none | maximum |
| **free_native_ala (planned)** | **native_ala Boltz** | **none** | **moderate** | **maximum** |

free_native_ala gives maximum design freedom on a more foldable backbone, while the native_ala backbone itself encodes the lyric sequence geometry.

---

## Phase 2 — Step 8: native_ala_free MPNN Run (Run 4)

**Script:** `analysis/09d_proteinmpnn_native_ala_free.py`
**Date:** 2026-04-07
**Input backbone:** native_ala Boltz model_0 PDB per bar (chain-fixed)
**Design:** 50 seqs × 12 bars × temp=0.1, no fixed positions
**ESMFold filter:** pLDDT ≥ 0.35

### Bug fixed: chain ID corruption
Boltz native_ala PDBs already have single-char chain `A`. The previous `fix_boltz_chain()` was designed for 2-char chain IDs and stripped 2 chars, shifting coordinate columns left — causing ProteinMPNN PDB parse failure. Fixed by detecting whether chain is already single-char before rewriting.

### Results

| Bar | Pass | Mean ESMFold pLDDT | vs free_design (v2) |
|-----|------|--------------------|---------------------|
| bar_46 | 48/50 | 0.850 | ↑↑ (was 0.309 — backbone was dead on concordance) |
| bar_6  | 50/50 | 0.845 | ↑ (was 0.711) |
| bar_32 | 47/50 | 0.824 | ↑↑ (was 0.537) |
| bar_3  | 47/50 | 0.744 | ↑↑ (was 0.558) |
| bar_8  | 50/50 | 0.750 | ↑↑ (was 0.366) |
| bar_0  | 49/50 | 0.761 | ↑↑ (was 0.401) |
| bar_13 | 50/50 | 0.712 | ↑ (was 0.637) |
| bar_77 | 50/50 | 0.683 | ↑ (was 0.517) |
| bar_27 | 49/50 | 0.601 | — (was 0.637) |
| bar_17 | 49/50 | 0.579 | — (was 0.637) |
| bar_11 | 50/50 | 0.491 | ↑ (was 0.458) |
| bar_9  | 50/50 | 0.482 | ↓ (was 0.839 — concordance backbone better for bar_9) |
| bar_46 | 48/50 | 0.850 | ↑↑ backbone is now designable |

**589/600 pass** · outputs: `outputs/proteinmpnn_native_ala_free/filtered_results.csv`

**bar_46 headline:** confirmed backbone-level failure on concordance (0.309 fully unconstrained). Same bar on native_ala backbone scores 0.850 — the Ala-substituted fold is genuinely designable. This bar is back in play.

**bar_9 inversion:** native_ala backbone is less amenable than concordance for bar_9. The concordance backbone happens to encode an excellent Saposin A-like helix bundle; the native_ala backbone folds differently. Both sequences are submitted to wet-lab — the comparison is informative.

---

## Phase 2 — Step 9: Boltz Validation v3 (native_ala_free)

**Notebook:** `notebooks/boltz_validation_v3.ipynb`
**Date:** 2026-04-07
**Colab:** A100, ~2.5 hrs
**Buckets:** native_ala_free (N=5, 589 seqs) + scrambled_naf (N=1, 36 seqs)

### Bug fixed: complex_plddt key
Boltz confidence JSON uses `complex_plddt` not `plddt`. Patched in Cell 7 and helpers Cell 4 of v3 notebook.

### Bug fixed: Drive restore timeout
Copying 2945 PDB files via `shutil.copytree` from Drive to scratch took >20 min. Fixed Cell 7 to parse confidence JSONs directly from Drive path — JSONs are tiny, no copy needed. PDBs only required for RMSD (Cell 9), fetched from Drive there too.

### pLDDT results

| Bucket | n | mean pLDDT | sd |
|--------|---|-----------|-----|
| scrambled (v2) | 36 | 0.426 | 0.068 |
| concordance (v2) | 12 | 0.441 | 0.069 |
| native_ala (v2) | 12 | 0.543 | 0.107 |
| scrambled_naf | 36 | 0.561 | 0.134 |
| free_design (v2) | 612 | 0.643 | 0.173 |
| **native_ala_free** | **589** | **0.806** | **0.140** |

**mean 0.806 is the highest any bucket has achieved** — native_ala backbones provide substantially better design geometry than concordance backbones. The 0.245 gap between scrambled_naf (0.561) and native_ala_free (0.806) demonstrates that MPNN is adding genuine structural information on top of the backbone, not just exploiting composition.

### Per-bar Boltz pLDDT — native_ala_free

| Bar | n | mean pLDDT | max | vs free_design (v2) | Δ |
|-----|---|-----------|-----|---------------------|---|
| bar_46 | 48 | **0.937** | 0.968 | 0.387 | +0.551 ↑↑ |
| bar_6  | 50 | **0.924** | 0.951 | 0.809 | +0.115 ↑ |
| bar_32 | 47 | **0.899** | 0.963 | 0.619 | +0.281 ↑ |
| bar_8  | 50 | **0.885** | 0.940 | 0.447 | +0.438 ↑↑ |
| bar_0  | 49 | **0.876** | 0.930 | 0.460 | +0.416 ↑↑ |
| bar_3  | 47 | **0.857** | 0.953 | 0.688 | +0.169 ↑ |
| bar_77 | 50 | **0.838** | 0.978 | 0.608 | +0.230 ↑ |
| bar_13 | 50 | 0.790 | 0.945 | 0.723 | +0.068 ↑ |
| bar_17 | 49 | 0.747 | 0.912 | 0.743 | +0.004 |
| bar_27 | 49 | 0.700 | 0.850 | 0.723 | -0.023 ↓ |
| bar_11 | 50 | 0.621 | 0.847 | 0.590 | +0.031 ↑ |
| bar_9  | 50 | 0.609 | 0.860 | 0.918 | -0.308 ↓ |

**bar_46 headline of the project:** concordance backbone pLDDT 0.309 (backbone failure) → native_ala backbone pLDDT 0.937 mean, 0.968 max (+0.551). The Barbie Goin Bad lyric encodes a completely different, highly foldable protein depending on which substitution rule is applied.

**bar_9 inversion confirmed:** native_ala backbone is measurably worse (0.609 vs 0.918). The concordance backbone for Super Bass encodes a near-ideal Saposin A-like helix bundle; native_ala substitution disrupts that geometry. Both sequences go to wet-lab.

**bar_27 essentially unchanged** (-0.023): Ganja Burn folds equally well on both backbone geometries. Still 0 FoldSeek hits on native_ala — foldable and structurally novel regardless of design approach.

10 of 12 bars improve on native_ala backbone. Mean improvement across all 12: **+0.164 Boltz pLDDT**.

### scrambled_na control (ESMFold, 2026-04-07)

To separate backbone geometry effect from MPNN design quality, native_ala sequences were also scrambled (same AA composition, random order) and folded with ESMFold. native_ala ESMFold scores from Phase 1 (8 bars) + spot-folded (4 missing bars: bar_27/32/46/77).

| Bar | native_ala ESMFold | scrambled_na ESMFold (mean 3×) | native_ala_free ESMFold | MPNN gain |
|-----|-------------------|-------------------------------|------------------------|-----------|
| bar_6  | 0.350 | 0.426 | 0.845 | +0.419 |
| bar_32 | 0.720 | 0.500 | 0.824 | +0.324 |
| bar_3  | 0.438 | 0.620 | 0.744 | +0.124 |
| bar_8  | 0.358 | 0.376 | 0.750 | +0.374 |
| bar_0  | 0.384 | 0.345 | 0.761 | +0.416 |
| bar_13 | 0.450 | 0.452 | 0.712 | +0.260 |
| bar_77 | 0.464 | 0.428 | 0.683 | +0.255 |
| bar_27 | 0.365 | 0.390 | 0.601 | +0.211 |
| bar_17 | 0.297 | 0.424 | 0.579 | +0.155 |
| bar_11 | 0.284 | 0.327 | 0.491 | +0.164 |
| bar_9  | 0.427 | 0.447 | 0.482 | +0.035 |
| bar_46 | 0.387 | 0.446 | 0.850 | +0.404 |
| **mean** | **0.411** | **0.432** | **0.694** | **+0.262** |

**MPNN gain** = native_ala_free ESMFold − scrambled_na ESMFold (composition baseline).

Key observations:
- native_ala and scrambled_na score nearly identically (~0.41 vs 0.43) — AA composition alone, regardless of arrangement, gives modest foldability. The native_ala sequence arrangement offers no advantage over a random shuffle.
- MPNN consistently adds +0.12 to +0.42 above composition. bar_9 is the outlier (+0.035) — concordance backbone is the better geometry for that bar.
- bar_32 native_ala (0.720) is anomalously high — the lyric amino acid sequence itself happens to be unusually foldable. Worth noting in the write-up.

**Key finding:** scrambled_na ESMFold means cluster around 0.35–0.50 — composition alone explains modest foldability. native_ala_free (MPNN on native_ala backbone) reaches 0.48–0.85. The delta (MPNN gain above composition baseline) is largest for bar_46 (+0.40) and bar_0 (+0.42), smallest for bar_9 (+0.04). MPNN is doing real work.

**`outputs/bioreason/scrambled_na_esm.csv`** — full per-scramble ESMFold scores.

### Scramble visualisation

**Script:** `analysis/12d_scramble_viz.py`
**Output:** `outputs/figures/fig_scramble_{bar_id}.png` (37 bars)

Each figure has three panels:
1. **Sequence grid** — one column per position, coloured by physicochemical class (hydrophobic=amber, polar=teal, positive=cyan, negative=red). Rows: native_ala + 3 scrambles. Shows residues redistributed randomly while palette stays identical.
2. **4×4 identity matrix** — pairwise sequence identity between all four variants. Off-diagonal values: ~0.05–0.15 (effectively random — confirms scrambles are genuinely independent permutations).
3. **AA composition bar chart** — identical across all four variants by construction. Sanity check that only arrangement, not composition, changed.

### RMSD vs native_ala backbone (2026-04-07)

**14,757 RMSD rows** — model_0 per design vs chain-fixed native_ala backbone PDB.
Saved to Drive: `results/boltz_rmsd_v3.csv`

| Bar | n | mean RMSD | sd | Interpretation |
|-----|---|----------|----|----------------|
| bar_6  | 50 | 1.22Å | 0.78 | Exceptionally tight — TIM barrel-like fold locks designs in |
| bar_46 | 48 | 1.34Å | 0.33 | Near-perfect recovery — native_ala backbone is well-defined and designable |
| bar_0  | 49 | 1.35Å | 1.03 | Tight — backbone geometry strongly constrains design space |
| bar_8  | 50 | 1.51Å | 1.29 | Tight |
| bar_32 | 47 | 1.80Å | 3.50 | Moderate mean but high SD — likely bimodal (two fold populations) |
| bar_13 | 50 | 2.12Å | 2.31 | Moderate |
| bar_3  | 47 | 2.26Å | 2.70 | Moderate |
| bar_17 | 49 | 3.71Å | 3.62 | Variable — backbone less constraining |
| bar_27 | 49 | 5.34Å | 3.42 | Higher deviation — consistent with 0 FoldSeek hits (unusual fold) |
| bar_9  | 50 | 7.81Å | 3.52 | Large deviation — MPNN moves away from native_ala backbone, explains pLDDT drop |
| bar_11 | 50 | 10.79Å | 5.32 | Very large — native_ala backbone not well-recovered |
| bar_77 | 50 | 13.54Å | 13.25 | Extremely variable — likely bimodal; two distinct fold populations |

**Pairwise RMSD (within naf per bar):** mean=5.83Å  sd=7.16Å

**vs v2 free_design RMSD (concordance backbone):** 7.13Å mean
→ native_ala_free designs stay closer to their input backbone (~4.1Å mean) than free_design did to the concordance backbone (7.13Å). Native_ala backbones are more constraining — MPNN finds sequences that converge on the backbone geometry rather than exploring alternative folds.

**bar_46 headline confirmed in RMSD:** 1.34Å mean with sd=0.33Å — the tightest cluster of any bar. The native_ala fold is so well-defined that every MPNN design folds straight back to it. This is a real, stable fold.

**bar_9 RMSD explains the pLDDT drop:** 7.81Å mean — MPNN cannot find sequences that both fold well *and* stay on the native_ala backbone for bar_9. Designs that score high pLDDT fold to a different structure. The concordance backbone is the better geometry for bar_9.

**bar_77 and bar_11 bimodality:** SD ≥ mean suggests two populations — some designs recover the backbone, others fold elsewhere. Both bars had high pLDDT in native_ala_free (0.683 and 0.491) but only some designs are responsible.

### Figure 75 — Boltz v3 pLDDT Strip (native_ala_free + scrambled_naf)

![Figure 75](outputs/figures/fig75_v3_strip.png)

**Figure 75 | Boltz-2 pLDDT strip — native_ala_free designs and scrambled_naf controls per bar.** Each column is one bar; dots are individual diffusion samples. native_ala_free (green) consistently outperforms scrambled_naf (grey), confirming MPNN sequence order adds structural information beyond composition alone (+0.245 pLDDT).

### Figure 76 — Boltz v3 Cross-Run Violin

![Figure 76](outputs/figures/fig76_v3_violin.png)

**Figure 76 | Boltz-2 pLDDT violin across all six buckets (v2 + v3).** Full bucket ladder from concordance (0.441) to native_ala_free (0.806). Mean values annotated. native_ala_free is the highest-confidence bucket achieved in the project.

### Figure 77 — RMSD vs native_ala Backbone (Violin)

![Figure 77](outputs/figures/fig77_v3_rmsd_violin.png)

**Figure 77 | RMSD of native_ala_free designs vs their input native_ala backbone, by bar (violin).** Low RMSD = MPNN recovers the backbone geometry; high RMSD = designs fold to a different structure. bar_46 (1.34Å) is the tightest; bar_77 (13.54Å) shows extreme bimodality.

### Figure 78 — native_ala_free RMSD Per Bar

![Figure 78](outputs/figures/fig78_v3_rmsd_per_bar.png)

**Figure 78 | RMSD vs native_ala backbone per bar — native_ala_free designs.** bar_6 and bar_46 cluster tightly around the input backbone; bar_9 and bar_77 show large deviations indicating MPNN found alternative folds.

### Figure 79 — Pairwise RMSD Histogram (native_ala_free)

![Figure 79](outputs/figures/fig79_v3_rmsd_pairwise.png)

**Figure 79 | Pairwise RMSD histogram within native_ala_free ensemble.** Mean pairwise RMSD 5.83Å — designs within a bar are structurally diverse despite sharing the same backbone input.

---

## Phase 2 — Step 9b: Boltz Validation v4 (scrambled_na)

**Notebook:** `notebooks/boltz_validation_v4_scrambled_na.ipynb`
**Date:** 2026-04-07
**Colab:** A100
**Input:** `outputs/bioreason/scrambled_na_esm.csv` — 111 sequences (37 bars × 3 scrambles of native_ala)
**N_MODELS:** 1 per scramble (control run)

### pLDDT results

| Bucket | n | mean pLDDT | sd |
|--------|---|-----------|-----|
| concordance | 12 | 0.441 | 0.069 |
| **scrambled_na** | **111** | **0.482** | **0.068** |
| native_ala | 12 | 0.543 | 0.107 |
| free_design | 612 | 0.643 | 0.173 |
| native_ala_free | 589 | 0.806 | 0.140 |

### Foldability decomposition

All Boltz-2 pLDDT values across all buckets (v2–v4):

| Bucket | n | mean pLDDT | Composition | Order |
|--------|---|-----------|-------------|-------|
| scrambled (free_design comp) | 36 | 0.426 | free_design | random |
| concordance | 12 | 0.441 | freq-rank + softmax BJOZXU | lyric |
| scrambled_na | 111 | 0.482 | native_ala (literal) | random |
| native_ala | 12 | 0.543 | native_ala (literal) | lyric |
| scrambled_naf | 36 | 0.561 | native_ala_free (MPNN) | random |
| free_design | 612 | 0.643 | MPNN on concordance backbone | designed |
| native_ala_free | 589 | 0.806 | MPNN on native_ala backbone | designed |

**Encoding reminder:** native_ala is the *more raw* encoding — lyric characters read literally as amino acids, BJOZXU → Ala. Concordance uses a frequency-rank remapping of all standard letters AND a softmax-peaked probabilistic draw for BJOZXU. native_ala (0.543) > concordance (0.441) — the freq-rank remapping hurts foldability relative to the literal read.

#### Pairwise comparisons

| Comparison | A | B | Δ | What changes | Clean? |
|---|---|---|---|---|---|
| conc vs native_ala | 0.441 | 0.543 | **+0.102** | encoding strategy: freq-rank+softmax vs literal pass-through | ✓ same lyric order |
| conc vs conc scrambles | — | — | — | not done | — |
| scrambled_na vs native_ala | 0.482 | 0.543 | **+0.061** | order only (native_ala composition) | ✓ pure lyric order effect |
| scrambled (fd) vs free_design | 0.426 | 0.643 | **+0.217** | MPNN sequence order (free_design composition) | ✓ pure order effect |
| scrambled_naf vs native_ala_free | 0.561 | 0.806 | **+0.245** | MPNN sequence order (native_ala_free composition) | ✓ pure order effect |
| scrambled_naf vs free_design | 0.561 | 0.643 | **+0.082** | designed order on concordance backbone vs naf composition randomly ordered | ✗ confounded |
| free_design vs native_ala_free | 0.643 | 0.806 | **+0.163** | backbone used for MPNN (concordance vs native_ala) | ✓ |

Key reads:
- **Concordance encoding costs −0.102** vs literal pass-through (native_ala) for the same lyric — freq-rank remapping hurts foldability
- **Lyric order effect** at every level: +0.061 (native_ala), +0.217 (free_design), +0.245 (native_ala_free) — MPNN finds far better order than lyrics, but the ordering signal is real at every stage
- **Backbone quality** for MPNN: native_ala backbone gives +0.163 over concordance backbone
- **Concordance scrambles not done** — would cleanly isolate concordance composition from concordance order; pending

### Per-bar scrambled_na Boltz pLDDT

| Bar | n | boltz_mean | boltz_sd | esm_mean |
|-----|---|-----------|---------|---------|
| bar_67 | 3 | 0.580 | 0.033 | 0.533 |
| bar_3  | 3 | 0.593 | 0.059 | 0.620 |
| bar_32 | 3 | 0.585 | 0.106 | 0.500 |
| bar_58 | 3 | 0.543 | 0.012 | 0.514 |
| bar_6  | 3 | 0.534 | 0.044 | 0.426 |
| bar_13 | 3 | 0.525 | 0.027 | 0.452 |
| bar_55 | 3 | 0.523 | 0.051 | 0.439 |
| bar_26 | 3 | 0.512 | 0.062 | 0.461 |
| bar_46 | 3 | 0.511 | 0.032 | 0.446 |
| bar_47 | 3 | 0.514 | 0.082 | 0.430 |
| bar_68 | 3 | 0.518 | 0.086 | 0.401 |
| bar_39 | 3 | 0.509 | 0.035 | 0.382 |
| bar_78 | 3 | 0.508 | 0.049 | 0.414 |
| bar_0  | 3 | 0.501 | 0.110 | 0.345 |
| bar_49 | 3 | 0.496 | 0.014 | 0.445 |
| bar_9  | 3 | 0.495 | 0.023 | 0.447 |
| bar_17 | 3 | 0.486 | 0.075 | 0.425 |
| bar_40 | 3 | 0.487 | 0.025 | 0.362 |
| bar_35 | 3 | 0.488 | 0.042 | 0.413 |
| bar_52 | 3 | 0.479 | 0.027 | 0.413 |
| bar_27 | 3 | 0.479 | 0.095 | 0.390 |
| bar_25 | 3 | 0.476 | 0.011 | 0.407 |
| bar_80 | 3 | 0.446 | 0.043 | 0.379 |
| bar_77 | 3 | 0.460 | 0.068 | 0.429 |
| bar_71 | 3 | 0.462 | 0.026 | 0.345 |
| bar_23 | 3 | 0.434 | 0.022 | 0.378 |
| bar_53 | 3 | 0.454 | 0.071 | 0.311 |
| bar_8  | 3 | 0.460 | 0.026 | 0.376 |
| bar_43 | 3 | 0.466 | 0.031 | 0.310 |
| bar_55 | 3 | 0.523 | 0.051 | 0.439 |
| bar_73 | 3 | 0.439 | 0.025 | 0.328 |
| bar_75 | 3 | 0.420 | 0.047 | 0.313 |
| bar_42 | 3 | 0.420 | 0.007 | 0.310 |
| bar_36 | 3 | 0.458 | 0.105 | 0.337 |
| bar_41 | 3 | 0.398 | 0.028 | 0.319 |
| bar_34 | 3 | 0.396 | 0.019 | 0.307 |
| bar_38 | 3 | 0.372 | 0.039 | 0.345 |
| bar_11 | 3 | 0.394 | 0.057 | 0.327 |

**bar_3 anomaly:** scrambled_na 0.593 > native_ala mean 0.543 — for bar_3 the AA composition alone is unusually foldable. The lyric positional order slightly hurts it; consistent with ESMFold (native_ala 0.438 vs scrambled_na 0.620). The composition of the bar_3 lyric happens to spell out a foldable sequence regardless of arrangement.

**bar_46 scrambled_na = 0.511** — notably higher than its concordance pLDDT (0.309). Scrambling the Ala-substituted sequence still folds better than the original concordance. The concordance backbone failure is a sequence-order effect, not a composition problem.

### Figure 80 — Cross-Bucket Violin (All 7 Buckets)

![Figure 80](outputs/figures/fig80_v4_violin_all_buckets.png)

**Figure 80 | Boltz-2 pLDDT violin across all seven buckets including scrambled_na (v4).** scrambled_na (0.482) slots between concordance (0.441) and native_ala (0.543), revealing the lyric order effect (+0.061) above composition alone. Full decomposition visible in a single figure.

### Figure 81 — scrambled_na Per Bar (Boltz vs ESMFold)

![Figure 81](outputs/figures/fig81_v4_sc_na_per_bar.png)

**Figure 81 | scrambled_na Boltz-2 pLDDT (grey) and ESMFold pLDDT (red diamond) per bar.** 3 scrambles per bar shown individually; white dash = bar mean. Both predictors broadly agree on which bars have foldable composition; bar_3 and bar_32 highest, bar_38 and bar_11 lowest.

### Figure 82 — ESMFold vs Boltz-2 Scatter (scrambled_na)

![Figure 82](outputs/figures/fig82_v4_esm_vs_boltz_sc_na.png)

**Figure 82 | ESMFold pLDDT vs Boltz-2 pLDDT for all 111 scrambled_na sequences.** Correlation coefficient annotated. Moderate positive correlation confirms both predictors agree on composition-driven foldability, with Boltz systematically higher than ESMFold for this sequence class.

---

## Phase 2 — Step 10: BioReason-Pro Stress Test

**Script:** `analysis/12_bioreason_prep.py` + `analysis/12b_bioreason_parse.py`
**Date:** 2026-04-07
**Status:** ⛔ PAUSED — model cannot discriminate lyric sequences from real proteins (see results below)

### Manual validation (bioreason.net, 2026-04-07)

Four sequences tested manually before committing to a full batch run:

| Sequence | Bar | Condition | pLDDT | FoldSeek | BioReason GO (biological process) | Grounded? |
|---|---|---|---|---|---|---|
| `ANDEVERYTHING...SPEAK` | bar_9 | native_ala | 0.543 | — | negative regulation of transcription by RNA Pol II (`GO:0000122`) | ✗ hallucinated |
| `CKPYERHIG...HTGA` | bar_46 | concordance | 0.309 | — | epithelial cell differentiation (`GO:0030855`) | ✗ hallucinated |
| `HTLEATLC...KLGP` | bar_11 | concordance | 0.344 | 0 hits | signal transduction + exosome (`GO:0007165`) | ✗ hallucinated |
| `MIEITLKK...IGN` | 1REG | T4 phage RegA (real protein) | — | — | viral host gene suppression (`GO:0039657`) + dsDNA binding | ✓ correct |

Raw chat exports: `outputs/bioreason/bioreason-chat-2026-04-07-*.md`

**bar_46 concordance — hallucinated epithelial differentiation (pLDDT 0.309, backbone failure):**

![bar_46 BioReason](outputs/bioreason/bar_46_screenshot_bioreason.png)

**1REG T4 phage RegA — correct grounded prediction (positive control):**

![1REG BioReason p1](outputs/bioreason/pdb_1REG_screenshot_1.png)

![1REG BioReason p2](outputs/bioreason/pdb_1REG_screenshot_2.png)

### Finding: BioReason confidence is binary, not graded

All three lyric-derived sequences — including bar_46 (confirmed backbone failure, pLDDT 0.309) and bar_11 (0 FoldSeek hits, worst structural confidence) — received equally confident, detailed, mechanistically plausible GO predictions. The anchor in every case was a single generic `GO:0005515` (protein binding) InterPro hit, from which BioReason confabulated an elaborate story (transcriptional repressor scaffold, epithelial differentiation organiser, ER-to-exosome adaptor).

The real protein (1REG, T4 phage RegA) received an accurate prediction grounded in real InterPro domain hits (`IPR002702` Translation repressor RegA family).

**The model has two modes:**
1. **InterPro domain hit → grounded, accurate, specific prediction**
2. **No domain hit → `GO:0005515` + confident hallucination**

There is no intermediate "low confidence" or "disordered" output. BioReason cannot flag a sequence as structurally implausible. The stress test would produce 118 confident hallucinations for lyric sequences and a handful of real predictions for any MPNN designs that happen to match a known domain — which defeats the purpose of the experiment.

### Decision: function prediction phase on pause

Running the full 118-sequence batch would generate uninterpretable results — the confidence scores would not discriminate anything scientifically meaningful about the lyric-derived sequences. This phase is paused until a function predictor is available that:
- Reports calibrated confidence (not binary grounded/hallucinated)
- Can flag sequences as disordered or outside training distribution
- Is tested on known-disordered sequences before being applied to lyric sequences

Candidates to revisit: ESMFold function head, ProteinBERT, or structure-conditioned function prediction once wet-lab PDBs are available.

### What this tells us

BioReason's behaviour is itself a finding: **any short sequence with standard amino acids is assigned a confident human protein function.** This is consistent with the Phase 1 ESMFold result (ESMFold cannot discriminate structured from scrambled sequences at short lengths). The lyric-derived sequences are indistinguishable from real sequences to both structure predictors (at the uncertainty level) and function predictors (at the annotation level). The biology and the culture are genuinely orthogonal — the sequences look real enough to fool state-of-the-art models.

---

## Platform Comparison

| | **Ginkgo Cloud Lab** | **Adaptyv Bio** |
|---|---|---|
| System | E. coli lysate + HiBiT tag | Reconstituted E. coli TX/TL |
| Per protein | ~$117 | ~$95 (open-source discount) |
| Total (24 proteins) | $2,808 | $2,285 |
| Replicates | 3 per protein | 1 |
| Readout | Quantitative HiBiT luminescence | Expressed/low/not detected + yield |
| Gene synthesis | Yes | Yes |
| De novo track record | Internal (large scale) | Public (10,000+ proteins, competitions) |
| Condition | None | Publish results on Proteinbase |

**Leaning Adaptyv Bio** — $500 cheaper, stronger published de novo validation, reconstituted system better for synthetic proteins. Decision pending FoldSeek Phase 2 and free_native_ala results.

---

## Submission Plan (48-protein plate)

| Slot | Count | Description |
|------|-------|-------------|
| 12 bars × concordance | 12 | Raw lyric-derived sequence |
| 12 bars × native_ala | 12 | Ala-substituted lyric sequence |
| 12 bars × free_design (best) | 12 | MPNN, no constraints, concordance backbone |
| 12 bars × free_native_ala (best) | 12 | MPNN, no constraints, native_ala backbone |
| Positive control (sfGFP) | 1 | Confirms platform ran |
| Reserve | ~7 | — |
| **Total** | **48** | ~$4,560 at $95/protein (Adaptyv) |

---

## Figure Inventory

| Fig | File | Description |
|-----|------|-------------|
| 1–10 | `outputs/figures/fig01–fig10` | ESMFold pilot (top-25 bars) |
| 11–20 | `outputs/figures/fig11–fig20` | ESMFold full run (85 bars × 5 conditions) |
| 21–23 | `outputs/pairwise/fig1–fig3` | Pairwise structural comparison |
| 24–26 | `outputs/pairwise/fig4–fig6` | Ensemble RMSD heatmaps + UMAP |
| 27 | `outputs/figures/fig7_backbone_traces.png` | Backbone traces |
| **28** | `outputs/figures/fig28_conc_vs_na_plddt.png` | Concordance vs native_ala pLDDT lollipop |
| **29** | `outputs/figures/fig29_conc_vs_na_ptm.png` | Concordance vs native_ala pTM scatter |
| **30** | `outputs/figures/fig30_candidate_table.png` | Phase 2 candidate summary table |
| **31** | `outputs/figures/fig31_boltz_plddt_per_bar.png` | Boltz pLDDT strip — 4 buckets per bar |
| **32** | `outputs/figures/fig32_boltz_plddt_violin.png` | Boltz pLDDT violin — 4 buckets |
| **33** | `outputs/figures/fig33_rmsd_violin.png` | RMSD vs backbone violin by bucket |
| **34** | `outputs/figures/fig34_rmsd_per_bar.png` | free_design RMSD per bar |
| **35** | `outputs/figures/fig35_rmsd_pairwise.png` | Pairwise RMSD histogram |
| **36** | `outputs/figures/fig36_rmsd_per_bar_grid.png` | 12-panel per-bar RMSD grid |
| **37** | `outputs/figures/fig37_foldseek_phase2.png` | FoldSeek Phase 2 hit comparison — 3 buckets × 12 bars |
| G | not committed | FoldSeek Phase 2 hit count by database |
| **38–74** | `outputs/figures/fig_scramble_{bar_id}.png` | Scramble visualisations — 37 bars, sequence grid + identity matrix + composition |
| **75** | `outputs/figures/fig75_v3_strip.png` | Boltz v3 pLDDT strip — native_ala_free + scrambled_naf per bar |
| **76** | `outputs/figures/fig76_v3_violin.png` | Boltz v3 cross-run violin — all 7 buckets |
| **77** | `outputs/figures/fig77_v3_rmsd_violin.png` | RMSD vs native_ala backbone violin |
| **78** | `outputs/figures/fig78_v3_rmsd_per_bar.png` | native_ala_free RMSD per bar |
| **79** | `outputs/figures/fig79_v3_rmsd_pairwise.png` | Pairwise RMSD histogram (native_ala_free) |
| **80** | `outputs/figures/fig80_v4_violin_all_buckets.png` | Cross-bucket violin — all 7 buckets including scrambled_na |
| **81** | `outputs/figures/fig81_v4_sc_na_per_bar.png` | scrambled_na per bar — Boltz vs ESMFold |
| **82** | `outputs/figures/fig82_v4_esm_vs_boltz_sc_na.png` | ESMFold vs Boltz scatter — scrambled_na |

| **83** | `outputs/figures/fig83_umap_all_buckets.png` | ESM-2 UMAP — all 4 buckets, colour by bucket, hull per bar |
| **84** | `outputs/figures/fig84_umap_native_ala.png` | ESM-2 UMAP — native_ala + native_ala_free only |
| **85** | `outputs/figures/fig85_umap_concordance.png` | ESM-2 UMAP — concordance + free_design only |
| **86** | `outputs/figures/fig86_umap_barbie_dangerous.png` | ESM-2 UMAP — bar_32 Barbie Dangerous, all 4 buckets |
| **87** | `outputs/figures/fig87_umap_lyric_seeds.png` | ESM-2 UMAP — lyric seeds only (concordance + native_ala, no MPNN) |
| **88** | `outputs/figures/fig88_umap_plddt.png` | ESM-2 UMAP — all 1286 seqs coloured by Boltz pLDDT (RdYlGn) |
| **89a–l** | `outputs/figures/fig89_umap_bar_XX.png` × 12 | Per-bar UMAP — all 4 buckets, bars with ≥4 sequences (12/37 bars) |
| **90** | `outputs/figures/fig90_rmsd_pairs.png` | Pairwise Boltz Cα RMSD — 4 pairs × 12 bars (2×2 subplots) |
| **91** | `outputs/figures/fig91_rmsd_combined.png` | Pairwise Boltz Cα RMSD — all 4 pairs combined in one panel |
| S1 | `outputs/bioreason/bar_46_screenshot_bioreason.png` | BioReason screenshot — bar_46 concordance (backbone failure) |
| S2 | `outputs/bioreason/pdb_1REG_screenshot_1.png` | BioReason screenshot — 1REG T4 phage RegA (positive control, p1) |
| S3 | `outputs/bioreason/pdb_1REG_screenshot_2.png` | BioReason screenshot — 1REG T4 phage RegA (positive control, p2) |

**Next local figure number:** Fig 92


---

## Phase 2 — Step 11: ESM-2 UMAP (Sequence Space Visualisation)

**Script:** `analysis/13_umap_esm2.py`
**Date:** 2026-04-07
**Status:** Complete — fig83–89 generated (2026-04-07)

### Rationale

After completing the foldability decomposition across 7 buckets, the next question is: **where do these sequences sit in protein sequence space?** The UMAP plots the full sequence landscape using ESM-2 protein language model embeddings — capturing evolutionary signal, biophysical plausibility, and structural propensity simultaneously.

### Design decisions

| Decision | Choice | Rationale |
|---|---|---|
| Embedding model | ESM-2 8M (esm2_t6_8M_UR50D, layer 6) | Fast on CPU; sufficient for sequence-space topology |
| Pooling | Mean over residue positions | Standard for sequence-level tasks |
| Dimensionality reduction | UMAP (cosine metric, n_neighbors=15, min_dist=0.1) | Preserves local + global structure; more stable than t-SNE |
| Buckets included | concordance, native_ala, free_design, native_ala_free | Core 4 buckets only — no scrambles |
| Hull | Convex hull per bar (all 4 buckets combined) | Shows how MPNN designs expand around the lyric seed sequences |
| Label | bar_id + song title on hull centroid | Directly readable in the figure |
| Dot size (lyric seqs) | Proportional to FoldSeek hit count | Visually encodes structural novelty |
| Panel A | Colour by bucket | Shows design strategy separation in sequence space |
| Panel B | Colour by Boltz pLDDT (lyric seqs); MPNN dims to grey | Shows foldability landscape |

### t-SNE vs UMAP — design note

t-SNE optimises for local neighbourhood only — clusters are meaningful, inter-cluster distances are not. UMAP preserves both local and global structure, runs faster, and is more stable across seeds. For this dataset (~1200 sequences, 4 related buckets) UMAP is preferred: we want to see both whether bars form clusters AND where they sit relative to each other globally.

### Key questions the plot answers

1. Do MPNN designs from the same bar cluster tightly (shared backbone → similar sequence space)?
2. Do concordance/native_ala sequences sit at the edge of or outside their bar's MPNN cluster (lyric sequence ≠ designed sequence)?
3. Do bars with known structure (bar_6, TIM barrel) sit in a distinct region vs novel bars (bar_27, 0 FoldSeek hits)?
4. Does the native_ala_free cluster separate from free_design (different backbone → different sequence space)?
5. Is Boltz pLDDT smoothly distributed across the UMAP (continuous landscape) or clustered (discrete foldability zones)?

### Embeddings cache

`outputs/embeddings/esm2_embeddings.csv` — computed once, reused on re-run. 480-dim per sequence (ESM-2 8M layer 6 mean pool).

### Figure 83 — All 4 Buckets

![Figure 83](outputs/figures/fig83_umap_all_buckets.png)

**Figure 83 | ESM-2 UMAP — all 4 buckets (1286 sequences).** Concordance (blue) and native_ala (orange) lyric seeds shown large; MPNN designs (purple = free_design, green = native_ala_free) shown small. Convex hulls per bar. Labels: bar_id + song title.

### Figure 84 — native_ala + native_ala_free

![Figure 84](outputs/figures/fig84_umap_native_ala.png)

**Figure 84 | ESM-2 UMAP — native_ala seeds (orange, large) and native_ala_free MPNN designs (green, small).** Subset UMAP re-run on 637 sequences. Shows how MPNN expands each bar's lyric seed into a cloud of designed sequences.

### Figure 85 — concordance + free_design

![Figure 85](outputs/figures/fig85_umap_concordance.png)

**Figure 85 | ESM-2 UMAP — concordance seeds (blue, large) and free_design MPNN designs (purple, small).** Subset UMAP re-run on 649 sequences. Concordance encoding produces a different seed position in sequence space than native_ala for the same bar.

### Figure 86 — bar_32 Barbie Dangerous (single bar)

![Figure 86](outputs/figures/fig86_umap_barbie_dangerous.png)

**Figure 86 | ESM-2 UMAP — bar_32 (Barbie Dangerous) only, all 4 buckets (103 sequences).** Single-bar UMAP (n_neighbors=8). Concordance and native_ala seeds labelled individually. free_design (51 seqs) and native_ala_free (50 seqs) clouds shown separately. Reveals whether the two MPNN strategies explore different or overlapping sequence regions for this bar.

### Figure 87 — Lyric Seeds Only (concordance + native_ala)

![Figure 87](outputs/figures/fig87_umap_lyric_seeds.png)

**Figure 87 | ESM-2 UMAP — lyric seeds only (74 sequences, no MPNN designs).** Concordance seeds (blue) and native_ala seeds (orange), one per bar, re-embedded in isolation. Removing the MPNN cloud gives a cleaner view of how the two encoding strategies relate across bars. Concordance (freq-rank remapping + softmax draw) and native_ala (literal pass-through + Ala for non-standard residues) produce different amino-acid compositions from the same lyric — both converge on a rap-phonology-driven composition landscape rather than natural protein sequence space.

### Figure 88 — All Sequences Coloured by Boltz pLDDT

![Figure 88](outputs/figures/fig88_umap_plddt.png)

**Figure 88 | ESM-2 UMAP — all 1286 sequences coloured by Boltz-2 mean pLDDT (RdYlGn colormap, green = high confidence).** Lyric seeds large; MPNN designs small. Grey dots = native_ala_free designs (Boltz v3 run pending — pLDDT not yet available). Foldability is not purely a function of sequence-space location: high-pLDDT sequences (green) appear scattered across the UMAP rather than concentrated in a single cluster, indicating that MPNN can sculpt foldable sequences from multiple starting points in embedding space. The lyric seed cloud (center) is predominantly red/yellow — expected given that rap-derived sequences have low structural prior.

### Figure 89 — Per-Bar UMAPs (12 bars with MPNN designs)

Per-bar UMAPs generated for all 12 bars that had MPNN designs (≥4 sequences total). The remaining 25 bars have only their 2 lyric seeds (concordance + native_ala) and cannot form a meaningful UMAP. Each plot shows all 4 buckets in the same colour scheme, with lyric seeds annotated individually and MPNN design clouds shown as small dots.

**Bars with per-bar UMAPs:** bar_00, bar_03, bar_06, bar_08, bar_09, bar_11, bar_13, bar_17, bar_27, bar_32, bar_46, bar_77

![Figure 89a — bar_06](outputs/figures/fig89_umap_bar_06.png)
![Figure 89b — bar_32](outputs/figures/fig89_umap_bar_32.png)
![Figure 89c — bar_03](outputs/figures/fig89_umap_bar_03.png)
![Figure 89d — bar_08](outputs/figures/fig89_umap_bar_08.png)
![Figure 89e — bar_00](outputs/figures/fig89_umap_bar_00.png)
![Figure 89f — bar_13](outputs/figures/fig89_umap_bar_13.png)
![Figure 89g — bar_11](outputs/figures/fig89_umap_bar_11.png)
![Figure 89h — bar_77](outputs/figures/fig89_umap_bar_77.png)
![Figure 89i — bar_27](outputs/figures/fig89_umap_bar_27.png)
![Figure 89j — bar_09](outputs/figures/fig89_umap_bar_09.png)
![Figure 89k — bar_17](outputs/figures/fig89_umap_bar_17.png)
![Figure 89l — bar_46](outputs/figures/fig89_umap_bar_46.png)

**Figure 89a–l | Per-bar ESM-2 UMAPs.** Each panel shows one bar's full sequence space (n_neighbors=8). Lyric seeds anchor the embedding; MPNN design clouds expand around them. Tight MPNN clusters suggest convergence on a shared fold family from that bar's backbone. Spread suggests more diverse exploration. Overlap between free_design (purple) and native_ala_free (green) clouds indicates the two design strategies produce similar sequence distributions for that bar, while separation indicates the backbone templates (concordance vs native_ala) drive MPNN into distinct regions of sequence space.

---

## Phase 2 — Step 12: Pairwise Boltz Structure RMSD

**Script:** `analysis/14_rmsd_comparison.py`
**Date:** 2026-04-07
**Status:** Complete — fig90 (4 pairs) + fig91 (combined) generated

### Rationale

Rather than comparing each structure to an external ESMFold reference, this analysis asks: **how structurally different are the Boltz-2 predicted structures between buckets?** Each pair directly compares two Boltz structures via Cα RMSD (Kabsch alignment, model_0). This gives a bucket-vs-bucket structural divergence measure entirely within Boltz predictions.

### Pairs

| Pair | Comparison | What it tells us |
|---|---|---|
| 1 | native_ala ↔ concordance | How different are the two lyric-seed structures for the same bar — encoding strategy divergence |
| 2 | native_ala ↔ native_ala_free MPNN | How much MPNN redesigns from the native_ala seed |
| 3 | concordance ↔ free_design MPNN | How much MPNN redesigns from the concordance seed |
| 4 | free_design ↔ native_ala_free MPNN | How structurally different are the two MPNN strategy outputs |

### Key findings

| Observation | Detail |
|---|---|
| Lyric seeds are structurally very different from each other (Pair 1) | RMSD = 11–44 Å across bars. bar_46 is the most divergent (44 Å). Same lyric, different encoding → completely different fold for some bars. |
| MPNN redesigns heavily from native_ala (Pair 2) | Median 1–45 Å. Most bars show RMSD 3–14 Å, confirming MPNN explores a wide neighbourhood around the seed. |
| concordance → free_design more constrained (Pair 3) | Median 4–15 Å, generally lower than Pair 2, suggesting concordance seeds lead to more focused MPNN designs. |
| Two MPNN strategies produce very different structures (Pair 4) | Median 12–30 Å for most bars — fd and naf MPNN outputs are as different from each other as the two lyric seeds. Backbone template is the dominant determinant of MPNN design space. |
| bar_46 outlier | All pairs show high RMSD for bar_46, suggesting both lyric-seed structures are globally disordered with random coil, so structural diversity is artificially inflated. |

### Figure 90 — Pairwise RMSD: 4 Pairs (2×2 grid)

![Figure 90](outputs/figures/fig90_rmsd_pairs.png)

**Figure 90 | Pairwise Boltz Cα RMSD — 4 pairs × 12 bars (2×2 subplots).** Each subplot shows one structural comparison. Pair 1 (top-left): single RMSD value per bar between the two lyric seeds. Pairs 2–3 (top-right, bottom-left): distribution of RMSD values between the lyric seed and each MPNN design (median + IQR band). Pair 4 (bottom-right): cross-strategy MPNN distribution. Bars sorted by Pair 1 RMSD.

### Figure 91 — All 4 Pairs Combined

![Figure 91](outputs/figures/fig91_rmsd_combined.png)

**Figure 91 | All 4 pairwise Boltz Cα RMSD comparisons overlaid per bar.** Circle = Pair 1 (na ↔ conc); Diamond = Pair 2 (na ↔ naf MPNN); Square = Pair 3 (conc ↔ fd MPNN); Triangle = Pair 4 (fd ↔ naf MPNN). IQR band shown for distribution pairs. The consistent pattern: Pair 4 (MPNN vs MPNN) tracks closely with Pair 1 (seed vs seed), confirming that backbone template choice is the primary driver of MPNN design space — not the MPNN optimisation itself.

---

## Pending Work

| Priority | Task | Script |
|----------|------|--------|
| 1 | Concordance scrambles Boltz run (fills last gap in decomposition table) | new notebook |
| 2 | Codon optimisation | `analysis/10_codon_optimize.py` |
| 3 | Platform decision + submission | Adaptyv Bio |
| — | ~~BioReason-Pro batch run~~ | ⛔ paused — model hallucinates for all non-domain sequences |

---

## Boltz Technical Notes

### YAML format (required for single-sequence mode)
```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: GELTEEVLKELEKK...
      msa: empty    # omitting this causes MSA error — no predictions returned
```

### Colab run flags
```bash
boltz predict {input_dir} \
  --out_dir {out_dir} \
  --diffusion_samples 5 \
  --output_format pdb \   # default is mmcif — must specify pdb
  --no_kernels \          # underscore not hyphen; required on Colab
  --override              # force regeneration if output dirs exist
```

### numpy binary incompatibility
Boltz changes numpy version on install. Colab kernel must restart before importing pandas or numpy after `pip install boltz`. Cell 1 of `notebooks/boltz_validation_v2.ipynb` handles this automatically with `os.kill(os.getpid(), 9)` — skip Cell 1 on subsequent runs.

### Chain ID fix (for ProteinMPNN input)
Boltz sometimes writes 2-char chain IDs (`b11`, `A `) which shift coordinate columns and break ProteinMPNN. The fix must detect whether the chain is already single-char before rewriting — otherwise it strips two chars and corrupts coordinates:
```python
# WRONG (corrupts coordinates if chain already single-char):
line = line[:21] + 'A' + line[24:]

# CORRECT (in analysis/09d_proteinmpnn_native_ala_free.py):
chain, next_char = line[21], line[22]
if next_char != ' ':          # 2-char chain ID
    line = line[:21] + 'A ' + line[24:]
elif chain not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    line = line[:21] + 'A' + line[22:]
```

### Boltz confidence JSON key
Boltz writes `complex_plddt` (not `plddt`) in confidence JSON files. Use:
```python
d.get('complex_plddt', float('nan'))   # pLDDT
d.get('ptm',           float('nan'))   # pTM
d.get('confidence_score', float('nan'))  # composite
```

### Parsing from Drive without copying PDBs
For parsing pLDDT only, read confidence JSONs directly from Drive — no `shutil.copytree` needed. Copying 2945 PDB files takes >20 min; JSONs parse in seconds:
```python
preds = list((DRIVE_RES / 'boltz_outputs_naf').rglob('predictions'))[0]
jsons = sorted((preds / name).glob('confidence_*.json'))
```

---

## References

- Dauparas et al. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. *Science* 378, 49–56.
- Kipnis et al. (2023). Improving Protein Expression, Stability, and Function with ProteinMPNN. *JACS* 145(29).
- van Kempen et al. (2024). Fast and accurate protein structure search with Foldseek. *Nature Biotechnology* 42, 243–246.
- Adaptyv Bio EGFR design competition. [bioRxiv 2025.04.17.648362](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2)
- Adaptyv Bio expression service. [start.adaptyvbio.com](https://start.adaptyvbio.com)
- Ginkgo Cloud Lab CFPS + HiBiT. [cloud.ginkgo.bio](https://cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit)

---

*Living document. Phase 2 active — native_ala_free MPNN + Boltz v3 complete (mean pLDDT 0.806), BioReason stress test pipeline ready, codon optimisation + platform decision next.*
