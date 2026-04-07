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

---

### Figure 28 — Concordance vs Native-Alanine pLDDT (Boltz-2)

![Figure 28](outputs/figures/fig28_conc_vs_na_plddt.png)

**Figure 28 | Boltz-2 pLDDT comparison between concordance and native-alanine conditions across Phase 2 candidate bars.** Each bar represents one of the 12 candidate sequences scored under two encoding conditions: *concordance* (full lyric-to-amino-acid mapping including BJOZXU non-standard characters) and *native_ala* (BJOZXU positions substituted with alanine, all other lyric-derived amino acids preserved). Bars are sorted by the concordance pLDDT score. Native-alanine outperforms concordance in the majority of cases, with a mean improvement of +0.10 pLDDT units. This suggests that the non-standard BJOZXU characters in the concordance sequence act as structural liabilities — their substitution with alanine, a structurally neutral and helix-favouring residue, consistently improves Boltz-2 confidence. The concordance condition is nonetheless retained in the submission panel as it is the purest expression of the lyric-to-protein mapping.

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

## Phase 2 — Step 7: FoldSeek Phase 2 (in progress)

**Script:** `analysis/10_foldseek_phase2.py`
**Date:** 2026-04-06
**Status:** Running (36 structures × 3 databases)

For each of the 12 bars, submitting the best PDB from each of 3 buckets (concordance, native_ala, free_design) to FoldSeek (pdb100, afdb-swissprot, mgnify_esm30). Results cached to Drive.

**Questions being answered:**
1. Do concordance, native_ala, and free_design structures hit the same protein families or different ones?
2. Does free_design find better-quality PDB hits (higher probability) than concordance?
3. Which bars remain structurally novel (zero hits) across all three buckets?

**Figure F** — `fig_foldseek_phase2.png` — FoldSeek hit comparison (pending, on Drive when complete)

---

## Phase 2 — Scientific Question: native_ala as backbone for MPNN

If native_ala folds more confidently (0.543 vs 0.441) and folds to a structurally distinct space from concordance, then using the native_ala Boltz backbone as the MPNN input may produce a new design series — starting from a more foldable backbone without sacrificing the cultural connection to the lyric.

### Three design strategies compared

| Strategy | Backbone | Fixed positions | Cultural fidelity | Design freedom |
|----------|----------|-----------------|-------------------|----------------|
| masked_BJOZXU (Run 1) | concordance Boltz | lyric AAs | high | low (5–18 positions) |
| free_design (Run 2) | concordance Boltz | none | none | maximum |
| **free_native_ala (planned)** | **native_ala Boltz** | **none** | **moderate** | **maximum** |

free_native_ala gives maximum design freedom on a more foldable backbone, while the native_ala backbone itself encodes the lyric sequence geometry. Planned as `analysis/09d_proteinmpnn_native_ala_free.py`.

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
| F | Drive: `fig_foldseek_phase2.png` | FoldSeek Phase 2 hits (pending) |

**Next figure number:** Fig 37

---

## Pending Work

| Priority | Task | Script |
|----------|------|--------|
| 1 | Parse FoldSeek Phase 2 results | `analysis/10_foldseek_phase2.py` |
| 2 | Build free_native_ala MPNN run | `analysis/09d_proteinmpnn_native_ala_free.py` |
| 3 | Boltz-2 validate free_native_ala designs | Update `notebooks/boltz_validation_v2.ipynb` |
| 4 | Codon optimisation | `analysis/10_codon_optimize.py` |
| 5 | Platform decision + submission | Adaptyv Bio |

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
Boltz writes 2-char chain IDs (`b11`, `A `) which break ProteinMPNN:
```python
line = line[:21] + 'A' + line[24:]
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

*Living document. Phase 2 active — Boltz-2 validation complete, FoldSeek Phase 2 running, free_native_ala MPNN run next.*
