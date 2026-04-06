# Rap Snacks v1 — Phase 2 Lab Notebook

**Repo:** `steamulater/rap-snacks-v1`
**Started:** March 2026
**Goal:** Select, design, and submit rap-lyric-derived proteins to a cell-free expression platform — closing the loop from lyric → fold → real protein.

---

## Concept

Phase 1 established that concordance-mapped rap lyrics produce structurally novel proteins with no homologs in PDB, afdb-swissprot, or MGnify. Phase 2 takes the best candidates from that structural analysis, optimizes their sequences for wet-lab foldability using ProteinMPNN, and submits codon-optimized DNA to a cell-free protein synthesis (CFPS) platform.

The core narrative: **the backbone encodes the bar. ProteinMPNN finds the most foldable sequence given that shape.** Every expressed protein traces directly back to a specific lyric.

---

## Target Platform — Platform Comparison (TBD)

### Why CFPS over cell-based expression

- No transformation, no cloning — linear DNA is sufficient
- Fast turnaround, plate-scale throughput (automated)
- De novo designed proteins are notoriously hard to express in cells; CFPS tolerates them better
- Both candidate platforms use E. coli cell-free TX/TL systems, suitable for our 80–155 AA range

### Sequence constraints (both platforms)

| Parameter | Target |
|---|---|
| Length | 50–300 AA (our bars: 80–155 AA — ideal) |
| DNA template | Linear or plasmid, under T7 promoter |
| Tag | N- or C-terminal tag for detection (His6 or HiBiT depending on platform) |
| Codon usage | Optimized for *E. coli* K12 |
| Avoid | Long hydrophobic stretches (aggregation), >5 consecutive charged residues (insolubility), rare codons |

---

## Platform Comparison — Ginkgo Cloud Lab vs Adaptyv Bio

Two full-service CFPS platforms are under consideration. Quotes obtained for 24 proteins (2026-03-25). Screenshots: `Screen_adaptyv.png`, `Sceen_ginkgo.png` (repo root).

### Ginkgo Cloud Lab — CFPS + HiBiT

**URL:** [cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit](https://cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit)

**System:** E. coli CFPS with HiBiT bioluminescent readout. HiBiT is an 11-AA tag (Promega NanoLuc) fused to the protein — luminescence signal is proportional to expression level, quantitative and highly sensitive. Fully automated on Ginkgo's robotic RAC platform.

**What's included:** CFPS expression + quantitative HiBiT luminescence readout. Autonomous lab execution — no hands-on lab work required from submitter.

**Pricing (actual quote, 24 proteins × 3 replicates = 72 samples):**

| Proteins | Replicates | Total Samples | Estimated Total | Per Protein |
|---|---|---|---|---|
| 24 | 3 | 72 | **$2,808** | ~$117 |

Screenshot: `Sceen_ginkgo.png`

**Scale:** 96-well plate scale. 3 replicates per protein included in the quoted price.

**De novo protein track record:** Validated internally at scale, but public benchmarking on de novo designed proteins is limited in published literature.

---

### Adaptyv Bio — Expression Service

**URL:** [start.adaptyvbio.com](https://start.adaptyvbio.com)

**System:** Fully reconstituted E. coli TX/TL cell-free system (individually purified components — lower background than lysate-based systems, better suited for sensitive/synthetic designs).

**What's included:** Gene synthesis + CFPS expression + readout. Submit AA sequence via UI or API; they handle DNA synthesis.

**Pricing (actual quote, 24 proteins, Standard Delivery):**

| Proteins | Delivery | Per Protein | Total | Condition |
|---|---|---|---|---|
| 24 | Standard (2–3 weeks) | **$95** | **$2,285** | Open-source discount applied (20%) — requires releasing results on Proteinbase |
| 24 | Economy (8+ weeks, ~15% off) | ~$81 | ~$1,940 | Economy tier + open-source discount |

Screenshot: `Screen_adaptyv.png`

> **Note:** Earlier estimate of $47/protein was incorrect — actual price is $95/protein with the open-source/Proteinbase discount applied. Without the discount the price is ~$119/protein. The open-source discount requires agreeing to publish results on [Proteinbase](https://proteinbase.io) — a condition that is compatible with this project's publication goals.

**Readout:** Expressed / low / not detected + relative yield + QC flags (insolubility, tag issues).

**De novo protein track record:** Extensive — 10,000+ designed proteins tested. Ran two public protein design competitions (EGFR binder challenge), validating 601 de novo designs including sub-100 nM binders. Over 30 companies and 10+ preprints using Adaptyv data in 2025.

---

### Side-by-side (24 proteins)

| | **Ginkgo Cloud Lab (HiBiT)** | **Adaptyv Bio** |
|---|---|---|
| Per protein | ~$117 | ~$95 (with open-source discount) |
| Total (24 proteins) | **$2,808** | **$2,285** |
| Replicates included | 3 per protein | 1 (single expression run) |
| Readout | Quantitative HiBiT luminescence | Expressed/low/not detected + yield |
| Gene synthesis included | Yes (DNA template upload) | Yes |
| Submission | EstiMate tool → upload DNA | Upload AA sequence via UI or API |
| Cell-free system | E. coli lysate + HiBiT tag | Reconstituted E. coli TX/TL |
| De novo protein validation | Internal (large scale) | Public (10,000+ proteins, competitions) |
| Turnaround | Not specified | 2–3 weeks (standard) |
| Open-source condition | None | Must publish on Proteinbase for discount |

**Summary:** Adaptyv is ~$500 cheaper for 24 proteins. Ginkgo includes 3 replicates, which provides technical reproducibility at the cost of price. Adaptyv's reconstituted TX/TL system has stronger published validation for de novo designs. **No platform decision made yet.**

---

## Step 1 — Candidate Selection

**Script (to build):** `analysis/08_candidate_selection.py`
**Input:** `outputs/boltz/boltz_summary.csv`, `outputs/pairwise/struct_rmsd.csv`, `outputs/pairwise/seq_identity.csv`, `data/bar_index_snapshot.json`
**Output:** `data/phase2_candidates.csv` — ranked shortlist of 8–12 bars

### Selection filters

| Filter | Threshold | Rationale |
|---|---|---|
| Boltz pTM (mean) | ≥ 0.35 | Model committed to a fold |
| Boltz pLDDT (mean) | ≥ 0.40 | Per-residue confidence |
| Sequence length | 80–150 AA | CFPS sweet spot |
| Structural novelty | 10th-pct RMSD ≥ 11 Å | No close structural neighbors — no point submitting known folds |
| Sequence novelty | Max pairwise identity ≤ 0.30 | All bars pass this (max = 30%) |
| Hydrophobic run length | ≤ 6 consecutive hydrophobic AA | Key solubility predictor for E. coli CFPS |
| No stop codons in lyric-derived sequence | — | Sanity check |

### Expected candidates

Based on Phase 1 results, the pool will likely be:
- **bar_27** (Ganja Burn) — only `confident_protein_like` bar, pLDDT=0.665, pTM=0.487. Top priority.
- The 12 `uncertain_protein_like` bars (pTM ≥ 0.4), filtered further by hydrophobicity and length
- Cross-reference with Fig 23 (novelty scatter) — bars in the top-right quadrant (novel in both sequence and structure axes)

Target shortlist: **8–12 bars**

---

## Step 2 — Soft Sequence Design with ProteinMPNN

**Tool:** ProteinMPNN (Dauparas et al., *Science* 2022)
**Script (to build):** `analysis/09_proteinmpnn_design.py`
**Input:** Boltz-2 concordance PDB (best model) for each candidate bar
**Output:** `outputs/proteinmpnn/{bar_id}/designed_seqs.fasta` — 50 sequences per bar

### Why ProteinMPNN over BoltzGen / BoltzDesign

| Tool | Use case | Fit for this project |
|---|---|---|
| **ProteinMPNN** | Fixed backbone → redesign sequence | ✓ Perfect — keep the lyric backbone, improve BJOZXU positions |
| BoltzGen | Unconditional sequence generation | ✗ Doesn't respect our backbone |
| BoltzDesign1 | All-atom binder design | ✗ Overkill for monomeric expression, less wet-lab validated |
| RFdiffusion | Backbone generation | ✗ We already have backbones |

ProteinMPNN has direct experimental validation for improving *E. coli* expression, solubility, and melting temperature (Kipnis et al., *JACS* 2023). It is the most reliable choice for a first wet-lab submission.

### Design strategy — partial sequence fixation

```
Boltz concordance PDB (best model)
        |
        ↓
ProteinMPNN
  - FIXED positions:   standard AA residues from the lyric mapping
                       (the lyric is preserved in the backbone)
  - DESIGNABLE positions: BJOZXU-derived positions (filled with
                           softmax draws in concordance condition)
                           → ProteinMPNN picks the most foldable AA here
        |
        ↓
50 candidate sequences per bar
        |
        ↓
Self-consistency filter (ESMFold each designed sequence)
  → Keep sequences where ESMFold TM-score vs original backbone ≥ 0.70
        |
        ↓
~5–10 sequences per bar pass to codon optimization
```

The BJOZXU positions are the design space — they were never real amino acids to begin with. Letting ProteinMPNN choose at those positions is biologically sound and preserves the lyric at every position that was already a canonical AA.

### Self-consistency filtering

After ProteinMPNN generates 50 sequences per bar:
1. Fold each with ESMFold (fast, free API)
2. Compute TM-score between designed fold and original Boltz backbone
3. Keep sequences with TM-score ≥ 0.70 (fold matches backbone)
4. Secondary filter: predicted solubility (CamSol or IUPred2 disorder score < 0.5)

### Running ProteinMPNN

ProteinMPNN runs locally from the cloned repo at `/Users/tamukamartin/Desktop/adaptyv_competiton/ProteinMPNN/`.
Script `analysis/09_proteinmpnn_design.py` handles the full pipeline end-to-end:

```bash
# Full pipeline: MPNN design + ESMFold self-consistency filter
python analysis/09_proteinmpnn_design.py --top-n 12 --n-seqs 50 --temp 0.1 --esm-plddt-min 0.35
```

The script:
1. Fixes Boltz 3-char chain IDs in PDBs before passing to MPNN
2. Builds `data/phase2_fixed_positions.jsonl` (lyric AA positions fixed per bar)
3. Runs `protein_mpnn_run.py` via subprocess for each bar
4. Folds each design with ESMFold API (pLDDT returned in 0–1 scale directly)
5. Keeps designs with ESMFold pLDDT ≥ 0.35 → `outputs/proteinmpnn/filtered_seqs.fasta`

> **ESMFold pLDDT note:** The ESMFold API returns pLDDT in 0–1 scale in the B-factor column (not 0–100). Do not divide by 100.

Fixed positions JSONL maps each bar's lyric-AA positions (non-BJOZXU) to be held constant.

---

### ProteinMPNN Run 1 — masked_BJOZXU, Temperature 0.1 (2026-04-01)

**Design strategy:** BJOZXU positions free · all lyric-AA positions fixed
**Run:** top-12 candidates × 50 sequences × temperature 0.1 · ESMFold self-consistency filter pLDDT ≥ 0.35
**Output:** `outputs/proteinmpnn/filtered_results.csv` (612 rows), `outputs/proteinmpnn/filtered_seqs.fasta` (393 sequences)
**Label:** `masked_BJOZXU` — lyric amino acids are locked; only the non-amino-acid lyric characters (B, J, O, Z, X, U) are redesigned

**Masking rationale:** BJOZXU characters have no valid amino acid translation. These positions are the natural design space — fixing lyric AAs preserves the semantic content of the bar while letting MPNN choose the most foldable residue at each non-AA position. This is the most conservative design strategy.

**Constraint analysis:**

| Bar | Seq len | BJOZXU (designable) | Lyric AA (fixed) | Unique designs |
|-----|---------|---------------------|------------------|----------------|
| bar_6  | 88 | 5  | 83 | 5/51  |
| bar_9  | 82 | 8  | 74 | 20/51 |
| bar_32 | 94 | 17 | 77 | 42/51 |
| bar_17 | 121| 18 | 103| 49/51 |

At temp=0.1 with few designable positions, MPNN converges to near-identical sequences (bar_6: only 5 unique out of 51). Mutation rate from concordance: 5–16% depending on BJOZXU density.

| Bar | Pass / Total | Mean pLDDT | Status |
|-----|-------------|------------|--------|
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

**Total passing:** 393 / 612 (64%)

---

### What "Dropout" Means

A bar **drops out** when 0 of 50 ProteinMPNN designs pass the self-consistency pLDDT threshold (≥ 0.35). There are two possible causes:

**Cause 1 — The backbone is disordered and unfixable.**
If the underlying Boltz prediction for a bar has low confidence (low pTM, no compact core), then *no* sequence can fold to that backbone geometry, because the backbone itself doesn't encode a stable fold. MPNN will generate plausible-looking sequences, but ESMFold folding them independently always returns disordered (pLDDT < 0.3). This is the bar failing on structural grounds — a valid scientific result meaning this lyric doesn't encode a foldable shape.

**Cause 2 — The backbone is foldable but temperature is too low.**
ProteinMPNN's sampling temperature controls sequence diversity. At `temp=0.1` (near-greedy), MPNN converges very tightly on the single most probable sequence — 50 designs are barely distinguishable from each other. If the optimal sequence is in a region of parameter space MPNN doesn't explore at low temperature, you miss it entirely. Raising to 0.2–0.3 explores a broader neighbourhood and can rescue bars that appear to drop out at 0.1.

**How to distinguish them:** Rerun the dropout bars at temperatures 0.2 and 0.3. If they recover (some designs pass), temperature was the issue. If all three temps give 0/50, the backbone is the problem and the bar is dropped from the submission.

---

### ProteinMPNN Run 2 — free_design (planned)

**Design strategy:** No fixed positions — MPNN optimizes the full sequence for the Boltz backbone
**Script:** `analysis/09_proteinmpnn_free_design.py`
**Output:** `outputs/proteinmpnn_free/`
**Label:** `free_design`

**Rationale:** The masked_BJOZXU run preserves lyric identity but is heavily constrained — bars with few BJOZXU positions have almost no design freedom (bar_6: only 5 designable positions). The free design run asks: *what is the most foldable sequence for this backbone, with no constraints?* This serves two purposes:

1. **Upper bound on foldability** — if free_design also scores low pLDDT, the backbone itself is the problem (not the masking constraint). Confirms backbone-level dropout for bar_0, bar_8, bar_46.
2. **Wet-lab comparison** — submitting both masked_BJOZXU and free_design for the same bar directly tests whether preserving the lyric AAs costs foldability. If free_design expresses and masked_BJOZXU doesn't, the lyric is the bottleneck.

**Run parameters:**
```bash
python analysis/09_proteinmpnn_free_design.py --top-n 12 --n-seqs 50 --temp 0.1
```

No `--fixed_positions_jsonl` — all positions designable.

**Expected output per bar:** 50 sequences with much higher sequence diversity than masked_BJOZXU (full 88–128 AA design space vs 5–18 positions). Mutation rate from concordance expected: 30–60%.

---

## Step 3 — Codon Optimization and Submission

**Script (to build):** `analysis/10_codon_optimize.py`
**Tool:** Benchling Codon Optimization API or Ginkgo's internal codon optimizer
**Input:** Selected designed sequences (FASTA)
**Output:** `outputs/submission/phase2_submission.csv` — codon-optimized DNA sequences + metadata

### Codon optimization for E. coli K12

Standard codon optimization:
- Replace rare codons (usage < 10% in *E. coli* K12) with common synonyms
- Avoid internal T7 terminator sequences in the DNA
- Avoid hairpin structures in the first 30 nt of the ORF (affects ribosome binding)
- Add: T7 promoter + RBS upstream, His6-tag (C-terminal preferred for de novo proteins), T7 terminator downstream

### Submission format

Each submission entry (platform-agnostic):
```
bar_id       | bar_27
song         | Ganja Burn
lyric        | "You gotta have real skill..."
design_id    | bar_27_mpnn_003
aa_sequence  | CPMND...  (ProteinMPNN output)
dna_sequence | ATGCCGATG...  (codon-optimized, only needed for Ginkgo)
length_aa    | 86
tag          | C-His6 (Adaptyv) or C-HiBiT (Ginkgo Cloud Lab)
notes        | BJOZXU positions redesigned; lyric positions fixed
```

### Readout plan

| Assay | What it tells us |
|---|---|
| SDS-PAGE (denaturing) | Expressed or not, correct molecular weight |
| Native PAGE / SEC | Folded monomer vs aggregate |
| Circular dichroism (CD) | Secondary structure content (α/β match to Boltz prediction?) |
| Thermal melt (nanoDSF) | Melting temperature — higher = better designed |

---

## Submission Scale

**Target: 48-protein plate** (Adaptyv Bio, Standard Delivery)

| Slot allocation | Count |
|---|---|
| 12 candidate bars × 4 sequences (concordance, native_ala, MPNN-1, MPNN-2) | 48 |
| 1 positive control (sfGFP or known-expressing protein) | 1 |
| Contingency / reserve | ~7 |
| **Total** | **48** |

**Cost: ~$4,560** (48 × $95, Adaptyv open-source discount — requires publishing results on Proteinbase)

---

## Controls

### Positive control (wet lab)
- 1 × known-expressing protein (sfGFP or equivalent) — confirms the platform ran correctly
- Occupies 1 well; essential

### Scrambled sequence (in silico only — no wet lab slot needed)

**Script:** `analysis/11_scrambled_control.py` — generates 3 shuffled variants of both concordance and native_ala sequences per bar, folds all with ESMFold API.

**Results (2026-04-01, full 12-bar candidate set, no 403 dropouts):**

| Bucket | n | mean pLDDT | sd | <0.3 |
|--------|---|------------|-----|------|
| concordance (original) | 37 | 0.324 | 0.041 | 12/37 |
| native_ala (original) | 9 | 0.370 | 0.060 | 2/9 |
| scrambled_concordance | 110 | 0.333 | 0.059 | 32/110 |
| scrambled_native_ala | 111 | 0.397 | 0.099 | 10/111 |

**Finding:** ESMFold does *not* discriminate between originals and scrambled sequences at this length. Scrambled concordance (0.333) ≈ concordance (0.324); scrambled native_ala (0.397) actually *exceeds* native_ala (0.370). This is consistent with known ESMFold limitations — it relies on PLM embeddings from ESM-2, which are partially order-insensitive for short sequences without strong long-range contacts.

**Implication:** The scrambled control does not work as an in silico negative control under ESMFold. The structural novelty argument rests on Boltz-2 + FoldSeek (zero homologs across PDB/afdb/MGnify) — these are the authoritative results for structural claims.

**Plan:** Re-run scrambled control through ColabFold (AlphaFold2 + MSA) using `notebooks/colabfold_validation.ipynb` Section B. AF2 is order-sensitive and MSA-aware — more likely to show a signal. If ColabFold shows scrambled < original, that strengthens the structural case. Either result is informative.

- If reviewers push back, wet-lab scrambled can be added in a follow-up run ($95/slot at Adaptyv).

### Native_ala (wet lab — already in main submission)
- The native_ala sequence for each bar is already one of the 4 submission slots per bar — serves as the deterministic alanine-substitution baseline

---

## Open Questions (Phase 2)

1. **Which positions to fix in ProteinMPNN?** The mutation mask from Phase 1 tracks exactly which positions came from standard AA vs BJOZXU draws. This is the exact input needed for `--fixed_positions_jsonl`.

2. **Which platform?** Adaptyv Bio favoured — $4,560 for 48 proteins vs Ginkgo $5,616. Ginkgo includes triplicates but that is more useful after you know something expresses. Decision pending native_ala fold results.

3. **IP / disclosure** — Check biosecurity review requirements for novel sequences on whichever platform is chosen. Adaptyv's terms publicly documented (Proteinbase release required for discount); Ginkgo may involve a custom services agreement.

4. **What constitutes "success"?** A band on SDS-PAGE at the correct MW = expression. CD spectrum with secondary structure features matching Boltz predictions = folding. Either result is publishable — expression failure of a lyric-derived sequence is also a finding.

---

## Repo Structure (Phase 2 additions)

```
analysis/
├── 08_candidate_selection.py       <- filter 85 bars → shortlist
├── 09_proteinmpnn_design.py        <- run ProteinMPNN, self-consistency filter
├── 10_codon_optimize.py            <- codon-optimize + build submission CSV
└── 11_scrambled_control.py         <- in silico scrambled control (ESMFold + pLDDT comparison)
data/
├── phase2_candidates.csv           <- ranked shortlist (output of step 1)
└── phase2_fixed_positions.jsonl    <- ProteinMPNN fixed positions per bar
outputs/
├── proteinmpnn/                    <- designed sequences per bar
├── scrambled/                      <- scrambled seqs + ESMFold results (in silico control)
└── submission/                     <- codon-optimized DNA, submission-ready
```

---

## References

- Dauparas et al. (2022). Robust deep learning–based protein sequence design using ProteinMPNN. *Science*, 378, 49–56. [doi:10.1126/science.add2187](https://www.science.org/doi/10.1126/science.add2187)
- Kipnis et al. (2023). Improving Protein Expression, Stability, and Function with ProteinMPNN. *JACS*, 145(29). [doi:10.1021/jacs.3c10941](https://pubs.acs.org/doi/10.1021/jacs.3c10941)
- Ginkgo Cloud Lab — CFPS + HiBiT protocol. [cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit](https://cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit)
- Ginkgo Bioworks + OpenAI GPT-5 autonomous lab (Feb 2026) — 36,000 CFPS experiments, 40% cost reduction. [PR Newswire](https://www.prnewswire.com/news-releases/ginkgo-bioworks-autonomous-laboratory-driven-by-openais-gpt-5-achieves-40-improvement-over-state-of-the-art-scientific-benchmark-302680619.html)
- Ginkgo HiBiT high-throughput screening application note. [ginkgo.bio](https://www.ginkgo.bio/resources/application-notes/high-throughput-screening-hibit-proteins)
- Adaptyv Bio expression service. [start.adaptyvbio.com](https://start.adaptyvbio.com)
- Adaptyv Bio cell-free expression technology. [docs.adaptyvbio.com](https://docs.adaptyvbio.com/wiki/technologies/cell-free-protein-expression)
- Adaptyv EGFR design competition — crowdsourced de novo protein validation. [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2)

---

*Living document. Phase 2 active — Boltz-2 native_ala complete, candidates selected, ProteinMPNN design next.*
