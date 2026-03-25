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

Two full-service CFPS platforms are under consideration. Neither choice has been made — a quote from Ginkgo Cloud Lab is needed before deciding.

### Ginkgo Cloud Lab — CFPS + HiBiT

**URL:** [cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit](https://cloud.ginkgo.bio/protocols/cell-free-protein-expression-hibit)

**System:** E. coli CFPS with HiBiT bioluminescent readout. HiBiT is an 11-AA tag (Promega NanoLuc) fused to the protein — luminescence signal is proportional to expression level, quantitative and highly sensitive. Fully automated on Ginkgo's robotic RAC platform.

**What's included:** CFPS expression + quantitative HiBiT luminescence readout. Autonomous lab execution — no hands-on lab work required from submitter.

**Pricing:** Not publicly listed. Quote-based via their EstiMate tool — requires account creation at cloud.ginkgo.bio. Ginkgo's autonomous lab infrastructure was used to run 36,000 CFPS conditions in their GPT-5 collaboration (Feb 2026), achieving $422/g for sfGFP after optimization. Per-protein pricing at small scale is unknown.

**Scale:** Platform is optimised for high-throughput campaigns; may be overengineered for 8–12 proteins but designed to scale to 96-well and beyond.

**De novo protein track record:** Validated internally at scale, but public benchmarking on de novo designed proteins is limited in published literature.

---

### Adaptyv Bio — Expression Service

**URL:** [start.adaptyvbio.com](https://start.adaptyvbio.com)

**System:** Fully reconstituted E. coli TX/TL cell-free system (individually purified components — lower background than lysate-based systems, better suited for sensitive/synthetic designs).

**What's included:** Gene synthesis + CFPS expression + readout. Submit AA sequence via UI or API; they handle DNA synthesis.

**Pricing:**

| Service | Price | Turnaround | Includes |
|---|---|---|---|
| Expression | **$47/protein** | 2–3 weeks | Gene synthesis + expression + readout |
| Thermostability | $86/protein | 2–3 weeks | Expression + Tm measurement |
| Binding (BLI/SPR) | $79/protein | 3–4 weeks | Expression + binding characterization |

**Readout:** Expressed / low / not detected + relative yield + QC flags (insolubility, tag issues).

**De novo protein track record:** Extensive — 10,000+ designed proteins tested. Ran two public protein design competitions (EGFR binder challenge), validating 601 de novo designs including sub-100 nM binders. Over 30 companies and 10+ preprints using Adaptyv data in 2025.

---

### Side-by-side

| | **Ginkgo Cloud Lab (HiBiT)** | **Adaptyv Bio** |
|---|---|---|
| Pricing | Quote required — not public | $47/protein, transparent |
| Readout | Quantitative HiBiT luminescence | Expressed/low/not detected + yield |
| Gene synthesis included | Unknown | Yes |
| Submission | Protocol via EstiMate tool | Upload AA sequence via UI or API |
| Cell-free system | E. coli lysate + HiBiT tag | Reconstituted E. coli TX/TL |
| De novo protein validation | Internal (large scale) | Public (10,000+ proteins, competitions) |
| Turnaround | Unknown without quote | 2–3 weeks |
| Scale fit for 8–12 proteins | Uncertain — get quote | Well-suited |
| Scale fit for 100+ proteins | Likely competitive | Competitive |

**To decide:** Get a Ginkgo Cloud Lab quote for 10 proteins at [cloud.ginkgo.bio](https://cloud.ginkgo.bio) using their EstiMate tool, then compare directly to Adaptyv's $470 total (10 × $47).

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

ProteinMPNN runs locally or on Colab:
```bash
# Install
pip install protein-mpnn

# Design with fixed positions
python protein_mpnn_run.py \
    --pdb_path outputs/boltz_outputs/.../bar_27_model_0_fixed.pdb \
    --out_folder outputs/proteinmpnn/bar_27/ \
    --num_seq_per_target 50 \
    --sampling_temp 0.1 \
    --fixed_positions_jsonl data/phase2_fixed_positions.jsonl
```

Fixed positions JSON maps each bar's lyric-AA positions (non-BJOZXU) to be held constant.

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

## Open Questions (Phase 2)

1. **Which positions to fix in ProteinMPNN?** The mutation mask from Phase 1 tracks exactly which positions came from standard AA vs BJOZXU draws. This is the exact input needed for `--fixed_positions_jsonl`.

2. **Which platform?** Get Ginkgo Cloud Lab quote via EstiMate before committing. Compare against Adaptyv's $47/protein. See Platform Comparison section above.

3. **How many designs to submit?** Both platforms support 96-well plate scale. Budget for ~48–96 sequences (6–12 bars × 5–10 designs each, plus controls).

4. **Controls to include:**
   - Native protein with known expression (positive control)
   - Scrambled sequence from same bar (negative control — same AA composition, random order)
   - Alanine condition sequence for same bar (compare BJOZXU strategy)

5. **IP / disclosure** — Check biosecurity review requirements for novel sequences on whichever platform is chosen. Adaptyv's terms are publicly documented; Ginkgo's may involve a custom services agreement.

5. **What constitutes "success"?** A band on SDS-PAGE at the correct MW = expression. CD spectrum with secondary structure features matching Boltz predictions = folding. Either result is publishable — expression failure of a lyric-derived sequence is also a finding.

---

## Repo Structure (Phase 2 additions)

```
analysis/
├── 08_candidate_selection.py       <- filter 85 bars → shortlist
├── 09_proteinmpnn_design.py        <- run ProteinMPNN, self-consistency filter
└── 10_codon_optimize.py            <- codon-optimize + build submission CSV
data/
├── phase2_candidates.csv           <- ranked shortlist (output of step 1)
└── phase2_fixed_positions.jsonl    <- ProteinMPNN fixed positions per bar
outputs/
├── proteinmpnn/                    <- designed sequences per bar
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

*Living document. Phase 2 begins after Boltz-2 native_alanine folds return from Colab.*
