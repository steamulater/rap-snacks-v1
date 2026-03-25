# Rap Snacks v1 — Phase 2 Lab Notebook

**Repo:** `steamulater/rap-snacks-v1`
**Started:** March 2026
**Goal:** Select, design, and submit rap-lyric-derived proteins to Ginkgo Bioworks cell-free expression — closing the loop from lyric → fold → real protein.

---

## Concept

Phase 1 established that concordance-mapped rap lyrics produce structurally novel proteins with no homologs in PDB, afdb-swissprot, or MGnify. Phase 2 takes the best candidates from that structural analysis, optimizes their sequences for wet-lab foldability using ProteinMPNN, and submits codon-optimized DNA to Ginkgo Bioworks' cell-free protein synthesis (CFPS) platform.

The core narrative: **the backbone encodes the bar. ProteinMPNN finds the most foldable sequence given that shape.** Every expressed protein traces directly back to a specific lyric.

---

## Target Platform — Ginkgo Bioworks CFPS

**System:** E. coli cell-free transcription/translation (T7 promoter)
**Input:** Codon-optimized linear DNA template
**Output:** Expressed protein from bacterial lysate, no cloning required
**Readout:** SDS-PAGE (expressed/not), optionally circular dichroism for fold confirmation

### Why CFPS over cell-based expression

- No transformation, no cloning — linear DNA is sufficient
- Fast turnaround, plate-scale throughput (automated)
- **Solubility-Enhanced Reagent Mix** — Ginkgo's kit specifically targets difficult/de novo proteins
- De novo designed proteins are notoriously hard to express in cells; CFPS tolerates them better
- Ginkgo's autonomous lab (in collaboration with OpenAI GPT-5, Feb 2026) runs exactly this kind of iterative design–express–learn loop at scale

### Sequence constraints for CFPS submission

| Parameter | Target |
|---|---|
| Length | 50–300 AA (our bars: 80–155 AA — ideal) |
| DNA template | Linear, under T7 promoter |
| Tag | N- or C-terminal His-tag for detection/purification |
| Codon usage | Optimized for *E. coli* K12 |
| Avoid | Long hydrophobic stretches (aggregation), >5 consecutive charged residues (insolubility), rare codons |

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

### Submission format to Ginkgo CFPS

Each submission entry:
```
bar_id       | bar_27
song         | Ganja Burn
lyric        | "You gotta have real skill..."
design_id    | bar_27_mpnn_003
aa_sequence  | CPMND...  (ProteinMPNN output)
dna_sequence | ATGCCGATG...  (codon-optimized)
length_aa    | 86
tag          | C-His6
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

2. **How many designs to submit?** Ginkgo CFPS is plate-based — 96-well is natural. Budget for ~48–96 sequences (6–12 bars × 5–10 designs each, plus controls).

3. **Controls to include:**
   - Native protein with known expression (positive control)
   - Scrambled sequence from same bar (negative control — same AA composition, random order)
   - Alanine condition sequence for same bar (compare BJOZXU strategy)

4. **IP / disclosure** — Ginkgo's research CFPS is a reagent kit (buyable). The actual submission pathway for a collaborative research project may involve their custom services team. Check biosecurity review requirements for novel sequences.

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
- Kipnis et al. (2023). Improving Protein Expression, Stability, and Function with ProteinMPNN. *JACS*, 145(29), 14audits–. [doi:10.1021/jacs.3c10941](https://pubs.acs.org/doi/10.1021/jacs.3c10941)
- Ginkgo Bioworks CFPS Kit (E. coli). [reagents.ginkgo.bio](https://reagents.ginkgo.bio/products/cell-free-protein-expression-e-coli)
- Ginkgo Bioworks + OpenAI GPT-5 autonomous lab (Feb 2026) — 36,000 CFPS experiments, 40% cost reduction. [PR Newswire](https://www.prnewswire.com/news-releases/ginkgo-bioworks-autonomous-laboratory-driven-by-openais-gpt-5-achieves-40-improvement-over-state-of-the-art-scientific-benchmark-302680619.html)
- Ginkgo Protein Engineering Services. [ginkgobioworks.com](https://www.ginkgobioworks.com/offerings/protein-services/)

---

*Living document. Phase 2 begins after Boltz-2 native_alanine folds return from Colab.*
