# Rap Snacks — Protein Bars

Convert culturally iconic rap lyrics into real proteins. Fold them. Design them. Express them in a cell-free system. Find out if Nicki Minaj's words exist in nature.

**Upstream input:** Iconicity score from organic social signal (#NeurodivergentRapLines, X/Twitter)
**Downstream output:** Wet lab–confirmed expressed proteins, each traceable to a specific lyric

---

## Status

| Phase | Status |
|-------|--------|
| Phase 1 — Data collection + concordance folding | ✅ Complete |
| Phase 2 — ProteinMPNN design + Boltz validation + selection | ✅ Complete |
| Phase 2 — Wet lab expression (Ginkgo Cloud Lab) | ⏳ In progress — results ~2026-04-23 |
| Phase 3 — Pipeline automation + Steamulater Games Season 5 | 🔜 Planning |

**24 proteins submitted to Ginkgo Cloud Lab** (Confirmation #DVZL22LN5, 2026-04-09)
Cell-free protein synthesis + HiBiT luminescence detection, $936 total.

---

## What This Project Does

```
Rap lyric
    ↓  concordance encoding — letter frequency → amino acid mapping
Protein sequence
    ↓  Boltz-2 — structure prediction (5 diffusion samples, A100)
3D structure
    ↓  ProteinMPNN — redesign on native_ala backbone (native_ala_free)
Optimised sequence
    ↓  Boltz-2 validation — confirm redesigned sequence folds
Validated design
    ↓  FoldSeek — structural homolog search (pdb100, afdb-swissprot, mgnify)
Biological context
    ↓  Ginkgo Cloud Lab — CFPS + HiBiT expression screen
Real protein
```

---

## Key Findings (Phase 1 + 2)

- **r(iconicity, structural confidence) = 0.051** — cultural resonance and protein foldability are orthogonal. A bar going viral tells you nothing about whether it folds
- **13 bars** have zero structural homologs in PDB, AlphaFold DB, and MGnify — structurally novel proteins that may not exist in nature
- **native_ala_free MPNN redesign** outperforms concordance encoding: mean Boltz pLDDT 0.806 vs 0.441
- **free_design RMSD vs concordance backbone = 7.13Å < concordance noise floor 9.58Å** — MPNN encodes backbone geometry better than re-folding the original sequence
- Notable hit: bar_96 "Big Foot but you still a small fry" → sulfur carrier protein from *Thermotoga maritima* (prob=0.956, extremophile living at 80°C hydrothermal vents)

---

## Selected Proteins (24 designs submitted for wet lab)

| Group | n | Description |
|-------|---|-------------|
| A | 11 | Best ProteinMPNN design per bar (native_ala_free) |
| B | 4 | Free design vs native_ala_free comparison pairs |
| C | 2 | Replicates — top performers with second MPNN seed |
| D | 5 | Raw lyric seeds — concordance + native_ala encodings |
| E | 1 | Negative control — bar_6_scrambled_02 (composition-matched) |
| + | 1 | Positive control — 6E5C (Baker lab de novo dsβH, NMR-verified) |

Full selection with sequences, metrics, and PyMOL renders: `outputs/selected_proteins.csv` and `LABNOTEBOOK_PHASE2.md`.

---

## Encoding: Five Modes

2×2 factorial design across standard letter mapping × BOJUXZ substitution:

| Mode | Standard 20 AA | BOJUXZ positions | Isolates |
|------|---------------|-----------------|---------|
| `concordance` | Freq-rank concordance | Softmax peaked draw | Canonical encoding |
| `alanine` | Freq-rank concordance | → A | Effect of BOJUXZ softmax |
| `random` | Freq-rank concordance | Uniform random | Is concordance better than chance? |
| `native` | AA letter pass-through | Softmax peaked draw | Effect of freq remapping |
| `native_alanine` | AA letter pass-through | → A | BOJUXZ strategy without remapping |

BOJUXZ = B, O, J, U, X, Z — letters with no standard single-letter amino acid code. Substituted probabilistically (softmax, λ=2.0, seed=42).

---

## Pipeline Scripts

```
pipeline/
├── collect.py              Stage 1 — tweet collection (#NeurodivergentRapLines)
├── process.py              Stage 2 — bar extraction + iconicity scoring
├── remerge.py              Stage 3 — fuzzy deduplication
├── verify.py               Stage 4 — Genius API attribution
├── convert.py              Stage 5 — lyric → FASTA (5 modes) + bar_index_snapshot.json
└── foldseek_lookup.py      FoldSeek REST API — structural homolog search

analysis/
├── 09d_proteinmpnn_native_ala_free.py   ProteinMPNN redesign on native_ala backbones
└── 10_foldseek_phase2.py                FoldSeek Phase 2 — 36 structures × 3 databases

notebooks/
├── boltz_validation_v4_scrambled_na.ipynb   Boltz-2 validation (Colab A100)
└── boltz_validation_v2.ipynb                Phase 1 Boltz run

outputs/
├── selected_proteins.csv       24 selected designs with sequences + metrics
├── pymol/                      PDB files + PyMOL session files per bar
│   └── screenshots/            25 PyMOL renders + Ginkgo order confirmation
├── boltz_validation/           Boltz-2 confidence scores + RMSD
└── figures/                    95+ analysis figures
```

---

## Data

| File | Description |
|------|-------------|
| `data/aggregated_lines_v2_frozen.csv` | Canonical input — 539 bars, 474 verified. **Immutable.** |
| `data/aggregated_lines_v2_enriched.csv` | Master CSV — all columns including sequences + structural metrics |
| `data/bar_index_snapshot.json` | `bar_N → lyric + all 5 sequences`. Authoritative attribution record |
| `data/phase2_candidates.csv` | 12 bars selected for Phase 2 MPNN design |
| `outputs/selected_proteins.csv` | 24 final designs submitted for wet lab |
| `Selected_Ginkgo_Cloud_Lab_Input_Template_CFPS (2).xlsx` | Exact file submitted to Ginkgo 2026-04-09 |

---

## Wet Lab Pipeline (Ginkgo Cloud Lab)

**Screen (submitted):** Cell-Free Protein Synthesis + HiBiT luminescence, $39/protein
**Scale-up (planned, for hit expressors):** E. coli BL21, His-tag purification, $85–160/protein
**Characterisation (planned):** DSF thermostability + intact mass spec + SEC-LC-MS, ~$92/sample
**Structural (if warranted):** CD, SAXS, NMR, X-ray, cryo-EM at university core facility

Full decision rationale and pricing in `LABNOTEBOOK_PHASE2.md`.

---

## Phase 3 — Steamulater Games: Season 5

The next phase automates this pipeline into a public platform where anyone can submit a lyric and receive their own protein. Participants compete on whose lyric produces the most structurally interesting, novel, or foldable design. Winner's protein gets expressed in a real cell-free system.

See `LABNOTEBOOK_PHASE3.md` for full vision, technical roadmap, and competition format.

---

## Lab Notebooks

| File | Contents |
|------|----------|
| `LABNOTEBOOK.md` | Phase 1 — data collection, encoding, Boltz pilot, FoldSeek, integrity audit |
| `LABNOTEBOOK_PHASE2.md` | Phase 2 — ProteinMPNN, Boltz validation, design selection, PyMOL gallery, wet lab |
| `LABNOTEBOOK_PHASE3.md` | Phase 3 — pipeline automation, Steamulater Games Season 5 |

---

## Environment

- **Local:** MacBook Air M3, Python 3.x, no GPU required for pipeline stages
- **Cloud:** Google Colab Pro+ (A100-SXM4-80GB) for Boltz-2 folding
- **APIs:** X API v2 (tweet collection), Genius API (attribution), FoldSeek REST, ESMFold API
- **Wet lab:** Ginkgo Cloud Lab (CFPS + HiBiT), Adaptyv Bio (purification, if needed)

```bash
pip install -r requirements.txt
```
