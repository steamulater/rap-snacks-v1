# Rap Snacks v1 — Phase 3 Lab Notebook

**Repo:** `steamulater/rap-snacks-v1`
**Phase 3 start:** 2026-04-09
**Last updated:** 2026-04-09
**Status:** Planning — pipeline automation + public platform design

---

## Phase 2 Completion Summary

Phase 2 closed on 2026-04-09 with the Ginkgo Cloud Lab order submission (Confirmation #DVZL22LN5, $936, 24 proteins). Full record in `LABNOTEBOOK_PHASE2.md`.

**What Phase 2 delivered:**
- 24 proteins selected from 12 candidate bars across 4 encoding conditions
- ProteinMPNN redesign on native_ala backbones (native_ala_free bucket) — best performers
- Boltz-2 validation: 5 diffusion models per design, pLDDT + pTM + RMSD characterised
- FoldSeek Phase 2: structural homolog search across pdb100, afdb-swissprot, mgnify
- PyMOL structure gallery: 25 renders, all 24 designs documented
- Platform decision: Ginkgo CFPS HiBiT at $39/protein vs Adaptyv $179/protein
- Full biophysical characterisation pathway mapped: screen → purification → DSF/MS/SEC → structural core
- Wet lab order submitted and confirmed

**AI session token usage (2026-04-09, single Claude Code session):**

This entire Phase 2 close-out — PyMOL diagnosis, design swap, gallery documentation, platform research, order submission — was completed in one Claude Code session.

| Source | Est. tokens |
|--------|------------|
| Images processed (28 PyMOL screenshots + desktop + order confirmation) | ~45,000 |
| Large file reads (lab notebooks, CSVs, PDB files, validation data) | ~35,000 |
| Bash/tool outputs (git, python, find, grep) | ~12,000 |
| Generated text (responses, notebook entries, commit messages) | ~20,000 |
| User messages + system context | ~8,000 |
| Subagents (Explore agent × 27 images, general-purpose Adaptyv research) | ~60,000 |
| **Total estimated** | **~180,000–210,000 tokens** |

Images were the single largest driver (~45K in main context). The 28 PyMOL screenshots alone account for ~30% of main-context token usage. Context compression occurred at least once during the session.

---

## Phase 3 Vision — Steamulater Games: Season 5

### Concept

**"Steamulater Games"** is a public-facing reality competition educational program built on the rap-snacks pipeline. Participants submit a lyric, phrase, or sentence of their choice. The pipeline converts it to a protein, folds it, designs the most foldable sequence, and returns a structural render with confidence metrics and biological context. Participants compete on whose lyric produces the most structurally interesting, novel, or foldable protein.

The name "Season 5" signals that this is part of a broader Steamulater Games series — each season uses a different scientific pipeline or cultural dataset as the competition substrate. This is Season 5.

**Why this works as a public program:**
- The lyric→protein pipeline is deterministic and reproducible — the same lyric always produces the same protein
- The structural output is visually striking and immediately shareable (PyMOL renders, pLDDT colour maps)
- The science is real — participants are doing actual protein design, not a simulation
- Cultural ownership: participants feel ownership over "their" protein because it came from their words
- Educational scaffold: each step (encoding, folding, design, validation) can be explained at varying depth

---

## Phase 3 Technical Roadmap

### Stage 1 — Pipeline Automation

Convert the current multi-step manual pipeline into a single callable function:

```
lyric (str)
    → concordance encoding       pipeline/convert.py
    → FASTA                      pipeline/convert.py
    → Boltz-2 fold               notebooks/boltz_predict.ipynb (Colab)
    → native_ala backbone        pipeline/convert.py (native_ala mode)
    → ProteinMPNN redesign       analysis/09d_proteinmpnn_native_ala_free.py
    → Boltz-2 validation         notebooks/boltz_validation_v4.ipynb
    → FoldSeek homolog search    pipeline/foldseek_lookup.py
    → PyMOL render               outputs/pymol/visualize.pml
    → result card                [new] pipeline/result_card.py
```

**Key automation targets:**
- End-to-end CLI: `python run_pipeline.py --lyric "your text here" --name "participant_id"`
- Colab notebook: single-cell execution for Boltz fold + MPNN + validation
- Result card generator: PNG card with structure render, pLDDT, pTM, FoldSeek result, lyric text, participant name
- Batch mode: process multiple lyrics in parallel (for competition rounds)

**Infrastructure decisions pending:**
- Where does Boltz run? (Colab Pro+ A100 currently — need to evaluate Modal, RunPod, or Replicate for public scale)
- MPNN: currently local — evaluate ESMFold API or Hugging Face Spaces for serverless
- FoldSeek: REST API rate limit (~27 submissions before 429) — implement queue + retry
- Storage: Drive currently — S3 or Supabase for public results

### Stage 2 — Public Platform

**Steamulater Games web app:**

| Component | Description |
|-----------|-------------|
| Submission form | Lyric/phrase input, participant name, optional song attribution |
| Pipeline runner | Async job queue — lyric submitted → result returned in ~15 min |
| Results gallery | Public gallery of all submitted proteins, searchable by participant/song/structure |
| Competition bracket | Voting + judging on structural novelty, confidence, FoldSeek result |
| Leaderboard | Ranked by composite score (pLDDT × pTM × FoldSeek novelty) |
| Result card | Shareable PNG per protein — designed for social media |

**Competition format (Season 5 rules — draft):**
- Open submission period: participants submit one lyric per round
- Each lyric produces one protein (concordance encoding, deterministic)
- Judging criteria: (1) Boltz confidence score, (2) FoldSeek novelty (zero hits = bonus), (3) community vote on cultural significance of lyric
- Rounds: weekly, 4 weeks → final
- Prize: winner's protein gets submitted to Ginkgo for wet lab expression

### Stage 3 — Educational Layer

Each result card and gallery page includes an educational scaffold:

- **Encoding explainer:** "The O in your lyric became Alanine because..."
- **Structure explainer:** "Your lyric folded into an alpha helix because it contains many..."
- **Novelty card:** "Your protein has no known homologs in the PDB — it may not exist in nature"
- **Biological context:** FoldSeek hit annotation — "Your lyric encodes a fold similar to [protein] found in [organism]"
- **Glossary:** pLDDT, pTM, alpha helix, FoldSeek — written for general audience

---

## Immediate Next Steps (Phase 3, Sprint 1)

| Priority | Task | Notes |
|----------|------|-------|
| 1 | Await Ginkgo HiBiT results | ~2 weeks |
| 2 | Build `run_pipeline.py` end-to-end CLI | Wraps existing scripts |
| 3 | Parameterise Colab notebook for any lyric input | Replace hardcoded bar list with `--lyric` arg |
| 4 | Build `result_card.py` — shareable PNG per protein | matplotlib + PyMOL render |
| 5 | Scoping doc: infrastructure for public scale | Modal vs RunPod vs Replicate for Boltz |
| 6 | Steamulater Games Season 5 — rules + brand doc | Separate file |

---

## Open Questions

1. **Encoding for non-English lyrics:** The concordance table is built on English letter frequencies. How does it handle Spanish, French, Yoruba, Mandarin romanisation? Does it need a language-specific frequency table per submission?
2. **BOJUXZ handling for arbitrary input:** The current softmax draw is seeded at 42 for reproducibility. For public submissions, seed = hash(lyric) to make each protein deterministic and unique per lyric
3. **Rate limits at scale:** Boltz is the bottleneck — ~15 min per protein on A100. At 100 submissions/week, need parallel Colab sessions or a cloud GPU provider
4. **IP / attribution:** Participants submit lyrics written by other artists. Need a clear terms of service on attribution and use
5. **Wet lab prize logistics:** Ginkgo CFPS for a winning lyric — who pays? Could be sponsorship, could be competition entry fee model

---

*Living document. Phase 3 begins 2026-04-09. Wet lab results expected ~2026-04-23.*
