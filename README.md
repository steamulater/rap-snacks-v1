# Rap Snacks v1

Convert culturally iconic rap lyrics into protein sequences, fold them, and find their structural homologs in nature.

**Upstream input:** Iconicity score from organic social signal (#NeurodivergentRapLines, X/Twitter)
**Downstream output:** Folded protein structures with lyric-origin mutation masks

---

## Quick Start

```bash
# 1. Audit the data — review 20 random rows before proceeding
make audit

# 2. Convert (all 5 modes, 80-300 AA filter)
make convert

# 3. Fold on Colab (see Makefile for Boltz-2 instructions)

# 4. ESMFold pilot — top 25 bars
make esm-pilot

# 5. Parse Boltz outputs (after downloading from Colab)
make parse-boltz

# 6. FoldSeek structural search
make foldseek

# 7. Enrich master CSV with all metrics
make enrich
```

---

## Pipeline

```
data/aggregated_lines_v2_frozen.csv    <- immutable input, never modified
         |
00_audit.py        validate + 20 random row spot-check
         |
01_convert.py      lyric -> FASTA (5 modes) + bar_index_snapshot.json
         |                   + enriched CSV with all sequences
     /fastas/
bars_v2_{mode}.fasta × 5 modes
         |
[Boltz-2 on Colab]  concordance condition, PDB output
         |
03_parse_boltz.py  extract pLDDT, pTM, confidence -> boltz_confidence_scores.csv
         |
02_esm_fold.py     ESMFold API, all 5 conditions, top 25 pilot
         |
04_foldseek.py     structural homolog search, pTM>=0.4 bars
         |
05_enrich_csv.py   join all metrics -> aggregated_lines_v2_enriched.csv
```

---

## V2 Improvements (vs v1 pilot)

| Issue | V1 | V2 fix |
|---|---|---|
| Lyric attribution drift | CSV re-sorted after conversion, bar_N mapping broke | `bar_index_snapshot.json` written at conversion time — never depends on CSV ordering |
| FASTA/CSV agreement | No check | `len(lyric_cleaned) == len(fasta_seq)` asserted for all bars, all modes |
| Frozen CSV modified | No guard | SHA-256 hash written on first run; subsequent runs refuse if hash changed |
| Length confound | 10-225 AA (uncontrolled) | 80-300 AA filter — tight cohort, minimal length variance |
| Metrics scattered | Separate output CSVs | All structural metrics joined back to one master enriched CSV |
| CIF output | Needed conversion for downstream tools | `--output_format pdb` flag in Boltz predict |
| Sequences not tracked | Only FASTA files | Per-condition sequence columns in master CSV |

---

## Data

| File | Description |
|---|---|
| `data/aggregated_lines_v2_frozen.csv` | Canonical input — 539 bars, 474 verified. **Immutable.** |
| `data/aggregated_lines_v2_enriched.csv` | Master working CSV — all columns including sequences + structural metrics |
| `data/bar_index_snapshot.json` | `bar_N -> lyric + all 5 sequences`. Authoritative attribution record. |
| `data/frozen_csv.sha256` | Hash guard for frozen CSV |

---

## Conversion: Five Modes

2×2 factorial design across standard letter mapping × BOJUXZ substitution:

| Mode | Standard 20 | BOJUXZ | Isolates |
|---|---|---|---|
| `concordance` | Freq-rank concordance | Softmax peaked draw | Canonical run |
| `alanine` | Freq-rank concordance | → A | Effect of softmax BOJUXZ |
| `random` | Freq-rank concordance | Uniform random | Is concordance better than chance? |
| `native` | AA pass-through | Softmax peaked draw | Effect of freq remapping |
| `native_alanine` | AA pass-through | → A | BOJUXZ strategy without remapping |

All runs: seed=42, λ=2.0, 80-300 AA filter.

---

## Length Filter Rationale

V1 pilot (225 bars, 28–136 AA) showed:
- `r(length, pLDDT) = −0.746` — length is the dominant structural predictor
- pTM > 0.4 (protein-like) peaks in the 40–90 AA range
- 80–300 AA chosen for v1: tight cohort eliminates residual length confound entirely

---

## Boltz-2 (Colab)

After `make convert`, upload `outputs/fastas/bars_v2_concordance.fasta` to Colab:

```python
!pip install boltz cuequivariance-torch
!boltz predict bars_v2_concordance.fasta \
    --use_msa_server \
    --output_format pdb \
    --diffusion_samples 1 \
    --out_dir boltz_outputs/
```

Download `boltz_outputs/` to `outputs/boltz/`, then `make parse-boltz`.

---

## Environment

- Local: Python 3.x, no GPU required for pipeline stages 0–2, 4–5
- Cloud: Google Colab Pro+ (A100-SXM4-80GB) for Boltz-2
- ESMFold: free API at `api.esmatlas.com` — no auth, no local install

```bash
pip install -r requirements.txt
```
