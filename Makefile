# rap-snacks-v1 pipeline
# ----------------------
# Run stages in order. Each stage depends on the previous.
# Boltz-2 (step 4) is a manual Colab step — see instructions below.
#
# Quick start:
#   make audit          <- run first, review output before proceeding
#   make convert        <- runs all 5 modes
#   make esm-pilot      <- ESMFold top 25 bars
#   make enrich         <- join available metrics to master CSV
#
# Full pipeline (after Boltz Colab run):
#   make all

PYTHON = python3
PIPELINE = pipeline

# ---------------------------------------------------------------------------
# Stage 0 — Audit
# ---------------------------------------------------------------------------

.PHONY: audit
audit:
	$(PYTHON) $(PIPELINE)/00_audit.py --sample 20 --seed 99

.PHONY: audit-tight
audit-tight:
	$(PYTHON) $(PIPELINE)/00_audit.py --min-length 80 --max-length 90 --sample 20

# ---------------------------------------------------------------------------
# Stage 1 — Convert (all 5 modes)
# ---------------------------------------------------------------------------

.PHONY: convert
convert:
	$(PYTHON) $(PIPELINE)/01_convert.py --min-length 80 --max-length 90 --seed 42

.PHONY: convert-dry
convert-dry:
	$(PYTHON) $(PIPELINE)/01_convert.py --min-length 80 --max-length 90 --dry-run

# ---------------------------------------------------------------------------
# Stage 2 — ESMFold ensemble
# ---------------------------------------------------------------------------

.PHONY: esm-pilot
esm-pilot:
	$(PYTHON) $(PIPELINE)/02_esm_fold.py --pilot --top-n 25 --delay 1.0

.PHONY: esm-pilot-dry
esm-pilot-dry:
	$(PYTHON) $(PIPELINE)/02_esm_fold.py --pilot --top-n 25 --dry-run

.PHONY: esm-full
esm-full:
	$(PYTHON) $(PIPELINE)/02_esm_fold.py --delay 1.0

.PHONY: esm-resume
esm-resume:
	$(PYTHON) $(PIPELINE)/02_esm_fold.py --pilot --top-n 25 --delay 1.0 --resume

# ---------------------------------------------------------------------------
# Stage 3 — Parse Boltz outputs
# (Run after downloading boltz outputs from Colab to outputs/boltz/)
# ---------------------------------------------------------------------------

.PHONY: parse-boltz
parse-boltz:
	$(PYTHON) $(PIPELINE)/03_parse_boltz.py

# ---------------------------------------------------------------------------
# Stage 4 — FoldSeek structural homolog search
# (Run after parse-boltz — needs confidence scores for pTM filter)
# ---------------------------------------------------------------------------

.PHONY: foldseek
foldseek:
	$(PYTHON) $(PIPELINE)/04_foldseek.py --min-ptm 0.4 --delay 5.0

.PHONY: foldseek-dry
foldseek-dry:
	$(PYTHON) $(PIPELINE)/04_foldseek.py --dry-run

.PHONY: foldseek-resume
foldseek-resume:
	$(PYTHON) $(PIPELINE)/04_foldseek.py --min-ptm 0.4 --delay 8.0 --resume

# ---------------------------------------------------------------------------
# Stage 5 — Enrich master CSV with all structural metrics
# ---------------------------------------------------------------------------

.PHONY: enrich
enrich:
	$(PYTHON) $(PIPELINE)/05_enrich_csv.py

# ---------------------------------------------------------------------------
# Convenience targets
# ---------------------------------------------------------------------------

.PHONY: all
all: convert esm-pilot parse-boltz foldseek enrich

.PHONY: status
status:
	@echo "=== Data files ==="
	@ls -lh data/ 2>/dev/null || echo "  (no data/ files yet)"
	@echo ""
	@echo "=== FASTA outputs ==="
	@ls -lh outputs/fastas/ 2>/dev/null || echo "  (not converted yet)"
	@echo ""
	@echo "=== ESMFold results ==="
	@ls -lh outputs/esm/ 2>/dev/null || echo "  (not run yet)"
	@echo ""
	@echo "=== Boltz outputs ==="
	@ls -lh outputs/boltz/ 2>/dev/null || echo "  (not downloaded yet)"
	@echo ""
	@echo "=== FoldSeek results ==="
	@ls -lh outputs/foldseek/ 2>/dev/null || echo "  (not run yet)"

.PHONY: clean-outputs
clean-outputs:
	@echo "This will delete all outputs/. Are you sure? Run: rm -rf outputs/"

# ---------------------------------------------------------------------------
# Boltz-2 on Colab — MANUAL STEP
# ---------------------------------------------------------------------------
# After running `make convert`, upload outputs/fastas/bars_v2_concordance.fasta
# to Google Colab Pro+ (A100-SXM4-80GB).
#
# Colab commands:
#   !pip install boltz cuequivariance-torch
#   !boltz predict bars_v2_concordance.fasta \
#       --use_msa_server \
#       --output_format pdb \
#       --diffusion_samples 1 \
#       --out_dir boltz_outputs/
#
# Download boltz_outputs/ and place at outputs/boltz/
# Then run: make parse-boltz
# ---------------------------------------------------------------------------
