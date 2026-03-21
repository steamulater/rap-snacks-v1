"""
01_convert.py
-------------
FASTA conversion v2. Runs all five modes in a single pass.

V2 improvements over v1:
  - All five modes converted in one run — sequences stored in master enriched CSV
  - bar_index_snapshot.json written at conversion time — lyric attribution is
    permanently locked and never depends on CSV row ordering again
  - Length filter: 80-90 AA by default (configurable --min-length / --max-length)
  - Pre/post conversion agreement check: len(lyric_cleaned) == len(fasta_seq_*)
    asserted for every row before writing output
  - SHA-256 hash of frozen CSV written on first run; subsequent runs verify it
  - Output format documented for Boltz: add --output_format pdb to boltz predict

Outputs:
  data/aggregated_lines_v2_enriched.csv   master CSV + all new columns
  data/bar_index_snapshot.json            bar_N -> lyric + sequences (authoritative)
  data/frozen_csv.sha256                  hash guard
  outputs/fastas/bars_v2_{mode}.fasta     one FASTA per mode
  outputs/masks/mask_v2_{mode}.json       mutation masks per mode

Run:
  python pipeline/01_convert.py
  python pipeline/01_convert.py --min-length 80 --max-length 90 --seed 42
  python pipeline/01_convert.py --dry-run   (show counts, write nothing)
"""

import argparse
import csv
import hashlib
import json
import math
import random
import re
import sys
from pathlib import Path

FROZEN_CSV = Path("data/aggregated_lines_v2_frozen.csv")
ENRICHED_CSV = Path("data/aggregated_lines_v2_enriched.csv")
SNAPSHOT_JSON = Path("data/bar_index_snapshot.json")
HASH_FILE = Path("data/frozen_csv.sha256")
FASTAS_DIR = Path("outputs/fastas")
MASKS_DIR = Path("outputs/masks")

MODES = ["concordance", "alanine", "random", "native", "native_alanine"]
SEEDS = {
    "concordance": 42,
    "alanine": 42,      # deterministic — seed ignored
    "random": 42,
    "native": 42,
    "native_alanine": 42,  # deterministic — seed ignored
}

FASTA_LINE_WIDTH = 60

# ---------------------------------------------------------------------------
# Frequency tables (identical to v1 — do not modify)
# ---------------------------------------------------------------------------

AA_BY_ABUNDANCE = [
    ("L", 9.66), ("A", 8.25), ("G", 7.07), ("V", 6.87), ("E", 6.75),
    ("S", 6.56), ("I", 5.96), ("K", 5.84), ("R", 5.53), ("D", 5.45),
    ("T", 5.34), ("P", 4.70), ("N", 4.06), ("Q", 3.93), ("F", 3.86),
    ("H", 2.27), ("Y", 2.92), ("M", 2.42), ("C", 1.37), ("W", 1.08),
]

AA_BY_RARITY = list(reversed(AA_BY_ABUNDANCE))  # W -> L
CANONICAL_AAS = set(aa for aa, _ in AA_BY_ABUNDANCE)
ALL_AAS = [aa for aa, _ in AA_BY_ABUNDANCE]

LETTER_TO_AA = {
    "E": "L", "T": "A", "A": "G", "I": "V", "N": "E", "S": "S",
    "R": "I", "H": "T", "L": "P", "D": "K", "C": "R", "M": "D",
    "F": "F", "P": "Q", "G": "N", "W": "H", "Y": "C", "V": "W",
    "K": "Y", "Q": "M",
}

NON_STANDARD = {"Z": 0.09, "J": 0.16, "X": 0.23, "B": 1.48, "U": 2.73, "O": 7.64}
ENGLISH_MAX_PCT = 12.49


# ---------------------------------------------------------------------------
# Softmax helpers (identical to v1)
# ---------------------------------------------------------------------------

def build_prob_vector(char_eng_pct: float, lambda_val: float) -> list:
    char_norm = 1.0 - (char_eng_pct / ENGLISH_MAX_PCT)
    n = len(AA_BY_RARITY)
    distances = []
    for i in range(n):
        aa_rarity_norm = 1.0 - (i / (n - 1))
        distances.append(-lambda_val * abs(char_norm - aa_rarity_norm))
    max_d = max(distances)
    exps = [math.exp(d - max_d) for d in distances]
    total = sum(exps)
    return [e / total for e in exps]


def sample_softmax(prob_vector: list, rng: random.Random) -> str:
    r = rng.random()
    cumulative = 0.0
    for i, p in enumerate(prob_vector):
        cumulative += p
        if r <= cumulative:
            return AA_BY_RARITY[i][0]
    return AA_BY_RARITY[-1][0]


# ---------------------------------------------------------------------------
# Per-mode converters — returns (sequence, mutation_mask)
# ---------------------------------------------------------------------------

def clean(bar: str) -> str:
    return re.sub(r"[^A-Za-z]", "", bar).upper()


def convert_concordance(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask, cache = [], [], {}
    for letter in cleaned:
        if letter in LETTER_TO_AA:
            seq.append(LETTER_TO_AA[letter])
        elif letter in NON_STANDARD:
            if letter not in cache:
                cache[letter] = build_prob_vector(NON_STANDARD[letter], lambda_val)
            mask.append(len(seq))
            seq.append(sample_softmax(cache[letter], rng))
    return "".join(seq), mask


def convert_alanine(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask = [], []
    for letter in cleaned:
        if letter in LETTER_TO_AA:
            seq.append(LETTER_TO_AA[letter])
        elif letter in NON_STANDARD:
            mask.append(len(seq))
            seq.append("A")
    return "".join(seq), mask


def convert_random(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask = [], []
    for letter in cleaned:
        if letter in LETTER_TO_AA:
            seq.append(LETTER_TO_AA[letter])
        elif letter in NON_STANDARD:
            mask.append(len(seq))
            seq.append(rng.choice(ALL_AAS))
    return "".join(seq), mask


def convert_native(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask, cache = [], [], {}
    for letter in cleaned:
        if letter in CANONICAL_AAS:
            seq.append(letter)
        elif letter in NON_STANDARD:
            if letter not in cache:
                cache[letter] = build_prob_vector(NON_STANDARD[letter], lambda_val)
            mask.append(len(seq))
            seq.append(sample_softmax(cache[letter], rng))
    return "".join(seq), mask


def convert_native_alanine(cleaned: str, lambda_val: float, rng: random.Random):
    seq, mask = [], []
    for letter in cleaned:
        if letter in CANONICAL_AAS:
            seq.append(letter)
        elif letter in NON_STANDARD:
            mask.append(len(seq))
            seq.append("A")
    return "".join(seq), mask


CONVERTERS = {
    "concordance": convert_concordance,
    "alanine": convert_alanine,
    "random": convert_random,
    "native": convert_native,
    "native_alanine": convert_native_alanine,
}


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def wrap(seq: str, width: int = FASTA_LINE_WIDTH) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def fasta_header(bar_n: int, mode: str, row: dict) -> str:
    iconicity = row.get("aggregate_iconicity", "n/a")
    try:
        iconicity = f"{float(iconicity):.4f}"
    except (ValueError, TypeError):
        iconicity = "n/a"
    song = (row.get("genius_song_title") or "").strip()
    attribution = row.get("attribution", "").strip()
    badge = row.get("divergence_badge", "").strip()
    return (
        f">bar_{bar_n} | mode={mode} | iconicity={iconicity} "
        f"| badge={badge} | {song} | {attribution}"
    )


# ---------------------------------------------------------------------------
# CSV hash guard
# ---------------------------------------------------------------------------

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def check_or_write_hash(csv_path: Path, hash_path: Path, dry_run: bool) -> bool:
    current = sha256_file(csv_path)
    if hash_path.exists():
        stored = hash_path.read_text().strip()
        if stored != current:
            print(
                f"\n[ERROR] Frozen CSV hash mismatch!\n"
                f"  Stored:  {stored}\n"
                f"  Current: {current}\n"
                f"  The frozen CSV has been modified since the last conversion run.\n"
                f"  If this is intentional, delete {hash_path} and re-run.\n"
                f"  If not, restore the original frozen CSV before proceeding.",
                file=sys.stderr,
            )
            return False
        print(f"[OK] Frozen CSV hash verified: {current[:12]}...")
    else:
        if not dry_run:
            hash_path.write_text(current)
            print(f"[OK] Frozen CSV hash written: {current[:12]}...")
        else:
            print(f"[DRY-RUN] Would write hash: {current[:12]}...")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Convert frozen CSV to FASTAs — all 5 modes")
    p.add_argument("--input", default=str(FROZEN_CSV))
    p.add_argument("--min-length", type=int, default=80)
    p.add_argument("--max-length", type=int, default=90)
    p.add_argument("--lambda-val", type=float, default=2.0)
    p.add_argument("--seed", type=int, default=42, help="Base seed (mode-specific offset applied)")
    p.add_argument("--dry-run", action="store_true", help="Print counts, write nothing")
    return p.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input)

    if not input_path.exists():
        print(f"[ERROR] Input not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Hash guard
    if not check_or_write_hash(input_path, HASH_FILE, args.dry_run):
        sys.exit(1)

    # Load CSV
    with open(input_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        all_rows = list(reader)
        original_cols = list(reader.fieldnames or [])

    print(f"\n[01_convert.py] Loaded {len(all_rows)} rows from {input_path}")
    print(f"  Length filter: {args.min_length}–{args.max_length} AA")
    print(f"  Lambda: {args.lambda_val}  |  Base seed: {args.seed}")

    # Filter eligible rows
    eligible = []
    for row in all_rows:
        fe = row.get("fasta_eligible", "").strip().lower()
        attr = row.get("attribution", "").strip().lower()
        bar = (row.get("canonical_bar") or "").strip()
        if fe == "true" and attr in ("primary", "featured") and bar:
            eligible.append(row)

    print(f"  Eligible (fasta_eligible=True, attributed): {len(eligible)}")

    # Clean and filter by length
    filtered = []
    too_short = too_long = 0
    for row in eligible:
        bar = (row.get("canonical_bar") or "").strip()
        cleaned = clean(bar)
        length = len(cleaned)
        if length < args.min_length:
            too_short += 1
        elif length > args.max_length:
            too_long += 1
        else:
            filtered.append((row, cleaned, length))

    print(f"  Filtered out (too short <{args.min_length}): {too_short}")
    print(f"  Filtered out (too long  >{args.max_length}): {too_long}")
    print(f"  Passing length filter: {len(filtered)}")

    if not filtered:
        print("\n[WARN] No bars passed the length filter. Adjust --min-length / --max-length.")
        sys.exit(0)

    if args.dry_run:
        print("\n[DRY-RUN] No files written.")
        return

    # Create output dirs
    FASTAS_DIR.mkdir(parents=True, exist_ok=True)
    MASKS_DIR.mkdir(parents=True, exist_ok=True)

    # Build per-mode RNGs with deterministic offsets
    rngs = {mode: random.Random(args.seed + i) for i, mode in enumerate(MODES)}

    # Convert all bars, all modes
    results = []  # list of dicts, one per bar
    for bar_n, (row, cleaned, length) in enumerate(filtered):
        entry = {
            "bar_id": f"bar_{bar_n}",
            "bar_n": bar_n,
            "canonical_bar": row.get("canonical_bar", "").strip(),
            "lyric_cleaned": cleaned,
            "lyric_cleaned_len": length,
            "genius_song_title": row.get("genius_song_title", "").strip(),
            "genius_url": row.get("genius_url", "").strip(),
            "attribution": row.get("attribution", "").strip(),
            "aggregate_iconicity": row.get("aggregate_iconicity", ""),
            "divergence_badge": row.get("divergence_badge", "").strip(),
            "sequences": {},
            "masks": {},
        }
        for mode in MODES:
            seq, mask = CONVERTERS[mode](cleaned, args.lambda_val, rngs[mode])
            entry["sequences"][mode] = seq
            entry["masks"][mode] = mask

        # --- Agreement check ---
        for mode in MODES:
            seq_len = len(entry["sequences"][mode])
            if seq_len != length:
                print(
                    f"\n[ERROR] Length mismatch at bar_{bar_n}, mode={mode}: "
                    f"lyric_cleaned_len={length}, fasta_seq_len={seq_len}",
                    file=sys.stderr,
                )
                sys.exit(1)

        results.append(entry)

    print(f"\n[OK] Agreement check passed: all {len(results)} bars, all 5 modes.")

    # Write FASTAs and masks
    fasta_buffers = {mode: [] for mode in MODES}
    mask_records = {mode: [] for mode in MODES}

    for entry in results:
        for mode in MODES:
            seq = entry["sequences"][mode]
            mask = entry["masks"][mode]
            header = fasta_header(entry["bar_n"], mode, {
                "genius_song_title": entry["genius_song_title"],
                "attribution": entry["attribution"],
                "aggregate_iconicity": entry["aggregate_iconicity"],
                "divergence_badge": entry["divergence_badge"],
            })
            fasta_buffers[mode].append(f"{header}\n{wrap(seq)}")
            mask_records[mode].append({
                "bar_id": entry["bar_id"],
                "canonical_bar": entry["canonical_bar"],
                "song": entry["genius_song_title"],
                "mode": mode,
                "lyric_cleaned": entry["lyric_cleaned"],
                "sequence": seq,
                "sequence_length": len(seq),
                "mutation_mask": mask,
                "mask_count": len(mask),
                "mask_pct": round(len(mask) / len(seq) * 100, 2) if seq else 0,
            })

    for mode in MODES:
        fasta_path = FASTAS_DIR / f"bars_v2_{mode}.fasta"
        with open(fasta_path, "w", encoding="utf-8") as f:
            f.write("\n\n".join(fasta_buffers[mode]) + "\n")
        print(f"  Wrote {fasta_path}  ({len(fasta_buffers[mode])} sequences)")

        mask_path = MASKS_DIR / f"mask_v2_{mode}.json"
        with open(mask_path, "w", encoding="utf-8") as f:
            json.dump(mask_records[mode], f, indent=2)
        print(f"  Wrote {mask_path}")

    # Write bar_index_snapshot.json — the authoritative lyric→sequence record
    snapshot = {}
    for entry in results:
        snapshot[entry["bar_id"]] = {
            "bar_id": entry["bar_id"],
            "bar_n": entry["bar_n"],
            "canonical_bar": entry["canonical_bar"],
            "lyric_cleaned": entry["lyric_cleaned"],
            "lyric_cleaned_len": entry["lyric_cleaned_len"],
            "genius_song_title": entry["genius_song_title"],
            "genius_url": entry["genius_url"],
            "attribution": entry["attribution"],
            "aggregate_iconicity": entry["aggregate_iconicity"],
            "divergence_badge": entry["divergence_badge"],
            "fasta_seq_concordance": entry["sequences"]["concordance"],
            "fasta_seq_alanine": entry["sequences"]["alanine"],
            "fasta_seq_random": entry["sequences"]["random"],
            "fasta_seq_native": entry["sequences"]["native"],
            "fasta_seq_native_alanine": entry["sequences"]["native_alanine"],
            "lambda_val": args.lambda_val,
            "seed": args.seed,
            "min_length": args.min_length,
            "max_length": args.max_length,
        }

    with open(SNAPSHOT_JSON, "w", encoding="utf-8") as f:
        json.dump(snapshot, f, indent=2)
    print(f"\n  Wrote {SNAPSHOT_JSON}  ({len(snapshot)} bars)")

    # Write enriched CSV — original columns + all new columns
    new_cols = [
        "bar_id", "lyric_cleaned", "lyric_cleaned_len",
        "fasta_seq_concordance", "fasta_seq_alanine", "fasta_seq_random",
        "fasta_seq_native", "fasta_seq_native_alanine", "fasta_seq_len",
    ]

    # Build lookup: canonical_bar -> enrichment (first match only)
    enrichment = {}
    for entry in results:
        key = entry["canonical_bar"]
        if key not in enrichment:
            enrichment[key] = {
                "bar_id": entry["bar_id"],
                "lyric_cleaned": entry["lyric_cleaned"],
                "lyric_cleaned_len": str(entry["lyric_cleaned_len"]),
                "fasta_seq_concordance": entry["sequences"]["concordance"],
                "fasta_seq_alanine": entry["sequences"]["alanine"],
                "fasta_seq_random": entry["sequences"]["random"],
                "fasta_seq_native": entry["sequences"]["native"],
                "fasta_seq_native_alanine": entry["sequences"]["native_alanine"],
                "fasta_seq_len": str(entry["lyric_cleaned_len"]),
            }

    out_cols = original_cols + [c for c in new_cols if c not in original_cols]

    with open(ENRICHED_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols, extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            bar = (row.get("canonical_bar") or "").strip()
            extra = enrichment.get(bar, {k: "" for k in new_cols})
            writer.writerow({**row, **extra})

    print(f"  Wrote {ENRICHED_CSV}")

    print(f"\n[01_convert.py] Done. {len(results)} bars converted across 5 modes.")
    print(f"\nNext steps:")
    print(f"  Boltz-2: upload outputs/fastas/bars_v2_concordance.fasta to Colab")
    print(f"  ESMFold: python pipeline/02_esm_fold.py --pilot --top-n 25")


if __name__ == "__main__":
    main()
