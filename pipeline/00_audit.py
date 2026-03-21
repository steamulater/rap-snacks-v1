"""
00_audit.py
-----------
Validates data/aggregated_lines_v2_frozen.csv before any conversion.

Outputs:
  - Column inventory, row counts, attribution breakdown
  - Sequence length distribution for fasta-eligible bars
  - How many bars pass the 80-90 AA length filter
  - 20 randomly sampled rows for human spot-check
  - Null/missing field report

Run:
  python pipeline/00_audit.py
  python pipeline/00_audit.py --sample 20 --seed 42
"""

import argparse
import csv
import random
import re
import sys
from collections import Counter
from pathlib import Path

FROZEN_CSV = Path("data/aggregated_lines_v2_frozen.csv")

REQUIRED_COLS = [
    "canonical_bar",
    "genius_song_title",
    "attribution",
    "fasta_eligible",
    "aggregate_iconicity",
    "divergence_badge",
    "genius_url",
]

# Chars not in the standard 20 AA single-letter codes
NON_STANDARD = set("BOJUXZ")


def clean_bar(text: str) -> str:
    return re.sub(r"[^A-Za-z]", "", text).upper()


def parse_args():
    p = argparse.ArgumentParser(description="Audit frozen CSV before conversion")
    p.add_argument("--input", default=str(FROZEN_CSV))
    p.add_argument("--sample", type=int, default=20, help="Rows to spot-check (default: 20)")
    p.add_argument("--seed", type=int, default=99, help="RNG seed for row sampling")
    p.add_argument(
        "--min-length", type=int, default=80, help="Lower AA length bound (default: 80)"
    )
    p.add_argument(
        "--max-length", type=int, default=90, help="Upper AA length bound (default: 90)"
    )
    return p.parse_args()


def hr(char="-", width=70):
    print(char * width)


def main():
    args = parse_args()
    path = Path(args.input)

    if not path.exists():
        print(f"[ERROR] File not found: {path}", file=sys.stderr)
        sys.exit(1)

    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        columns = reader.fieldnames or []

    hr("=")
    print(f"AUDIT REPORT — {path}")
    hr("=")

    # --- Column inventory ---
    print(f"\nColumns ({len(columns)}):")
    for col in columns:
        print(f"  {col}")

    missing_cols = [c for c in REQUIRED_COLS if c not in columns]
    if missing_cols:
        print(f"\n[WARN] Missing required columns: {missing_cols}")
    else:
        print("\n[OK] All required columns present.")

    # --- Row counts ---
    hr()
    total = len(rows)
    eligible = [r for r in rows if r.get("fasta_eligible", "").strip().lower() == "true"]
    attributed = [
        r for r in eligible
        if r.get("attribution", "").strip().lower() in ("primary", "featured")
    ]

    print(f"\nRow counts:")
    print(f"  Total rows:                {total}")
    print(f"  fasta_eligible=True:       {len(eligible)}")
    print(f"  eligible + attributed:     {len(attributed)}")

    # --- Attribution breakdown ---
    attr_counts = Counter(r.get("attribution", "").strip().lower() for r in rows)
    print(f"\nAttribution breakdown (all rows):")
    for k, v in sorted(attr_counts.items(), key=lambda x: -x[1]):
        print(f"  {k or '(blank)':20s}  {v}")

    # --- Badge breakdown ---
    badge_counts = Counter(r.get("divergence_badge", "").strip() for r in attributed)
    print(f"\nDivergence badge breakdown (eligible+attributed):")
    for k, v in sorted(badge_counts.items(), key=lambda x: -x[1]):
        print(f"  {k or '(blank)':20s}  {v}")

    # --- Sequence length analysis ---
    hr()
    print(f"\nSequence length analysis (eligible + attributed rows):")
    lengths = []
    for r in attributed:
        bar = (r.get("canonical_bar") or "").strip()
        cleaned = clean_bar(bar)
        lengths.append(len(cleaned))

    if lengths:
        lengths.sort()
        buckets = {
            "<40": sum(1 for l in lengths if l < 40),
            "40-79": sum(1 for l in lengths if 40 <= l < 80),
            f"{args.min_length}-{args.max_length} (TARGET)": sum(
                1 for l in lengths if args.min_length <= l <= args.max_length
            ),
            "91-120": sum(1 for l in lengths if 91 <= l <= 120),
            ">120": sum(1 for l in lengths if l > 120),
        }
        print(f"  Min: {min(lengths)}  Max: {max(lengths)}  Median: {lengths[len(lengths)//2]}")
        print(f"\n  Length buckets:")
        for bucket, count in buckets.items():
            bar_fill = "#" * count
            print(f"    {bucket:30s}  {count:4d}  {bar_fill}")

        in_window = [
            r for r, l in zip(attributed, [len(clean_bar(r.get("canonical_bar","").strip())) for r in attributed])
            if args.min_length <= l <= args.max_length
        ]
        print(f"\n  [RESULT] Bars passing {args.min_length}-{args.max_length} AA filter: {len(in_window)}")
    else:
        print("  No eligible attributed rows found.")

    # --- Null / missing field report ---
    hr()
    print(f"\nNull/missing field report (eligible + attributed rows):")
    null_counts = Counter()
    for r in attributed:
        for col in REQUIRED_COLS:
            val = r.get(col, "").strip()
            if not val or val.lower() in ("none", "nan", ""):
                null_counts[col] += 1

    if null_counts:
        for col, count in sorted(null_counts.items(), key=lambda x: -x[1]):
            print(f"  {col:30s}  {count} missing")
    else:
        print("  [OK] No missing values in required columns.")

    # --- BOJUXZ density ---
    hr()
    print(f"\nBOJUXZ character density (eligible + attributed rows):")
    bojuxz_counts = []
    for r in attributed:
        bar = (r.get("canonical_bar") or "").strip()
        cleaned = clean_bar(bar)
        count = sum(1 for c in cleaned if c in NON_STANDARD)
        if cleaned:
            bojuxz_counts.append(count / len(cleaned))

    if bojuxz_counts:
        mean_density = sum(bojuxz_counts) / len(bojuxz_counts)
        nonzero = sum(1 for d in bojuxz_counts if d > 0)
        print(f"  Bars with any BOJUXZ:  {nonzero} / {len(bojuxz_counts)}")
        print(f"  Mean BOJUXZ density:   {mean_density:.3%}")

    # --- 20 random row spot-check ---
    hr("=")
    print(f"\nSPOT-CHECK: {args.sample} randomly sampled rows (seed={args.seed})")
    hr("=")

    rng = random.Random(args.seed)
    sample = rng.sample(attributed, min(args.sample, len(attributed)))

    for idx, row in enumerate(sample, 1):
        bar = (row.get("canonical_bar") or "").strip()
        cleaned = clean_bar(bar)
        length = len(cleaned)
        in_window = args.min_length <= length <= args.max_length
        flag = "  <-- IN WINDOW" if in_window else ""

        print(f"\n  [{idx:02d}]  bar: {bar[:80]}")
        print(f"        song: {row.get('genius_song_title','?')}")
        print(f"        attribution: {row.get('attribution','?')}  |  badge: {row.get('divergence_badge','?')}")
        print(f"        iconicity: {row.get('aggregate_iconicity','?')}")
        print(f"        cleaned len: {length}{flag}")
        print(f"        url: {row.get('genius_url','?')}")

    hr("=")
    print(f"\nReview the spot-check rows above before running 01_convert.py.")
    print(f"If everything looks correct:")
    print(f"  python pipeline/01_convert.py")
    hr("=")


if __name__ == "__main__":
    main()
