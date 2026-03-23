"""
02_esm_fold.py
--------------
ESMFold API ensemble runner for all 5 conditions.

Reads sequences from data/bar_index_snapshot.json (post-conversion).
Re-draws stochastic sequences at each seed using the same converter logic.
Appends pLDDT to outputs/esm/plddt_scores.csv after every call (crash-safe).

Nstruct plan:
  concordance:    15 seeds  (peaked softmax)
  alanine:         1 seed   (deterministic)
  random:         30 seeds  (uniform distribution)
  native:         15 seeds  (same softmax as concordance)
  native_alanine:  1 seed   (deterministic)

Pilot mode (--pilot): top N bars by iconicity only (default: 25).

Usage:
  python pipeline/02_esm_fold.py --pilot
  python pipeline/02_esm_fold.py --pilot --top-n 25 --delay 1.5
  python pipeline/02_esm_fold.py                           # full run all bars
  python pipeline/02_esm_fold.py --condition concordance   # single condition
  python pipeline/02_esm_fold.py --resume                  # skip completed
  python pipeline/02_esm_fold.py --dry-run                 # show plan, no calls
"""

import argparse
import csv
import json
import random
import ssl
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

# Ensure pipeline/ is on path for shared utils
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import CONVERTERS, clean as clean_bar

SNAPSHOT_JSON = Path("data/bar_index_snapshot.json")
ESM_DIR = Path("outputs/esm")
PLDDT_CSV = ESM_DIR / "plddt_scores.csv"

NSTRUCT = {
    "concordance": 15,
    "alanine": 1,
    "random": 30,
    "native": 15,
    "native_alanine": 1,
}

ESM_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
SSL_CONTEXT = ssl.create_default_context()
SSL_CONTEXT.check_hostname = False
SSL_CONTEXT.verify_mode = ssl.CERT_NONE

DETERMINISTIC = {"alanine", "native_alanine"}


def mean_plddt(pdb_text: str) -> float | None:
    """Extract mean pLDDT from PDB B-factor column.
    ESMFold's API stores pLDDT as 0-1 directly in the B-factor field."""
    values = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            try:
                values.append(float(line[60:66].strip()))
            except ValueError:
                pass
    if not values:
        return None
    return sum(values) / len(values)  # already 0-1


def call_esm(sequence: str, retries: int = 3, base_delay: float = 30.0) -> str | None:
    data = sequence.encode("utf-8")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                ESM_URL, data=data, method="POST",
                headers={"Content-Type": "application/x-www-form-urlencoded"},
            )
            with urllib.request.urlopen(req, context=SSL_CONTEXT, timeout=60) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            wait = base_delay * (attempt + 1)
            if e.code == 429:
                print(f"\n      [429] Rate limited — waiting {wait:.0f}s...")
                time.sleep(wait)
            elif e.code >= 500:
                print(f"\n      [HTTP {e.code}] {e.reason} — retry {attempt+1}/{retries} in {wait:.0f}s...")
                time.sleep(wait)
            else:
                print(f"\n      [HTTP {e.code}] {e.reason}")
                return None
        except Exception as e:
            print(f"\n      [ERROR] {e}")
            time.sleep(5)
    return None


def load_completed(path: Path) -> set:
    """Return set of (bar_id, condition, seed) that succeeded. Failed entries
    are intentionally excluded so --resume will retry them."""
    if not path.exists():
        return set()
    completed = set()
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("status", "ok") == "ok":
                completed.add((row["bar_id"], row["condition"], str(row["seed"])))
    return completed


def append_result(path: Path, row: dict):
    is_new = not path.exists()
    with open(path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        if is_new:
            writer.writeheader()
        writer.writerow(row)


def parse_args():
    p = argparse.ArgumentParser(description="ESMFold API ensemble runner")
    p.add_argument("--pilot", action="store_true", help="Top-N bars by iconicity")
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--condition", choices=list(NSTRUCT.keys()), default=None)
    p.add_argument("--delay", type=float, default=1.0, help="Seconds between calls")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--lambda-val", type=float, default=2.0)
    p.add_argument("--seed", type=int, default=42, help="Base seed")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()

    if not SNAPSHOT_JSON.exists():
        print(f"[ERROR] {SNAPSHOT_JSON} not found. Run 01_convert.py first.", file=sys.stderr)
        sys.exit(1)

    with open(SNAPSHOT_JSON, encoding="utf-8") as f:
        snapshot = json.load(f)

    bars = sorted(
        [{"bar_id": k, **v} for k, v in snapshot.items()],
        key=lambda b: float(b.get("aggregate_iconicity") or 0),
        reverse=True,
    )

    if args.pilot:
        bars = bars[: args.top_n]
        print(f"[02_esm_fold.py] Pilot: top {len(bars)} bars by iconicity")
    else:
        print(f"[02_esm_fold.py] Full run: {len(bars)} bars")

    conditions = [args.condition] if args.condition else list(NSTRUCT.keys())
    total_calls = sum(NSTRUCT[c] for c in conditions) * len(bars)

    print(f"  Conditions: {conditions}")
    print(f"  Planned calls: {total_calls}")
    for c in conditions:
        print(f"    {c}: {len(bars)} bars × {NSTRUCT[c]} seeds = {len(bars) * NSTRUCT[c]}")

    if args.dry_run:
        print("[DRY-RUN] No calls made.")
        return

    ESM_DIR.mkdir(parents=True, exist_ok=True)
    completed = load_completed(PLDDT_CSV) if args.resume else set()
    if completed:
        print(f"  Resuming — {len(completed)} already done, skipping.")

    call_n = 0
    for condition in conditions:
        pdb_dir = ESM_DIR / condition / "pdbs"
        pdb_dir.mkdir(parents=True, exist_ok=True)
        nstruct = NSTRUCT[condition]

        for bar in bars:
            bar_id = bar["bar_id"]
            cleaned = bar.get("lyric_cleaned", "")
            if not cleaned:
                continue

            for seed_offset in range(nstruct):
                seed = args.seed + seed_offset

                # Deterministic conditions: only one run
                if condition in DETERMINISTIC and seed_offset > 0:
                    break

                key = (bar_id, condition, str(seed))
                if key in completed:
                    continue

                # Draw sequence
                rng = random.Random(seed)
                sequence, _ = CONVERTERS[condition](cleaned, args.lambda_val, rng)

                call_n += 1
                pdb_path = pdb_dir / f"{bar_id}_seed{seed}.pdb"

                print(
                    f"  [{call_n}/{total_calls}] {bar_id} | {condition} | seed={seed} "
                    f"| len={len(sequence)}",
                    end="", flush=True,
                )

                result_row = {
                    "bar_id": bar_id,
                    "condition": condition,
                    "seed": seed,
                    "sequence_len": len(sequence),
                    "plddt": "",
                    "status": "failed",
                    "song": bar.get("genius_song_title", ""),
                    "iconicity": bar.get("aggregate_iconicity", ""),
                }

                # Reuse saved PDB if it exists — skip API call
                if pdb_path.exists():
                    pdb_text = pdb_path.read_text()
                    print(" [cached]", end="", flush=True)
                else:
                    pdb_text = call_esm(sequence)
                    if pdb_text is not None:
                        pdb_path.write_text(pdb_text)

                if pdb_text is not None:
                    plddt = mean_plddt(pdb_text)
                    result_row["plddt"] = f"{plddt:.6f}" if plddt is not None else ""
                    result_row["status"] = "ok"
                    print(f"  -> pLDDT={plddt:.4f}" if plddt is not None else "  -> pLDDT=?")
                else:
                    print("  -> FAILED")

                append_result(PLDDT_CSV, result_row)
                time.sleep(args.delay)

    print(f"\n[02_esm_fold.py] Done. {call_n} calls. Results: {PLDDT_CSV}")
    print(f"Next: python pipeline/05_enrich_csv.py")


if __name__ == "__main__":
    main()
