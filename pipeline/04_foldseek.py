"""
04_foldseek.py
--------------
FoldSeek structural homolog search via REST API.

Posts CIF/PDB files from outputs/boltz/predictions to search.foldseek.com.
Searches three databases: pdb100, afdb-swissprot, mgnify_esm30.
Writes results to outputs/foldseek/foldseek_hits.csv.

Key API notes (learned from v1):
  - target field format: "ACCESSION Description" — split on first space
  - No TM-score in web API — use `prob` (0-1) as hit quality metric
  - Rate limit: HTTP 429 after ~27 submissions -> exponential backoff
  - Effective length floor: ~35 AA — shorter bars return zero hits regardless
  - MGnify hits have empty target_name (metagenomic, no annotation)
  - Use --resume to skip bars with existing raw JSON in outputs/foldseek/raw/

Usage:
  python pipeline/04_foldseek.py --dry-run
  python pipeline/04_foldseek.py
  python pipeline/04_foldseek.py --min-ptm 0.4 --delay 8.0
  python pipeline/04_foldseek.py --resume
  python pipeline/04_foldseek.py --skip bar_5 bar_12
"""

import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
import urllib.error
from pathlib import Path

BOLTZ_PREDICTIONS = Path("outputs/boltz/predictions")
BOLTZ_CSV = Path("outputs/boltz/boltz_confidence_scores.csv")
FOLDSEEK_DIR = Path("outputs/foldseek")
RAW_DIR = FOLDSEEK_DIR / "raw"
HITS_CSV = FOLDSEEK_DIR / "foldseek_hits.csv"
SNAPSHOT_JSON = Path("data/bar_index_snapshot.json")

FOLDSEEK_SUBMIT = "https://search.foldseek.com/api/ticket"
FOLDSEEK_RESULT = "https://search.foldseek.com/api/result/{ticket}/0"
DATABASES = ["pdb100", "afdb-swissprot", "mgnify_esm30"]
FOLDSEEK_FLOOR_AA = 35  # bars shorter than this return zero hits


def submit_foldseek(struct_path: Path, databases: list[str]) -> str | None:
    """Submit a CIF/PDB file to FoldSeek. Returns ticket ID or None."""
    with open(struct_path, "rb") as f:
        file_data = f.read()

    boundary = "----FoldSeekBoundary"
    body = b""
    body += f"--{boundary}\r\n".encode()
    body += (
        f'Content-Disposition: form-data; name="q"; filename="{struct_path.name}"\r\n'
        f"Content-Type: application/octet-stream\r\n\r\n"
    ).encode()
    body += file_data + b"\r\n"

    for db in databases:
        body += f"--{boundary}\r\n".encode()
        body += f'Content-Disposition: form-data; name="database[]"\r\n\r\n{db}\r\n'.encode()

    body += f"--{boundary}--\r\n".encode()

    req = urllib.request.Request(
        FOLDSEEK_SUBMIT,
        data=body,
        method="POST",
        headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
    )
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            result = json.loads(resp.read())
            return result.get("id")
    except urllib.error.HTTPError as e:
        if e.code == 429:
            raise
        print(f"    [HTTP {e.code}] submit failed")
        return None
    except Exception as e:
        print(f"    [ERROR] submit: {e}")
        return None


def poll_foldseek(ticket: str, max_wait: int = 120, interval: int = 5) -> dict | None:
    """Poll FoldSeek until job complete. Returns raw result dict or None."""
    result_url = FOLDSEEK_RESULT.format(ticket=ticket)
    status_url = f"https://search.foldseek.com/api/ticket/{ticket}"
    waited = 0
    while waited < max_wait:
        try:
            with urllib.request.urlopen(status_url, timeout=15) as resp:
                status = json.loads(resp.read())
            if status.get("status") == "COMPLETE":
                with urllib.request.urlopen(result_url, timeout=30) as resp:
                    return json.loads(resp.read())
            elif status.get("status") == "ERROR":
                print(f"    [FOLDSEEK ERROR] ticket={ticket}")
                return None
        except Exception as e:
            print(f"    [POLL ERROR] {e}")
        time.sleep(interval)
        waited += interval
    print(f"    [TIMEOUT] ticket={ticket} after {max_wait}s")
    return None


def parse_hits(raw: dict, bar_id: str, snap: dict) -> list[dict]:
    """Parse raw FoldSeek result into flat rows."""
    rows = []
    results = raw.get("results", [])
    for db_result in results:
        db = db_result.get("db", "")
        alignments = db_result.get("alignments", [])
        if not alignments:
            continue
        for aln in alignments[0] if isinstance(alignments[0], list) else [alignments[0]]:
            target = aln.get("target", "")
            parts = target.split(" ", 1)
            accession = parts[0]
            target_name = parts[1] if len(parts) > 1 else ""
            rows.append({
                "bar_id": bar_id,
                "canonical_bar": snap.get("canonical_bar", ""),
                "song": snap.get("genius_song_title", ""),
                "attribution": snap.get("attribution", ""),
                "iconicity": snap.get("aggregate_iconicity", ""),
                "seq_len": snap.get("lyric_cleaned_len", ""),
                "db": db,
                "accession": accession,
                "target_name": target_name,
                "prob": aln.get("prob", ""),
                "evalue": aln.get("evalue", ""),
                "qlen": aln.get("qLen", ""),
                "tlen": aln.get("tLen", ""),
                "foldseek_result": "hit",
            })
    if not rows:
        rows.append({
            "bar_id": bar_id,
            "canonical_bar": snap.get("canonical_bar", ""),
            "song": snap.get("genius_song_title", ""),
            "attribution": snap.get("attribution", ""),
            "iconicity": snap.get("aggregate_iconicity", ""),
            "seq_len": snap.get("lyric_cleaned_len", ""),
            "db": "", "accession": "", "target_name": "",
            "prob": "", "evalue": "", "qlen": "", "tlen": "",
            "foldseek_result": "no_hits_novel",
        })
    return rows


def parse_args():
    p = argparse.ArgumentParser(description="FoldSeek structural homolog search")
    p.add_argument("--boltz-dir", default=str(BOLTZ_PREDICTIONS))
    p.add_argument("--min-ptm", type=float, default=0.4,
                   help="Min pTM to search (default: 0.4)")
    p.add_argument("--min-plddt", type=float, default=0.0,
                   help="Min pLDDT to search (default: 0.0 = no floor)")
    p.add_argument("--min-seq-len", type=int, default=FOLDSEEK_FLOOR_AA,
                   help=f"Min sequence length (default: {FOLDSEEK_FLOOR_AA})")
    p.add_argument("--delay", type=float, default=5.0,
                   help="Seconds between submissions (use 8.0 after rate limiting)")
    p.add_argument("--resume", action="store_true",
                   help="Skip bars with existing raw JSON in outputs/foldseek/raw/")
    p.add_argument("--skip", nargs="*", default=[],
                   help="bar IDs to skip (e.g. --skip bar_5 bar_12)")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    boltz_dir = Path(args.boltz_dir)

    if not SNAPSHOT_JSON.exists():
        print(f"[ERROR] {SNAPSHOT_JSON} not found. Run 01_convert.py first.", file=sys.stderr)
        sys.exit(1)

    with open(SNAPSHOT_JSON, encoding="utf-8") as f:
        snapshot = json.load(f)

    # Load Boltz confidence scores for filtering
    boltz_scores = {}
    if BOLTZ_CSV.exists():
        with open(BOLTZ_CSV, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                boltz_scores[row["bar_id"]] = {
                    "ptm": float(row.get("ptm") or 0),
                    "plddt": float(row.get("plddt") or 0),
                    "seq_len": int(row.get("boltz_seq_len") or 0),
                }
    else:
        print(f"[WARN] {BOLTZ_CSV} not found — running on all bars with structure files")

    # Select bars to search
    bar_dirs = sorted(boltz_dir.glob("bar_*"))
    to_search = []

    for bar_dir in bar_dirs:
        bar_id = bar_dir.name
        if bar_id in args.skip:
            continue

        snap = snapshot.get(bar_id, {})
        seq_len = snap.get("lyric_cleaned_len", 0) or 0

        scores = boltz_scores.get(bar_id, {})
        ptm = scores.get("ptm", 0)
        plddt = scores.get("plddt", 0)

        # Apply filters
        if boltz_scores and ptm < args.min_ptm:
            continue
        if boltz_scores and plddt < args.min_plddt:
            continue
        if seq_len < args.min_seq_len:
            print(f"  [SKIP] {bar_id} — seq_len={seq_len} < {args.min_seq_len} (FoldSeek floor)")
            continue

        # Find structure file
        structs = list(bar_dir.glob("*_model_0.cif")) + list(bar_dir.glob("*_model_0.pdb"))
        if not structs:
            continue

        to_search.append((bar_id, structs[0], snap, ptm, plddt))

    print(f"[04_foldseek.py] Bars to search: {len(to_search)}")
    print(f"  Filters: pTM≥{args.min_ptm}, pLDDT≥{args.min_plddt}, seq_len≥{args.min_seq_len}")

    if args.dry_run:
        for bar_id, struct, snap, ptm, plddt in to_search:
            seq_len = snap.get("lyric_cleaned_len", "?")
            print(f"  {bar_id:10s}  ptm={ptm:.3f}  plddt={plddt:.3f}  len={seq_len}  {struct.name}")
        print("[DRY-RUN] No submissions made.")
        return

    FOLDSEEK_DIR.mkdir(parents=True, exist_ok=True)
    RAW_DIR.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for idx, (bar_id, struct_path, snap, ptm, plddt) in enumerate(to_search, 1):
        raw_path = RAW_DIR / f"{bar_id}.json"

        if args.resume and raw_path.exists():
            print(f"  [{idx}/{len(to_search)}] {bar_id} — resuming from cache")
            with open(raw_path, encoding="utf-8") as f:
                raw = json.load(f)
            all_rows.extend(parse_hits(raw, bar_id, snap))
            continue

        print(f"  [{idx}/{len(to_search)}] {bar_id} (pTM={ptm:.3f}, pLDDT={plddt:.3f})... ", end="", flush=True)

        attempt = 0
        ticket = None
        while ticket is None and attempt < 5:
            try:
                ticket = submit_foldseek(struct_path, DATABASES)
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    wait = 30 * (attempt + 1)
                    print(f"\n    [429] Rate limited — waiting {wait}s...")
                    time.sleep(wait)
                    attempt += 1
                else:
                    print(f"\n    [HTTP {e.code}] giving up")
                    break

        if ticket is None:
            print("FAILED (no ticket)")
            continue

        raw = poll_foldseek(ticket)
        if raw is None:
            print("FAILED (no result)")
            continue

        with open(raw_path, "w", encoding="utf-8") as f:
            json.dump(raw, f)

        hits = parse_hits(raw, bar_id, snap)
        hit_count = sum(1 for h in hits if h["foldseek_result"] == "hit")
        print(f"{hit_count} hits")
        all_rows.extend(hits)

        time.sleep(args.delay)

    if all_rows:
        with open(HITS_CSV, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"\n  Wrote {HITS_CSV}  ({len(all_rows)} rows)")

    hit_bars = sum(1 for r in all_rows if r["foldseek_result"] == "hit")
    novel_bars = sum(1 for r in all_rows if r["foldseek_result"] == "no_hits_novel")
    print(f"\n  Bars with hits:       {hit_bars}")
    print(f"  Bars novel (no hits): {novel_bars}")
    print(f"\nNext: python pipeline/05_enrich_csv.py")


if __name__ == "__main__":
    main()
