"""
05_enrich_csv.py
----------------
Joins all structural metrics back into the master enriched CSV.
Run after each new data source is available.

Sources joined:
  outputs/boltz/boltz_confidence_scores.csv  -> boltz_plddt, boltz_ptm, etc.
  outputs/esm/plddt_scores.csv               -> esm_plddt_* columns
  outputs/foldseek/foldseek_hits.csv         -> foldseek_result, best_hit, prob

Output: data/aggregated_lines_v2_enriched.csv (updated in-place)
        data/aggregated_lines_v2_enriched_backup.csv (backup before overwrite)

Run at any point — missing sources are skipped with a note.

Usage:
  python pipeline/05_enrich_csv.py
  python pipeline/05_enrich_csv.py --no-backup
"""

import argparse
import csv
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

ENRICHED_CSV = Path("data/aggregated_lines_v2_enriched.csv")
SNAPSHOT_JSON = Path("data/bar_index_snapshot.json")
BOLTZ_CSV = Path("outputs/boltz/boltz_summary.csv")
ESM_CSV = Path("outputs/esm/plddt_scores.csv")
FOLDSEEK_CSV = Path("outputs/foldseek/foldseek_hits.csv")


def load_boltz(path: Path) -> dict:
    """bar_id -> boltz metric dict"""
    data = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            bar_id = row["bar_id"]
            data[bar_id] = {
                "boltz_plddt":      row.get("plddt_mean") or row.get("plddt", ""),
                "boltz_plddt_best": row.get("plddt_best", ""),
                "boltz_plddt_sd":   row.get("plddt_sd", ""),
                "boltz_ptm":        row.get("ptm_mean") or row.get("ptm", ""),
                "boltz_ptm_best":   row.get("ptm_best", ""),
                "boltz_confidence": row.get("confidence_mean") or row.get("confidence", ""),
                "boltz_structural_class": row.get("structural_class", ""),
                "boltz_n_models":   row.get("n_models", ""),
            }
    return data


def load_esm(path: Path) -> dict:
    """bar_id -> {condition -> list of plddt values}"""
    data = defaultdict(lambda: defaultdict(list))
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            bar_id = row.get("bar_id", "")
            condition = row.get("condition", "")
            plddt = row.get("plddt", "")
            if bar_id and condition and plddt:
                try:
                    data[bar_id][condition].append(float(plddt))
                except ValueError:
                    pass
    # Compute per-bar per-condition mean and sd
    result = {}
    for bar_id, conditions in data.items():
        entry = {}
        for cond, values in conditions.items():
            if values:
                mean = sum(values) / len(values)
                variance = sum((v - mean) ** 2 for v in values) / len(values) if len(values) > 1 else 0
                sd = variance ** 0.5
                entry[f"esm_plddt_{cond}_mean"] = f"{mean:.6f}"
                entry[f"esm_plddt_{cond}_sd"] = f"{sd:.6f}"
                entry[f"esm_plddt_{cond}_n"] = str(len(values))
        result[bar_id] = entry
    return result


def load_foldseek(path: Path) -> dict:
    """bar_id -> best hit info"""
    best = {}
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            bar_id = row.get("bar_id", "")
            result_type = row.get("foldseek_result", "")
            prob_str = row.get("prob", "")

            if bar_id not in best:
                best[bar_id] = {
                    "foldseek_result": result_type,
                    "foldseek_best_hit": "",
                    "foldseek_best_prob": "",
                    "foldseek_best_db": "",
                    "foldseek_organism": "",
                }

            if result_type == "hit" and prob_str:
                try:
                    prob = float(prob_str)
                    existing = best[bar_id].get("foldseek_best_prob", "")
                    existing_prob = float(existing) if existing else -1
                    if prob > existing_prob:
                        best[bar_id].update({
                            "foldseek_result": "hit",
                            "foldseek_best_hit": row.get("target_name", ""),
                            "foldseek_best_prob": prob_str,
                            "foldseek_best_db": row.get("db", ""),
                            "foldseek_organism": "",  # could join taxon if available
                        })
                except ValueError:
                    pass
    return best


def structural_class(plddt_str: str, ptm_str: str) -> str:
    """Assign 2D quadrant label."""
    try:
        plddt = float(plddt_str)
        ptm = float(ptm_str)
    except (ValueError, TypeError):
        return ""
    high_plddt = plddt >= 0.6
    high_ptm = ptm >= 0.4
    if high_plddt and high_ptm:
        return "confident_protein_like"
    elif high_plddt and not high_ptm:
        return "confident_alien"
    elif not high_plddt and high_ptm:
        return "uncertain_protein_like"
    else:
        return "disordered"


def parse_args():
    p = argparse.ArgumentParser(description="Join structural metrics into master CSV")
    p.add_argument("--no-backup", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()

    if not ENRICHED_CSV.exists():
        print(f"[ERROR] {ENRICHED_CSV} not found. Run 01_convert.py first.", file=sys.stderr)
        sys.exit(1)

    with open(ENRICHED_CSV, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_cols = list(reader.fieldnames or [])

    print(f"[05_enrich_csv.py] Loaded {len(rows)} rows from {ENRICHED_CSV}")

    # Load available data sources
    boltz_data = {}
    if BOLTZ_CSV.exists():
        boltz_data = load_boltz(BOLTZ_CSV)
        print(f"  Boltz-2 scores loaded: {len(boltz_data)} bars")
    else:
        print(f"  [SKIP] {BOLTZ_CSV} not found")

    esm_data = {}
    if ESM_CSV.exists():
        esm_data = load_esm(ESM_CSV)
        print(f"  ESMFold scores loaded: {len(esm_data)} bars")
    else:
        print(f"  [SKIP] {ESM_CSV} not found")

    foldseek_data = {}
    if FOLDSEEK_CSV.exists():
        foldseek_data = load_foldseek(FOLDSEEK_CSV)
        print(f"  FoldSeek hits loaded: {len(foldseek_data)} bars")
    else:
        print(f"  [SKIP] {FOLDSEEK_CSV} not found")

    # Determine new columns to add
    new_cols = []

    if boltz_data:
        new_cols += [
            "boltz_plddt", "boltz_plddt_best", "boltz_plddt_sd",
            "boltz_ptm", "boltz_ptm_best",
            "boltz_confidence", "boltz_structural_class", "boltz_n_models",
        ]

    if esm_data:
        sample = next(iter(esm_data.values()))
        new_cols += sorted(sample.keys())

    if foldseek_data:
        new_cols += [
            "foldseek_result", "foldseek_best_hit", "foldseek_best_prob",
            "foldseek_best_db",
        ]

    # Deduplicate while preserving order
    seen = set(original_cols)
    added_cols = []
    for c in new_cols:
        if c not in seen:
            added_cols.append(c)
            seen.add(c)

    out_cols = original_cols + added_cols

    # Enrich rows
    updated = 0
    for row in rows:
        bar_id = row.get("bar_id", "")
        if not bar_id:
            continue

        if boltz_data and bar_id in boltz_data:
            row.update(boltz_data[bar_id])
            updated += 1

        if esm_data and bar_id in esm_data:
            row.update(esm_data[bar_id])

        if foldseek_data and bar_id in foldseek_data:
            row.update(foldseek_data[bar_id])

    # Backup then overwrite
    if not args.no_backup:
        backup = ENRICHED_CSV.with_name("aggregated_lines_v2_enriched_backup.csv")
        shutil.copy2(ENRICHED_CSV, backup)
        print(f"  Backup written: {backup}")

    with open(ENRICHED_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"  Updated {ENRICHED_CSV}  ({updated} bars enriched, {len(added_cols)} new columns)")
    print(f"\nNew columns added: {added_cols}")


if __name__ == "__main__":
    main()
