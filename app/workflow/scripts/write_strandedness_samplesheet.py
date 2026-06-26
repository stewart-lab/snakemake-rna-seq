#!/usr/bin/env python3
"""Write a copy of samples.csv with the 'strandedness' column fully populated.

Samples that already give a strandedness in the input samples.csv keep it;
samples left blank are filled in from their inferred call file
(<strand-dir>/<sample>.strandedness.txt, written by the infer_strandedness rule).
All values are normalized to forward / reverse / unstranded.
"""
import argparse
import csv
import os

VALUES = {"forward", "reverse", "unstranded"}


def normalize(value):
    value = (value or "").strip().lower()
    if value == "none":
        value = "unstranded"
    if value and value not in VALUES:
        raise SystemExit(
            f"[samplesheet] invalid strandedness '{value}'; expected one of "
            f"{sorted(VALUES)} or 'none'."
        )
    return value


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", required=True, help="input samples.csv")
    ap.add_argument("--strand-dir", required=True,
                    help="directory holding <sample>.strandedness.txt call files")
    ap.add_argument("--out", required=True, help="output samples.csv")
    args = ap.parse_args()

    with open(args.samples, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    if "strandedness" not in fieldnames:
        fieldnames.append("strandedness")

    for row in rows:
        sample = row["sample"]
        strandedness = normalize(row.get("strandedness"))
        if not strandedness:
            call_path = os.path.join(args.strand_dir, f"{sample}.strandedness.txt")
            with open(call_path) as cf:
                strandedness = normalize(cf.read())
        row["strandedness"] = strandedness

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"[samplesheet] wrote {args.out} with strandedness for {len(rows)} samples")


if __name__ == "__main__":
    main()
