#!/usr/bin/env python3
"""Turn RSeQC `infer_experiment.py` output into a strandedness call.

`infer_experiment.py` reports, for a subsample of aligned reads, the fraction
explained by each read-orientation pattern. This script reads that report,
decides whether the library is forward-stranded, reverse-stranded, or
unstranded, and writes a small shell-sourceable file with the matching RSEM and
Picard settings so the RSEM and Picard rules can just `source` it.

Orientation patterns (from infer_experiment.py):
  paired-end  forward:  "1++,1--,2+-,2-+"   reverse:  "1+-,1-+,2++,2--"
  single-end  forward:  "++,--"             reverse:  "+-,-+"

"forward" means the read (or read 1) is on the same strand as the transcript
(fr-secondstrand / ISF). "reverse" is the dUTP/Illumina-stranded case
(fr-firststrand / ISR), which is the most common stranded protocol.

Mapping applied:
  forward    -> RSEM --strandedness forward (--forward-prob 1),
                Picard STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND
  reverse    -> RSEM --strandedness reverse (--forward-prob 0),
                Picard STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND
  unstranded -> RSEM --strandedness none    (--forward-prob 0.5),
                Picard STRAND_SPECIFICITY=NONE
"""
import argparse
import os
import re

# (label, regex capturing the trailing fraction) for the forward and reverse
# orientation lines, covering both paired-end and single-end report wording.
FORWARD_PATTERNS = (
    r'"1\+\+,1--,2\+-,2-\+":\s*([0-9.]+)',  # paired-end
    r'"\+\+,--":\s*([0-9.]+)',              # single-end
)
REVERSE_PATTERNS = (
    r'"1\+-,1-\+,2\+\+,2--":\s*([0-9.]+)',  # paired-end
    r'"\+-,-\+":\s*([0-9.]+)',              # single-end
)

# strandedness call -> (rsem --strandedness, rsem --forward-prob, Picard STRAND_SPECIFICITY)
SETTINGS = {
    "forward":    ("forward", "1",   "SECOND_READ_TRANSCRIPTION_STRAND"),
    "reverse":    ("reverse", "0",   "FIRST_READ_TRANSCRIPTION_STRAND"),
    "unstranded": ("none",    "0.5", "NONE"),
}


def _first_fraction(text, patterns):
    """First trailing fraction matched by any of `patterns`, or None."""
    for pattern in patterns:
        m = re.search(pattern, text)
        if m:
            return float(m.group(1))
    return None


def call_strandedness(forward, reverse, threshold):
    """Classify the library from the forward/reverse fractions.

    The decision is made over the *explained* reads only (undetermined reads are
    ignored): if one orientation accounts for at least `threshold` of the
    explained reads, the library is stranded that way; otherwise unstranded.
    """
    forward = forward or 0.0
    reverse = reverse or 0.0
    explained = forward + reverse
    if explained <= 0:
        return "unstranded", 0.0
    forward_frac = forward / explained
    if forward_frac >= threshold:
        return "forward", forward_frac
    if forward_frac <= 1 - threshold:
        return "reverse", forward_frac
    return "unstranded", forward_frac


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--infer-experiment", required=True,
                    help="infer_experiment.py output text file")
    ap.add_argument("--threshold", type=float, default=0.8,
                    help="min fraction of explained reads in one orientation to "
                         "call the library stranded (default: 0.8)")
    ap.add_argument("--out", required=True,
                    help="output shell-sourceable strandedness file")
    args = ap.parse_args()

    with open(args.infer_experiment) as fh:
        text = fh.read()

    forward = _first_fraction(text, FORWARD_PATTERNS)
    reverse = _first_fraction(text, REVERSE_PATTERNS)
    if forward is None and reverse is None:
        raise SystemExit(
            f"[strandedness] could not find orientation fractions in "
            f"{args.infer_experiment}; is this infer_experiment.py output?"
        )

    call, forward_frac = call_strandedness(forward, reverse, args.threshold)
    rsem_strandedness, rsem_forward_prob, picard_strand = SETTINGS[call]

    undetermined = max(0.0, 1.0 - (forward or 0.0) - (reverse or 0.0))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(f"STRANDEDNESS={call}\n")
        fh.write(f"RSEM_STRANDEDNESS={rsem_strandedness}\n")
        fh.write(f"RSEM_FORWARD_PROB={rsem_forward_prob}\n")
        fh.write(f"PICARD_STRAND={picard_strand}\n")
        fh.write(f"FORWARD_FRACTION={forward or 0.0:.4f}\n")
        fh.write(f"REVERSE_FRACTION={reverse or 0.0:.4f}\n")
        fh.write(f"UNDETERMINED_FRACTION={undetermined:.4f}\n")

    print(f"[strandedness] call={call} "
          f"(forward={forward or 0.0:.4f}, reverse={reverse or 0.0:.4f}, "
          f"threshold={args.threshold}) -> "
          f"RSEM --strandedness {rsem_strandedness}, "
          f"Picard STRAND_SPECIFICITY={picard_strand}")


if __name__ == "__main__":
    main()
