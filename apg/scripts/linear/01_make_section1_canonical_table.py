#!/usr/bin/env python3
import sys
import pandas as pd

INFILE = "apg/results/01_linear_behavior/tables/bwa_vs_mm2_contig_comparison.with_quality.tsv"
OUTFILE = "apg/results/01_linear_behavior/tables/section1_linear_canonical.tsv"

# Expected quality tier labels
VALID = {"NP", "RG", "PP", "UM"}

def normalize_q(x):
    if pd.isna(x):
        return "UM"
    s = str(x).strip()
    # Some pipelines use lowercase or words; normalize lightly
    s_up = s.upper()
    # Keep only known tiers if present; else pass through (we'll still classify conservatively)
    return s_up

def behavior_class(bwa_q, mm2_q):
    bwa_q = normalize_q(bwa_q)
    mm2_q = normalize_q(mm2_q)

    # Your definitions
    if (bwa_q == "NP") and (mm2_q == "NP"):
        return "L1"
    if (bwa_q in {"NP", "RG"}) and (mm2_q == "RG"):
        return "L2"
    if ((bwa_q == "NP" and mm2_q == "RG") or (bwa_q == "RG" and mm2_q == "NP")):
        return "L3"
    if (bwa_q in {"NP", "RG"}) and (mm2_q == "PP"):
        return "G1"
    if (bwa_q in {"NP", "RG", "PP"}) and (mm2_q == "UM"):
        return "G2"
    if (bwa_q == "PP") and (mm2_q == "PP"):
        return "G3"
    if (bwa_q == "PP") and (mm2_q in {"RG", "NP"}):
        return "X"
    # Anything else falls into X for inspection
    return "X"

def main():
    df = pd.read_csv(INFILE, sep="\t", dtype=str)

    # Cast numeric columns where useful (keep strings elsewhere)
    num_cols = ["bwa_qlen","bwa_coverage","bwa_identity","bwa_mapq",
                "mm2_qlen","mm2_coverage","mm2_identity","mm2_mapq"]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Consistent contig length
    df["contig_len"] = df["bwa_qlen"].fillna(df["mm2_qlen"])

    # Normalize quality tiers
    df["bwa_tier"] = df["bwa_quality"].apply(normalize_q)
    df["mm2_tier"] = df["mm2_quality"].apply(normalize_q)

    # Assign behavior class
    df["behavior_class"] = [behavior_class(b, m) for b, m in zip(df["bwa_tier"], df["mm2_tier"])]

    # Keep only canonical columns for the manuscript + plotting
    keep = [
        "qname", "contig_len",
        "bwa_coverage", "bwa_identity", "bwa_mapq", "bwa_tier",
        "mm2_coverage", "mm2_identity", "mm2_mapq", "mm2_tier",
        "behavior_class",
        "bwa_mapped", "mm2_mapped", "same_chr",
        "bwa_good", "mm2_good",
        "bwa_nearly_perfect", "mm2_nearly_perfect",
        "bwa_tname", "mm2_tname"
    ]
    keep = [c for c in keep if c in df.columns]
    out = df[keep].copy()

    # Rename qname to contig_id for clarity
    out = out.rename(columns={"qname": "contig_id"})

    out.to_csv(OUTFILE, sep="\t", index=False)
    print(f"[OK] Wrote {OUTFILE} with {out.shape[0]} contigs and {out.shape[1]} columns.", file=sys.stderr)

if __name__ == "__main__":
    main()
