#!/usr/bin/env python3
import sys
import pandas as pd

INFILE = "apg/results/01_linear_behavior/tables/bwa_vs_mm2_contig_comparison.with_quality.tsv"
OUTFILE = "apg/results/01_linear_behavior/tables/section1_linear_canonical.tsv"

def normalize_q(x):
    """Map various pipeline labels to canonical tiers: NP/RG/PP/UM."""
    if pd.isna(x):
        return "UM"
    s = str(x).strip().upper()

    # Your current table uses labels like "NEARLY PERFECT"
    if s in {"NEARLY PERFECT", "NEARLY_PERFECT", "NP"}:
        return "NP"
    if s in {"REASONABLY GOOD", "REASONABLY_GOOD", "RG"}:
        return "RG"
    if s in {"PARTIAL", "POOR", "PARTIALLY MAPPED", "PARTIAL/POOR", "PP"}:
        return "PP"
    if s in {"UNMAPPED", "UM", "NA", "NONE", ""}:
        return "UM"

    # If already canonical, keep
    if s in {"NP", "RG", "PP", "UM"}:
        return s

    # Conservative fallback
    return "UM"

def behavior_class(bwa_q, mm2_q):
    """Assign behavior class based on your tier-combination definitions."""
    bwa_q = normalize_q(bwa_q)
    mm2_q = normalize_q(mm2_q)

    # L classes (stable / aligner-sensitive)
    if (bwa_q == "NP") and (mm2_q == "NP"):
        return "L1"
    if (bwa_q in {"NP", "RG"}) and (mm2_q == "RG"):
        return "L2"
    if ((bwa_q == "NP" and mm2_q == "RG") or (bwa_q == "RG" and mm2_q == "NP")):
        return "L3"

    # Graph-candidate classes
    if (bwa_q in {"NP", "RG"}) and (mm2_q == "PP"):
        return "G1"
    if (bwa_q in {"NP", "RG", "PP"}) and (mm2_q == "UM"):
        return "G2"
    if (bwa_q == "PP") and (mm2_q == "PP"):
        return "G3"

    # Rare/odd
    if (bwa_q == "PP") and (mm2_q in {"RG", "NP"}):
        return "X"

    # Anything else -> X
    return "X"

def main():
    df = pd.read_csv(INFILE, sep="\t", dtype=str)

    # Convert key numeric columns where helpful
    num_cols = ["bwa_qlen", "bwa_coverage", "bwa_identity", "bwa_mapq",
                "mm2_qlen", "mm2_coverage", "mm2_identity", "mm2_mapq"]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    # Choose a consistent contig length
    df["contig_len"] = df.get("bwa_qlen", pd.Series([pd.NA]*len(df))).fillna(df.get("mm2_qlen", pd.Series([pd.NA]*len(df))))

    # Normalize tiers
    df["bwa_tier"] = df["bwa_quality"].apply(normalize_q) if "bwa_quality" in df.columns else "UM"
    df["mm2_tier"] = df["mm2_quality"].apply(normalize_q) if "mm2_quality" in df.columns else "UM"

    # Assign behavior class
    df["behavior_class"] = [behavior_class(b, m) for b, m in zip(df["bwa_tier"], df["mm2_tier"])]

    # Canonical output columns
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
    out = df[keep].copy().rename(columns={"qname": "contig_id"})

    out.to_csv(OUTFILE, sep="\t", index=False)
    print(f"[OK] Wrote {OUTFILE} with {out.shape[0]} contigs and {out.shape[1]} columns.", file=sys.stderr)

if __name__ == "__main__":
    main()
