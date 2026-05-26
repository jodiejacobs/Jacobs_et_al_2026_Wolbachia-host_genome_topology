#!/usr/bin/env python3
"""
insulation_to_tads.py

Convert a cooltools insulation table (with boundary calls) into a TAD list that
diffDomain's `dvsd multiple` can consume.

cooltools insulation gives BOUNDARIES, not domains. A TAD is defined here as the
interval between two consecutive boundaries on the same chromosome. The first
domain on a chromosome runs from the chromosome start (or first bin) to the first
boundary; the last runs from the last boundary to the chromosome end.

diffDomain expects a tab-separated BED-like file:
    - column 1: chrom
    - column 2: start
    - column 3: end
    - MUST have a header line (any names; only first 3 cols are used)
    - chrom names must match the cooler exactly (here: bare '2L', 'X', ... no 'chr')

Usage:
    python insulation_to_tads.py \
        --insulation chromatin_structure_DOX_cooltools_insulation.tsv \
        --boundary-col is_boundary_150000 \
        --output DOX_tads.bed \
        [--min-bins 2] [--reso 64000] [--chroms 2L 2R 3L 3R 4 X Y] \
        [--drop-bad-bins]
"""
import argparse
import sys
import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--insulation", required=True,
                   help="cooltools insulation TSV (one condition).")
    p.add_argument("--output", required=True,
                   help="Output TAD bed (tab-separated, with header).")
    p.add_argument("--boundary-col", default=None,
                   help="Name of the boolean boundary column "
                        "(e.g. is_boundary_150000). If omitted, the script "
                        "auto-detects the first column starting with 'is_boundary'.")
    p.add_argument("--strength-col", default=None,
                   help="Optional boundary_strength_* column. If given, only "
                        "boundaries with non-NaN strength are used (cooltools "
                        "writes strength only at called boundaries anyway).")
    p.add_argument("--min-bins", type=int, default=2,
                   help="Minimum number of bins a TAD must span to be kept. "
                        "diffDomain skips TADs shorter than --min_nbin (default 10) "
                        "at test time, but pre-filtering keeps the list clean. "
                        "[default: 2]")
    p.add_argument("--reso", type=int, default=None,
                   help="Bin size in bp. If omitted, inferred from the table "
                        "(end-start of the first row).")
    p.add_argument("--chroms", nargs="+", default=None,
                   help="Restrict to these chromosomes, in this order. "
                        "If omitted, uses all chroms present, in file order.")
    p.add_argument("--drop-bad-bins", action="store_true",
                   help="If set, boundaries flagged is_bad_bin==True are ignored.")
    return p.parse_args()


def main():
    args = parse_args()

    df = pd.read_csv(args.insulation, sep="\t")
    required = {"chrom", "start", "end"}
    if not required.issubset(df.columns):
        sys.exit(f"ERROR: insulation table missing one of {required}; "
                 f"found columns: {list(df.columns)}")

    # --- resolve boundary column ---
    bcol = args.boundary_col
    if bcol is None:
        cands = [c for c in df.columns if c.startswith("is_boundary")]
        if not cands:
            sys.exit("ERROR: no 'is_boundary*' column found and --boundary-col "
                     "not given.")
        bcol = cands[0]
        sys.stderr.write(f"[insulation_to_tads] auto-detected boundary column: {bcol}\n")
    if bcol not in df.columns:
        sys.exit(f"ERROR: boundary column '{bcol}' not in table.")

    # --- infer resolution ---
    reso = args.reso
    if reso is None:
        reso = int(df["end"].iloc[0] - df["start"].iloc[0])
        sys.stderr.write(f"[insulation_to_tads] inferred resolution: {reso} bp\n")

    # cooltools writes booleans as True/False; coerce robustly
    is_bound = df[bcol].astype(str).str.lower().isin(["true", "1", "1.0"])
    if args.strength_col and args.strength_col in df.columns:
        is_bound &= df[args.strength_col].notna()
    if args.drop_bad_bins and "is_bad_bin" in df.columns:
        bad = df["is_bad_bin"].astype(str).str.lower().isin(["true", "1", "1.0"])
        is_bound &= ~bad

    df = df.assign(_is_bound=is_bound)

    chroms = args.chroms if args.chroms else list(dict.fromkeys(df["chrom"]))

    tads = []
    for chrom in chroms:
        sub = df[df["chrom"] == chrom].sort_values("start").reset_index(drop=True)
        if sub.empty:
            sys.stderr.write(f"[insulation_to_tads] WARNING: no rows for {chrom}\n")
            continue
        # Boundary positions = the 'start' of bins flagged as boundaries.
        bsites = sub.loc[sub["_is_bound"], "start"].tolist()
        chrom_start = int(sub["start"].iloc[0])
        chrom_end = int(sub["end"].iloc[-1])

        # Build edge list: chrom_start, each boundary, chrom_end
        edges = [chrom_start] + [int(b) for b in bsites] + [chrom_end]
        edges = sorted(set(edges))
        for s, e in zip(edges[:-1], edges[1:]):
            span_bins = (e - s) / reso
            if span_bins >= args.min_bins:
                tads.append((chrom, s, e))

    if not tads:
        sys.exit("ERROR: no TADs produced. Check boundary column / chrom names.")

    out = pd.DataFrame(tads, columns=["chrom", "start", "end"])
    out.to_csv(args.output, sep="\t", index=False)
    sys.stderr.write(
        f"[insulation_to_tads] wrote {len(out)} TADs across "
        f"{out['chrom'].nunique()} chroms to {args.output}\n")


if __name__ == "__main__":
    main()