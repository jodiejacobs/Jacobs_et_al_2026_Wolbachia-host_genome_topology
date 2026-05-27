#!/usr/bin/env python
"""Annotate diffHic differential interactions with overlapping genes.

Uses FlyBase GFF3 (dmel-all-r6.XX.gff.gz) which already has 2L/2R/etc.
chromosome names and FBgn IDs as feature IDs.
"""

import argparse
import gzip
import re
import pandas as pd
import pybedtools

ID_RE = re.compile(r'ID=([^;]+)')
NAME_RE = re.compile(r'Name=([^;]+)')
BIOTYPE_RE = re.compile(r'(?:gene_biotype|biotype|so_term_name)=([^;]+)')


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gff_to_gene_bed(gff_path, out_bed):
    """Extract gene features from FlyBase GFF3, write BED6+."""
    rows = []
    with _open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = f[0]
            start = int(f[3]) - 1   # GFF 1-based inclusive -> BED 0-based
            end = int(f[4])
            strand = f[6] if f[6] in ("+", "-") else "."
            attrs = f[8]
            fbgn = ID_RE.search(attrs)
            name = NAME_RE.search(attrs)
            biotype = BIOTYPE_RE.search(attrs)
            label = "|".join([
                fbgn.group(1) if fbgn else ".",
                name.group(1) if name else ".",
                biotype.group(1) if biotype else ".",
            ])
            rows.append((chrom, start, end, label, ".", strand))
    bed = pybedtools.BedTool(rows).sort()
    bed.saveas(out_bed)
    return bed


def anchors_to_bed(df, side):
    rows = [
        (str(r[f"chr{side}"]), int(r[f"start{side}"]), int(r[f"end{side}"]), str(i))
        for i, r in df.iterrows()
    ]
    return pybedtools.BedTool(rows).sort()


def overlap_genes(anchor_bed, gene_bed):
    """anchor_idx -> list of 'FBgn|name|biotype'."""
    out = {}
    for iv in anchor_bed.intersect(gene_bed, wa=True, wb=True):
        idx = iv.fields[3]
        gene = iv.fields[7]
        out.setdefault(idx, []).append(gene)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--interactions", required=True)
    ap.add_argument("--gff", required=True, help="FlyBase GFF3 (.gff or .gff.gz)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--gene-bed", default="flybase_genes.bed")
    args = ap.parse_args()

    df = pd.read_csv(args.interactions)
    df["chr1"] = df["chr1"].astype(str)
    df["chr2"] = df["chr2"].astype(str)
    print(f"Loaded {len(df)} interactions")

    gene_bed = gff_to_gene_bed(args.gff, args.gene_bed)
    print(f"Extracted {gene_bed.count()} gene features -> {args.gene_bed}")

    hits1 = overlap_genes(anchors_to_bed(df, 1), gene_bed)
    hits2 = overlap_genes(anchors_to_bed(df, 2), gene_bed)

    df["anchor1_genes"] = [";".join(hits1.get(str(i), [])) for i in df.index]
    df["anchor2_genes"] = [";".join(hits2.get(str(i), [])) for i in df.index]
    df["n_anchor1_genes"] = df["anchor1_genes"].apply(
        lambda s: 0 if s == "" else s.count(";") + 1)
    df["n_anchor2_genes"] = df["anchor2_genes"].apply(
        lambda s: 0 if s == "" else s.count(";") + 1)

    for col in ["anchor1_genes", "anchor2_genes", "intervening_genes"]:
        if col in df.columns:
            df[col] = df[col].replace("", "NA")

    df.to_csv(args.out, index=False)
    n_either = ((df["n_anchor1_genes"] > 0) | (df["n_anchor2_genes"] > 0)).sum()
    n_both = ((df["n_anchor1_genes"] > 0) & (df["n_anchor2_genes"] > 0)).sum()
    print(f"≥1 gene at either anchor: {n_either} / {len(df)}")
    print(f"genes at both anchors:    {n_both} / {len(df)}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
