#!/usr/bin/env python3
"""
Replicate concordance for Micro-C contact matrices.

For every condition, compares replicate 1 vs replicate 2 (and, as a reference,
all other sample pairs) at each resolution using:
  * Pearson r of raw cis contact counts (zeros included, upper triangle,
    distance <= max_dist), pooled across chromosomes and per chromosome
  * Pearson r of log2(count + 1)
  * Stratum-adjusted correlation coefficient (SCC; HiCRep, Yang et al. 2017),
    re-implemented here (mean-filter smoothing, per-diagonal Pearson, weighted by
    n_d * sqrt(var_x * var_y)) so it does not depend on the hicrep package version.
    SCC removes the distance-decay inflation of plain Pearson. Dense matrices are
    used, so SCC is skipped for resolutions < 4 kb (NaN).

Inputs are the per-replicate mcool files {mcool_dir}/{sample}-{rep}.matrix.mcool.

Outputs (prefix = --output_prefix):
  _pairwise.tsv         one row per sample pair x resolution x chromosome
  _replicate_summary.tsv  rep1 vs rep2 per condition x resolution (pooled + SCC)
  _concordance.pdf      supplemental figure
"""
import argparse
import itertools
import sys

import cooler
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter

DEFAULT_H = {1000: 10, 8000: 3, 32000: 1, 128000: 1}


def cis_sums(clr_a, clr_b, chrom, max_dist, res):
    """Sufficient statistics for Pearson r between two cis matrices of one chromosome.

    Works on sparse pixels so 1 kb matrices never get densified. Zeros are
    counted via n (number of upper-triangle bins within max_dist).
    """
    size = clr_a.chromsizes[chrom]
    nb = int(np.ceil(size / res))
    kmax = min(int(max_dist // res), nb - 1)
    # number of upper-triangle bins with 0 <= j - i <= kmax
    n = nb * (kmax + 1) - kmax * (kmax + 1) // 2

    def load(clr):
        m = clr.matrix(balance=False, sparse=True).fetch(chrom).tocoo()
        keep = (m.col >= m.row) & ((m.col - m.row) <= kmax)
        return pd.Series(m.data[keep].astype(np.float64),
                         index=pd.MultiIndex.from_arrays([m.row[keep], m.col[keep]]))

    a, b = load(clr_a), load(clr_b)
    df = pd.concat([a.rename("a"), b.rename("b")], axis=1).fillna(0.0)
    out = {"n": n}
    for tag, f in (("raw", lambda x: x), ("log", lambda x: np.log2(x + 1.0))):
        # log2(0+1)=0 so zero bins still contribute nothing to the sums
        x, y = f(df["a"].values), f(df["b"].values)
        out[tag] = (x.sum(), y.sum(), (x * x).sum(), (y * y).sum(), (x * y).sum())
    return out


def scc_chrom(clr_a, clr_b, chrom, h, max_dist, res, rng):
    """HiCRep-style SCC for one chromosome (cis, raw counts, depth-matched)."""
    nb = int(np.ceil(clr_a.chromsizes[chrom] / res))
    if res < 4000 or nb > 8000:
        return np.nan
    kmax = min(int(max_dist // res), nb - 1)
    A = clr_a.matrix(balance=False).fetch(chrom).astype(np.float64)
    B = clr_b.matrix(balance=False).fetch(chrom).astype(np.float64)
    A, B = np.triu(A), np.triu(B)  # cis upper triangle incl. diagonal
    ta, tb = A.sum(), B.sum()
    # depth match: binomially thin the deeper replicate to the shallower one
    if ta > tb:
        A = rng.binomial(A.astype(np.int64), tb / ta).astype(np.float64)
    elif tb > ta:
        B = rng.binomial(B.astype(np.int64), ta / tb).astype(np.float64)
    # smooth on the symmetric matrix so the filter sees both sides of the diagonal
    A = uniform_filter(A + np.triu(A, 1).T, size=2 * h + 1, mode="constant")
    B = uniform_filter(B + np.triu(B, 1).T, size=2 * h + 1, mode="constant")
    num = den = 0.0
    for d in range(1, kmax + 1):
        x, y = np.diagonal(A, d), np.diagonal(B, d)
        keep = (x > 0) | (y > 0)
        n = int(keep.sum())
        if n < 3:
            continue
        x, y = x[keep], y[keep]
        sx, sy = x.std(), y.std()
        if sx == 0 or sy == 0:
            continue
        rho = np.corrcoef(x, y)[0, 1]
        w = n * sx * sy
        num += w * rho
        den += w
    return num / den if den > 0 else np.nan


def pearson_from_sums(n, s):
    sx, sy, sxx, syy, sxy = s
    cov = sxy - sx * sy / n
    vx = sxx - sx * sx / n
    vy = syy - sy * sy / n
    return cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mcool_dir", required=True)
    p.add_argument("--samples", nargs="+", required=True, help="e.g. JW18DOX JW18wMel")
    p.add_argument("--replicates", nargs="+", default=["1", "2"])
    p.add_argument("--resolutions", nargs="+", type=int, default=[8000, 32000])
    p.add_argument("--chromosomes", nargs="+", default=["2L", "2R", "3L", "3R", "4", "X"])
    p.add_argument("--max_dist", type=int, default=2_000_000,
                   help="max genomic distance (bp) for Pearson and SCC")
    p.add_argument("--h", nargs="*", default=None,
                   help="HiCRep smoothing as RES:H pairs, e.g. 8000:3 32000:1")
    p.add_argument("--output_prefix", required=True)
    args = p.parse_args()

    h_map = dict(DEFAULT_H)
    for kv in args.h or []:
        r, h = kv.split(":")
        h_map[int(r)] = int(h)

    names = [f"{s}-{r}" for s in args.samples for r in args.replicates]
    rows = []
    for res in args.resolutions:
        clrs = {nm: cooler.Cooler(f"{args.mcool_dir}/{nm}.matrix.mcool::resolutions/{res}")
                for nm in names}
        for a, b in itertools.combinations(names, 2):
            print(f"[{res}] {a} vs {b}", file=sys.stderr, flush=True)
            tot = {"n": 0, "raw": np.zeros(5), "log": np.zeros(5)}
            per_chrom = {}
            for chrom in args.chromosomes:
                st = cis_sums(clrs[a], clrs[b], chrom, args.max_dist, res)
                per_chrom[chrom] = {"pearson_raw": pearson_from_sums(st["n"], st["raw"]),
                                    "pearson_log2": pearson_from_sums(st["n"], st["log"])}
                tot["n"] += st["n"]
                tot["raw"] += np.array(st["raw"])
                tot["log"] += np.array(st["log"])
            # SCC is computed per chromosome; downsample so depth differences
            # between replicates do not drive the score
            rng = np.random.default_rng(0)
            scc = [scc_chrom(clrs[a], clrs[b], c, h_map[res], args.max_dist, res, rng)
                   for c in args.chromosomes]
            for chrom, s in zip(args.chromosomes, scc):
                per_chrom[chrom]["scc"] = s
            for chrom, d in per_chrom.items():
                rows.append(dict(res=res, a=a, b=b, chrom=chrom, **d))
            # "all" row: pooled Pearson over all chromosomes, mean SCC
            rows.append(dict(res=res, a=a, b=b, chrom="all",
                             pearson_raw=pearson_from_sums(tot["n"], tot["raw"]),
                             pearson_log2=pearson_from_sums(tot["n"], tot["log"]),
                             scc=float(np.nanmean(scc))))

    df = pd.DataFrame(rows)
    df["cond_a"] = df["a"].str.rsplit("-", n=1).str[0]
    df["cond_b"] = df["b"].str.rsplit("-", n=1).str[0]
    df["same_condition"] = df["cond_a"] == df["cond_b"]
    df["h"] = df["res"].map(h_map)
    df["max_dist"] = args.max_dist
    df.to_csv(f"{args.output_prefix}_pairwise.tsv", sep="\t", index=False)

    reps = df[df.same_condition & (df.chrom == "all")]
    summ = reps[["res", "cond_a", "pearson_raw", "pearson_log2", "scc", "h", "max_dist"]] \
        .rename(columns={"cond_a": "condition"})
    summ.to_csv(f"{args.output_prefix}_replicate_summary.tsv", sep="\t", index=False)
    print(summ.to_string(index=False))

    make_figure(df, args, f"{args.output_prefix}_concordance.pdf")


def make_figure(df, args, out_pdf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res_list = args.resolutions
    metrics = [("pearson_raw", "Pearson r (raw counts)"), ("scc", "HiCRep SCC")]
    fig, axes = plt.subplots(len(metrics), len(res_list),
                             figsize=(4.2 * len(res_list), 3.8 * len(metrics)), squeeze=False)
    names = [f"{s}-{r}" for s in args.samples for r in args.replicates]
    for i, (m, label) in enumerate(metrics):
        for j, res in enumerate(res_list):
            ax = axes[i][j]
            sub = df[(df.res == res) & (df.chrom == "all")]
            mat = pd.DataFrame(np.nan, index=names, columns=names)
            for _, r in sub.iterrows():
                mat.loc[r.a, r.b] = mat.loc[r.b, r.a] = r[m]
            arr = mat.to_numpy(copy=True)
            np.fill_diagonal(arr, 1.0)
            vals = arr[~np.eye(len(names), dtype=bool)]
            im = ax.imshow(arr, cmap="viridis",
                           vmin=np.nanmin(vals) - 0.01, vmax=1.0)
            ax.set_xticks(range(len(names)))
            ax.set_yticks(range(len(names)))
            ax.set_xticklabels(names, rotation=90, fontsize=6)
            ax.set_yticklabels(names, fontsize=6)
            for x in range(len(names)):
                for y in range(len(names)):
                    ax.text(x, y, f"{arr[y, x]:.2f}", ha="center", va="center",
                            fontsize=5, color="w")
            ax.set_title(f"{label}\n{res // 1000} kb", fontsize=8)
            fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(out_pdf)


if __name__ == "__main__":
    main()
