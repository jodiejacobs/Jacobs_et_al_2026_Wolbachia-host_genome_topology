#!/usr/bin/env python3
"""
Supplemental figure: Micro-C replicate concordance.

Reads the table written by replicate_concordance.py plus the mcool files.

  A, B  Replicate 1 vs replicate 2 cis contacts per condition (hexbin, depth-
        scaled log2(count+1), pixels nonzero in either replicate, <= max_dist)
  C, D  Per-chromosome SCC and Pearson r, replicate 1 vs replicate 2
  E, F  Pairwise SCC and Pearson r among all replicates
Annotated r / SCC values come from the pairwise table (computed on all pixels,
zeros included), not from the plotted subset.
"""
import argparse

import cooler
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Project colour convention: wMel dark green, uninfected (DOX) light green
COLORS = {"JW18DOX": "#8fcb84", "JW18wMel": "#09aa4b"}
LABELS = {"JW18DOX": "DOX (uninfected)", "JW18wMel": "wMel"}
EDGE = "#2b2b2b"
MUTED = "#6b6b6b"

plt.rcParams.update({
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
})


def panel_letter(ax, s):
    ax.text(-0.18, 1.08, s, transform=ax.transAxes, fontsize=10, fontweight="bold")


def load_pixels(clr, max_dist, res):
    """cis upper-triangle pixels within max_dist, all chromosomes, as {(chrom,i,j): count}."""
    out = {}
    for chrom in clr.chromnames:
        if chrom == "Y":
            continue
        m = clr.matrix(balance=False, sparse=True).fetch(chrom).tocoo()
        keep = (m.col >= m.row) & ((m.col - m.row) * res <= max_dist)
        out[chrom] = pd.Series(m.data[keep].astype(float),
                               index=pd.MultiIndex.from_arrays([m.row[keep], m.col[keep]]))
    return pd.concat(out)


def scatter_panel(ax, mcool_dir, cond, res, max_dist, stats):
    c1 = cooler.Cooler(f"{mcool_dir}/{cond}-1.matrix.mcool::resolutions/{res}")
    c2 = cooler.Cooler(f"{mcool_dir}/{cond}-2.matrix.mcool::resolutions/{res}")
    a, b = load_pixels(c1, max_dist, res), load_pixels(c2, max_dist, res)
    df = pd.concat([a.rename("a"), b.rename("b")], axis=1).fillna(0.0)
    # scale both replicates to their mean depth so depth differences do not shift the cloud
    ta, tb = df["a"].sum(), df["b"].sum()
    mean_t = (ta + tb) / 2
    x = np.log2(df["a"].values * mean_t / ta + 1)
    y = np.log2(df["b"].values * mean_t / tb + 1)
    hb = ax.hexbin(x, y, gridsize=60, bins="log", mincnt=1, cmap="Greys", linewidths=0, rasterized=True)
    lim = max(x.max(), y.max()) * 1.02
    ax.plot([0, lim], [0, lim], color=COLORS[cond], lw=1.0, zorder=3)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Replicate 1, log$_2$(count + 1)")
    ax.set_ylabel("Replicate 2, log$_2$(count + 1)")
    ax.set_title(f"{LABELS[cond]}, {res // 1000} kb", loc="left", fontsize=7.5)
    s = stats[(stats.res == res) & (stats.cond_a == cond) & (stats.same_condition)
              & (stats.chrom == "all")].iloc[0]
    ax.text(0.04, 0.96, f"Pearson r = {s.pearson_raw:.3f}\nSCC = {s.scc:.3f}",
            transform=ax.transAxes, va="top", fontsize=6.5)
    return hb


def chrom_panel(ax, df, metric, ylabel, conds, chroms):
    off = {conds[0]: -0.14, conds[1]: 0.14} if len(conds) == 2 else {c: 0 for c in conds}
    for res, filled in ((8000, True), (32000, False)):
        for cond in conds:
            sub = df[(df.res == res) & (df.cond_a == cond) & df.same_condition
                     & (df.chrom != "all")].set_index("chrom").reindex(chroms)
            xs = np.arange(len(chroms)) + off[cond] + (0.0 if filled else 0.0)
            ax.scatter(xs, sub[metric], s=18, marker="o" if cond == conds[0] else "s",
                       facecolor=COLORS[cond] if filled else "white",
                       edgecolor=COLORS[cond] if not filled else EDGE,
                       linewidth=1.0 if not filled else 0.5, zorder=3)
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms)
    ax.set_ylabel(ylabel)
    vals = df[(df.chrom != "all") & df.same_condition][metric]
    lo, hi = np.nanmin(vals), min(np.nanmax(vals), 1.0)
    pad = max((hi - lo) * 0.15, 0.002)
    ax.set_ylim(lo - pad, min(hi + pad, 1.0005))
    ax.grid(axis="y", color="#e6e6e6", lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def heatmap_panel(ax, df, metric, res, names, title):
    sub = df[(df.res == res) & (df.chrom == "all")]
    mat = pd.DataFrame(np.nan, index=names, columns=names)
    for _, r in sub.iterrows():
        mat.loc[r.a, r.b] = mat.loc[r.b, r.a] = r[metric]
    arr = mat.to_numpy(copy=True)
    np.fill_diagonal(arr, 1.0)
    off = arr[~np.eye(len(names), dtype=bool)]
    vmin = max(0.0, np.nanmin(off) - 0.02)
    im = ax.imshow(arr, cmap="Blues", vmin=vmin, vmax=1.0)
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    short = [n.replace("JW18", "").replace("-", " rep") for n in names]
    ax.set_xticklabels(short, rotation=45, ha="right")
    ax.set_yticklabels(short)
    for i in range(len(names)):
        for j in range(len(names)):
            v = arr[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                    color="white" if (v - vmin) / (1 - vmin) > 0.6 else EDGE)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(title, loc="left", fontsize=7.5)
    return im


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mcool_dir", required=True)
    p.add_argument("--pairwise", required=True, help="*_pairwise.tsv from replicate_concordance.py")
    p.add_argument("--conditions", nargs="+", default=["JW18DOX", "JW18wMel"])
    p.add_argument("--scatter_res", type=int, default=32000)
    p.add_argument("--heatmap_res", type=int, default=8000)
    p.add_argument("--max_dist", type=int, default=2_000_000)
    p.add_argument("--output_prefix", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.pairwise, sep="\t")
    conds = args.conditions
    chroms = [c for c in df.chrom.unique() if c != "all"]
    names = [f"{c}-{r}" for c in conds for r in (1, 2)]

    fig, axes = plt.subplots(2, 3, figsize=(7.4, 5.2), constrained_layout=True)
    for ax, cond, letter in zip(axes[0][:2], conds, "AB"):
        scatter_panel(ax, args.mcool_dir, cond, args.scatter_res, args.max_dist, df)
        panel_letter(ax, letter)
    chrom_panel(axes[1][0], df, "scc", "SCC, replicate 1 vs 2", conds, chroms)
    chrom_panel(axes[1][1], df, "pearson_raw", "Pearson r, replicate 1 vs 2", conds, chroms)
    panel_letter(axes[1][0], "C")
    panel_letter(axes[1][1], "D")
    axes[1][0].set_title("Per chromosome", loc="left", fontsize=7.5)
    axes[1][1].set_title("Per chromosome", loc="left", fontsize=7.5)

    # legend for C/D: colour/shape = condition, fill = resolution
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", ls="", markerfacecolor=COLORS[conds[0]], markeredgecolor=EDGE,
               markersize=5, label=LABELS[conds[0]]),
        Line2D([], [], marker="s", ls="", markerfacecolor=COLORS[conds[1]], markeredgecolor=EDGE,
               markersize=5, label=LABELS[conds[1]]),
        Line2D([], [], marker="o", ls="", markerfacecolor=MUTED, markeredgecolor=MUTED,
               markersize=5, label="8 kb (filled)"),
        Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor=MUTED,
               markersize=5, label="32 kb (open)"),
    ]
    axes[1][1].legend(handles=handles, frameon=False, fontsize=6, loc="upper center",
                      bbox_to_anchor=(0.5, -0.12), ncol=2, handletextpad=0.3,
                      columnspacing=1.0)

    im1 = heatmap_panel(axes[0][2], df, "scc", args.heatmap_res, names,
                        f"SCC, all replicates, {args.heatmap_res // 1000} kb")
    im2 = heatmap_panel(axes[1][2], df, "pearson_raw", args.heatmap_res, names,
                        f"Pearson r, all replicates, {args.heatmap_res // 1000} kb")
    panel_letter(axes[0][2], "E")
    panel_letter(axes[1][2], "F")
    fig.colorbar(im1, ax=axes[0][2], fraction=0.046, pad=0.02, shrink=0.8)
    fig.colorbar(im2, ax=axes[1][2], fraction=0.046, pad=0.02, shrink=0.8)

    fig.savefig(f"{args.output_prefix}.pdf")
    fig.savefig(f"{args.output_prefix}.png", dpi=300)


if __name__ == "__main__":
    main()
