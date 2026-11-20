import cooler
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# =============================================================================
# Parameters
# =============================================================================

base_path = "/private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/processed_files"
save_path  = (
    "/private/groups/russelllab/jodie/wolbachia_induced_DE/"
    "Jacobs_et_al_2026_Wolbachia-host_genome_topology/dc_by_infection/results/contact_ratio"
)

RESOLUTION  = 1_000     # cooler bin resolution (bp)
WINDOW_SIZE = 100_000   # window size (bp); must be >> RESOLUTION; 100 bins @ 1 kb
STEP_SIZE   = 100_000   # step between windows (bp); set < WINDOW_SIZE for sliding overlap
SMOOTH_BINS = 5         # rolling-mean window (number of windows); 1 = no smoothing

X_CHROM   = 'X'
AUTOSOMES = ['2L', '2R', '3L', '3R']

files = {
    'JW18-DOX-1':  f"{base_path}/JW18-DOX-1.matrix_1kb.mcool::/resolutions/{RESOLUTION}",
    'JW18-DOX-2':  f"{base_path}/JW18-DOX-2.matrix_1kb.mcool::/resolutions/{RESOLUTION}",
    'JW18-wMel-1': f"{base_path}/JW18-wMel-1.matrix_1kb.mcool::/resolutions/{RESOLUTION}",
    'JW18-wMel-2': f"{base_path}/JW18-wMel-2.matrix_1kb.mcool::/resolutions/{RESOLUTION}",
}

CONDITIONS = {
    'JW18-DOX-1':  'Uninfected',
    'JW18-DOX-2':  'Uninfected',
    'JW18-wMel-1': 'Infected',
    'JW18-wMel-2': 'Infected',
}

COLORS = {'Uninfected': '#d62728', 'Infected': '#1f77b4'}

# Optional: mark loci of interest on the X chromosome (name -> coordinate in bp)
# e.g. Ndf locus.  Set to {} to skip.
LOCI_OF_INTEREST = {
    # 'Ndf': 8_400_000,   # update with correct dm6 coordinate
}

# =============================================================================
# Functions
# =============================================================================

def windowed_intra_density(clr, chrom, window_size, step_size):
    """
    Slide a window across `chrom` and compute intra-window contact density.

    Density for each window [start, end):
        density = sum(sub-matrix contacts) / window_length * 1e6   [contacts / Mb]

    Uses raw (unbalanced) counts for consistency with the whole-chromosome script.
    Windows shorter than half the target size (trailing edge) are dropped.
    Returns a DataFrame: chrom | start | end | mid | contacts | density
    """
    chrom_len = int(clr.chromsizes[chrom])
    records = []

    for start in range(0, chrom_len, step_size):
        end = min(start + window_size, chrom_len)
        if (end - start) < window_size // 2:          # skip tiny trailing window
            continue
        region = f'{chrom}:{start}-{end}'
        try:
            mat = clr.matrix(balance=False, sparse=True).fetch(region)
            contacts = float(mat.sum())
        except Exception as e:
            print(f"    Warning: could not fetch {region}: {e}")
            contacts = np.nan
        density = contacts / (end - start) * 1e6 if not np.isnan(contacts) else np.nan
        records.append({
            'chrom': chrom, 'start': start, 'end': end,
            'mid': (start + end) // 2,
            'contacts': contacts, 'density': density,
        })

    return pd.DataFrame(records)


def autosome_baseline(clr, autosomes, window_size, step_size):
    """
    Compute the genome-wide autosome contact density baseline for a sample.
    Mirrors the per-sample normalization in the original whole-chromosome script:
        baseline = mean(windowed densities across all autosome windows)
    Returns (mean, std).
    """
    all_densities = []
    for chrom in autosomes:
        df = windowed_intra_density(clr, chrom, window_size, step_size)
        all_densities.extend(df['density'].dropna().tolist())
    arr = np.array(all_densities)
    return np.nanmean(arr), np.nanstd(arr)


def rolling_smooth(series, n):
    """Center-aligned rolling mean; n=1 is a no-op."""
    return series.rolling(window=n, center=True, min_periods=1).mean()


# =============================================================================
# Main: compute windowed X / autosome ratios for each sample
# =============================================================================

all_results = []

for sample, filepath in files.items():
    print(f"\n{'='*60}")
    print(f"Processing {sample}  ({CONDITIONS[sample]})")
    clr = cooler.Cooler(filepath)
    print(f"  Chromosomes: {clr.chromnames}")

    # --- X chromosome windowed density ---
    print(f"  Computing X windows ({WINDOW_SIZE//1_000} kb, step {STEP_SIZE//1_000} kb)...")
    x_df = windowed_intra_density(clr, X_CHROM, WINDOW_SIZE, STEP_SIZE)
    print(f"  {len(x_df)} windows on X  |  "
          f"mean density: {x_df['density'].mean():.2f} contacts/Mb")

    # --- Autosome baseline ---
    print(f"  Computing autosome baseline across {AUTOSOMES}...")
    auto_mean, auto_sd = autosome_baseline(clr, AUTOSOMES, WINDOW_SIZE, STEP_SIZE)
    print(f"  Autosome baseline: {auto_mean:.2f} ± {auto_sd:.2f} contacts/Mb")

    # --- Ratio ---
    x_df['auto_baseline']  = auto_mean
    x_df['ratio']          = x_df['density'] / auto_mean
    x_df['sample']         = sample
    x_df['condition']      = CONDITIONS[sample]
    all_results.append(x_df)

# Combine and save raw per-window data
full_df = pd.concat(all_results, ignore_index=True)
out_csv = f"{save_path}/X_autosome_windowed_ratios.csv"
full_df.to_csv(out_csv, index=False)
print(f"\nPer-window data saved: {out_csv}")

# =============================================================================
# Aggregate across replicates per condition
# =============================================================================

agg_df = (
    full_df
    .groupby(['mid', 'condition'], as_index=False)
    .agg(
        mean_ratio   = ('ratio',   'mean'),
        sd_ratio     = ('ratio',   'std'),
        mean_density = ('density', 'mean'),
        sd_density   = ('density', 'std'),
        n_reps       = ('ratio',   'count'),
    )
    .sort_values(['condition', 'mid'])
)

# Fill NaN sd (single-rep windows) with 0
agg_df['sd_ratio']   = agg_df['sd_ratio'].fillna(0)
agg_df['sd_density'] = agg_df['sd_density'].fillna(0)

# Apply smoothing per condition
for cond in agg_df['condition'].unique():
    mask = agg_df['condition'] == cond
    agg_df.loc[mask, 'smooth_ratio']   = rolling_smooth(agg_df.loc[mask, 'mean_ratio'],   SMOOTH_BINS).values
    agg_df.loc[mask, 'smooth_density'] = rolling_smooth(agg_df.loc[mask, 'mean_density'], SMOOTH_BINS).values

# Fold change: pivot on condition
pivot = (
    agg_df
    .pivot(index='mid', columns='condition', values='smooth_ratio')
    .reset_index()
)
pivot.columns.name = None
pivot['log2_FC'] = np.log2(pivot['Infected'] / pivot['Uninfected'])
pivot = pivot.sort_values('mid')

# =============================================================================
# Global stats (mirrors original script output)
# =============================================================================

print(f"\n{'='*60}")
print("SUMMARY STATISTICS (per-window ratios pooled)")
print(f"{'='*60}")

for cond in ['Uninfected', 'Infected']:
    vals = full_df.loc[full_df['condition'] == cond, 'ratio'].dropna().values
    print(f"{cond}: n={len(vals)} windows  |  "
          f"mean={vals.mean():.4f}  |  sd={vals.std():.4f}  |  "
          f"median={np.median(vals):.4f}")

dox_vals  = full_df.loc[full_df['condition'] == 'Uninfected', 'ratio'].dropna().values
wmel_vals = full_df.loc[full_df['condition'] == 'Infected',   'ratio'].dropna().values
t_stat, t_pval = stats.ttest_ind(dox_vals, wmel_vals)
u_stat, u_pval = stats.mannwhitneyu(dox_vals, wmel_vals, alternative='two-sided')
print(f"\nFold change (Infected / Uninfected): {wmel_vals.mean() / dox_vals.mean():.4f}")
print(f"Percent change: {(wmel_vals.mean() - dox_vals.mean()) / dox_vals.mean() * 100:.2f}%")
print(f"T-test:          t = {t_stat:.4f}, P = {t_pval:.4e}")
print(f"Mann-Whitney U:  U = {u_stat:.4f}, P = {u_pval:.4e}")

# =============================================================================
# Plotting
# =============================================================================

fig, axes = plt.subplots(
    3, 1, figsize=(14, 12), sharex=True,
    gridspec_kw={'height_ratios': [2.5, 2, 1.5]},
)

# ── Panel 1: X / Auto ratio per condition ────────────────────────────────────
ax1 = axes[0]
for cond, grp in agg_df.groupby('condition'):
    grp = grp.sort_values('mid')
    x_pos = grp['mid'].values / 1e6
    ax1.plot(x_pos, grp['smooth_ratio'],  label=cond, color=COLORS[cond], linewidth=1.5, zorder=3)
    ax1.fill_between(
        x_pos,
        grp['mean_ratio'] - grp['sd_ratio'],
        grp['mean_ratio'] + grp['sd_ratio'],
        color=COLORS[cond], alpha=0.20, zorder=2,
    )

ax1.axhline(1.0, color='#888888', linestyle='--', linewidth=0.8, label='Autosome baseline')
ax1.set_ylabel('X / Autosome Contact Ratio', fontsize=11)
ax1.set_title(
    f'X Chromosome Contact Density Relative to Autosome Baseline\n'
    f'({WINDOW_SIZE//1_000} kb windows, {STEP_SIZE//1_000} kb step, '
    f'smoothed {SMOOTH_BINS}×)',
    fontsize=11,
)
ax1.legend(frameon=False, fontsize=9)

# ── Panel 2: log₂ fold change ────────────────────────────────────────────────
ax2 = axes[1]
x_fc = pivot['mid'].values / 1e6
fc   = pivot['log2_FC'].values

ax2.plot(x_fc, fc, color='#333333', linewidth=1.2, zorder=3)
ax2.axhline(0, color='#888888', linestyle='--', linewidth=0.8)
ax2.fill_between(x_fc, fc, 0, where=(fc >= 0), alpha=0.35, color=COLORS['Infected'],   label='Increased in wMel')
ax2.fill_between(x_fc, fc, 0, where=(fc  < 0), alpha=0.35, color=COLORS['Uninfected'], label='Decreased in wMel')
ax2.set_ylabel('log₂ FC (Infected / Uninfected)', fontsize=11)
ax2.set_title('Wolbachia Effect on X Chromosome Contact Density', fontsize=11)
ax2.legend(frameon=False, fontsize=9)

# ── Panel 3: raw contact density ─────────────────────────────────────────────
ax3 = axes[2]
for cond, grp in agg_df.groupby('condition'):
    grp = grp.sort_values('mid')
    ax3.plot(
        grp['mid'].values / 1e6, grp['smooth_density'],
        label=cond, color=COLORS[cond], linewidth=1.5,
    )
ax3.set_ylabel('Contact Density\n(contacts / Mb)', fontsize=11)
ax3.set_xlabel('X Chromosome Position (Mb)', fontsize=11)
ax3.set_title('Raw Contact Density along X', fontsize=11)
ax3.legend(frameon=False, fontsize=9)

# ── Mark loci of interest on all panels ──────────────────────────────────────
for locus_name, locus_coord in LOCI_OF_INTEREST.items():
    locus_mb = locus_coord / 1e6
    for ax in axes:
        ax.axvline(locus_mb, color='darkorange', linestyle=':', linewidth=1.2, alpha=0.8)
    axes[0].text(locus_mb, axes[0].get_ylim()[1], locus_name,
                 rotation=90, va='top', ha='right', fontsize=8, color='darkorange')

plt.tight_layout()
out_fig = f"{save_path}/X_autosome_windowed_ratio.pdf"
plt.savefig(out_fig, dpi=300, bbox_inches='tight')
print(f"\nFigure saved: {out_fig}")
