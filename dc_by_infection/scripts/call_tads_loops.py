#!/usr/bin/env python3
"""
Compare chromatin structure across conditions using existing diffHic results.
Uses null_model_results.csv for null model and all_results_combined.csv for interaction data.

MODIFIED: Now accepts variable number of conditions (not hardcoded to 4)
CORRECTED: Fixed compartment significance testing to use proper statistical tests
FIXED (2026-05): Eigenvectors are now phased against GC content so E1 > 0 == A
                 (gene-rich/active) consistently across chromosomes and conditions.
                 This makes A->B / B->A switches biologically interpretable rather
                 than dependent on the arbitrary sign of each eigendecomposition.
UPDATED: Renamed "TAD boundary" differential interaction logic to "hotspots".
UPDATED: Switched from binomial test to label-permutation null for switch-rate tests.
ADDED: Implemented direct cooltools.insulation derived TAD analysis with plotting.
ADDED: Compartment Saddle plots and compartmentalization strength (AA*BB)/(AB*BA).
       This metric quantifies global compartmentalization strength without per-bin testing.
       Method popularized for differential work by:
         - Schwarzer et al. 2017 (Nature 551:51–56)
         - Nora et al. 2017 (Cell 169(5):930–944)
         - cooltools (Venev et al.)
ADDED: save_tad_boundary_data() exports the raw per-boundary table that underlies
       Panel 1 of the TAD insulation summary PDF as a tab-delimited .txt file
       with the same name stem as the PDF.
"""
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import cooler
import cooltools
from cooltools import expected_cis
from cooltools import dots
from cooltools import eigs_cis

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns
from multiprocessing import Pool
import bioframe


def get_closest_resolution(mcool_file, target_resolution):
    """Get the closest available resolution in the mcool file."""
    try:
        f = cooler.fileops.list_coolers(mcool_file)
        available_resolutions = [int(path.split('/')[-1]) for path in f]
        available_resolutions.sort()
        closest = min(available_resolutions, key=lambda x: abs(x - target_resolution))
        print(f"Target resolution: {target_resolution}, using: {closest}")
        return closest
    except Exception as e:
        print(f"Error getting resolutions from {mcool_file}: {e}")
        return None


def load_diffhic_results(results_file):
    """Load the combined diffHic results."""
    print(f"Loading diffHic results from {results_file}...")
    results_df = pd.read_csv(results_file)

    significant_results = results_df[
        (results_df['FDR'] < 0.05) &
        (abs(results_df['logFC']) > 1)
    ].copy()

    print(f"Loaded {len(results_df)} total interactions, {len(significant_results)} significant")
    return results_df, significant_results


def load_null_model(null_file):
    """Load null model results from diffHic analysis"""
    print(f"Loading null model from {null_file}...")
    null_df = pd.read_csv(null_file)
    print(f"Loaded null model with {len(null_df)} entries")
    return null_df


def identify_hotspots_from_interactions(significant_interactions, conditions, window_size=50000):
    """
    Identify differential-interaction hotspots.
    Hotspots are regions with enriched differential interactions.
    """
    print("Identifying differential-interaction hotspots from interactions...")

    all_hotspots = {}

    condition_col = None
    for possible_col in ['infection', 'condition', 'comparison', 'contrast', 'group']:
        if possible_col in significant_interactions.columns:
            condition_col = possible_col
            print(f"  Found condition column: {condition_col}")
            break

    if condition_col is not None:
        groups_to_process = significant_interactions[condition_col].unique()
        print(f"  Processing {len(groups_to_process)} groups: {groups_to_process}")

        for group_name in groups_to_process:
            if pd.isna(group_name):
                continue

            print(f"  Processing {group_name}...")
            group_data = significant_interactions[
                significant_interactions[condition_col] == group_name
            ].copy()

            hotspots = process_interactions_for_hotspots(group_data, window_size)

            if hotspots:
                hotspots_df = pd.DataFrame(hotspots)
                clean_name = str(group_name).replace('infectionJW18', '').replace('infection', '')
                hotspots_df['condition'] = clean_name
                all_hotspots[clean_name] = hotspots_df
                print(f"  Found {len(hotspots_df)} hotspots for {clean_name}")
    else:
        print(f"  No condition column found. Processing all interactions together.")
        print(f"  Will assign hotspots to conditions based on logFC direction")

        for condition in conditions:
            print(f"  Creating hotspots for {condition}...")
            if condition in ['DOX', 'uninfected']:
                condition_data = significant_interactions.copy()
            else:
                condition_data = significant_interactions[
                    significant_interactions['logFC'] > 0
                ].copy()

            hotspots = process_interactions_for_hotspots(condition_data, window_size)

            if hotspots:
                hotspots_df = pd.DataFrame(hotspots)
                hotspots_df['condition'] = condition
                all_hotspots[condition] = hotspots_df
                print(f"  Found {len(hotspots_df)} hotspots for {condition}")

    return all_hotspots


def process_interactions_for_hotspots(interaction_data, window_size=50000):
    """Helper function to process interactions and identify hotspots."""
    hotspots = []

    for chrom in interaction_data['chr1'].unique():
        chrom_interactions = interaction_data[
            (interaction_data['chr1'] == chrom) &
            (interaction_data['chr2'] == chrom)
        ].copy()

        if len(chrom_interactions) == 0:
            continue

        max_pos = max(
            chrom_interactions['end1'].max(),
            chrom_interactions['end2'].max()
        )
        bins = np.arange(0, max_pos + window_size, window_size)
        bin_counts = np.zeros(len(bins) - 1)

        for _, interaction in chrom_interactions.iterrows():
            start_bin = np.searchsorted(bins, interaction['start1'], side='right') - 1
            end_bin = np.searchsorted(bins, interaction['end2'], side='right') - 1
            for b in range(max(0, start_bin), min(len(bin_counts), end_bin + 1)):
                bin_counts[b] += abs(interaction['logFC'])

        from scipy.signal import find_peaks
        peaks, properties = find_peaks(
            bin_counts,
            height=np.percentile(bin_counts[bin_counts > 0], 75) if np.any(bin_counts > 0) else 0,
            distance=2
        )

        for peak in peaks:
            if peak < len(bins) - 1:
                hotspot_start = bins[peak]
                hotspot_end = bins[peak + 1]

                hotspots.append({
                    'chrom': chrom,
                    'start': hotspot_start,
                    'end': hotspot_end,
                    'hotspot_strength': bin_counts[peak],
                    'n_interactions': np.sum(
                        (chrom_interactions['start1'] >= hotspot_start) &
                        (chrom_interactions['end2'] <= hotspot_end)
                    )
                })

    return hotspots


def compare_hotspots(hotspots_data, null_model, fdr_threshold=0.05):
    """Compare hotspots between conditions using the null model."""
    print("Comparing hotspots between conditions...")

    if not hotspots_data:
        print("No hotspots to compare")
        return pd.DataFrame()

    ref_condition = None
    for condition in ['DOX', 'uninfected']:
        if condition in hotspots_data:
            ref_condition = condition
            break

    if ref_condition is None:
        print("Warning: No reference condition found")
        return pd.DataFrame()

    ref_hotspots = hotspots_data[ref_condition]
    comparison_results = []

    for infected_condition in hotspots_data.keys():
        if infected_condition == ref_condition:
            continue

        infected_hotspots = hotspots_data[infected_condition]
        print(f"  Comparing {ref_condition} ({len(ref_hotspots)}) vs {infected_condition} ({len(infected_hotspots)})")

        overlaps = []
        for _, ref_hotspot in ref_hotspots.iterrows():
            chrom_hotspots = infected_hotspots[
                infected_hotspots['chrom'] == ref_hotspot['chrom']
            ]

            for _, inf_hotspot in chrom_hotspots.iterrows():
                overlap_start = max(ref_hotspot['start'], inf_hotspot['start'])
                overlap_end = min(ref_hotspot['end'], inf_hotspot['end'])

                if overlap_start < overlap_end:
                    overlap_length = overlap_end - overlap_start
                    ref_length = ref_hotspot['end'] - ref_hotspot['start']
                    inf_length = inf_hotspot['end'] - inf_hotspot['start']

                    if (overlap_length / ref_length >= 0.5 and
                            overlap_length / inf_length >= 0.5):

                        strength_change = (inf_hotspot['hotspot_strength'] -
                                           ref_hotspot['hotspot_strength'])

                        overlaps.append({
                            'chrom': ref_hotspot['chrom'],
                            'ref_start': ref_hotspot['start'],
                            'ref_end': ref_hotspot['end'],
                            'inf_start': inf_hotspot['start'],
                            'inf_end': inf_hotspot['end'],
                            'ref_strength': ref_hotspot['hotspot_strength'],
                            'inf_strength': inf_hotspot['hotspot_strength'],
                            'strength_change': strength_change,
                            'log2_fold_change': np.log2(
                                (inf_hotspot['hotspot_strength'] + 1) /
                                (ref_hotspot['hotspot_strength'] + 1)
                            )
                        })

        if not overlaps:
            print(f"    No overlapping hotspots found")
            continue

        overlaps_df = pd.DataFrame(overlaps)

        n_permutations = 1000
        null_changes = []
        for _ in range(n_permutations):
            perm_strengths = np.random.permutation(overlaps_df['inf_strength'])
            null_change = np.mean(perm_strengths - overlaps_df['ref_strength'])
            null_changes.append(null_change)

        null_changes = np.array(null_changes)
        observed_change = overlaps_df['strength_change'].mean()
        p_value = np.sum(np.abs(null_changes) >= np.abs(observed_change)) / n_permutations

        hotspot_p_values = []
        for _, row in overlaps_df.iterrows():
            null_dist = np.random.normal(0, null_changes.std(), 1000)
            p_val = np.sum(np.abs(null_dist) >= np.abs(row['strength_change'])) / 1000
            hotspot_p_values.append(p_val)

        overlaps_df['p_value'] = hotspot_p_values

        if len(hotspot_p_values) > 0:
            _, overlaps_df['fdr'], _, _ = multipletests(hotspot_p_values, method='fdr_bh')
            overlaps_df['significant'] = overlaps_df['fdr'] < fdr_threshold
        else:
            overlaps_df['fdr'] = []
            overlaps_df['significant'] = []

        overlaps_df['comparison'] = f"{ref_condition}_vs_{infected_condition}"
        overlaps_df['global_p_value'] = p_value

        comparison_results.append(overlaps_df)

        print(f"    Found {len(overlaps_df)} overlapping hotspots, "
              f"{overlaps_df['significant'].sum() if len(overlaps_df) > 0 else 0} significant")

    if comparison_results:
        return pd.concat(comparison_results, ignore_index=True)
    else:
        return pd.DataFrame()


def calculate_cooltools_tads(mcool_files, conditions, resolution=50000, window=150000):
    """
    Calculate TADs (insulation scores and boundaries) using cooltools.
    Auto-corrects window sizes to ensure they are exact multiples of the matrix resolution.
    """
    print("\nCalculating cooltools insulation scores and boundaries...")
    tads_data = {}

    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"  Processing {condition}...")
        actual_res = get_closest_resolution(mcool_file, resolution)
        if actual_res is None:
            continue

        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_res}")

        # CRITICAL FIX: cooltools requires the window to be an EXACT multiple of the bin size.
        bin_multiplier = max(3, round(window / actual_res))  # Force at least a 3-bin window
        actual_window = actual_res * bin_multiplier

        try:
            print(f"    Using resolution: {actual_res} bp, actual window: {actual_window} bp ({bin_multiplier} bins)")
            ins_df = cooltools.insulation(clr, [actual_window], nproc=4)

            # Standardize column names back to what the plotting/saving functions expect
            ins_df = ins_df.rename(columns={
                f'is_boundary_{actual_window}': f'is_boundary_{window}',
                f'boundary_strength_{actual_window}': f'boundary_strength_{window}',
                f'log2_insulation_score_{actual_window}': f'log2_insulation_score_{window}'
            })

            tads_data[condition] = ins_df
            boundary_col = f'is_boundary_{window}'
            n_boundaries = ins_df[boundary_col].sum() if boundary_col in ins_df.columns else 0
            print(f"    Found {n_boundaries} TAD boundaries for {condition}")

        except Exception as e:
            print(f"    Warning: Failed to calculate insulation for {condition}: {e}")
            import traceback
            traceback.print_exc()

    return tads_data


def save_tad_boundary_data(tads_data, window, output_prefix):
    """
    Write the raw per-boundary table that underlies Panel 1 of the TAD
    insulation summary PDF.

    Column definitions
    ------------------
    condition               sample label (e.g. wMel, DOX)
    chrom / start / end     genomic coordinates of the boundary bin
    boundary_strength       delta-insulation score at this boundary
                            (how sharply the score dips; higher = sharper)
    log2_insulation_score   raw log2 insulation score for this bin
    boundary_class          coarse classification for quick filtering —
                              shared   : boundary present in both this condition
                                         and the reference (within window/2 bp)
                              unique   : boundary present in only one condition
                                         (gained or lost relative to reference)
    direction               finer detail on the relationship to the reference —
                              present_in_{ref}       reference-condition boundary
                              {cond}_only            gained (not in reference)
                              shared_stronger_{X}    shared; stronger in condition X
                              shared_equal           shared; equal strength
    reference_condition     which sample is used as the reference (DOX or
                            uninfected if available, otherwise first condition)
    ref_boundary_strength   boundary strength at this locus in the reference
                            (NaN for unique boundaries)
    strength_diff_vs_ref    boundary_strength − ref_boundary_strength
                            (positive = stronger than reference; NaN for unique)
    p_value                 two-sided z-score p-value of this boundary's strength
                            vs. the all-bin boundary_strength distribution of the
                            reference condition (genomic background null).
                            Answers: "is this a significant boundary?"
                            Consistent across all rows including reference.
    diff_p_value            two-sided p-value testing whether strength_diff_vs_ref
                            is significantly different from zero, using an empirical
                            null: std of all shared-boundary strength differences
                            genome-wide (null centred at 0, spread = observed noise).
                            Answers: "does this boundary change significantly between
                            conditions?"
                            NaN for unique boundaries (categorically differential
                            by definition) and for reference-condition rows.
    p_value_note            plain-English reminder of the statistical approach.

    Output file: {output_prefix}_tad_insulation_summary.txt  (TSV)
    """
    if not tads_data:
        print("  save_tad_boundary_data: no TAD data — skipping.")
        return

    boundary_col   = f'is_boundary_{window}'
    strength_col   = f'boundary_strength_{window}'
    insulation_col = f'log2_insulation_score_{window}'

    # ── choose reference condition ────────────────────────────────────────────
    ref_condition = None
    for cand in ['DOX', 'uninfected']:
        if cand in tads_data:
            ref_condition = cand
            break
    if ref_condition is None:
        ref_condition = next(iter(tads_data))
    print(f"\n  save_tad_boundary_data: reference = '{ref_condition}'")

    # ── extract flagged boundaries per condition ──────────────────────────────
    condition_boundaries = {}
    for cond, df in tads_data.items():
        if boundary_col in df.columns:
            condition_boundaries[cond] = df[df[boundary_col]].copy().reset_index(drop=True)
        else:
            condition_boundaries[cond] = pd.DataFrame()
            print(f"  Warning: '{boundary_col}' not found for {cond}; no boundaries extracted.")

    ref_df = condition_boundaries.get(ref_condition, pd.DataFrame())

    # ── reference background (ALL bins) for p_value z-score null ─────────────
    # All bins (not just flagged boundaries) so the null represents the full
    # genomic distribution of boundary-strength values.
    ref_all_df = tads_data.get(ref_condition, pd.DataFrame())
    if strength_col in ref_all_df.columns:
        ref_bg = ref_all_df[strength_col].dropna().values
    else:
        ref_bg = np.array([])
    ref_bg_mean = float(np.nanmean(ref_bg)) if len(ref_bg) > 0 else 0.0
    ref_bg_std  = float(np.nanstd(ref_bg))  if len(ref_bg) > 0 else 1.0
    if ref_bg_std == 0.0:
        ref_bg_std = 1.0

    # Two boundaries overlap when their centres are within half a window.
    overlap_bp = window // 2

    records = []

    for cond, boundaries in condition_boundaries.items():
        if boundaries.empty:
            continue

        for _, row in boundaries.iterrows():
            chrom    = row['chrom']
            start    = int(row['start'])
            end      = int(row['end'])
            strength = (float(row[strength_col])
                        if strength_col in row.index and not pd.isna(row[strength_col])
                        else np.nan)
            insul    = (float(row[insulation_col])
                        if insulation_col in row.index and not pd.isna(row[insulation_col])
                        else np.nan)

            # p_value: z-score of this boundary's strength vs. reference background
            if not np.isnan(strength):
                z_bg    = (strength - ref_bg_mean) / ref_bg_std
                p_value = float(2.0 * (1.0 - stats.norm.cdf(abs(z_bg))))
            else:
                p_value = np.nan

            # ── reference-condition rows ──────────────────────────────────────
            if cond == ref_condition:
                records.append(dict(
                    condition             = cond,
                    chrom                 = chrom,
                    start                 = start,
                    end                   = end,
                    boundary_strength     = round(strength, 6) if not np.isnan(strength) else np.nan,
                    log2_insulation_score = round(insul, 6)    if not np.isnan(insul)    else np.nan,
                    boundary_class        = 'reference',
                    direction             = f'present_in_{ref_condition}',
                    reference_condition   = ref_condition,
                    ref_boundary_strength = np.nan,
                    strength_diff_vs_ref  = np.nan,
                    p_value               = round(p_value, 6) if not np.isnan(p_value) else np.nan,
                    diff_p_value          = np.nan,
                    p_value_note          = ('p_value: two-sided z-score vs. all-bin background '
                                             'of reference (self); diff_p_value: N/A for '
                                             'reference condition'),
                ))
                continue

            # ── non-reference rows: find overlapping reference boundary ────────
            if not ref_df.empty and 'chrom' in ref_df.columns:
                ref_chrom_df = ref_df[ref_df['chrom'] == chrom]
            else:
                ref_chrom_df = pd.DataFrame()

            overlapping = pd.DataFrame()
            if not ref_chrom_df.empty and strength_col in ref_chrom_df.columns:
                mask = (
                    (ref_chrom_df['start'] <= end   + overlap_bp) &
                    (ref_chrom_df['end']   >= start - overlap_bp)
                )
                overlapping = ref_chrom_df[mask]

            if overlapping.empty:
                # Boundary unique to this condition (gained relative to reference)
                boundary_class = 'unique'
                direction      = f'{cond}_only'
                ref_strength   = np.nan
                strength_diff  = np.nan

            else:
                # Closest overlapping reference boundary by centre distance
                centres  = (overlapping['start'] + overlapping['end']) / 2.0
                this_ctr = (start + end) / 2.0
                best_idx = (centres - this_ctr).abs().idxmin()
                best_row = overlapping.loc[best_idx]

                ref_strength = (float(best_row[strength_col])
                                if strength_col in best_row.index
                                and not pd.isna(best_row[strength_col])
                                else np.nan)

                boundary_class = 'shared'

                if np.isnan(strength) or np.isnan(ref_strength):
                    strength_diff = np.nan
                    direction     = 'shared_equal'
                else:
                    strength_diff = strength - ref_strength
                    if   strength_diff >  1e-9:
                        direction = f'shared_stronger_{cond}'
                    elif strength_diff < -1e-9:
                        direction = f'shared_stronger_{ref_condition}'
                    else:
                        direction = 'shared_equal'

            records.append(dict(
                condition             = cond,
                chrom                 = chrom,
                start                 = start,
                end                   = end,
                boundary_strength     = round(strength, 6)      if not np.isnan(strength)      else np.nan,
                log2_insulation_score = round(insul, 6)         if not np.isnan(insul)          else np.nan,
                boundary_class        = boundary_class,
                direction             = direction,
                reference_condition   = ref_condition,
                ref_boundary_strength = round(ref_strength, 6)  if not np.isnan(ref_strength)  else np.nan,
                strength_diff_vs_ref  = round(strength_diff, 6) if not np.isnan(strength_diff) else np.nan,
                p_value               = round(p_value, 6)       if not np.isnan(p_value)       else np.nan,
                diff_p_value          = np.nan,   # filled in below for shared boundaries
                p_value_note          = ('p_value: two-sided z-score vs. all-bin background of '
                                         'reference; diff_p_value: two-sided z-score of '
                                         'strength_diff vs. empirical null (std of all shared '
                                         'diffs, centred at 0); both descriptive — no replicates'),
            ))

    if not records:
        print("  save_tad_boundary_data: no boundary records to write.")
        return

    out_df = (pd.DataFrame(records)
              .sort_values(['condition', 'chrom', 'start'])
              .reset_index(drop=True))

    # ── diff_p_value: empirical null on strength_diff_vs_ref ─────────────────
    # Null: strength differences at shared boundaries are centred at 0;
    # spread estimated from the observed distribution of all shared diffs.
    # This tests "does this boundary change more than typical shared boundaries?"
    # rather than assuming a theoretical distribution.
    shared_mask  = out_df['boundary_class'] == 'shared'
    shared_diffs = out_df.loc[shared_mask, 'strength_diff_vs_ref'].dropna()

    if len(shared_diffs) > 1:
        diff_null_std = float(shared_diffs.std())
        if diff_null_std > 0:
            z_diff = out_df.loc[shared_mask, 'strength_diff_vs_ref'] / diff_null_std
            out_df.loc[shared_mask, 'diff_p_value'] = (
                2.0 * (1.0 - stats.norm.cdf(z_diff.abs()))
            ).round(6)
            print(f"  diff_p_value null: std of {len(shared_diffs)} shared boundary "
                  f"differences = {diff_null_std:.4f}")
        else:
            print("  Warning: std of shared diffs is 0; diff_p_value not computed.")
    else:
        print(f"  Warning: only {len(shared_diffs)} shared boundaries — "
              f"diff_p_value empirical null unreliable, leaving as NaN.")

    # ── FDR correction on diff_p_value across shared boundaries ──────────────
    diff_pvals = out_df.loc[shared_mask, 'diff_p_value'].dropna()
    if len(diff_pvals) > 0:
        valid_idx = diff_pvals.index
        _, fdr_vals, _, _ = multipletests(diff_pvals.values, method='fdr_bh')
        out_df.loc[valid_idx, 'diff_fdr'] = fdr_vals.round(6)
    else:
        out_df['diff_fdr'] = np.nan

    out_file = f"{output_prefix}_tad_insulation_summary.txt"
    out_df.to_csv(out_file, sep='\t', index=False)
    print(f"  Saved TAD boundary table → {out_file}")
    print(f"  Rows: {len(out_df)} | Columns: {list(out_df.columns)}")
    for cond, grp in out_df.groupby('condition'):
        class_counts = grp['boundary_class'].value_counts().to_dict()
        dir_counts   = grp['direction'].value_counts().to_dict()
        print(f"    {cond}: {len(grp)} boundaries | class: {class_counts} | direction: {dir_counts}")


def create_tad_venn_diagrams(tads_data, window, chromosomes, output_prefix):
    """
    Per-chromosome Venn diagrams of TAD boundary overlap between the reference
    condition and each other condition.

    Layout: summary panel (all chromosomes) + one panel per chromosome,
    arranged in a 2 × 4 grid.  Ellipse patches in normalised [0,1] axes
    coordinates — no equal-aspect requirement, no external venn library needed.

    Outputs
    -------
    {output_prefix}_tad_venn_{ref}_vs_{cond}.pdf
    {output_prefix}_tad_venn_{ref}_vs_{cond}_stats.tsv
        Columns: chrom, comparison, {ref}_only, shared, {cond}_only,
                 total_{ref}, total_{cond}, pct_shared_of_ref, pct_shared_of_cond
        Row 'ALL' = genome-wide totals.
    """
    if not tads_data or len(tads_data) < 2:
        print("  create_tad_venn_diagrams: need >=2 conditions -- skipping.")
        return

    boundary_col = f'is_boundary_{window}'
    overlap_bp   = window // 2

    # -- reference condition --------------------------------------------------
    ref_condition = None
    for cand in ['DOX', 'uninfected']:
        if cand in tads_data:
            ref_condition = cand
            break
    if ref_condition is None:
        ref_condition = next(iter(tads_data))

    non_ref = [c for c in tads_data if c != ref_condition]
    if not non_ref:
        print("  create_tad_venn_diagrams: no non-reference conditions -- skipping.")
        return

    # -- chromosome order -----------------------------------------------------
    chrom_order = ['2L', '2R', '3L', '3R', 'X', '4']
    all_chroms  = set()
    for df in tads_data.values():
        if 'chrom' in df.columns:
            all_chroms.update(df['chrom'].unique())
    plot_chroms = [c for c in chrom_order if c in all_chroms]
    plot_chroms += sorted(c for c in all_chroms
                          if c not in chrom_order and c in chromosomes)

    # -- helpers --------------------------------------------------------------
    def get_boundary_centres(df, chrom=None):
        if boundary_col not in df.columns:
            return []
        # fillna(False): cooltools sets is_boundary NaN for low-coverage bins;
        # boolean indexing raises ValueError on NaN without this guard.
        mask = df[boundary_col].fillna(False).astype(bool)
        sub = df[mask].copy()
        if chrom is not None:
            sub = sub[sub['chrom'] == chrom]
        return list(zip(sub['chrom'],
                        ((sub['start'] + sub['end']) / 2).astype(int)))

    def count_overlap(ref_centres, cond_centres):
        """Greedy nearest-neighbour matching within overlap_bp. O(n.m) -- fine at TAD scale."""
        matched_ref  = set()
        matched_cond = set()
        for ri, (r_chr, r_ctr) in enumerate(ref_centres):
            best_dist, best_ci = overlap_bp + 1, None
            for ci, (c_chr, c_ctr) in enumerate(cond_centres):
                if c_chr != r_chr or ci in matched_cond:
                    continue
                dist = abs(r_ctr - c_ctr)
                if dist <= overlap_bp and dist < best_dist:
                    best_dist, best_ci = dist, ci
            if best_ci is not None:
                matched_ref.add(ri)
                matched_cond.add(best_ci)
        shared    = len(matched_ref)
        ref_only  = len(ref_centres)  - shared
        cond_only = len(cond_centres) - len(matched_cond)
        return shared, ref_only, cond_only

    def draw_venn2(ax, left_n, shared_n, right_n,
                   left_label, right_label, title,
                   left_color='#3498db', right_color='#e74c3c'):
        """
        Two-circle Venn on ax.
        Uses Ellipse patches in normalised [0,1] data coordinates.
        No set_aspect('equal') -- safe inside any subplot grid.
        """
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')

        for cx, fc in [(0.35, left_color), (0.65, right_color)]:
            ax.add_patch(Ellipse(
                xy=(cx, 0.47), width=0.52, height=0.60,
                facecolor=fc, edgecolor=fc, alpha=0.38, linewidth=1.5
            ))

        # counts
        for x, n in [(0.15, left_n), (0.50, shared_n), (0.85, right_n)]:
            ax.text(x, 0.47, str(n),
                    ha='center', va='center', fontsize=13,
                    fontweight='bold', color='#1a252f')

        # condition labels
        ax.text(0.24, 0.87, left_label,  ha='center', va='center',
                fontsize=8, fontweight='bold', color=left_color)
        ax.text(0.76, 0.87, right_label, ha='center', va='center',
                fontsize=8, fontweight='bold', color=right_color)

        # role labels under counts
        ax.text(0.15, 0.10, 'lost',   ha='center', va='center',
                fontsize=7, color=left_color,  style='italic')
        ax.text(0.50, 0.10, 'shared', ha='center', va='center',
                fontsize=7, color='#555555',   style='italic')
        ax.text(0.85, 0.10, 'gained', ha='center', va='center',
                fontsize=7, color=right_color, style='italic')

        total = left_n + shared_n + right_n
        pct   = f'({shared_n / total * 100:.0f}% shared)' if total > 0 else ''
        ax.set_title(f'{title}\n{pct}', fontsize=9, pad=3)

    error_log = f"{output_prefix}_tad_venn_errors.txt"
    print(f"\n  create_tad_venn_diagrams: ref='{ref_condition}', "
          f"non_ref={non_ref}, plot_chroms={plot_chroms}")
    print(f"  boundary_col='{boundary_col}', overlap_bp={overlap_bp}")
    for cond_check, df_check in tads_data.items():
        has_col = boundary_col in df_check.columns
        if has_col:
            n_true = df_check[boundary_col].fillna(False).astype(bool).sum()
            print(f"  {cond_check}: {boundary_col} present, {n_true} boundaries called")
        else:
            print(f"  {cond_check}: {boundary_col} MISSING -- cols: {list(df_check.columns)}")

    # -- one figure per non-reference condition --------------------------------
    for cond in non_ref:
        print(f"\n  TAD Venn: {ref_condition} vs {cond}")
        try:
            ref_df  = tads_data[ref_condition]
            cond_df = tads_data[cond]

            n_panels = 1 + len(plot_chroms)
            n_cols   = 4
            n_rows   = -(-n_panels // n_cols)  # ceiling division

            fig, axes = plt.subplots(n_rows, n_cols,
                                      figsize=(n_cols * 3.5, n_rows * 3.5))
            # reshape to 1-D regardless of n_rows
            axes_flat = np.array(axes).reshape(-1)

            # Panel 0: all chromosomes combined
            ref_all  = get_boundary_centres(ref_df)
            cond_all = get_boundary_centres(cond_df)
            sh_all, ref_only_all, cond_only_all = count_overlap(ref_all, cond_all)

            draw_venn2(axes_flat[0],
                       ref_only_all, sh_all, cond_only_all,
                       ref_condition, cond, 'All Chromosomes')

            print(f"    ALL : {ref_only_all} lost | {sh_all} shared | {cond_only_all} gained")

            stats_rows = [{
                'chrom':                    'ALL',
                'comparison':               f'{ref_condition}_vs_{cond}',
                f'{ref_condition}_only':    ref_only_all,
                'shared':                   sh_all,
                f'{cond}_only':             cond_only_all,
                f'total_{ref_condition}':   len(ref_all),
                f'total_{cond}':            len(cond_all),
                'pct_shared_of_ref':
                    round(sh_all / len(ref_all)  * 100, 1) if ref_all  else 0.0,
                'pct_shared_of_cond':
                    round(sh_all / len(cond_all) * 100, 1) if cond_all else 0.0,
            }]

            # Per-chromosome panels
            for i, chrom in enumerate(plot_chroms):
                ax_idx = i + 1
                if ax_idx >= len(axes_flat):
                    break

                ref_c  = get_boundary_centres(ref_df,  chrom)
                cond_c = get_boundary_centres(cond_df, chrom)
                shared, ref_only, cond_only = count_overlap(ref_c, cond_c)

                draw_venn2(axes_flat[ax_idx],
                           ref_only, shared, cond_only,
                           ref_condition, cond, f'Chr {chrom}')

                print(f"    Chr {chrom}: {ref_only} lost | {shared} shared | {cond_only} gained")

                stats_rows.append({
                    'chrom':                    chrom,
                    'comparison':               f'{ref_condition}_vs_{cond}',
                    f'{ref_condition}_only':    ref_only,
                    'shared':                   shared,
                    f'{cond}_only':             cond_only,
                    f'total_{ref_condition}':   len(ref_c),
                    f'total_{cond}':            len(cond_c),
                    'pct_shared_of_ref':
                        round(shared / len(ref_c)  * 100, 1) if ref_c  else 0.0,
                    'pct_shared_of_cond':
                        round(shared / len(cond_c) * 100, 1) if cond_c else 0.0,
                })

            # hide unused axes
            for j in range(len(plot_chroms) + 1, len(axes_flat)):
                axes_flat[j].axis('off')

            plt.suptitle(f'TAD Boundary Overlap:  {ref_condition}  vs  {cond}',
                         fontsize=13, fontweight='bold')
            plt.tight_layout()

            out_pdf = f"{output_prefix}_tad_venn_{ref_condition}_vs_{cond}.pdf"
            plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"  Saved: {out_pdf}")

            # stats TSV
            stats_df = pd.DataFrame(stats_rows)
            out_tsv  = f"{output_prefix}_tad_venn_{ref_condition}_vs_{cond}_stats.tsv"
            stats_df.to_csv(out_tsv, sep='\t', index=False)
            print(f"  Saved: {out_tsv}")

            # sentence fill-ins
            print(f"\n  === Sentence fill-ins: {ref_condition} vs {cond} ===")
            if ref_all:
                pct_ref  = sh_all / len(ref_all)  * 100
                pct_cond = sh_all / len(cond_all) * 100 if cond_all else 0
                print(f"  Overlap: {sh_all} shared of {len(ref_all)} {ref_condition} "
                      f"boundaries ({pct_ref:.0f}% of ref, {pct_cond:.0f}% of {cond})")
            for label, key in [('Lost  (ref-only)', f'{ref_condition}_only'),
                                ('Gained (cond-only)', f'{cond}_only')]:
                hits = [(r['chrom'], r[key]) for r in stats_rows[1:] if r.get(key, 0) > 0]
                if hits:
                    print(f"  {label}: " + ", ".join(f"Chr {c}: {n}" for c, n in hits))

        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            msg = f"ERROR in create_tad_venn_diagrams for {cond}: {e}\n{tb}"
            print(msg)
            # also write to a file so errors are visible in SLURM runs
            with open(error_log, 'a') as _ef:
                _ef.write(msg + "\n")
            plt.close('all')

def create_tad_insulation_plots(tads_data, window, output_prefix):
    """
    Create summary visualizations for Cooltools TADs/Insulation.
    Generates three panels:
    1. Total TAD boundaries.
    2. KDE of boundary strengths.
    3. KDE of global insulation scores.

    Also writes a tab-delimited raw-data file backing Panel 1:
        {output_prefix}_tad_insulation_summary.txt
    """
    if not tads_data:
        print("No TAD insulation data to plot")
        return

    print("\nCreating TAD insulation plots...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    boundary_col   = f'is_boundary_{window}'
    strength_col   = f'boundary_strength_{window}'
    insulation_col = f'log2_insulation_score_{window}'

    conditions      = list(tads_data.keys())
    boundary_counts = []

    df_list_str = []
    df_list_ins = []

    for condition in conditions:
        df = tads_data[condition]

        # 1. Boundary counts
        if boundary_col in df.columns:
            count = df[boundary_col].sum()
            boundary_counts.append(count)
        else:
            boundary_counts.append(0)

        # 2. Boundary strengths setup for KDE
        if strength_col in df.columns and boundary_col in df.columns:
            tmp_str = df[df[boundary_col]][[strength_col]].copy()
            tmp_str.rename(columns={strength_col: 'value'}, inplace=True)
            tmp_str['condition'] = condition
            df_list_str.append(tmp_str)

        # 3. Global insulation setup for KDE
        if insulation_col in df.columns:
            tmp_ins = df[[insulation_col]].dropna().copy()
            if len(tmp_ins) > 50000:
                tmp_ins = tmp_ins.sample(50000, random_state=42)
            tmp_ins.rename(columns={insulation_col: 'value'}, inplace=True)
            tmp_ins['condition'] = condition
            df_list_ins.append(tmp_ins)

    # Panel 1: Total Boundaries
    ax = axes[0]
    ax.bar(conditions, boundary_counts, color='#3498db', edgecolor='black', alpha=0.8)
    ax.set_title('Total TAD Boundaries')
    ax.set_ylabel('Count')
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=45)
    for i, count in enumerate(boundary_counts):
        ax.text(i, count, str(count), ha='center', va='bottom', fontweight='bold')

    # Panel 2: Boundary Strengths KDE
    ax = axes[1]
    if df_list_str:
        df_str_all = pd.concat(df_list_str, ignore_index=True)
        sns.kdeplot(data=df_str_all, x='value', hue='condition', ax=ax, common_norm=False, fill=True, alpha=0.3)
        ax.set_title('Distribution of Boundary Strengths')
        ax.set_xlabel('Boundary Strength')
        ax.set_ylabel('Density')
    else:
        ax.text(0.5, 0.5, "No boundary strength data", ha='center', va='center')

    # Panel 3: Global Insulation KDE
    ax = axes[2]
    if df_list_ins:
        df_ins_all = pd.concat(df_list_ins, ignore_index=True)
        sns.kdeplot(data=df_ins_all, x='value', hue='condition', ax=ax, common_norm=False, fill=True, alpha=0.3)
        ax.set_title('Global Insulation Scores')
        ax.set_xlabel('Log2 Insulation Score')
        ax.set_ylabel('Density')
    else:
        ax.text(0.5, 0.5, "No insulation data", ha='center', va='center')

    plt.tight_layout()
    out_file = f"{output_prefix}_tad_insulation_summary.pdf"
    plt.savefig(out_file, dpi=300)
    print(f"  Saved TAD insulation plots to {out_file}")
    plt.close()

    # ── export the raw boundary table backing Panel 1 ────────────────────────
    save_tad_boundary_data(tads_data, window, output_prefix)


def calculate_compartments_all_conditions(mcool_files, conditions, chromosomes,
                                          genome_fasta, resolution=64000, n_eigs=3):
    """
    Calculate A/B compartments for all conditions using cooltools.
    """
    all_compartments = {}

    print(f"Loading genome FASTA for GC phasing: {genome_fasta}")
    fasta_records = bioframe.load_fasta(genome_fasta)

    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"\nCalculating compartments for {condition}...")

        actual_resolution = get_closest_resolution(mcool_file, resolution)
        if actual_resolution is None:
            print(f"Error: Could not determine resolution for {mcool_file}")
            continue

        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_resolution}")

        # Build viewframe restricted to chromosomes of interest (drop Y, Mt)
        from cooltools.lib.common import make_cooler_view
        view_df = make_cooler_view(clr)
        filtered_chroms = [c for c in chromosomes if c not in ['Y', 'Mt']]
        view_df = view_df[view_df['chrom'].isin(filtered_chroms)].reset_index(drop=True)

        if view_df.empty:
            print(f"Warning: No valid chromosomes found for {condition}")
            continue

        # Build GC phasing track at this cooler's binning
        bins = clr.bins()[:]
        bins = bins[bins['chrom'].isin(view_df['chrom'].unique())].reset_index(drop=True)
        try:
            gc_track = bioframe.frac_gc(
                bins[['chrom', 'start', 'end']],
                fasta_records,
                mapped_only=True,
            )
            print(f"  Built GC phasing track: {len(gc_track)} bins, "
                  f"GC range [{gc_track['GC'].min():.3f}, {gc_track['GC'].max():.3f}]")
        except Exception as e:
            print(f"Warning: Could not build GC phasing track for {condition}: {e}")
            print(f"  Falling back to UNphased eigendecomposition — A/B labels may flip "
                  f"between chromosomes/conditions.")
            gc_track = None

        try:
            print(f"  Running eigs_cis for {condition}...")
            eigvals, eigvecs = eigs_cis(
                clr,
                phasing_track=gc_track,
                view_df=view_df,
                n_eigs=n_eigs,
                clr_weight_name='weight',
                ignore_diags=2,
                clip_percentile=99.9,
            )

            if eigvecs.empty:
                print(f"Warning: empty eigvecs for {condition}")
                continue

            eigvecs['condition'] = condition
            eigvecs['resolution'] = actual_resolution

            if 'E1' not in eigvecs.columns:
                print(f"Warning: No E1 column for {condition}")
                continue

            # With GC phasing, E1>0 == high GC == A compartment
            eigvecs['compartment'] = np.where(eigvecs['E1'] > 0, 'A', 'B')

            # QC: per-chromosome corr(E1, GC) should be positive after phasing
            if gc_track is not None:
                merged_qc = eigvecs.merge(gc_track, on=['chrom', 'start', 'end'])
                merged_qc = merged_qc.dropna(subset=['E1', 'GC'])
                if len(merged_qc) > 10:
                    r_global = merged_qc[['E1', 'GC']].corr().iloc[0, 1]
                    print(f"  QC global corr(E1, GC) = {r_global:+.3f} "
                          f"({'OK' if r_global > 0 else 'NEGATIVE — phasing failed!'})")
                    for chrom in sorted(merged_qc['chrom'].unique()):
                        sub = merged_qc[merged_qc['chrom'] == chrom]
                        if len(sub) > 10:
                            r = sub[['E1', 'GC']].corr().iloc[0, 1]
                            flag = '' if r > 0 else '  <-- sign-flipped, check!'
                            print(f"    {chrom}: r = {r:+.3f}, n={len(sub)}{flag}")

            all_compartments[condition] = eigvecs

        except Exception as e:
            print(f"Warning: Could not calculate compartments for {condition}: {e}")
            import traceback
            traceback.print_exc()
            continue

    return all_compartments


def compare_compartments_to_null(compartment_data, null_model, conditions, fdr_threshold=0.05):
    """
    Compare compartment changes between conditions.
    """
    print("\nComparing compartment changes...")

    ref_condition = None
    for cond in ['DOX', 'uninfected']:
        if cond in compartment_data:
            ref_condition = cond
            break

    if ref_condition is None:
        print("Warning: No reference condition found in compartment data")
        return pd.DataFrame()

    uninfected = compartment_data[ref_condition]
    comparison_results = []

    for infected_condition in conditions:
        if infected_condition == ref_condition or infected_condition not in compartment_data:
            continue

        infected = compartment_data[infected_condition]

        merged = pd.merge(
            uninfected[['chrom', 'start', 'end', 'E1', 'compartment']],
            infected[['chrom', 'start', 'end', 'E1', 'compartment']],
            on=['chrom', 'start', 'end'],
            suffixes=('_uninf', '_inf')
        )

        if len(merged) == 0:
            print(f"Warning: No matching bins between {ref_condition} and {infected_condition}")
            continue

        merged['switch'] = merged['compartment_uninf'] != merged['compartment_inf']
        merged['E1_diff'] = merged['E1_inf'] - merged['E1_uninf']
        merged['abs_E1_diff'] = np.abs(merged['E1_diff'])

        n_before = len(merged)
        merged = merged.dropna(subset=['E1_uninf', 'E1_inf', 'E1_diff'])
        n_after = len(merged)
        if n_before != n_after:
            print(f"    Removed {n_before - n_after} bins with missing E1 values "
                  f"({(n_before-n_after)/n_before*100:.1f}%)")

        if len(merged) == 0:
            print(f"    Warning: No valid bins remaining after filtering NaNs")
            continue

        print(f"  Testing {len(merged)} bins for {infected_condition}...")

        # Global paired test
        t_stat, global_p_value = stats.ttest_rel(merged['E1_inf'], merged['E1_uninf'])
        print(f"    Global E1 change: t={t_stat:.3f}, p={global_p_value:.3e}")

        # Per-bin z-score test on |E1_diff|
        mean_abs_diff = merged['abs_E1_diff'].mean()
        std_abs_diff = merged['abs_E1_diff'].std()
        print(f"    Mean |E1_diff|: {mean_abs_diff:.4f}, SD: {std_abs_diff:.4f}")

        merged['z_score'] = (merged['abs_E1_diff'] - mean_abs_diff) / std_abs_diff
        merged['p_value'] = 2 * (1 - stats.norm.cdf(np.abs(merged['z_score'])))
        _, merged['fdr'], _, _ = multipletests(merged['p_value'], method='fdr_bh')
        merged['significant'] = merged['fdr'] < fdr_threshold

        # Label permutation null for compartment switch rate
        observed_switches = merged['switch'].sum()
        total_bins = len(merged)
        switch_rate = observed_switches / total_bins

        n_perms = 1000
        uninf_labels = merged['compartment_uninf'].values
        inf_labels = merged['compartment_inf'].values
        null_switches = np.zeros(n_perms)

        for i in range(n_perms):
            shuffled = np.random.permutation(inf_labels)
            null_switches[i] = np.sum(uninf_labels != shuffled)

        mean_null = np.mean(null_switches)
        # two-sided permutation p-value
        switch_p_value = np.sum(np.abs(null_switches - mean_null) >= np.abs(observed_switches - mean_null)) / n_perms

        merged['comparison'] = f"{ref_condition}_vs_{infected_condition}"
        merged['switch_rate'] = switch_rate
        merged['switch_p_value'] = switch_p_value
        merged['global_e1_p_value'] = global_p_value

        comparison_results.append(merged)

        n_sig = merged['significant'].sum()
        print(f"    {n_sig} bins significant at FDR < {fdr_threshold} "
              f"({n_sig/total_bins*100:.1f}%)")
        print(f"    P-value range: [{merged['p_value'].min():.4f}, {merged['p_value'].max():.4f}]")
        print(f"    Top 5 most significant bins have |E1_diff| >= "
              f"{merged.nsmallest(5, 'p_value')['abs_E1_diff'].min():.3f}")
        print(f"    E1_diff range: [{merged['E1_diff'].min():.3f}, {merged['E1_diff'].max():.3f}]")
        print(f"    Switch rate: {switch_rate:.2%}, perm_p={switch_p_value:.3e}")

    if comparison_results:
        return pd.concat(comparison_results, ignore_index=True)
    else:
        return pd.DataFrame()


def calculate_saddle_plots(mcool_files, conditions, compartment_data, chromosomes, resolution=64000, n_bins=30, output_prefix="out"):
    """
    Calculate and plot compartment saddle plots and compartmentalization strength.

    This method quantifies global compartmentalization strength without per-bin testing.
    Popularized for differential Hi-C/Micro-C work by:
      - Schwarzer et al. 2017 (Nature 551:51-56)
      - Nora et al. 2017 (Cell 169(5):930-944)
    and implemented here utilizing the cooltools API.
    """
    print("\n" + "=" * 60)
    print("CALCULATING COMPARTMENT SADDLE PLOTS")
    print("=" * 60)

    saddle_strengths = []

    # Setup plotting grid
    fig, axes = plt.subplots(1, len(conditions), figsize=(5 * len(conditions), 4))
    if len(conditions) == 1:
        axes = [axes]

    for idx, (condition, mcool_file) in enumerate(zip(conditions, mcool_files)):
        print(f"  Processing saddle for {condition}...")
        actual_res = get_closest_resolution(mcool_file, resolution)
        if actual_res is None or condition not in compartment_data:
            print(f"  Skipping {condition}: missing mcool resolution or E1 tracks")
            continue

        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_res}")

        from cooltools.lib.common import make_cooler_view
        view_df = make_cooler_view(clr)
        filtered_chroms = [c for c in chromosomes if c not in ['Y', 'Mt']]
        view_df = view_df[view_df['chrom'].isin(filtered_chroms)].reset_index(drop=True)

        # Get E1 track and ensure proper format
        eigvecs = compartment_data[condition]
        track = eigvecs[['chrom', 'start', 'end', 'E1']].copy()
        track = track.dropna()

        try:
            print(f"    Calculating expected cis contacts...")
            expected = expected_cis(clr, view_df=view_df, nproc=4)

            print(f"    Digitizing E1 track into {n_bins} bins (numpy quantile method)...")
            # Use numpy instead of cooltools.digitize to avoid version-specific
            # API breakage (v_name was removed in cooltools 0.7.x).
            E1_vals   = track['E1'].values
            fin_mask  = np.isfinite(E1_vals)
            pct_edges = np.linspace(0, 100, n_bins + 1)
            bin_edges = np.nanpercentile(E1_vals[fin_mask], pct_edges)
            bin_edges = np.unique(bin_edges)        # collapse duplicates
            states    = np.zeros(len(E1_vals), dtype=int)
            states[fin_mask] = (
                np.digitize(E1_vals[fin_mask], bin_edges[1:], right=True) + 1
            )
            q_track          = track.copy()
            q_track['state'] = states

            print(f"    Computing saddle map...")
            interaction_sum, interaction_count = cooltools.saddle(
                clr,
                expected,
                q_track,
                'cis',
                view_df=view_df,
                q_name='state',
                nproc=4
            )

            # Generate observed/expected saddle map
            saddle_map = interaction_sum / interaction_count

            # Calculate global compartmentalization strength from corners
            # Utilizing top and bottom 20% bins as the active/inactive corners
            n_corner = max(1, n_bins // 5)
            BB = np.nanmean(saddle_map[:n_corner, :n_corner])
            AA = np.nanmean(saddle_map[-n_corner:, -n_corner:])
            BA = np.nanmean(saddle_map[:n_corner, -n_corner:])
            AB = np.nanmean(saddle_map[-n_corner:, :n_corner])

            # Compartmentalization strength: (AA * BB) / (AB * BA)
            strength = (AA * BB) / (AB * BA)

            saddle_strengths.append({
                'condition': condition,
                'strength': strength,
                'AA_corner': AA,
                'BB_corner': BB,
                'AB_corner': AB,
                'BA_corner': BA
            })

            # Plot heatmap
            ax = axes[idx]
            im = ax.imshow(
                np.log2(saddle_map),
                cmap='coolwarm',
                vmin=-1, vmax=1
            )
            ax.set_title(f"{condition}\nStrength: {strength:.2f}")
            ax.set_xticks([])
            ax.set_yticks([])
            if idx == 0:
                ax.set_ylabel("B -> A")
                ax.set_xlabel("B -> A")

        except Exception as e:
            print(f"  Warning: Failed to calculate saddle for {condition}: {e}")
            import traceback
            traceback.print_exc()

    if saddle_strengths:
        plt.colorbar(im, ax=axes, label="log2(obs/exp)")
        plt.savefig(f"{output_prefix}_saddle_plots.pdf", bbox_inches='tight')
        plt.close()

        df_strength = pd.DataFrame(saddle_strengths)
        df_strength.to_csv(f"{output_prefix}_saddle_strengths.tsv", sep='\t', index=False)
        print(f"  Saved saddle plots to {output_prefix}_saddle_plots.pdf")
        print(f"  Saved saddle strengths to {output_prefix}_saddle_strengths.tsv")
        return df_strength
    else:
        plt.close()
        return pd.DataFrame()


def call_loops_all_conditions(mcool_files, conditions, chromosomes,
                              resolution=8000, fdr_threshold=0.1):
    """Call loops for all conditions using cooltools dots."""
    all_loops = {}

    for condition, mcool_file in zip(conditions, mcool_files):
        print(f"\nCalling loops for {condition}...")

        actual_resolution = get_closest_resolution(mcool_file, resolution)
        if actual_resolution is None:
            print(f"Error: Could not determine resolution for {mcool_file}")
            continue

        clr = cooler.Cooler(f"{mcool_file}::resolutions/{actual_resolution}")

        try:
            print(f"Calculating expected for {condition}...")
            expected_df = expected_cis(clr, nproc=4)

            print(f"Calling dots for {condition}...")
            dot_calls = dots(
                clr,
                expected_df,
                expected_value_col="balanced.avg",
                clr_weight_name="weight",
                max_loci_separation=2000000,
                lambda_bin_fdr=fdr_threshold,
                clustering_radius=20000,
                nproc=4
            )

            dot_calls['condition'] = condition
            dot_calls['resolution'] = actual_resolution
            all_loops[condition] = dot_calls

        except Exception as e:
            print(f"Warning: Could not call loops for {condition}: {e}")
            try:
                print(f"  Trying simplified loop calling for {condition}...")
                dot_calls = dots(
                    clr,
                    expected_df,
                    max_loci_separation=2000000,
                    lambda_bin_fdr=fdr_threshold,
                    clustering_radius=20000,
                    nproc=1
                )
                dot_calls['condition'] = condition
                dot_calls['resolution'] = actual_resolution
                all_loops[condition] = dot_calls
            except Exception as e2:
                print(f"  Simplified approach also failed: {e2}")
                continue

    return all_loops


def compare_loops_to_null(loop_data, null_model, conditions, fdr_threshold=0.05):
    """Compare loop calls between conditions and test for significance."""
    print("\nComparing loops between conditions...")

    ref_condition = None
    for cond in ['DOX', 'uninfected']:
        if cond in loop_data:
            ref_condition = cond
            break

    if ref_condition is None:
        print("Warning: No reference condition found in loop data")
        print(f"Available conditions: {list(loop_data.keys())}")
        return pd.DataFrame()

    uninfected_loops = loop_data[ref_condition]
    comparison_results = []

    for infected_condition in conditions:
        if infected_condition == ref_condition or infected_condition not in loop_data:
            continue

        infected_loops = loop_data[infected_condition]

        def create_loop_set(df):
            loop_set = set()
            for _, row in df.iterrows():
                chrom1 = row['chrom1']
                start1 = round(row['start1'] / 10000) * 10000
                chrom2 = row['chrom2']
                start2 = round(row['start2'] / 10000) * 10000
                loop_set.add((chrom1, start1, chrom2, start2))
            return loop_set

        print(f"  Comparing {len(uninfected_loops)} uninfected loops vs "
              f"{len(infected_loops)} infected loops")

        uninf_set = create_loop_set(uninfected_loops)
        inf_set = create_loop_set(infected_loops)

        gained = inf_set - uninf_set
        lost = uninf_set - inf_set

        print(f"  Found {len(gained)} gained loops and {len(lost)} lost loops")

        n_permutations = 1000
        null_gained = []
        null_lost = []

        all_loops_list = list(uninf_set | inf_set)
        n_uninf = len(uninf_set)
        n_total = len(all_loops_list)

        if n_total == 0:
            print(f"  No loops found for comparison in {infected_condition}")
            continue

        for _ in range(n_permutations):
            perm_indices = np.random.permutation(n_total)
            perm_uninf_indices = perm_indices[:min(n_uninf, n_total)]
            perm_uninf = set(all_loops_list[i] for i in perm_uninf_indices)
            perm_inf = set(all_loops_list) - perm_uninf
            null_gained.append(len(perm_inf - perm_uninf))
            null_lost.append(len(perm_uninf - perm_inf))

        if len(null_gained) > 0:
            p_gained = np.sum(np.array(null_gained) >= len(gained)) / n_permutations
            p_lost = np.sum(np.array(null_lost) >= len(lost)) / n_permutations
        else:
            p_gained = 1.0
            p_lost = 1.0

        n_comparisons = len(conditions) - 1
        fdr_gained = min(p_gained * n_comparisons, 1.0)
        fdr_lost = min(p_lost * n_comparisons, 1.0)

        comparison_results.append({
            'comparison': f"{ref_condition}_vs_{infected_condition}",
            'gained_loops': len(gained),
            'lost_loops': len(lost),
            'total_uninf_loops': len(uninf_set),
            'total_inf_loops': len(inf_set),
            'p_value_gained': p_gained,
            'p_value_lost': p_lost,
            'fdr_gained': fdr_gained,
            'fdr_lost': fdr_lost,
            'significant_gained': fdr_gained < fdr_threshold,
            'significant_lost': fdr_lost < fdr_threshold
        })

        print(f"  {infected_condition}: gained={len(gained)} (p={p_gained:.3f}), "
              f"lost={len(lost)} (p={p_lost:.3f})")

    return pd.DataFrame(comparison_results)


def create_compartment_switch_plots(compartment_comp, output_prefix):
    """
    Detailed plots for compartment switches.
    """
    if compartment_comp.empty:
        print("No compartment comparison data to plot")
        return

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    for comparison in compartment_comp['comparison'].unique():
        data = compartment_comp[compartment_comp['comparison'] == comparison]
        sig_data = data[data['significant']]
        not_sig_data = data[~data['significant']]

        comp_name = comparison.replace('_vs_', ' vs ')

        # Plot 1: Significant A->B vs B->A switches
        ax = axes[0, 0]
        if len(sig_data) > 0:
            a_to_b = sig_data[(sig_data['compartment_uninf'] == 'A') &
                              (sig_data['compartment_inf'] == 'B')].shape[0]
            b_to_a = sig_data[(sig_data['compartment_uninf'] == 'B') &
                              (sig_data['compartment_inf'] == 'A')].shape[0]
            no_switch_sig = sig_data[~sig_data['switch']].shape[0]

            categories = ['A→B', 'B→A', 'No Switch']
            counts = [a_to_b, b_to_a, no_switch_sig]
            colors = ['#e74c3c', '#3498db', '#95a5a6']

            ax.bar(categories, counts, color=colors, alpha=0.8, edgecolor='black')
            ax.set_ylabel('Number of Significant Bins')
            ax.set_title(f'Significant Compartment Changes\n{comp_name}')
            for i, (cat, count) in enumerate(zip(categories, counts)):
                ax.text(i, count, str(count), ha='center', va='bottom', fontweight='bold')

        # Plot 2: % Changed vs Unchanged (significant only)
        ax = axes[0, 1]
        if len(sig_data) > 0:
            switched = sig_data[sig_data['switch']].shape[0]
            not_switched = sig_data[~sig_data['switch']].shape[0]

            total_sig = switched + not_switched
            pct_switched = (switched / total_sig * 100) if total_sig > 0 else 0
            pct_not_switched = (not_switched / total_sig * 100) if total_sig > 0 else 0

            categories = ['Switched', 'Unchanged']
            percentages = [pct_switched, pct_not_switched]
            colors = ['#e67e22', '#2ecc71']

            bars = ax.bar(categories, percentages, color=colors, alpha=0.8, edgecolor='black')
            ax.set_ylabel('Percentage (%)')
            ax.set_title(f'Significant Bins: Changed vs Unchanged\n{comp_name}')
            ax.set_ylim([0, 100])
            for bar, pct in zip(bars, percentages):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2., height,
                        f'{pct:.1f}%', ha='center', va='bottom', fontweight='bold')

        # Plot 3: Chromosome-level switches by direction
        ax = axes[0, 2]
        if len(sig_data) > 0:
            chrom_order = ['2L', '2R', '3L', '3R', '4', 'X']
            chroms_present = [c for c in chrom_order if c in sig_data['chrom'].unique()]

            a_to_b_by_chrom = []
            b_to_a_by_chrom = []
            for chrom in chroms_present:
                chrom_data = sig_data[sig_data['chrom'] == chrom]
                a_to_b = chrom_data[(chrom_data['compartment_uninf'] == 'A') &
                                    (chrom_data['compartment_inf'] == 'B')].shape[0]
                b_to_a = chrom_data[(chrom_data['compartment_uninf'] == 'B') &
                                    (chrom_data['compartment_inf'] == 'A')].shape[0]
                a_to_b_by_chrom.append(a_to_b)
                b_to_a_by_chrom.append(b_to_a)

            x = np.arange(len(chroms_present))
            width = 0.35
            ax.bar(x - width / 2, a_to_b_by_chrom, width, label='A→B',
                   color='#e74c3c', alpha=0.8, edgecolor='black')
            ax.bar(x + width / 2, b_to_a_by_chrom, width, label='B→A',
                   color='#3498db', alpha=0.8, edgecolor='black')
            ax.set_xlabel('Chromosome')
            ax.set_ylabel('Number of Significant Switches')
            ax.set_title(f'Switches by Chromosome\n{comp_name}')
            ax.set_xticks(x)
            ax.set_xticklabels(chroms_present)
            ax.legend()

        # Plot 4: E1 difference distributions for switched vs non-switched
        ax = axes[1, 0]
        if len(sig_data) > 0:
            switched_e1 = sig_data[sig_data['switch']]['abs_E1_diff']
            not_switched_e1 = sig_data[~sig_data['switch']]['abs_E1_diff']

            ax.hist(switched_e1, bins=30, alpha=0.6, label='Switched',
                    color='#e67e22', edgecolor='black')
            ax.hist(not_switched_e1, bins=30, alpha=0.6, label='Unchanged',
                    color='#2ecc71', edgecolor='black')
            ax.set_xlabel('|E1 Difference|')
            ax.set_ylabel('Count')
            ax.set_title(f'E1 Changes in Significant Bins\n{comp_name}')
            ax.legend()

        # Plot 5: Switch rate, significant vs non-significant
        ax = axes[1, 1]
        sig_switched = sig_data[sig_data['switch']].shape[0]
        sig_invariant = sig_data[~sig_data['switch']].shape[0]
        not_sig_switched = not_sig_data[not_sig_data['switch']].shape[0]
        not_sig_invariant = not_sig_data[~not_sig_data['switch']].shape[0]

        sig_switch_rate = sig_switched / (sig_switched + sig_invariant) if (sig_switched + sig_invariant) > 0 else 0
        not_sig_switch_rate = not_sig_switched / (not_sig_switched + not_sig_invariant) if (not_sig_switched + not_sig_invariant) > 0 else 0

        categories = ['Significant Bins', 'Non-Significant Bins']
        switch_rates = [sig_switch_rate * 100, not_sig_switch_rate * 100]
        colors = ['#e67e22', '#95a5a6']

        bars = ax.bar(categories, switch_rates, color=colors, alpha=0.8, edgecolor='black')
        ax.set_ylabel('Switch Rate (%)')
        ax.set_title(f'Switch Rate: Significant vs Non-Significant\n{comp_name}')
        ax.set_ylim([0, 100])
        for bar, rate in zip(bars, switch_rates):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2., height,
                    f'{rate:.1f}%', ha='center', va='bottom', fontweight='bold')

        # Plot 6: Statistical summary table
        ax = axes[1, 2]
        ax.axis('off')

        if len(sig_data) > 0:
            total_bins = len(data)
            n_significant = len(sig_data)

            switch_rate_sig = sig_switched / sig_invariant if sig_invariant > 0 else float('inf')
            switch_rate_not_sig = not_sig_switched / not_sig_invariant if not_sig_invariant > 0 else float('inf')

            from scipy.stats import fisher_exact
            contingency_table = [[sig_switched, sig_invariant],
                                 [not_sig_switched, not_sig_invariant]]
            odds_ratio, fisher_p = fisher_exact(contingency_table)

            # Label permutation for significant bins
            uninf_labels = sig_data['compartment_uninf'].values
            inf_labels = sig_data['compartment_inf'].values
            n_perms = 1000
            null_switches = np.zeros(n_perms)
            for i in range(n_perms):
                shuffled = np.random.permutation(inf_labels)
                null_switches[i] = np.sum(uninf_labels != shuffled)

            mean_null = np.mean(null_switches)
            perm_p_sig = np.sum(np.abs(null_switches - mean_null) >= np.abs(sig_switched - mean_null)) / n_perms

            summary_text = [
                f"Statistical Summary - {comp_name}",
                "=" * 50,
                f"Total bins: {total_bins}",
                f"Significant bins: {n_significant} ({n_significant/total_bins*100:.1f}%)",
                "",
                "Significant Bins:",
                f"  Switched: {sig_switched}",
                f"  Invariant: {sig_invariant}",
                f"  Switch rate: {sig_switch_rate*100:.1f}%",
                f"  Ratio (sw/inv): {switch_rate_sig:.3f}" if sig_invariant > 0 else "  Ratio: undefined",
                "",
                "Non-Significant Bins:",
                f"  Switched: {not_sig_switched}",
                f"  Invariant: {not_sig_invariant}",
                f"  Switch rate: {not_sig_switch_rate*100:.1f}%",
                f"  Ratio (sw/inv): {switch_rate_not_sig:.3f}" if not_sig_invariant > 0 else "  Ratio: undefined",
                "",
                "Statistical Tests:",
                f"  Fisher's exact test:",
                f"    Odds ratio: {odds_ratio:.3f}",
                f"    p-value: {fisher_p:.3e}",
                f"    → {'NOT significant' if fisher_p > 0.05 else 'SIGNIFICANT'}",
                "",
                f"  Sig bins switch vs permutation null:",
                f"    p-value: {perm_p_sig:.3e}",
                f"    → {'NOT significant' if perm_p_sig > 0.05 else 'SIGNIFICANT'}",
            ]

            ax.text(0.05, 0.95, '\n'.join(summary_text),
                    transform=ax.transAxes, fontsize=8, verticalalignment='top',
                    fontfamily='monospace',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            print("\n" + "=" * 60)
            for line in summary_text:
                print(line)
            print("=" * 60)

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_compartment_switches.pdf", dpi=300, bbox_inches='tight')
    print(f"\nSaved compartment switch plots: {output_prefix}_compartment_switches.pdf")
    plt.close()


def create_summary_plots(hotspot_comparison, compartment_comp, loop_comp, output_prefix):
    """Create summary visualizations for all comparisons."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Plot 1: Hotspot changes
    ax = axes[0, 0]
    if not hotspot_comparison.empty and 'comparison' in hotspot_comparison.columns:
        try:
            boundary_stats = hotspot_comparison.groupby('comparison')['significant'].sum()
            if len(boundary_stats) > 0:
                ax.bar(range(len(boundary_stats)), boundary_stats.values)
                ax.set_xticks(range(len(boundary_stats)))
                ax.set_xticklabels([x.replace('_vs_', ' vs ') for x in boundary_stats.index],
                                   rotation=45)
        except Exception as e:
            print(f"Warning: Could not plot Hotspot changes: {e}")
    ax.set_ylabel('Number of Significant Changes')
    ax.set_title('Hotspot Changes')

    # Plot 2: Hotspot strength changes
    ax = axes[0, 1]
    if not hotspot_comparison.empty and 'log2_fold_change' in hotspot_comparison.columns:
        try:
            for comparison in hotspot_comparison['comparison'].unique():
                data = hotspot_comparison[hotspot_comparison['comparison'] == comparison]
                sig_data = data[data['significant']]
                if len(sig_data) > 0:
                    ax.hist(sig_data['log2_fold_change'], alpha=0.5,
                            label=comparison.replace('_vs_', ' vs '), bins=20)
            ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot hotspot strength changes: {e}")
    ax.set_xlabel('Log2 Fold Change in Strength')
    ax.set_ylabel('Count')
    ax.set_title('Hotspot Strength Changes')

    # Plot 3: Compartment switches
    ax = axes[0, 2]
    if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
        try:
            comparisons = compartment_comp['comparison'].unique()
            switch_rates = []
            sig_switch_rates = []
            for comparison in comparisons:
                data = compartment_comp[compartment_comp['comparison'] == comparison]
                switch_rate = data['switch'].sum() / len(data) * 100 if len(data) > 0 else 0
                sig_switch_rate = data[data['significant']]['switch'].sum() / len(data) * 100 if len(data) > 0 else 0
                switch_rates.append(switch_rate)
                sig_switch_rates.append(sig_switch_rate)

            if len(comparisons) > 0:
                x = np.arange(len(comparisons))
                width = 0.35
                ax.bar(x - width / 2, switch_rates, width, label='All', alpha=0.7)
                ax.bar(x + width / 2, sig_switch_rates, width, label='Significant', alpha=0.7)
                ax.set_xticks(x)
                ax.set_xticklabels([x.replace('_vs_', ' vs ') for x in comparisons], rotation=45)
                ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot compartment switches: {e}")
    ax.set_ylabel('Switch Rate (%)')
    ax.set_title('Compartment Switch Rates')

    # Plot 4: E1 value changes
    ax = axes[1, 0]
    if not compartment_comp.empty and 'E1_uninf' in compartment_comp.columns:
        try:
            for comparison in compartment_comp['comparison'].unique():
                data = compartment_comp[compartment_comp['comparison'] == comparison]
                sig_data = data[data['significant']]
                if len(sig_data) > 0:
                    valid_data = sig_data.dropna(subset=['E1_uninf', 'E1_inf'])
                    if len(valid_data) > 0:
                        ax.scatter(valid_data['E1_uninf'], valid_data['E1_inf'],
                                   alpha=0.5, s=10, label=comparison.replace('_vs_', ' vs '))
            ax.plot([-1, 1], [-1, 1], 'k--', alpha=0.5)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot E1 changes: {e}")
    ax.set_xlabel('E1 Uninfected')
    ax.set_ylabel('E1 Infected')
    ax.set_title('Compartment Strength Changes')

    # Plot 5: Loop changes
    ax = axes[1, 1]
    if not loop_comp.empty and 'gained_loops' in loop_comp.columns:
        try:
            x = np.arange(len(loop_comp))
            width = 0.35
            ax.bar(x - width / 2, loop_comp['gained_loops'], width, label='Gained')
            ax.bar(x + width / 2, loop_comp['lost_loops'], width, label='Lost')
            ax.set_xticks(x)
            ax.set_xticklabels([x.replace('_vs_', ' vs ') for x in loop_comp['comparison']],
                               rotation=45)
            ax.legend()
        except Exception as e:
            print(f"Warning: Could not plot loop changes: {e}")
    ax.set_ylabel('Number of Loops')
    ax.set_title('Loop Changes')

    # Plot 6: Summary statistics
    ax = axes[1, 2]
    summary_data = []

    all_comparisons = set()
    if not hotspot_comparison.empty and 'comparison' in hotspot_comparison.columns:
        all_comparisons.update(hotspot_comparison['comparison'].unique())
    if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
        all_comparisons.update(compartment_comp['comparison'].unique())
    if not loop_comp.empty and 'comparison' in loop_comp.columns:
        all_comparisons.update(loop_comp['comparison'].unique())

    for comp in all_comparisons:
        hotspot_sig = 0
        comp_sig = 0
        loop_count = 0
        try:
            if not hotspot_comparison.empty and 'comparison' in hotspot_comparison.columns:
                hotspot_sig = hotspot_comparison[hotspot_comparison['comparison'] == comp]['significant'].sum()
            if not compartment_comp.empty and 'comparison' in compartment_comp.columns:
                comp_sig = compartment_comp[compartment_comp['comparison'] == comp]['significant'].sum()
            if not loop_comp.empty and 'comparison' in loop_comp.columns:
                loop_data = loop_comp[loop_comp['comparison'] == comp]
                if len(loop_data) > 0:
                    loop_count = loop_data['gained_loops'].values[0] + loop_data['lost_loops'].values[0]
        except Exception as e:
            print(f"Warning: Could not calculate summary for {comp}: {e}")

        comp_name = comp.replace('DOX_vs_', '').replace('uninfected_vs_', '')
        summary_data.append({
            'Comparison': comp_name,
            'Hotspots': hotspot_sig,
            'Compartments': comp_sig,
            'Loops': loop_count
        })

    summary_df = pd.DataFrame(summary_data)
    if not summary_df.empty:
        try:
            summary_df.set_index('Comparison').plot(kind='bar', ax=ax)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        except Exception as e:
            print(f"Warning: Could not plot summary: {e}")
    ax.set_ylabel('Count')
    ax.set_title('Summary of Significant Changes')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_chromatin_structure_comparison.pdf", dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Compare chromatin structure using existing diffHic results')
    parser.add_argument('--mcool_files', nargs='+', required=True,
                        help='Micro-C mcool files (must match number of conditions)')
    parser.add_argument('--conditions', nargs='+', required=True,
                        help='Condition names (must match number of mcool files)')
    parser.add_argument('--diffhic_results', required=True,
                        help='Combined diffHic results file (all_results_combined.csv)')
    parser.add_argument('--null_model', required=True,
                        help='Null model results file (null_model_results.csv)')
    parser.add_argument('--genome_fasta', required=True,
                        help='Genome FASTA (e.g. dm6.fa) used to phase compartment '
                             'eigenvectors by GC content. Without phasing, E1 sign '
                             'is arbitrary and A/B labels can flip between '
                             'chromosomes/conditions.')
    parser.add_argument('--chromosomes', nargs='+',
                        default=['2L', '2R', '3L', '3R', '4', 'X'],
                        help='Chromosomes to analyze')
    parser.add_argument('--resolution_compartment', type=int, default=50000,
                        help='Resolution for compartment analysis')
    parser.add_argument('--n_saddle_bins', type=int, default=30,
                        help='Number of quantiles for saddle plot compartmentalization analysis')
    parser.add_argument('--resolution_loop', type=int, default=5000,
                        help='Resolution for loop calling')
    parser.add_argument('--fdr_threshold', type=float, default=0.1,
                        help='FDR threshold for significance')
    parser.add_argument('--hotspot_window_size', type=int, default=50000,
                        help='Window size for differential-interaction hotspots identification')
    parser.add_argument('--resolution_insulation', type=int, default=50000,
                        help='Resolution for cooltools TAD insulation calculation')
    parser.add_argument('--window_insulation', type=int, default=150000,
                        help='Window size for cooltools TAD insulation calculation')
    parser.add_argument('--skip_loops', action='store_true', default=False,
                        help='Skip loop calling (cooltools dots). Useful when the '
                             'job is memory/time-constrained — all other analyses '
                             'still run.')
    parser.add_argument('--output_prefix', required=True,
                        help='Output file prefix')

    args = parser.parse_args()

    if len(args.mcool_files) != len(args.conditions):
        raise ValueError(f"Number of mcool files ({len(args.mcool_files)}) must match "
                         f"number of conditions ({len(args.conditions)})")

    print(f"Analyzing {len(args.conditions)} conditions: {', '.join(args.conditions)}")
    print(f"Using cooltools version: {cooltools.__version__}")
    print(f"Genome FASTA for GC phasing: {args.genome_fasta}")

    # Load diffHic results and null model
    all_results, significant_results = load_diffhic_results(args.diffhic_results)
    null_model = load_null_model(args.null_model)

    # Hotspots from differential interactions
    print("\n" + "=" * 60)
    print("IDENTIFYING DIFFERENTIAL-INTERACTION HOTSPOTS")
    print("=" * 60)

    hotspots = identify_hotspots_from_interactions(
        significant_results, args.conditions, window_size=args.hotspot_window_size
    )

    hotspot_comparison = compare_hotspots(
        hotspots, null_model, args.fdr_threshold
    )

    # Cooltools Insulation / TADs
    print("\n" + "=" * 60)
    print("CALCULATING COOLTOOLS TADS (INSULATION)")
    print("=" * 60)

    cooltools_tads = calculate_cooltools_tads(
        args.mcool_files, args.conditions,
        resolution=args.resolution_insulation,
        window=args.window_insulation
    )

    # TAD visualizations run immediately so they are saved even if loops kill the job
    print("\n" + "=" * 60)
    print("TAD INSULATION PLOTS + VENN DIAGRAMS (saved before loop calling)")
    print("=" * 60)
    create_tad_insulation_plots(cooltools_tads, args.window_insulation, args.output_prefix)
    create_tad_venn_diagrams(cooltools_tads, args.window_insulation,
                             args.chromosomes, args.output_prefix)

    # Compartments (GC-phased)
    print("\n" + "=" * 60)
    print("CALCULATING COMPARTMENTS (GC-phased)")
    print("=" * 60)

    compartment_data = calculate_compartments_all_conditions(
        args.mcool_files, args.conditions, args.chromosomes,
        genome_fasta=args.genome_fasta,
        resolution=args.resolution_compartment,
    )

    compartment_comparison = compare_compartments_to_null(
        compartment_data, null_model, args.conditions, args.fdr_threshold
    )

    # Saddle Plots & Compartmentalization Strength
    saddle_strengths_df = calculate_saddle_plots(
        args.mcool_files, args.conditions, compartment_data, args.chromosomes,
        resolution=args.resolution_compartment,
        n_bins=args.n_saddle_bins,
        output_prefix=args.output_prefix
    )

    # Loops (skippable with --skip_loops)
    print("\n" + "=" * 60)
    print("CALLING LOOPS")
    print("=" * 60)

    if args.skip_loops:
        print("  --skip_loops set: skipping loop calling.")
        loop_data       = {}
        loop_comparison = pd.DataFrame()
    else:
        loop_data = call_loops_all_conditions(
            args.mcool_files, args.conditions, args.chromosomes,
            resolution=args.resolution_loop
        )
        loop_comparison = compare_loops_to_null(
            loop_data, null_model, args.conditions, args.fdr_threshold
        )

    # Visualizations
    print("\n" + "=" * 60)
    print("CREATING VISUALIZATIONS")
    print("=" * 60)

    create_compartment_switch_plots(compartment_comparison, args.output_prefix)
    create_summary_plots(
        hotspot_comparison, compartment_comparison, loop_comparison,
        args.output_prefix
    )

    # Save results
    print("\n" + "=" * 60)
    print("SAVING RESULTS")
    print("=" * 60)

    if not hotspot_comparison.empty:
        hotspot_comparison.to_csv(
            f"{args.output_prefix}_hotspot_comparison.tsv", sep='\t', index=False
        )
        print(f"Saved hotspot comparison: {args.output_prefix}_hotspot_comparison.tsv")

    if not compartment_comparison.empty:
        compartment_comparison.to_csv(
            f"{args.output_prefix}_compartment_comparison.tsv", sep='\t', index=False
        )
        print(f"Saved compartment comparison: {args.output_prefix}_compartment_comparison.tsv")

    if not loop_comparison.empty:
        loop_comparison.to_csv(
            f"{args.output_prefix}_loop_comparison.tsv", sep='\t', index=False
        )
        print(f"Saved loop comparison: {args.output_prefix}_loop_comparison.tsv")

    # Save raw compartment calls per condition
    for condition, eigvecs in compartment_data.items():
        out_file = f"{args.output_prefix}_{condition}_compartments.tsv"
        eigvecs.to_csv(out_file, sep='\t', index=False)
        print(f"Saved {condition} compartments: {out_file}")

    # Save Cooltools TAD insulation per condition
    for condition, ins_df in cooltools_tads.items():
        out_file = f"{args.output_prefix}_{condition}_cooltools_insulation.tsv"
        ins_df.to_csv(out_file, sep='\t', index=False)
        print(f"Saved {condition} cooltools TAD insulation: {out_file}")

    # Save Hotspots as BED
    for condition, hotspots_df in hotspots.items():
        bed_file = f"{args.output_prefix}_{condition}_hotspots.bed"
        hotspots_df[['chrom', 'start', 'end']].to_csv(
            bed_file, sep='\t', index=False, header=False
        )
        print(f"Saved {condition} hotspots: {bed_file}")

    # Summary
    summary_stats = {
        'total_interactions_analyzed': len(all_results),
        'significant_interactions': len(significant_results),
        'hotspots_identified': sum(len(df) for df in hotspots.values()),
        'significant_hotspot_changes': hotspot_comparison['significant'].sum() if not hotspot_comparison.empty else 0,
        'significant_compartment_changes': compartment_comparison['significant'].sum() if not compartment_comparison.empty else 0,
        'total_loop_changes': loop_comparison[['gained_loops', 'lost_loops']].sum().sum() if not loop_comparison.empty else 0
    }

    with open(f"{args.output_prefix}_analysis_summary.txt", 'w') as f:
        f.write("CHROMATIN STRUCTURE ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Conditions analyzed: {', '.join(args.conditions)}\n")
        f.write(f"Genome FASTA for GC phasing: {args.genome_fasta}\n\n")
        f.write(f"Total interactions analyzed: {summary_stats['total_interactions_analyzed']}\n")
        f.write(f"Significant differential interactions: {summary_stats['significant_interactions']}\n")
        f.write(f"Hotspots identified: {summary_stats['hotspots_identified']}\n")
        f.write(f"Significant hotspot changes: {summary_stats['significant_hotspot_changes']}\n")
        f.write(f"Significant compartment changes: {summary_stats['significant_compartment_changes']}\n")
        f.write(f"Total loop changes: {summary_stats['total_loop_changes']}\n\n")

        f.write("HOTSPOTS PER CONDITION:\n")
        for condition, hs_df in hotspots.items():
            f.write(f"  {condition}: {len(hs_df)} hotspots\n")

        f.write("\nCOMPARTMENTALIZATION STRENGTH:\n")
        if not saddle_strengths_df.empty:
            for _, row in saddle_strengths_df.iterrows():
                f.write(f"  {row['condition']}: {row['strength']:.3f}\n")
        else:
            f.write("  No saddle strength data calculated.\n")

        f.write("\nSIGNIFICANT CHANGES PER COMPARISON:\n")
        if not hotspot_comparison.empty:
            for comparison in hotspot_comparison['comparison'].unique():
                sig_count = hotspot_comparison[
                    hotspot_comparison['comparison'] == comparison
                ]['significant'].sum()
                f.write(f"  Hotspots {comparison}: {sig_count}\n")

        if not compartment_comparison.empty:
            for comparison in compartment_comparison['comparison'].unique():
                sub = compartment_comparison[compartment_comparison['comparison'] == comparison]
                sig_count = sub['significant'].sum()
                switch_rate = sub['switch_rate'].iloc[0] if len(sub) > 0 else 0
                f.write(f"  Compartments {comparison}: {sig_count} significant changes, "
                        f"{switch_rate:.2%} switch rate\n")

        if not loop_comparison.empty:
            for _, row in loop_comparison.iterrows():
                f.write(f"  Loops {row['comparison']}: +{row['gained_loops']} gained, "
                        f"-{row['lost_loops']} lost\n")

    print(f"\nAnalysis complete! Results saved to {args.output_prefix}_*")
    print(f"Summary saved to: {args.output_prefix}_analysis_summary.txt")


if __name__ == '__main__':
    main()