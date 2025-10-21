#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pybedtools
from scipy import stats
from statsmodels.stats.multitest import multipletests
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict

def load_interactions(interactions_file):
    """Load differential interactions (wMel vs DOX comparison)."""
    print(f"Loading interactions from {interactions_file}")
    interactions = pd.read_csv(interactions_file)
    
    print(f"Loaded {len(interactions)} interactions")
    print("Available columns:", interactions.columns.tolist())
    
    # Validate required columns
    required_cols = ['logFC', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
    missing_cols = [col for col in required_cols if col not in interactions.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    return interactions


def load_null_model(null_file):
    """Load and characterize the null model."""
    print(f"\nLoading null model from {null_file}")
    null_model = pd.read_csv(null_file)
    
    print(f"Loaded {len(null_model)} null model entries")
    
    # Extract null distribution parameters
    null_logfc = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
    
    null_stats = {
        'n': len(null_logfc),
        'mean': null_logfc.mean(),
        'std': null_logfc.std(),
        'median': null_logfc.median(),
        'min': null_logfc.min(),
        'max': null_logfc.max()
    }
    
    print("\nNull Model Characterization:")
    print(f"  Size: {null_stats['n']} interactions")
    print(f"  Mean logFC: {null_stats['mean']:.4f}")
    print(f"  Std logFC: {null_stats['std']:.4f}")
    print(f"  Median: {null_stats['median']:.4f}")
    print(f"  Range: [{null_stats['min']:.4f}, {null_stats['max']:.4f}]")
    
    if abs(null_stats['mean']) > 0.5:
        print(f"  WARNING: Null mean ({null_stats['mean']:.4f}) is far from 0!")
    
    return null_model, null_stats


def classify_enhancers(enhancer_file, classification_file=None):
    """Load and classify enhancers."""
    print("\nLoading enhancer annotations...")
    
    try:
        enhancers = pd.read_csv(enhancer_file, sep='\t', header=None,
                               names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'class'])
        print("Loaded enhancer file with 7 columns")
    except:
        try:
            enhancers = pd.read_csv(enhancer_file, sep='\t')
            print(f"Loaded enhancer file with header: {enhancers.columns.tolist()}")
        except Exception as e:
            print(f"Error reading enhancer file: {e}")
            return pd.DataFrame()
    
    # Ensure required columns
    required_cols = ['chrom', 'start', 'end']
    for col in required_cols:
        if col not in enhancers.columns:
            raise ValueError(f"Required column '{col}' not found")
    
    # Convert coordinates
    enhancers['start'] = pd.to_numeric(enhancers['start'], errors='coerce')
    enhancers['end'] = pd.to_numeric(enhancers['end'], errors='coerce')
    enhancers = enhancers.dropna(subset=['start', 'end'])
    enhancers['start'] = enhancers['start'].astype(int)
    enhancers['end'] = enhancers['end'].astype(int)
    
    # Remove invalid ranges
    enhancers = enhancers[enhancers['start'] < enhancers['end']]
    enhancers = enhancers[(enhancers['start'] >= 0) & (enhancers['end'] >= 0)]
    
    print(f"Valid enhancers: {len(enhancers)}")
    
    # Create name column if needed
    if 'name' not in enhancers.columns:
        enhancers['name'] = enhancers.apply(lambda x: f"{x['chrom']}:{x['start']}-{x['end']}", axis=1)
    
    # Handle classification
    if 'class' in enhancers.columns and not enhancers['class'].isna().all():
        print("Using classification from BED file")
        enhancers['class'] = enhancers['class'].str.lower().str.strip()
        enhancers['class'] = enhancers['class'].replace({
            'hk': 'housekeeping',
            'dev': 'developmental',
            'tissue_specific': 'developmental',
            'tissue-specific': 'developmental'
        })
    elif classification_file:
        print(f"Loading classification from {classification_file}")
        classification = pd.read_csv(classification_file, sep='\t', header=None, 
                                   names=['name', 'class'])
        classification['class'] = classification['class'].str.lower().str.strip()
        enhancers = enhancers.merge(classification, on='name', how='left')
    else:
        print("Using heuristic classification")
        housekeeping_markers = ['ubiq', 'house', 'const', 'rp', 'ef1', 'gapdh', 'actb', 'tubulin', 'actin']
        enhancers['class'] = enhancers['name'].apply(
            lambda x: 'housekeeping' if any(marker in str(x).lower() for marker in housekeeping_markers)
            else 'developmental'
        )
    
    enhancers['class'] = enhancers['class'].fillna('developmental')
    
    print(f"Classification counts:")
    print(f"  Housekeeping: {sum(enhancers['class'] == 'housekeeping')}")
    print(f"  Developmental: {sum(enhancers['class'] == 'developmental')}")
    
    return enhancers


def find_enhancer_interactions(interactions, enhancers):
    """Find enhancer-enhancer and enhancer-TSS interactions."""
    print("\nFinding enhancer interactions...")
    
    # Clean coordinates
    interactions = interactions.copy()
    for col in ['start1', 'end1', 'start2', 'end2']:
        interactions[col] = pd.to_numeric(interactions[col], errors='coerce')
    interactions = interactions.dropna(subset=['start1', 'end1', 'start2', 'end2'])
    for col in ['start1', 'end1', 'start2', 'end2']:
        interactions[col] = interactions[col].astype(int)
    
    enhancers = enhancers.copy()
    enhancers['start'] = enhancers['start'].astype(int)
    enhancers['end'] = enhancers['end'].astype(int)
    
    print(f"Processing {len(interactions)} interactions...")
    
    # Index enhancers by chromosome
    enh_by_chr = defaultdict(list)
    for idx, enh in enhancers.iterrows():
        enh_by_chr[enh['chrom']].append({
            'start': enh['start'],
            'end': enh['end'],
            'name': enh['name'],
            'class': enh['class']
        })
    
    def find_overlapping_enhancers(chrom, start, end):
        if chrom not in enh_by_chr:
            return []
        overlapping = []
        for enh in enh_by_chr[chrom]:
            if not (end <= enh['start'] or start >= enh['end']):
                overlapping.append(enh)
        return overlapping
    
    enh_interactions = []
    progress_step = max(len(interactions) // 20, 1)
    
    for idx, interaction in interactions.iterrows():
        if idx % progress_step == 0:
            print(f"  Progress: {idx}/{len(interactions)} ({100*idx/len(interactions):.1f}%)")
        
        anchor1_enhancers = find_overlapping_enhancers(
            interaction['chr1'], interaction['start1'], interaction['end1']
        )
        anchor2_enhancers = find_overlapping_enhancers(
            interaction['chr2'], interaction['start2'], interaction['end2']
        )
        
        # E-E interactions
        if anchor1_enhancers and anchor2_enhancers:
            for enh1 in anchor1_enhancers:
                for enh2 in anchor2_enhancers:
                    if enh1['class'] == enh2['class']:
                        contact_class = enh1['class']
                    else:
                        contact_class = 'cross_class'
                    
                    interaction_data = interaction.copy()
                    interaction_data['interaction_type'] = 'E-E'
                    interaction_data['contact_class'] = contact_class
                    interaction_data['enh1_name'] = enh1['name']
                    interaction_data['enh1_class'] = enh1['class']
                    interaction_data['enh2_name'] = enh2['name']
                    interaction_data['enh2_class'] = enh2['class']
                    enh_interactions.append(interaction_data)
        
        # E-TSS interactions
        elif anchor1_enhancers or anchor2_enhancers:
            for enh in (anchor1_enhancers + anchor2_enhancers):
                interaction_data = interaction.copy()
                interaction_data['interaction_type'] = 'E-TSS'
                interaction_data['contact_class'] = enh['class']
                interaction_data['enh_name'] = enh['name']
                interaction_data['enh_class'] = enh['class']
                enh_interactions.append(interaction_data)
    
    if enh_interactions:
        enh_df = pd.DataFrame(enh_interactions)
        enh_df = enh_df.drop_duplicates()
        
        print(f"\nFound {len(enh_df)} enhancer interactions:")
        print(f"  E-E interactions: {sum(enh_df['interaction_type'] == 'E-E')}")
        print(f"  E-TSS interactions: {sum(enh_df['interaction_type'] == 'E-TSS')}")
        print("  By class:")
        print(enh_df['contact_class'].value_counts())
        
        return enh_df
    else:
        print("No enhancer interactions found")
        return pd.DataFrame()


def compare_to_null_model(enh_interactions, null_stats, fdr_threshold=0.05):
    """
    Compare enhancer interactions to null model.
    
    Tests whether observed logFC values for different enhancer classes
    are significantly different from the null distribution.
    """
    print("\n" + "="*70)
    print("STATISTICAL COMPARISON TO NULL MODEL")
    print("="*70)
    
    if enh_interactions.empty:
        print("No enhancer interactions to analyze")
        return pd.DataFrame()
    
    null_mean = null_stats['mean']
    null_std = null_stats['std']
    
    print(f"\nNull expectation: mean={null_mean:.4f}, std={null_std:.4f}")
    print(f"\nTesting if enhancer interactions deviate from null model...\n")
    
    results = []
    
    # Compare each combination of interaction type and contact class
    for int_type in ['E-E', 'E-TSS']:
        type_data = enh_interactions[enh_interactions['interaction_type'] == int_type]
        
        if type_data.empty:
            continue
        
        for contact_class in type_data['contact_class'].unique():
            class_data = type_data[type_data['contact_class'] == contact_class]
            
            if len(class_data) < 3:
                continue
            
            # Get logFC values
            logfc_values = pd.to_numeric(class_data['logFC'], errors='coerce').dropna()
            
            if len(logfc_values) == 0:
                continue
            
            print(f"{int_type} - {contact_class}:")
            print(f"  n={len(logfc_values)}")
            print(f"  mean logFC={logfc_values.mean():.3f}")
            print(f"  median logFC={logfc_values.median():.3f}")
            
            try:
                # Method 1: One-sample t-test against null mean
                t_stat, p_t_test = stats.ttest_1samp(logfc_values, null_mean)
                print(f"  T-test vs null: t={t_stat:.3f}, p={p_t_test:.4e}")
                
                # Method 2: Z-test against null distribution
                if null_std > 0:
                    z_score = (logfc_values.mean() - null_mean) / (null_std / np.sqrt(len(logfc_values)))
                    p_z_test = 2 * (1 - stats.norm.cdf(abs(z_score)))
                    print(f"  Z-test vs null: z={z_score:.3f}, p={p_z_test:.4e}")
                else:
                    z_score = 0
                    p_z_test = 1.0
                
                # Method 3: Permutation test
                n_permutations = 1000
                null_samples = np.random.normal(null_mean, null_std, 
                                               (n_permutations, len(logfc_values)))
                null_means = null_samples.mean(axis=1)
                p_perm = np.sum(np.abs(null_means - null_mean) >= 
                               np.abs(logfc_values.mean() - null_mean)) / n_permutations
                print(f"  Permutation test: p={p_perm:.4f}")
                
                # Bootstrap confidence interval
                n_bootstrap = 1000
                bootstrap_means = []
                for _ in range(n_bootstrap):
                    boot_sample = np.random.choice(logfc_values, len(logfc_values), replace=True)
                    bootstrap_means.append(boot_sample.mean())
                
                ci_lower = np.percentile(bootstrap_means, 2.5)
                ci_upper = np.percentile(bootstrap_means, 97.5)
                ci_excludes_null = not (ci_lower <= null_mean <= ci_upper)
                print(f"  95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
                print(f"  CI excludes null: {ci_excludes_null}")
                
                # Effect size (Cohen's d)
                cohens_d = (logfc_values.mean() - null_mean) / logfc_values.std()
                print(f"  Cohen's d: {cohens_d:.3f}")
                
                # Standardized effect size
                effect_size = (logfc_values.mean() - null_mean) / null_std
                print(f"  Effect size (std from null): {effect_size:.3f}\n")
                
                # Biological significance
                biologically_significant = abs(logfc_values.mean()) > 0.5
                
                results.append({
                    'interaction_type': int_type,
                    'contact_class': contact_class,
                    'n_interactions': len(logfc_values),
                    'mean_logfc': logfc_values.mean(),
                    'median_logfc': logfc_values.median(),
                    'std_logfc': logfc_values.std(),
                    'null_mean': null_mean,
                    'null_std': null_std,
                    'deviation_from_null': logfc_values.mean() - null_mean,
                    't_statistic': t_stat,
                    'p_value_t_test': p_t_test,
                    'z_score': z_score,
                    'p_value_z_test': p_z_test,
                    'p_value_permutation': p_perm,
                    'ci_lower': ci_lower,
                    'ci_upper': ci_upper,
                    'ci_excludes_null': ci_excludes_null,
                    'cohens_d': cohens_d,
                    'effect_size': effect_size,
                    'biologically_significant': biologically_significant
                })
                
            except Exception as e:
                print(f"  Error in analysis: {e}\n")
                continue
    
    if not results:
        print("No comparisons could be performed")
        return pd.DataFrame()
    
    # Convert to DataFrame and apply FDR correction
    results_df = pd.DataFrame(results)
    
    print("="*70)
    print("FDR CORRECTION")
    print("="*70)
    
    for p_col in ['p_value_t_test', 'p_value_z_test', 'p_value_permutation']:
        if p_col in results_df.columns:
            _, fdr_values, _, _ = multipletests(results_df[p_col], method='fdr_bh')
            results_df[f'{p_col}_fdr'] = fdr_values
            results_df[f'{p_col}_significant'] = fdr_values < fdr_threshold
            
            n_sig = sum(fdr_values < fdr_threshold)
            print(f"{p_col}: {n_sig}/{len(results_df)} significant (FDR < {fdr_threshold})")
    
    # Combined significance
    sig_cols = [col for col in results_df.columns if col.endswith('_significant')]
    if sig_cols:
        results_df['any_significant'] = results_df[sig_cols].any(axis=1)
        print(f"\nOverall: {sum(results_df['any_significant'])}/{len(results_df)} significant by any method")
    
    # Add interpretation
    def interpret_result(row):
        interp = []
        if row.get('any_significant', False):
            interp.append("statistically_significant")
        if row.get('biologically_significant', False):
            interp.append("biologically_significant")
        
        effect = abs(row.get('effect_size', 0))
        if effect > 1:
            interp.append("large_effect")
        elif effect > 0.5:
            interp.append("medium_effect")
        elif effect > 0.2:
            interp.append("small_effect")
        else:
            interp.append("negligible_effect")
        
        if row.get('mean_logfc', 0) > null_mean:
            interp.append("increased_in_wMel")
        elif row.get('mean_logfc', 0) < null_mean:
            interp.append("decreased_in_wMel")
        
        return "; ".join(interp)
    
    results_df['interpretation'] = results_df.apply(interpret_result, axis=1)
    
    return results_df


def create_visualizations(enh_interactions, comparison_results, null_model, output_prefix):
    """Create comprehensive visualizations."""
    
    if enh_interactions.empty:
        print("No data to visualize")
        return
    
    print("\nCreating visualizations...")
    
    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.35)
    
    # Plot 1: LogFC distribution vs null
    ax1 = fig.add_subplot(gs[0, 0])
    null_logfc = pd.to_numeric(null_model['logFC'], errors='coerce').dropna()
    obs_logfc = pd.to_numeric(enh_interactions['logFC'], errors='coerce').dropna()
    
    ax1.hist(null_logfc, bins=50, alpha=0.5, label='Null model', color='lightblue', density=True)
    ax1.hist(obs_logfc, bins=50, alpha=0.5, label='Enhancer interactions', color='salmon', density=True)
    ax1.axvline(x=null_logfc.mean(), color='blue', linestyle='--', label=f'Null mean: {null_logfc.mean():.3f}')
    ax1.axvline(x=obs_logfc.mean(), color='red', linestyle='--', label=f'Observed mean: {obs_logfc.mean():.3f}')
    ax1.set_xlabel('LogFC')
    ax1.set_ylabel('Density')
    ax1.set_title('LogFC Distribution: Enhancers vs Null Model')
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)
    
    # Plot 2: LogFC by enhancer class
    ax2 = fig.add_subplot(gs[0, 1])
    for contact_class in enh_interactions['contact_class'].unique():
        class_data = enh_interactions[enh_interactions['contact_class'] == contact_class]
        class_logfc = pd.to_numeric(class_data['logFC'], errors='coerce').dropna()
        if len(class_logfc) > 0:
            ax2.hist(class_logfc, bins=30, alpha=0.6, label=contact_class)
    ax2.axvline(x=null_logfc.mean(), color='black', linestyle='--', alpha=0.5, label='Null mean')
    ax2.set_xlabel('LogFC')
    ax2.set_ylabel('Frequency')
    ax2.set_title('LogFC Distribution by Enhancer Class')
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    
    # Plot 3: LogFC by interaction type
    ax3 = fig.add_subplot(gs[0, 2])
    for int_type in enh_interactions['interaction_type'].unique():
        type_data = enh_interactions[enh_interactions['interaction_type'] == int_type]
        type_logfc = pd.to_numeric(type_data['logFC'], errors='coerce').dropna()
        if len(type_logfc) > 0:
            ax3.hist(type_logfc, bins=30, alpha=0.6, label=int_type)
    ax3.axvline(x=null_logfc.mean(), color='black', linestyle='--', alpha=0.5, label='Null mean')
    ax3.set_xlabel('LogFC')
    ax3.set_ylabel('Frequency')
    ax3.set_title('LogFC Distribution by Interaction Type')
    ax3.legend(fontsize=8)
    ax3.grid(alpha=0.3)
    
    # Plot 4: Interaction counts
    ax4 = fig.add_subplot(gs[0, 3])
    count_data = enh_interactions.groupby(['interaction_type', 'contact_class']).size().unstack(fill_value=0)
    count_data.plot(kind='bar', ax=ax4, width=0.8)
    ax4.set_title('Enhancer Interaction Counts')
    ax4.set_ylabel('Count')
    ax4.set_xlabel('Interaction Type')
    ax4.legend(title='Contact Class', fontsize=8)
    ax4.grid(alpha=0.3, axis='y')
    plt.setp(ax4.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    if not comparison_results.empty:
        # Plot 5: Effect sizes
        ax5 = fig.add_subplot(gs[1, 0])
        sig_data = comparison_results[comparison_results.get('any_significant', False)]
        nonsig_data = comparison_results[~comparison_results.get('any_significant', True)]
        
        if not sig_data.empty:
            ax5.scatter(sig_data['deviation_from_null'], sig_data['effect_size'],
                       c='red', s=100, alpha=0.7, label='Significant', edgecolors='black')
        if not nonsig_data.empty:
            ax5.scatter(nonsig_data['deviation_from_null'], nonsig_data['effect_size'],
                       c='gray', s=60, alpha=0.4, label='Not significant')
        
        ax5.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5, label='Medium effect')
        ax5.axhline(y=1.0, color='red', linestyle=':', alpha=0.5, label='Large effect')
        ax5.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        ax5.set_xlabel('Deviation from Null Mean')
        ax5.set_ylabel('Effect Size (standardized)')
        ax5.set_title('Effect Sizes vs Null Model')
        ax5.legend(fontsize=8)
        ax5.grid(alpha=0.3)
        
        # Plot 6: P-values
        ax6 = fig.add_subplot(gs[1, 1])
        if 'p_value_t_test' in comparison_results.columns:
            ax6.scatter(comparison_results['p_value_t_test'],
                       comparison_results.get('p_value_z_test', comparison_results['p_value_t_test']),
                       alpha=0.6, s=80, edgecolors='black', linewidths=0.5)
            ax6.plot([1e-10, 1], [1e-10, 1], 'k--', alpha=0.5)
            ax6.axhline(y=0.05, color='red', linestyle=':', alpha=0.5)
            ax6.axvline(x=0.05, color='red', linestyle=':', alpha=0.5)
            ax6.set_xlabel('P-value (T-test)')
            ax6.set_ylabel('P-value (Z-test)')
            ax6.set_title('P-value Comparison')
            ax6.set_xscale('log')
            ax6.set_yscale('log')
            ax6.grid(alpha=0.3)
        
        # Plot 7: Mean logFC by group
        ax7 = fig.add_subplot(gs[1, 2])
        comparison_results['group'] = comparison_results['interaction_type'] + '\n' + comparison_results['contact_class']
        colors = ['red' if sig else 'gray' for sig in comparison_results.get('any_significant', [False]*len(comparison_results))]
        
        x_pos = range(len(comparison_results))
        ax7.bar(x_pos, comparison_results['mean_logfc'], color=colors, alpha=0.7, edgecolor='black')
        ax7.axhline(y=comparison_results['null_mean'].iloc[0], color='black', 
                   linestyle='--', linewidth=2, label='Null mean')
        ax7.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5)
        ax7.axhline(y=-0.5, color='orange', linestyle=':', alpha=0.5)
        ax7.set_xticks(x_pos)
        ax7.set_xticklabels(comparison_results['group'], rotation=45, ha='right', fontsize=8)
        ax7.set_ylabel('Mean LogFC')
        ax7.set_title('Mean LogFC by Group')
        ax7.legend(fontsize=8)
        ax7.grid(alpha=0.3, axis='y')
        
        # Plot 8: Confidence intervals
        ax8 = fig.add_subplot(gs[1, 3])
        y_pos = range(len(comparison_results))
        colors_ci = ['red' if sig else 'gray' for sig in comparison_results.get('any_significant', [False]*len(comparison_results))]
        
        for i, (idx, row) in enumerate(comparison_results.iterrows()):
            ax8.plot([row['ci_lower'], row['ci_upper']], [i, i], 
                    color=colors_ci[i], linewidth=3, alpha=0.7)
            ax8.scatter(row['mean_logfc'], i, color=colors_ci[i], s=80, 
                       edgecolor='black', linewidths=1, zorder=3)
        
        ax8.axvline(x=comparison_results['null_mean'].iloc[0], color='black', 
                   linestyle='--', linewidth=2, label='Null mean')
        ax8.set_yticks(y_pos)
        ax8.set_yticklabels(comparison_results['group'], fontsize=8)
        ax8.set_xlabel('LogFC')
        ax8.set_title('95% Confidence Intervals')
        ax8.legend(fontsize=8)
        ax8.grid(alpha=0.3, axis='x')
        
        # Plot 9: Effect size comparison
        ax9 = fig.add_subplot(gs[2, :2])
        groups = comparison_results['group']
        x_pos = range(len(groups))
        
        ax9.bar(x_pos, comparison_results['effect_size'], color=colors, alpha=0.7, edgecolor='black')
        ax9.axhline(y=0.2, color='blue', linestyle=':', alpha=0.5, label='Small effect')
        ax9.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5, label='Medium effect')
        ax9.axhline(y=1.0, color='red', linestyle=':', alpha=0.5, label='Large effect')
        ax9.set_xticks(x_pos)
        ax9.set_xticklabels(groups, rotation=45, ha='right', fontsize=8)
        ax9.set_ylabel('Effect Size (std from null)')
        ax9.set_title('Effect Sizes Across Groups')
        ax9.legend(fontsize=8)
        ax9.grid(alpha=0.3, axis='y')
        
        # Plot 10: Summary barplot
        ax10 = fig.add_subplot(gs[2, 2:])
        summary_data = {
            'Total\nInteractions': len(enh_interactions),
            'E-E': sum(enh_interactions['interaction_type'] == 'E-E'),
            'E-TSS': sum(enh_interactions['interaction_type'] == 'E-TSS'),
            'Housekeeping': sum(enh_interactions['contact_class'] == 'housekeeping'),
            'Developmental': sum(enh_interactions['contact_class'] == 'developmental'),
            'Significant\nGroups': sum(comparison_results.get('any_significant', [False]*len(comparison_results)))
        }
        
        ax10.bar(range(len(summary_data)), list(summary_data.values()), 
                color=['steelblue', 'lightcoral', 'lightgreen', 'gold', 'plum', 'tomato'],
                alpha=0.8, edgecolor='black')
        ax10.set_xticks(range(len(summary_data)))
        ax10.set_xticklabels(list(summary_data.keys()), rotation=45, ha='right')
        ax10.set_ylabel('Count')
        ax10.set_title('Summary Statistics')
        ax10.grid(alpha=0.3, axis='y')
        
        # Add values on bars
        for i, (key, val) in enumerate(summary_data.items()):
            ax10.text(i, val, str(val), ha='center', va='bottom', fontweight='bold')
    
    plt.savefig(f"{output_prefix}_comprehensive_analysis.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_prefix}_comprehensive_analysis.pdf")


def create_report(enh_interactions, comparison_results, null_model, null_stats, output_prefix):
    """Create detailed text report."""
    
    with open(f"{output_prefix}_analysis_report.txt", 'w') as f:
        f.write("="*80 + "\n")
        f.write("ENHANCER CLASS ANALYSIS REPORT\n")
        f.write("wMel vs DOX Comparison\n")
        f.write("="*80 + "\n\n")
        
        # Null model
        f.write("NULL MODEL SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write(f"Null model size: {null_stats['n']} interactions\n")
        f.write(f"Mean logFC: {null_stats['mean']:.4f}\n")
        f.write(f"Std logFC: {null_stats['std']:.4f}\n")
        f.write(f"Median: {null_stats['median']:.4f}\n")
        f.write(f"Range: [{null_stats['min']:.4f}, {null_stats['max']:.4f}]\n\n")
        
        # Enhancer interactions
        f.write("ENHANCER INTERACTIONS SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write(f"Total enhancer interactions: {len(enh_interactions)}\n")
        f.write(f"E-E interactions: {sum(enh_interactions['interaction_type'] == 'E-E')}\n")
        f.write(f"E-TSS interactions: {sum(enh_interactions['interaction_type'] == 'E-TSS')}\n\n")
        
        f.write("By contact class:\n")
        for cls in enh_interactions['contact_class'].unique():
            n = sum(enh_interactions['contact_class'] == cls)
            pct = 100 * n / len(enh_interactions)
            f.write(f"  {cls}: {n} ({pct:.1f}%)\n")
        
        obs_logfc = pd.to_numeric(enh_interactions['logFC'], errors='coerce').dropna()
        f.write(f"\nOverall logFC statistics:\n")
        f.write(f"  Mean: {obs_logfc.mean():.3f}\n")
        f.write(f"  Median: {obs_logfc.median():.3f}\n")
        f.write(f"  Std: {obs_logfc.std():.3f}\n")
        f.write(f"  Range: [{obs_logfc.min():.3f}, {obs_logfc.max():.3f}]\n\n")
        
        # Comparison results
        if not comparison_results.empty:
            f.write("STATISTICAL COMPARISONS TO NULL MODEL\n")
            f.write("-"*80 + "\n")
            f.write(f"Total comparisons: {len(comparison_results)}\n")
            
            if 'any_significant' in comparison_results.columns:
                n_sig = sum(comparison_results['any_significant'])
                f.write(f"Significant comparisons: {n_sig}/{len(comparison_results)}\n\n")
            
            f.write("Effect size distribution:\n")
            if 'effect_size' in comparison_results.columns:
                large = sum(comparison_results['effect_size'].abs() > 1)
                medium = sum((comparison_results['effect_size'].abs() > 0.5) & 
                           (comparison_results['effect_size'].abs() <= 1))
                small = sum((comparison_results['effect_size'].abs() > 0.2) & 
                          (comparison_results['effect_size'].abs() <= 0.5))
                f.write(f"  Large effects (>1 std): {large}\n")
                f.write(f"  Medium effects (0.5-1 std): {medium}\n")
                f.write(f"  Small effects (0.2-0.5 std): {small}\n\n")
            
            # Detailed significant results
            if 'any_significant' in comparison_results.columns:
                sig_results = comparison_results[comparison_results['any_significant']]
                if not sig_results.empty:
                    f.write("\n" + "="*80 + "\n")
                    f.write("SIGNIFICANT RESULTS DETAILS\n")
                    f.write("="*80 + "\n\n")
                    
                    for idx, row in sig_results.iterrows():
                        f.write(f"\n{row['interaction_type']} - {row['contact_class']}\n")
                        f.write("-"*40 + "\n")
                        f.write(f"  N interactions: {row['n_interactions']}\n")
                        f.write(f"  Mean logFC: {row['mean_logfc']:.3f}\n")
                        f.write(f"  Null expectation: {row['null_mean']:.3f}\n")
                        f.write(f"  Deviation: {row['deviation_from_null']:.3f}\n")
                        f.write(f"  Effect size: {row['effect_size']:.3f} std\n")
                        f.write(f"  Cohen's d: {row['cohens_d']:.3f}\n")
                        f.write(f"  P-value (t-test): {row['p_value_t_test']:.4e}\n")
                        f.write(f"  P-value (z-test): {row['p_value_z_test']:.4e}\n")
                        f.write(f"  P-value (permutation): {row['p_value_permutation']:.4f}\n")
                        f.write(f"  95% CI: [{row['ci_lower']:.3f}, {row['ci_upper']:.3f}]\n")
                        f.write(f"  CI excludes null: {row['ci_excludes_null']}\n")
                        f.write(f"  Interpretation: {row['interpretation']}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"Report saved to {output_prefix}_analysis_report.txt")


def main():
    parser = argparse.ArgumentParser(
        description='Enhancer class analysis for wMel vs DOX comparison',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--enhancers', required=True,
                       help='Enhancer BED file')
    parser.add_argument('--interactions', required=True,
                       help='Differential interactions CSV (wMel vs DOX)')
    parser.add_argument('--null_model', required=True,
                       help='Null model CSV from diffHic')
    parser.add_argument('--classification',
                       help='Optional enhancer classification file')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold (default: 0.05)')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    print("\n" + "="*80)
    print("ENHANCER CLASS ANALYSIS: wMel vs DOX")
    print("="*80 + "\n")
    
    # Load data
    interactions = load_interactions(args.interactions)
    null_model, null_stats = load_null_model(args.null_model)
    enhancers = classify_enhancers(args.enhancers, args.classification)
    
    if enhancers.empty:
        print("ERROR: No valid enhancers!")
        return
    
    # Find enhancer interactions
    enh_interactions = find_enhancer_interactions(interactions, enhancers)
    
    if enh_interactions.empty:
        print("WARNING: No enhancer interactions found!")
        return
    
    # Statistical comparison
    comparison_results = compare_to_null_model(enh_interactions, null_stats, args.fdr_threshold)
    
    # Create outputs
    create_visualizations(enh_interactions, comparison_results, null_model, args.output_prefix)
    create_report(enh_interactions, comparison_results, null_model, null_stats, args.output_prefix)
    
    # Save results
    print("\nSaving results...")
    
    enh_interactions.to_csv(f"{args.output_prefix}_enhancer_interactions.tsv", 
                           sep='\t', index=False)
    print(f"  Saved: {args.output_prefix}_enhancer_interactions.tsv")
    
    if not comparison_results.empty:
        comparison_results.to_csv(f"{args.output_prefix}_statistical_comparisons.tsv",
                                sep='\t', index=False)
        print(f"  Saved: {args.output_prefix}_statistical_comparisons.tsv")
        
        if 'any_significant' in comparison_results.columns:
            sig_results = comparison_results[comparison_results['any_significant']]
            if not sig_results.empty:
                sig_results.to_csv(f"{args.output_prefix}_significant_results.tsv",
                                 sep='\t', index=False)
                print(f"  Saved: {args.output_prefix}_significant_results.tsv")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nResults written to: {args.output_prefix}*")


if __name__ == '__main__':
    main()
