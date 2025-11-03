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
import itertools

# Color scheme
WMEL_COLOR = '#09aa4b'      # wMel (infected) - negative logFC
UNINF_COLOR = '#8fcb84'     # Uninfected (DOX) - positive logFC
HK_COLOR = '#e74c3c'        # Housekeeping
DEV_COLOR = '#3498db'       # Developmental
CROSS_COLOR = '#9b59b6'     # Cross-class
SIG_COLOR = '#09aa4b'       # Significant results
NONSIG_COLOR = '#999999'    # Non-significant results

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


def compare_enhancer_classes(enh_interactions, fdr_threshold=0.05):
    """
    Compare enhancer classes to each other.
    
    Tests whether different enhancer classes show significantly different
    logFC distributions in response to Wolbachia infection.
    """
    print("\n" + "="*70)
    print("STATISTICAL COMPARISON BETWEEN ENHANCER CLASSES")
    print("="*70)
    
    if enh_interactions.empty:
        print("No enhancer interactions to analyze")
        return pd.DataFrame()
    
    obs_logfc = pd.to_numeric(enh_interactions['logFC'], errors='coerce').dropna()
    print(f"\nOverall distribution:")
    print(f"  Total interactions: {len(enh_interactions)}")
    print(f"  Mean logFC: {obs_logfc.mean():.3f}")
    print(f"  Median logFC: {obs_logfc.median():.3f}")
    print(f"  Std logFC: {obs_logfc.std():.3f}")
    print(f"  Range: [{obs_logfc.min():.3f}, {obs_logfc.max():.3f}]")
    
    results = []
    
    # Compare within each interaction type
    for int_type in ['E-E', 'E-TSS']:
        type_data = enh_interactions[enh_interactions['interaction_type'] == int_type]
        
        if type_data.empty:
            continue
        
        print(f"\n{'='*70}")
        print(f"ANALYZING {int_type} INTERACTIONS")
        print('='*70)
        
        # Get available classes for this interaction type
        available_classes = type_data['contact_class'].unique()
        
        # Compare all pairs of classes
        for class1, class2 in itertools.combinations(available_classes, 2):
            data1 = pd.to_numeric(type_data[type_data['contact_class'] == class1]['logFC'], 
                                 errors='coerce').dropna()
            data2 = pd.to_numeric(type_data[type_data['contact_class'] == class2]['logFC'], 
                                 errors='coerce').dropna()
            
            if len(data1) < 10 or len(data2) < 10:
                continue
            
            print(f"\n{class1} vs {class2}:")
            print(f"  {class1}: n={len(data1)}, mean={data1.mean():.3f}, median={data1.median():.3f}")
            print(f"  {class2}: n={len(data2)}, mean={data2.mean():.3f}, median={data2.median():.3f}")
            
            # Calculate statistics
            mean_diff = data1.mean() - data2.mean()
            print(f"  Mean difference: {mean_diff:.4f}")
            
            # T-test
            t_stat, p_t_test = stats.ttest_ind(data1, data2)
            print(f"  T-test: t={t_stat:.3f}, p={p_t_test:.4e}")
            
            # Mann-Whitney U test (non-parametric)
            u_stat, p_mann = stats.mannwhitneyu(data1, data2, alternative='two-sided')
            print(f"  Mann-Whitney U: U={u_stat:.0f}, p={p_mann:.4e}")
            
            # Effect size (Cohen's d)
            pooled_std = np.sqrt((data1.std()**2 + data2.std()**2) / 2)
            cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
            print(f"  Cohen's d: {cohens_d:.4f}")
            
            # Interpret effect size
            if abs(cohens_d) < 0.2:
                effect_interp = "negligible"
            elif abs(cohens_d) < 0.5:
                effect_interp = "small"
            elif abs(cohens_d) < 0.8:
                effect_interp = "medium"
            else:
                effect_interp = "large"
            print(f"  Effect size: {effect_interp}")
            
            # Bootstrap confidence interval for mean difference
            n_bootstrap = 1000
            bootstrap_diffs = []
            for _ in range(n_bootstrap):
                boot1 = np.random.choice(data1, len(data1), replace=True)
                boot2 = np.random.choice(data2, len(data2), replace=True)
                bootstrap_diffs.append(boot1.mean() - boot2.mean())
            
            ci_lower = np.percentile(bootstrap_diffs, 2.5)
            ci_upper = np.percentile(bootstrap_diffs, 97.5)
            ci_excludes_zero = not (ci_lower <= 0 <= ci_upper)
            print(f"  95% CI for difference: [{ci_lower:.4f}, {ci_upper:.4f}]")
            print(f"  CI excludes zero: {ci_excludes_zero}")
            
            # Biological significance
            biologically_significant = abs(mean_diff) > 0.5
            print(f"  Biologically significant (|diff| > 0.5): {biologically_significant}")
            
            results.append({
                'interaction_type': int_type,
                'class1': class1,
                'class2': class2,
                'n1': len(data1),
                'n2': len(data2),
                'mean1': data1.mean(),
                'mean2': data2.mean(),
                'median1': data1.median(),
                'median2': data2.median(),
                'std1': data1.std(),
                'std2': data2.std(),
                'mean_difference': mean_diff,
                't_statistic': t_stat,
                'p_value_ttest': p_t_test,
                'u_statistic': u_stat,
                'p_value_mannwhitney': p_mann,
                'cohens_d': cohens_d,
                'effect_size_category': effect_interp,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'ci_excludes_zero': ci_excludes_zero,
                'biologically_significant': biologically_significant
            })
    
    if not results:
        print("\nNo comparisons could be performed")
        return pd.DataFrame()
    
    # Convert to DataFrame and apply FDR correction
    results_df = pd.DataFrame(results)
    
    print("\n" + "="*70)
    print("FDR CORRECTION")
    print("="*70)
    
    for p_col in ['p_value_ttest', 'p_value_mannwhitney']:
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
        else:
            interp.append("not_statistically_significant")
        
        if row.get('biologically_significant', False):
            interp.append("biologically_significant")
        else:
            interp.append("not_biologically_significant")
        
        interp.append(f"{row.get('effect_size_category', 'unknown')}_effect")
        
        if row.get('mean_difference', 0) > 0:
            interp.append(f"{row['class1']}_higher_than_{row['class2']}")
        else:
            interp.append(f"{row['class2']}_higher_than_{row['class1']}")
        
        return "; ".join(interp)
    
    results_df['interpretation'] = results_df.apply(interpret_result, axis=1)
    
    return results_df


def format_pvalue(pval):
    """Format p-value for display"""
    if pval < 0.001:
        return f"p < 0.001"
    elif pval < 0.01:
        return f"p = {pval:.3f}"
    else:
        return f"p = {pval:.2f}"


def create_visualizations(enh_interactions, comparison_results, output_prefix):
    """Create comprehensive visualizations as separate PDF files."""
    
    if enh_interactions.empty:
        print("No data to visualize")
        return
    
    print("\nCreating visualizations...")
    
    obs_logfc = pd.to_numeric(enh_interactions['logFC'], errors='coerce').dropna()
    
    # Get wMel-enriched and uninf-enriched interactions
    wmel_enriched = enh_interactions[enh_interactions['logFC'] < 0].copy()
    uninf_enriched = enh_interactions[enh_interactions['logFC'] > 0].copy()
    
    # PLOT 1: Overall LogFC distribution
    print("  Creating plot 1: Overall LogFC distribution...")
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(obs_logfc, bins=100, alpha=0.7, color=WMEL_COLOR, edgecolor='black', linewidth=0.5)
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='No change (logFC=0)')
    ax.axvline(x=obs_logfc.mean(), color='black', linestyle='--', linewidth=2, 
              label=f'Mean: {obs_logfc.mean():.2f}')
    ax.axvline(x=obs_logfc.median(), color='gray', linestyle=':', linewidth=2,
              label=f'Median: {obs_logfc.median():.2f}')
    
    # Add stats box
    stats_text = f'n = {len(obs_logfc):,}\nStd = {obs_logfc.std():.2f}\nRange = [{obs_logfc.min():.2f}, {obs_logfc.max():.2f}]'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('Log2 Fold Change (wMel vs uninf)\nNegative = wMel-enriched | Positive = Uninf-enriched', 
                  fontsize=12, fontweight='bold')
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_title('Overall LogFC Distribution: Enhancer Interactions', 
                fontsize=14, fontweight='bold', pad=20)
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(alpha=0.3, linestyle=':', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_plot1_overall_distribution.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # PLOT 2: LogFC distribution by enhancer class
    print("  Creating plot 2: LogFC by enhancer class...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    for idx, int_type in enumerate(['E-E', 'E-TSS']):
        ax = axes[idx]
        type_data = enh_interactions[enh_interactions['interaction_type'] == int_type]
        
        if type_data.empty:
            continue
        
        classes = sorted(type_data['contact_class'].unique())
        colors_dict = {'housekeeping': HK_COLOR, 'developmental': DEV_COLOR, 'cross_class': CROSS_COLOR}
        
        for cls in classes:
            class_data = pd.to_numeric(type_data[type_data['contact_class'] == cls]['logFC'], 
                                      errors='coerce').dropna()
            color = colors_dict.get(cls, '#888888')
            ax.hist(class_data, bins=50, alpha=0.5, label=f'{cls} (n={len(class_data):,})', 
                   color=color, edgecolor='black', linewidth=0.5)
        
        ax.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.5)
        ax.set_xlabel('Log2 Fold Change', fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_title(f'{int_type} Interactions', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, framealpha=0.9)
        ax.grid(alpha=0.3, linestyle=':', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_plot2_distribution_by_class.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # PLOT 3: Box plots comparing classes
    print("  Creating plot 3: Box plots by class...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for idx, int_type in enumerate(['E-E', 'E-TSS']):
        ax = axes[idx]
        type_data = enh_interactions[enh_interactions['interaction_type'] == int_type].copy()
        
        if type_data.empty:
            continue
        
        # Prepare data for boxplot
        classes = sorted(type_data['contact_class'].unique())
        data_by_class = [pd.to_numeric(type_data[type_data['contact_class'] == cls]['logFC'],
                                       errors='coerce').dropna() for cls in classes]
        
        colors_list = [{'housekeeping': HK_COLOR, 'developmental': DEV_COLOR, 
                       'cross_class': CROSS_COLOR}.get(cls, '#888888') for cls in classes]
        
        bp = ax.boxplot(data_by_class, labels=classes, patch_artist=True,
                       showfliers=False, widths=0.6)
        
        for patch, color in zip(bp['boxes'], colors_list):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            patch.set_edgecolor('black')
            patch.set_linewidth(1.5)
        
        for element in ['whiskers', 'fliers', 'means', 'medians', 'caps']:
            plt.setp(bp[element], color='black', linewidth=1.5)
        
        ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.5, label='No change')
        
        # Add sample sizes
        for i, (cls, data) in enumerate(zip(classes, data_by_class)):
            ax.text(i+1, ax.get_ylim()[0], f'n={len(data):,}', 
                   ha='center', va='top', fontsize=8)
        
        ax.set_ylabel('Log2 Fold Change', fontsize=11, fontweight='bold')
        ax.set_xlabel('Enhancer Class', fontsize=11, fontweight='bold')
        ax.set_title(f'{int_type} Interactions', fontsize=12, fontweight='bold')
        ax.grid(alpha=0.3, axis='y', linestyle=':', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_plot3_boxplots_by_class.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # PLOT 4: Mean logFC comparison with error bars
    print("  Creating plot 4: Mean logFC by class...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for idx, int_type in enumerate(['E-E', 'E-TSS']):
        ax = axes[idx]
        type_data = enh_interactions[enh_interactions['interaction_type'] == int_type]
        
        if type_data.empty:
            continue
        
        classes = sorted(type_data['contact_class'].unique())
        means = []
        sems = []
        colors_list = []
        
        for cls in classes:
            class_data = pd.to_numeric(type_data[type_data['contact_class'] == cls]['logFC'],
                                      errors='coerce').dropna()
            means.append(class_data.mean())
            sems.append(class_data.sem())
            colors_list.append({'housekeeping': HK_COLOR, 'developmental': DEV_COLOR,
                              'cross_class': CROSS_COLOR}.get(cls, '#888888'))
        
        x_pos = np.arange(len(classes))
        bars = ax.bar(x_pos, means, yerr=sems, capsize=10, alpha=0.8,
                     color=colors_list, edgecolor='black', linewidth=1.5,
                     error_kw={'linewidth': 2, 'ecolor': 'black'})
        
        # Add value labels
        for i, (m, s) in enumerate(zip(means, sems)):
            ax.text(i, m + s + 0.1, f'{m:.3f}', ha='center', va='bottom',
                   fontsize=10, fontweight='bold')
        
        ax.axhline(y=0, color='red', linestyle='--', linewidth=2, alpha=0.5)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(classes, fontsize=10)
        ax.set_ylabel('Mean Log2 Fold Change ± SEM', fontsize=11, fontweight='bold')
        ax.set_xlabel('Enhancer Class', fontsize=11, fontweight='bold')
        ax.set_title(f'{int_type} Interactions', fontsize=12, fontweight='bold')
        ax.grid(alpha=0.3, axis='y', linestyle=':', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_plot4_mean_logfc_by_class.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    # PLOT 5: Interaction counts - wMel vs uninf by class
    print("  Creating plot 5: Interaction counts by direction...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for idx, int_type in enumerate(['E-E', 'E-TSS']):
        ax = axes[idx]
        type_data = enh_interactions[enh_interactions['interaction_type'] == int_type]
        
        if type_data.empty:
            continue
        
        classes = sorted(type_data['contact_class'].unique())
        wmel_counts = []
        uninf_counts = []
        
        for cls in classes:
            class_data = type_data[type_data['contact_class'] == cls]
            wmel_counts.append(len(class_data[class_data['logFC'] < 0]))
            uninf_counts.append(len(class_data[class_data['logFC'] > 0]))
        
        x_pos = np.arange(len(classes))
        width = 0.35
        
        bars1 = ax.bar(x_pos - width/2, wmel_counts, width, label='wMel-enriched',
                      color=WMEL_COLOR, alpha=0.8, edgecolor='black', linewidth=1.5)
        bars2 = ax.bar(x_pos + width/2, uninf_counts, width, label='Uninf-enriched',
                      color=UNINF_COLOR, alpha=0.8, edgecolor='black', linewidth=1.5)
        
        # Add value labels
        for bars in [bars1, bars2]:
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height):,}', ha='center', va='bottom', fontsize=9)
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(classes, fontsize=10)
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        ax.set_xlabel('Enhancer Class', fontsize=11, fontweight='bold')
        ax.set_title(f'{int_type} Interactions', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10, framealpha=0.9)
        ax.grid(alpha=0.3, axis='y', linestyle=':', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_plot5_counts_by_direction.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    if not comparison_results.empty:
        # PLOT 6: Effect sizes from comparisons
        print("  Creating plot 6: Effect sizes...")
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create comparison labels
        comparison_results['comparison'] = (comparison_results['interaction_type'] + '\n' +
                                           comparison_results['class1'] + ' vs ' +
                                           comparison_results['class2'])
        
        colors = [SIG_COLOR if sig else NONSIG_COLOR 
                 for sig in comparison_results.get('any_significant', [False]*len(comparison_results))]
        
        y_pos = range(len(comparison_results))
        
        bars = ax.barh(y_pos, comparison_results['cohens_d'], color=colors, alpha=0.8,
                      edgecolor='black', linewidth=1.5)
        
        # Add effect size thresholds
        ax.axvline(x=0.2, color='lightblue', linestyle=':', alpha=0.7, linewidth=2, label='Small (0.2)')
        ax.axvline(x=0.5, color='orange', linestyle=':', alpha=0.7, linewidth=2, label='Medium (0.5)')
        ax.axvline(x=0.8, color='red', linestyle=':', alpha=0.7, linewidth=2, label='Large (0.8)')
        ax.axvline(x=-0.2, color='lightblue', linestyle=':', alpha=0.7, linewidth=2)
        ax.axvline(x=-0.5, color='orange', linestyle=':', alpha=0.7, linewidth=2)
        ax.axvline(x=-0.8, color='red', linestyle=':', alpha=0.7, linewidth=2)
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.5, linewidth=1.5)
        
        # Add p-values
        for i, (idx, row) in enumerate(comparison_results.iterrows()):
            d_val = row['cohens_d']
            p_val = row['p_value_ttest']
            ax.text(d_val + 0.01 if d_val > 0 else d_val - 0.01, i,
                   format_pvalue(p_val), va='center',
                   ha='left' if d_val > 0 else 'right',
                   fontsize=8, fontweight='bold' if row.get('any_significant', False) else 'normal')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(comparison_results['comparison'], fontsize=9)
        ax.set_xlabel("Cohen's d (Effect Size)", fontsize=12, fontweight='bold')
        ax.set_title("Effect Sizes: Class Comparisons", fontsize=14, fontweight='bold', pad=20)
        
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=SIG_COLOR, edgecolor='black', label='Significant'),
            Patch(facecolor=NONSIG_COLOR, edgecolor='black', label='Not significant'),
        ]
        ax.legend(handles=legend_elements, fontsize=10, framealpha=0.9, loc='best')
        ax.grid(alpha=0.3, axis='x', linestyle=':', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_plot6_effect_sizes.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        
        # PLOT 7: Mean differences with confidence intervals
        print("  Creating plot 7: Mean differences with CIs...")
        fig, ax = plt.subplots(figsize=(12, 8))
        
        colors_ci = [SIG_COLOR if sig else NONSIG_COLOR 
                    for sig in comparison_results.get('any_significant', [False]*len(comparison_results))]
        
        y_pos = range(len(comparison_results))
        
        for i, (idx, row) in enumerate(comparison_results.iterrows()):
            # Draw confidence interval line
            ax.plot([row['ci_lower'], row['ci_upper']], [i, i],
                   color=colors_ci[i], linewidth=4, alpha=0.7, solid_capstyle='round')
            # Draw mean point
            ax.scatter(row['mean_difference'], i, color=colors_ci[i], s=120,
                      edgecolor='black', linewidths=2, zorder=3)
            # Add p-value
            if row.get('any_significant', False):
                ax.text(row['ci_upper'] + 0.05, i, format_pvalue(row['p_value_ttest']),
                       va='center', fontsize=8, fontweight='bold')
        
        ax.axvline(x=0, color='black', linestyle='-', alpha=0.5, linewidth=1.5)
        ax.axvline(x=0.5, color='orange', linestyle=':', alpha=0.5, linewidth=2, label='|diff| = 0.5')
        ax.axvline(x=-0.5, color='orange', linestyle=':', alpha=0.5, linewidth=2)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(comparison_results['comparison'], fontsize=9)
        ax.set_xlabel('Mean LogFC Difference (95% CI)', fontsize=12, fontweight='bold')
        ax.set_title('Mean Differences Between Enhancer Classes', fontsize=14, fontweight='bold', pad=20)
        
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=SIG_COLOR, edgecolor='black', label='Significant'),
            Patch(facecolor=NONSIG_COLOR, edgecolor='black', label='Not significant')
        ]
        ax.legend(handles=legend_elements, fontsize=10, framealpha=0.9, loc='best')
        ax.grid(alpha=0.3, axis='x', linestyle=':', linewidth=0.5)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_plot7_mean_differences.pdf", dpi=300, bbox_inches='tight')
        plt.close()
    
    # PLOT 8: Summary statistics
    print("  Creating plot 8: Summary statistics...")
    fig, ax = plt.subplots(figsize=(12, 6))
    
    summary_data = {
        'Total\nInteractions': len(enh_interactions),
        'wMel-enriched\n(logFC < 0)': len(wmel_enriched),
        'Uninf-enriched\n(logFC > 0)': len(uninf_enriched),
        'E-E': sum(enh_interactions['interaction_type'] == 'E-E'),
        'E-TSS': sum(enh_interactions['interaction_type'] == 'E-TSS'),
        'Housekeeping': sum(enh_interactions['contact_class'] == 'housekeeping'),
        'Developmental': sum(enh_interactions['contact_class'] == 'developmental'),
        'Cross-class': sum(enh_interactions['contact_class'] == 'cross_class')
    }
    
    if not comparison_results.empty:
        summary_data['Significant\nComparisons'] = sum(comparison_results.get('any_significant', 
                                                                              [False]*len(comparison_results)))
    
    colors_summary = []
    for key in summary_data.keys():
        if 'wMel' in key:
            colors_summary.append(WMEL_COLOR)
        elif 'Uninf' in key:
            colors_summary.append(UNINF_COLOR)
        elif 'Housekeeping' in key:
            colors_summary.append(HK_COLOR)
        elif 'Developmental' in key:
            colors_summary.append(DEV_COLOR)
        elif 'Cross-class' in key:
            colors_summary.append(CROSS_COLOR)
        elif 'Significant' in key:
            colors_summary.append(SIG_COLOR)
        else:
            colors_summary.append('#7fa8d1')
    
    bars = ax.bar(range(len(summary_data)), list(summary_data.values()),
            color=colors_summary, alpha=0.8, edgecolor='black', linewidth=1.5)
    
    ax.set_xticks(range(len(summary_data)))
    ax.set_xticklabels(list(summary_data.keys()), rotation=45, ha='right', fontsize=10)
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_title('Summary Statistics: Enhancer Interactions', fontsize=14, fontweight='bold', pad=20)
    ax.grid(alpha=0.3, axis='y', linestyle=':', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add values on bars
    for i, (key, val) in enumerate(summary_data.items()):
        ax.text(i, val, f'{val:,}', ha='center', va='bottom',
               fontweight='bold', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_plot8_summary_statistics.pdf", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"\nAll plots saved with prefix: {output_prefix}_plot*.pdf")


def create_report(enh_interactions, comparison_results, output_prefix):
    """Create detailed text report."""
    
    with open(f"{output_prefix}_analysis_report.txt", 'w') as f:
        f.write("="*80 + "\n")
        f.write("ENHANCER CLASS COMPARISON REPORT\n")
        f.write("wMel vs Uninf (DOX) Comparison\n")
        f.write("="*80 + "\n\n")
        
        # Enhancer interactions summary
        f.write("ENHANCER INTERACTIONS SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write(f"Total enhancer interactions: {len(enh_interactions)}\n")
        
        wmel_enriched = len(enh_interactions[enh_interactions['logFC'] < 0])
        uninf_enriched = len(enh_interactions[enh_interactions['logFC'] > 0])
        
        f.write(f"wMel-enriched (logFC < 0): {wmel_enriched} ({100*wmel_enriched/len(enh_interactions):.1f}%)\n")
        f.write(f"Uninf-enriched (logFC > 0): {uninf_enriched} ({100*uninf_enriched/len(enh_interactions):.1f}%)\n\n")
        
        f.write(f"E-E interactions: {sum(enh_interactions['interaction_type'] == 'E-E')}\n")
        f.write(f"E-TSS interactions: {sum(enh_interactions['interaction_type'] == 'E-TSS')}\n\n")
        
        f.write("By contact class:\n")
        for cls in sorted(enh_interactions['contact_class'].unique()):
            n = sum(enh_interactions['contact_class'] == cls)
            pct = 100 * n / len(enh_interactions)
            f.write(f"  {cls}: {n} ({pct:.1f}%)\n")
        
        obs_logfc = pd.to_numeric(enh_interactions['logFC'], errors='coerce').dropna()
        f.write(f"\nOverall logFC statistics:\n")
        f.write(f"  Mean: {obs_logfc.mean():.3f}\n")
        f.write(f"  Median: {obs_logfc.median():.3f}\n")
        f.write(f"  Std: {obs_logfc.std():.3f}\n")
        f.write(f"  Range: [{obs_logfc.min():.3f}, {obs_logfc.max():.3f}]\n\n")
        
        # Class-specific statistics
        f.write("BY ENHANCER CLASS:\n")
        f.write("-"*80 + "\n")
        for int_type in ['E-E', 'E-TSS']:
            type_data = enh_interactions[enh_interactions['interaction_type'] == int_type]
            f.write(f"\n{int_type} interactions:\n")
            
            for cls in sorted(type_data['contact_class'].unique()):
                class_data = pd.to_numeric(type_data[type_data['contact_class'] == cls]['logFC'],
                                          errors='coerce').dropna()
                f.write(f"  {cls}:\n")
                f.write(f"    n = {len(class_data)}\n")
                f.write(f"    Mean logFC = {class_data.mean():.3f}\n")
                f.write(f"    Median logFC = {class_data.median():.3f}\n")
                f.write(f"    Std = {class_data.std():.3f}\n\n")
        
        # Comparison results
        if not comparison_results.empty:
            f.write("\n" + "="*80 + "\n")
            f.write("CLASS COMPARISON RESULTS\n")
            f.write("="*80 + "\n")
            f.write(f"Total comparisons: {len(comparison_results)}\n")
            
            if 'any_significant' in comparison_results.columns:
                n_sig = sum(comparison_results['any_significant'])
                f.write(f"Significant comparisons (FDR < 0.05): {n_sig}/{len(comparison_results)}\n\n")
                
                if n_sig > 0:
                    f.write("SIGNIFICANT COMPARISONS:\n")
                    f.write("-"*80 + "\n")
                    
                    sig_results = comparison_results[comparison_results['any_significant']]
                    for idx, row in sig_results.iterrows():
                        f.write(f"\n{row['interaction_type']}: {row['class1']} vs {row['class2']}\n")
                        f.write("-"*40 + "\n")
                        f.write(f"  {row['class1']}: n={row['n1']}, mean={row['mean1']:.3f}\n")
                        f.write(f"  {row['class2']}: n={row['n2']}, mean={row['mean2']:.3f}\n")
                        f.write(f"  Mean difference: {row['mean_difference']:.4f}\n")
                        f.write(f"  Cohen's d: {row['cohens_d']:.4f} ({row['effect_size_category']} effect)\n")
                        f.write(f"  P-value (t-test): {row['p_value_ttest']:.4e}\n")
                        f.write(f"  FDR: {row['p_value_ttest_fdr']:.4e}\n")
                        f.write(f"  95% CI: [{row['ci_lower']:.4f}, {row['ci_upper']:.4f}]\n")
                        f.write(f"  CI excludes zero: {row['ci_excludes_zero']}\n")
                        f.write(f"  Biologically significant (|diff| > 0.5): {row['biologically_significant']}\n")
                        f.write(f"  Interpretation: {row['interpretation']}\n")
                
                # Non-significant results
                nonsig_results = comparison_results[~comparison_results['any_significant']]
                if not nonsig_results.empty:
                    f.write("\n\nNON-SIGNIFICANT COMPARISONS:\n")
                    f.write("-"*80 + "\n")
                    
                    for idx, row in nonsig_results.iterrows():
                        f.write(f"\n{row['interaction_type']}: {row['class1']} vs {row['class2']}\n")
                        f.write(f"  Mean difference: {row['mean_difference']:.4f}\n")
                        f.write(f"  Cohen's d: {row['cohens_d']:.4f} ({row['effect_size_category']} effect)\n")
                        f.write(f"  P-value: {row['p_value_ttest']:.4e}, FDR: {row['p_value_ttest_fdr']:.4e}\n")
            
            # Summary interpretation
            f.write("\n\n" + "="*80 + "\n")
            f.write("OVERALL INTERPRETATION\n")
            f.write("="*80 + "\n")
            
            if 'any_significant' in comparison_results.columns:
                n_sig = sum(comparison_results['any_significant'])
                n_bio_sig = sum(comparison_results.get('biologically_significant', [False]*len(comparison_results)))
                
                if n_sig == 0:
                    f.write("\nNo statistically significant differences were found between enhancer classes.\n")
                    f.write("This suggests that Wolbachia infection affects all enhancer types similarly.\n")
                elif n_bio_sig == 0:
                    f.write(f"\n{n_sig} comparison(s) were statistically significant but none showed\n")
                    f.write("biologically meaningful differences (|mean difference| > 0.5 logFC).\n")
                    f.write("This suggests minimal biological differences between enhancer classes\n")
                    f.write("in their response to Wolbachia infection.\n")
                else:
                    f.write(f"\n{n_sig} comparison(s) were statistically significant, with {n_bio_sig}\n")
                    f.write("showing biologically meaningful differences (|mean difference| > 0.5 logFC).\n")
                    f.write("\nEnhancer classes show differential responses to Wolbachia infection.\n")
                
                # Effect size summary
                large_effects = sum(comparison_results['effect_size_category'] == 'large')
                medium_effects = sum(comparison_results['effect_size_category'] == 'medium')
                small_effects = sum(comparison_results['effect_size_category'] == 'small')
                
                f.write(f"\nEffect size distribution:\n")
                f.write(f"  Large effects (|d| > 0.8): {large_effects}\n")
                f.write(f"  Medium effects (0.5 < |d| ≤ 0.8): {medium_effects}\n")
                f.write(f"  Small effects (0.2 < |d| ≤ 0.5): {small_effects}\n")
                f.write(f"  Negligible effects (|d| ≤ 0.2): {len(comparison_results) - large_effects - medium_effects - small_effects}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF REPORT\n")
        f.write("="*80 + "\n")
    
    print(f"Report saved to {output_prefix}_analysis_report.txt")


def main():
    parser = argparse.ArgumentParser(
        description='Enhancer class comparison for wMel vs uninf (DOX)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script compares how different enhancer classes respond to Wolbachia infection
by analyzing differential chromatin interactions between wMel-infected and uninfected cells.

It tests whether housekeeping vs developmental enhancers show significantly different
logFC distributions, using proper statistical comparisons between groups.
        """
    )
    
    parser.add_argument('--enhancers', required=True,
                       help='Enhancer BED file')
    parser.add_argument('--interactions', required=True,
                       help='Differential interactions CSV (wMel vs uninf)')
    parser.add_argument('--classification',
                       help='Optional enhancer classification file')
    parser.add_argument('--fdr_threshold', type=float, default=0.05,
                       help='FDR threshold (default: 0.05)')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    print("\n" + "="*80)
    print("ENHANCER CLASS COMPARISON: wMel vs Uninf (DOX)")
    print("="*80 + "\n")
    
    # Load data
    interactions = load_interactions(args.interactions)
    enhancers = classify_enhancers(args.enhancers, args.classification)
    
    if enhancers.empty:
        print("ERROR: No valid enhancers!")
        return
    
    # Find enhancer interactions
    enh_interactions = find_enhancer_interactions(interactions, enhancers)
    
    if enh_interactions.empty:
        print("WARNING: No enhancer interactions found!")
        return
    
    # Statistical comparison between classes
    comparison_results = compare_enhancer_classes(enh_interactions, args.fdr_threshold)
    
    # Create outputs
    create_visualizations(enh_interactions, comparison_results, args.output_prefix)
    create_report(enh_interactions, comparison_results, args.output_prefix)
    
    # Save results
    print("\nSaving results...")
    
    enh_interactions.to_csv(f"{args.output_prefix}_enhancer_interactions.tsv",
                           sep='\t', index=False)
    print(f"  Saved: {args.output_prefix}_enhancer_interactions.tsv")
    
    if not comparison_results.empty:
        comparison_results.to_csv(f"{args.output_prefix}_class_comparisons.tsv",
                                sep='\t', index=False)
        print(f"  Saved: {args.output_prefix}_class_comparisons.tsv")
        
        if 'any_significant' in comparison_results.columns:
            sig_results = comparison_results[comparison_results['any_significant']]
            if not sig_results.empty:
                sig_results.to_csv(f"{args.output_prefix}_significant_comparisons.tsv",
                                 sep='\t', index=False)
                print(f"  Saved: {args.output_prefix}_significant_comparisons.tsv")
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nResults written to: {args.output_prefix}*")
    print(f"Individual plots saved as: {args.output_prefix}_plot1.pdf through plot8.pdf")
    
    # Print final summary
    if not comparison_results.empty and 'any_significant' in comparison_results.columns:
        n_sig = sum(comparison_results['any_significant'])
        n_bio = sum(comparison_results.get('biologically_significant', [False]*len(comparison_results)))
        
        print("\n" + "="*80)
        print("SUMMARY")
        print("="*80)
        print(f"Comparisons performed: {len(comparison_results)}")
        print(f"Statistically significant: {n_sig}")
        print(f"Biologically significant (|diff| > 0.5): {n_bio}")
        
        if n_sig == 0:
            print("\n>>> No significant differences between enhancer classes.")
            print(">>> Wolbachia appears to affect all enhancer types similarly.")
        elif n_bio == 0:
            print("\n>>> Statistically significant but not biologically meaningful.")
            print(">>> Differences are very small - likely not biologically relevant.")
        else:
            print("\n>>> Enhancer classes show differential responses to infection!")
            print(">>> Check the detailed report for biological interpretation.")


if __name__ == '__main__':
    main()