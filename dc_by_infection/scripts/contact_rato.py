import cooler
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths
base_path = "/private/groups/russelllab/jodie/wolbachia_induced_DE/micro-c/processed_files"
save_path = "/private/groups/russelllab/jodie/wolbachia_induced_DE/Jacobs_et_al_2026_Wolbachia-host_genome_topology/dc_by_infection/results/contact_ratio"
resolution = 1000

files = {
    'JW18-DOX-1': f"{base_path}/JW18-DOX-1.matrix_1kb.mcool::/resolutions/{resolution}",
    'JW18-DOX-2': f"{base_path}/JW18-DOX-2.matrix_1kb.mcool::/resolutions/{resolution}",
    'JW18-wMel-1': f"{base_path}/JW18-wMel-1.matrix_1kb.mcool::/resolutions/{resolution}",
    'JW18-wMel-2': f"{base_path}/JW18-wMel-2.matrix_1kb.mcool::/resolutions/{resolution}",
}

def calculate_x_vs_individual_autosomes(clr, x_chrom='X', autosomes=['2L', '2R', '3L', '3R']):
    """
    Calculate ratio of X chromosome vs each individual autosome intra-chromosomal contact density.
    Uses RAW counts only.
    Returns a dictionary with one ratio per autosome.
    """
    
    # Get RAW X chromosome matrix
    x_matrix = clr.matrix(balance=False, sparse=True).fetch(x_chrom)
    x_contacts = float(x_matrix.sum())
    x_length = clr.chromsizes[x_chrom]
    x_density = (x_contacts / x_length) * 1e6
    
    print(f"    X chromosome: {x_contacts:,.0f} total contacts, density: {x_density:.2f} contacts/Mb")
    
    # Calculate ratio for each autosome
    ratios = {}
    for chrom in autosomes:
        auto_matrix = clr.matrix(balance=False, sparse=True).fetch(chrom)
        auto_contacts = float(auto_matrix.sum())
        auto_length = clr.chromsizes[chrom]
        auto_density = (auto_contacts / auto_length) * 1e6
        
        ratio = x_density / auto_density
        ratios[chrom] = {
            'ratio': ratio,
            'x_density': x_density,
            'auto_density': auto_density,
            'x_contacts': x_contacts,
            'auto_contacts': auto_contacts
        }
        
        print(f"    {chrom}: {auto_contacts:,.0f} contacts, density: {auto_density:.2f} contacts/Mb, X/{chrom} ratio: {ratio:.4f}")
    
    return ratios

# Calculate ratios for all samples
results = []

for sample_name, file_path in files.items():
    print(f"\nProcessing {sample_name}...")
    try:
        clr = cooler.Cooler(file_path)
        
        # Check available chromosomes
        print(f"  Available chromosomes: {clr.chromnames}")
        
        ratios_dict = calculate_x_vs_individual_autosomes(clr)
        
        condition = 'Uninfected' if 'DOX' in sample_name else 'Infected'
        replicate = sample_name.split('-')[-1]
        
        # Add one row per autosome comparison
        for chrom, data in ratios_dict.items():
            results.append({
                'Sample': sample_name,
                'Condition': condition,
                'Replicate': replicate,
                'Autosome': chrom,
                'Comparison': f'X/{chrom}',
                'X/Auto_Ratio': data['ratio'],
                'X_Density': data['x_density'],
                'Auto_Density': data['auto_density'],
                'X_Total_Contacts': data['x_contacts'],
                'Auto_Total_Contacts': data['auto_contacts']
            })
    except Exception as e:
        print(f"  ERROR processing {sample_name}: {e}")
        import traceback
        traceback.print_exc()

# Create DataFrame
df = pd.DataFrame(results)
print(f"\n{'='*60}")
print("Using RAW contact counts (not balanced)")
print("Individual X/Autosome ratios for each chromosome arm")
print(f"{'='*60}")
print("\n=== Results Summary ===")
print(df.to_string(index=False))

# Check for NaN values
if df['X/Auto_Ratio'].isna().any():
    print("\nERROR: Still getting NaN values!")
    exit(1)

# Statistical testing
dox_ratios = df[df['Condition'] == 'Uninfected']['X/Auto_Ratio'].values
wmel_ratios = df[df['Condition'] == 'Infected']['X/Auto_Ratio'].values

# Two-sample t-test
t_stat, t_pval = stats.ttest_ind(dox_ratios, wmel_ratios)

# Mann-Whitney U test (non-parametric alternative)
u_stat, u_pval = stats.mannwhitneyu(dox_ratios, wmel_ratios, alternative='two-sided')

print("\n=== Statistical Tests ===")
print(f"Uninfected (JW18-DOX) ratios: n={len(dox_ratios)}")
print(f"  Mean ± SD: {dox_ratios.mean():.4f} ± {dox_ratios.std():.4f}")
print(f"Infected (JW18-wMel) ratios: n={len(wmel_ratios)}")
print(f"  Mean ± SD: {wmel_ratios.mean():.4f} ± {wmel_ratios.std():.4f}")
print(f"\nFold change: {wmel_ratios.mean() / dox_ratios.mean():.4f}")
print(f"Percent change: {((wmel_ratios.mean() - dox_ratios.mean()) / dox_ratios.mean() * 100):.2f}%")
print(f"\nT-test: t = {t_stat:.4f}, P = {t_pval:.4e}")
print(f"Mann-Whitney U test: U = {u_stat:.4f}, P = {u_pval:.4e}")

# Visualization
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Bar plot with individual points
ax1 = axes[0]
sns.barplot(data=df, x='Condition', y='X/Auto_Ratio', hue='Condition', ax=ax1, 
            palette=['#d62728', '#1f77b4'], alpha=0.7, errorbar='sd', legend=False)
sns.stripplot(data=df, x='Condition', y='X/Auto_Ratio', ax=ax1, 
              color='black', size=8, alpha=0.8)
ax1.set_ylabel('X/Autosome Contact Ratio')
ax1.set_title('X Chromosome Contact Density\nRelative to Individual Autosomes')
ax1.set_ylim(bottom=0)

# Add significance annotation
y_max = df['X/Auto_Ratio'].max()
if u_pval < 0.001:
    sig_text = '***'
elif u_pval < 0.01:
    sig_text = '**'
elif u_pval < 0.05:
    sig_text = '*'
else:
    sig_text = 'ns'

ax1.plot([0, 1], [y_max * 1.1, y_max * 1.1], 'k-', linewidth=1)
ax1.text(0.5, y_max * 1.15, sig_text, ha='center', fontsize=14)

# Box plot alternative
ax2 = axes[1]
sns.boxplot(data=df, x='Condition', y='X/Auto_Ratio', hue='Condition', ax=ax2,
            palette=['#d62728', '#1f77b4'], legend=False)
sns.stripplot(data=df, x='Condition', y='X/Auto_Ratio', ax=ax2,
              color='black', size=8, alpha=0.8)
ax2.set_ylabel('X/Autosome Contact Ratio')
ax2.set_title('Distribution of All X/Auto Comparisons')

plt.tight_layout()
plt.savefig(f"{save_path}/X_autosome_contact_ratio_analysis.pdf", dpi=300)
print(f"\nFigure saved to: {save_path}/X_autosome_contact_ratio_analysis.pdf")

# Save results to CSV
df.to_csv(f"{save_path}/X_autosome_contact_ratios.csv", index=False)
print(f"Data saved to: {save_path}/X_autosome_contact_ratios.csv")