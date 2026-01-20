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

def calculate_intrachrom_ratio(clr, x_chrom='X', autosomes=['2L', '2R', '3L', '3R']):
    """
    Calculate ratio of X chromosome vs autosome intra-chromosomal contact density.
    Uses RAW counts only.
    """
    
    # Get RAW matrices only
    x_matrix = clr.matrix(balance=False, sparse=True).fetch(x_chrom)
    x_contacts = float(x_matrix.sum())
    
    print(f"    X chromosome: {x_contacts:,.0f} total contacts")
    
    # Get total intra-chromosomal contacts for autosomes
    auto_contacts = 0
    for chrom in autosomes:
        auto_matrix = clr.matrix(balance=False, sparse=True).fetch(chrom)
        chrom_contacts = float(auto_matrix.sum())
        print(f"    {chrom}: {chrom_contacts:,.0f} contacts")
        auto_contacts += chrom_contacts
    
    print(f"    Total autosome contacts: {auto_contacts:,.0f}")
    
    # Get chromosome lengths
    x_length = clr.chromsizes[x_chrom]
    auto_length = sum(clr.chromsizes[chrom] for chrom in autosomes)
    
    print(f"    X length: {x_length:,} bp")
    print(f"    Autosome total length: {auto_length:,} bp")
    
    # Calculate contact density (contacts per Mb)
    x_density = (x_contacts / x_length) * 1e6
    auto_density = (auto_contacts / auto_length) * 1e6
    
    print(f"    X density: {x_density:.2f} contacts/Mb")
    print(f"    Auto density: {auto_density:.2f} contacts/Mb")
    
    # Return ratio
    ratio = x_density / auto_density
    print(f"    Ratio: {ratio:.4f}")
    
    return ratio, x_density, auto_density, x_contacts, auto_contacts

# Calculate ratios for all samples
results = []

for sample_name, file_path in files.items():
    print(f"\nProcessing {sample_name}...")
    try:
        clr = cooler.Cooler(file_path)
        
        # Check available chromosomes
        print(f"  Available chromosomes: {clr.chromnames}")
        
        ratio, x_dens, auto_dens, x_contacts, auto_contacts = calculate_intrachrom_ratio(clr)
        
        condition = 'Uninfected' if 'DOX' in sample_name else 'Infected'
        replicate = sample_name.split('-')[-1]
        
        results.append({
            'Sample': sample_name,
            'Condition': condition,
            'Replicate': replicate,
            'X/Auto_Ratio': ratio,
            'X_Density': x_dens,
            'Auto_Density': auto_dens,
            'X_Total_Contacts': x_contacts,
            'Auto_Total_Contacts': auto_contacts
        })
    except Exception as e:
        print(f"  ERROR processing {sample_name}: {e}")
        import traceback
        traceback.print_exc()

# Create DataFrame
df = pd.DataFrame(results)
print(f"\n{'='*60}")
print("Using RAW contact counts (not balanced)")
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
print(f"Uninfected (JW18-DOX) ratios: {dox_ratios}")
print(f"Infected (JW18-wMel) ratios: {wmel_ratios}")
print(f"\nMean ± SD:")
print(f"  Uninfected: {dox_ratios.mean():.4f} ± {dox_ratios.std():.4f}")
print(f"  Infected: {wmel_ratios.mean():.4f} ± {wmel_ratios.std():.4f}")
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
ax1.set_title('X Chromosome Contact Density\nRelative to Autosomes')
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
ax2.set_title('Distribution of Replicates')

plt.tight_layout()
plt.savefig(f"{save_path}/X_autosome_contact_ratio_analysis.pdf", dpi=300)
print(f"\nFigure saved to: {save_path}/X_autosome_contact_ratio_analysis.pdf")

# Save results to CSV
df.to_csv(f"{save_path}/X_autosome_contact_ratios.csv", index=False)
print(f"Data saved to: {save_path}/X_autosome_contact_ratios.csv")