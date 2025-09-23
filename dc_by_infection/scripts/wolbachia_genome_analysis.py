#!/usr/bin/env python3
"""
Analyze Wolbachia genome contacts across different strains.
Compare chromatin organization between wMel, wRi, and wWil strains.
"""

import pandas as pd
import numpy as np
import cooler
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path

def load_wolbachia_genomes():
    """Load Wolbachia genome information for each strain"""
    # These are approximate genome sizes - adjust based on your actual references
    wolbachia_genomes = {
        'wMel': {'size': 1267782, 'contigs': ['wMel']},
        'wRi': {'size': 1445873, 'contigs': ['wRi']},
        'wWil': {'size': 1482455, 'contigs': ['wWil']}
    }
    return wolbachia_genomes

def extract_wolbachia_contacts(mcool_file, strain, resolution=1000):
    """Extract Wolbachia-specific contacts from mcool file"""
    print(f"Extracting Wolbachia contacts for {strain} at {resolution}bp...")
    
    try:
        clr = cooler.Cooler(f"{mcool_file}::resolutions/{resolution}")
        
        # Check if Wolbachia contigs are in the cooler
        all_chroms = clr.chromnames
        wolbachia_chroms = [c for c in all_chroms if strain.lower() in c.lower()]
        
        if not wolbachia_chroms:
            print(f"Warning: No Wolbachia chromosomes found for {strain}")
            return None
        
        # Extract Wolbachia-Wolbachia contacts
        wol_contacts = []
        for chrom in wolbachia_chroms:
            try:
                matrix = clr.matrix(balance=True).fetch(chrom)
                if matrix.size > 0:
                    wol_contacts.append({
                        'strain': strain,
                        'chrom': chrom,
                        'matrix': matrix,
                        'size': matrix.shape[0]
                    })
            except:
                continue
        
        # Extract Wolbachia-host contacts
        host_chroms = ['2L', '2R', '3L', '3R', '4', 'X']
        wol_host_contacts = []
        
        for wol_chrom in wolbachia_chroms:
            for host_chrom in host_chroms:
                try:
                    matrix = clr.matrix(balance=True).fetch(wol_chrom, host_chrom)
                    if matrix.size > 0:
                        wol_host_contacts.append({
                            'strain': strain,
                            'wol_chrom': wol_chrom,
                            'host_chrom': host_chrom,
                            'matrix': matrix
                        })
                except:
                    continue
        
        return {
            'intra_wolbachia': wol_contacts,
            'wolbachia_host': wol_host_contacts
        }
    
    except Exception as e:
        print(f"Error processing {strain}: {e}")
        return None

def analyze_wolbachia_compaction(contacts_dict):
    """Analyze chromatin compaction of Wolbachia genomes"""
    results = {}
    
    for strain, contacts in contacts_dict.items():
        if not contacts or not contacts['intra_wolbachia']:
            continue
            
        # Analyze intra-Wolbachia contacts
        for wol_data in contacts['intra_wolbachia']:
            matrix = wol_data['matrix']
            
            # Calculate contact decay
            distances = []
            mean_contacts = []
            
            for d in range(1, min(matrix.shape[0], 100)):
                diagonal = np.diagonal(matrix, d)
                valid_diagonal = diagonal[~np.isnan(diagonal)]
                if len(valid_diagonal) > 0:
                    distances.append(d)
                    mean_contacts.append(np.mean(valid_diagonal))
            
            # Fit power law decay
            if len(distances) > 10:
                log_dist = np.log10(distances)
                log_contacts = np.log10(mean_contacts)
                
                # Linear fit in log space
                slope, intercept = np.polyfit(log_dist, log_contacts, 1)
                
                results[strain] = {
                    'decay_exponent': -slope,
                    'compaction_score': intercept,
                    'matrix_size': matrix.shape[0],
                    'total_contacts': np.nansum(matrix),
                    'distances': distances,
                    'mean_contacts': mean_contacts
                }
    
    return results

def analyze_wolbachia_host_interactions(contacts_dict):
    """Analyze interactions between Wolbachia and host genomes"""
    results = {}
    
    for strain, contacts in contacts_dict.items():
        if not contacts or not contacts['wolbachia_host']:
            continue
            
        strain_results = {
            'total_contacts': 0,
            'by_chromosome': {},
            'contact_distribution': []
        }
        
        for interaction in contacts['wolbachia_host']:
            matrix = interaction['matrix']
            host_chrom = interaction['host_chrom']
            
            # Sum contacts for this chromosome pair
            total = np.nansum(matrix)
            strain_results['total_contacts'] += total
            strain_results['by_chromosome'][host_chrom] = total
            
            # Store flattened contact values
            flat_contacts = matrix.flatten()
            valid_contacts = flat_contacts[~np.isnan(flat_contacts) & (flat_contacts > 0)]
            strain_results['contact_distribution'].extend(valid_contacts)
        
        results[strain] = strain_results
    
    return results

def compare_strain_differences(compaction_results, interaction_results):
    """Compare chromatin organization between Wolbachia strains"""
    comparisons = []
    
    # Compare compaction between strains
    strains = list(compaction_results.keys())
    for i in range(len(strains)):
        for j in range(i+1, len(strains)):
            strain1, strain2 = strains[i], strains[j]
            
            comparison = {
                'strain1': strain1,
                'strain2': strain2,
                'decay_diff': (compaction_results[strain1]['decay_exponent'] - 
                              compaction_results[strain2]['decay_exponent']),
                'compaction_diff': (compaction_results[strain1]['compaction_score'] - 
                                   compaction_results[strain2]['compaction_score'])
            }
            
            # Compare host interactions
            if strain1 in interaction_results and strain2 in interaction_results:
                comparison['host_contact_ratio'] = (
                    interaction_results[strain1]['total_contacts'] / 
                    interaction_results[strain2]['total_contacts']
                )
                
                # Statistical test on contact distributions
                if (len(interaction_results[strain1]['contact_distribution']) > 0 and 
                    len(interaction_results[strain2]['contact_distribution']) > 0):
                    stat, pval = stats.mannwhitneyu(
                        interaction_results[strain1]['contact_distribution'],
                        interaction_results[strain2]['contact_distribution']
                    )
                    comparison['interaction_pvalue'] = pval
            
            comparisons.append(comparison)
    
    return comparisons

def create_wolbachia_plots(compaction_results, interaction_results, output_prefix):
    """Create visualizations for Wolbachia genome analysis"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: Contact decay curves
    ax = axes[0, 0]
    for strain, results in compaction_results.items():
        if 'distances' in results:
            ax.loglog(results['distances'], results['mean_contacts'], 
                     marker='o', label=f"{strain} (α={results['decay_exponent']:.2f})")
    
    ax.set_xlabel('Genomic Distance (bins)')
    ax.set_ylabel('Mean Contact Frequency')
    ax.set_title('Wolbachia Genome Contact Decay')
    ax.legend()
    
    # Plot 2: Compaction comparison
    ax = axes[0, 1]
    strains = list(compaction_results.keys())
    decay_exponents = [results['decay_exponent'] for results in compaction_results.values()]
    
    ax.bar(strains, decay_exponents)
    ax.set_ylabel('Decay Exponent')
    ax.set_title('Genome Compaction Comparison')
    
    # Plot 3: Host interaction levels
    ax = axes[1, 0]
    if interaction_results:
        strains = list(interaction_results.keys())
        total_contacts = [results['total_contacts'] for results in interaction_results.values()]
        
        ax.bar(strains, total_contacts)
        ax.set_ylabel('Total Wolbachia-Host Contacts')
        ax.set_title('Inter-genomic Interactions')
    
    # Plot 4: Host chromosome preference
    ax = axes[1, 1]
    if interaction_results:
        host_chroms = ['2L', '2R', '3L', '3R', '4', 'X']
        
        for strain, results in interaction_results.items():
            if 'by_chromosome' in results:
                values = [results['by_chromosome'].get(chrom, 0) for chrom in host_chroms]
                ax.plot(host_chroms, values, marker='o', label=strain)
        
        ax.set_xlabel('Host Chromosome')
        ax.set_ylabel('Contact Frequency')
        ax.set_title('Wolbachia-Host Contact Distribution')
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_wolbachia_genome_analysis.pdf", dpi=300)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze Wolbachia genome contacts')
    parser.add_argument('--mcool_files', nargs='+', required=True, 
                       help='Micro-C mcool files for infected samples')
    parser.add_argument('--strains', nargs='+', required=True,
                       help='Wolbachia strain names corresponding to mcool files')
    parser.add_argument('--resolution', type=int, default=1000,
                       help='Resolution for analysis')
    parser.add_argument('--output_prefix', required=True,
                       help='Output file prefix')
    
    args = parser.parse_args()
    
    # Validate inputs
    if len(args.mcool_files) != len(args.strains):
        print("Error: Number of mcool files must match number of strains")
        return
    
    # Extract contacts for each strain
    all_contacts = {}
    for mcool_file, strain in zip(args.mcool_files, args.strains):
        contacts = extract_wolbachia_contacts(mcool_file, strain, args.resolution)
        if contacts:
            all_contacts[strain] = contacts
    
    if not all_contacts:
        print("Error: No Wolbachia contacts found in any sample")
        return
    
    # Analyze genome compaction
    print("\nAnalyzing Wolbachia genome compaction...")
    compaction_results = analyze_wolbachia_compaction(all_contacts)
    
    # Analyze host interactions
    print("Analyzing Wolbachia-host interactions...")
    interaction_results = analyze_wolbachia_host_interactions(all_contacts)
    
    # Compare strains
    print("Comparing between strains...")
    comparisons = compare_strain_differences(compaction_results, interaction_results)
    
    # Create visualizations
    create_wolbachia_plots(compaction_results, interaction_results, args.output_prefix)
    
    # Save results
    # Compaction results
    compaction_df = pd.DataFrame([
        {
            'strain': strain,
            'decay_exponent': results['decay_exponent'],
            'compaction_score': results['compaction_score'],
            'total_contacts': results['total_contacts']
        }
        for strain, results in compaction_results.items()
    ])
    compaction_df.to_csv(f"{args.output_prefix}_compaction_results.tsv", sep='\t', index=False)
    
    # Interaction results
    interaction_df = pd.DataFrame([
        {
            'strain': strain,
            'total_host_contacts': results['total_contacts'],
            **{f'contacts_{chrom}': results['by_chromosome'].get(chrom, 0) 
               for chrom in ['2L', '2R', '3L', '3R', '4', 'X']}
        }
        for strain, results in interaction_results.items()
    ])
    interaction_df.to_csv(f"{args.output_prefix}_host_interactions.tsv", sep='\t', index=False)
    
    # Strain comparisons
    comparison_df = pd.DataFrame(comparisons)
    comparison_df.to_csv(f"{args.output_prefix}_strain_comparisons.tsv", sep='\t', index=False)
    
    print(f"\nAnalysis complete. Results saved to {args.output_prefix}_*")

if __name__ == '__main__':
    main()