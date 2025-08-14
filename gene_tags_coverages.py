#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def parse_site_metadata(metadata_file):
    """Parse the site metadata file"""
    site_data = {}
    with open(metadata_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sample = parts[0]
                site = parts[1]
                site_data[sample] = site
    return site_data

def parse_coverage_data(csv_file):
    """Parse the bcftools output CSV"""
    variant_data = []
    sample_data = []
    
    with open(csv_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 7:
                continue
                
            chrom, pos, ref, alt, maf, ns, f_missing = parts[:7]
            pos = int(pos)
            
            # Add variant data only once per variant
            variant_data.append({
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'maf': float(maf) if maf != '.' else 0,
                'ns': int(ns) if ns != '.' else 0,
                'f_missing': float(f_missing) if f_missing != '.' else 0
            })
            
            # Replace spaces with colons in the remaining part, then parse sample:DP pairs
            remaining = ','.join(parts[7:])  # Rejoin in case there were commas in sample names
            remaining = remaining.replace(' ', ':')  # Replace spaces with colons
            samples_dp = remaining.split(',')  # Split by commas
            
            for sample_dp in samples_dp:
                if ':' in sample_dp:
                    sample, dp = sample_dp.split(':')
                    dp = int(dp) if dp != '.' else 0
                    
                    sample_data.append({
                        'sample': sample,
                        'pos': pos,
                        'dp': dp
                    })
    
    variant_df = pd.DataFrame(variant_data)
    sample_df = pd.DataFrame(sample_data)
    
    print(f"Variant data shape: {variant_df.shape}")
    print(f"Sample data shape: {sample_df.shape}")
    print(f"Sample columns: {sample_df.columns.tolist()}")
    print(f"First few rows of sample data:")
    print(sample_df.head())
    
    return variant_df, sample_df

def create_coverage_heatmap(sample_df, output_file, site_data=None):
    """Create coverage heatmap with optional site grouping"""
    print(f"Total rows: {len(sample_df)}")
    print(f"Unique samples: {sample_df['sample'].nunique()}")
    print(f"Unique positions: {sample_df['pos'].nunique()}")
    print(f"Expected rows: {sample_df['sample'].nunique() * sample_df['pos'].nunique()}")
    
    # Check for actual duplicates in sample:position combinations
    duplicates = sample_df.groupby(['sample', 'pos']).size()
    if duplicates.max() > 1:
        print(f"Found duplicate sample:position combinations, max count: {duplicates.max()}")
        print("First few duplicates:")
        print(duplicates[duplicates > 1].head())
        
        # Show the actual duplicate rows to understand whats happening
        print("\nExamining duplicate rows:")
        duplicate_pairs = duplicates[duplicates > 1].index.tolist()[:5]  # First 5 duplicate pairs
        for sample, pos in duplicate_pairs:
            dup_rows = sample_df[(sample_df['sample'] == sample) & (sample_df['pos'] == pos)]
            print(f"\nSample: {sample}, Position: {pos}")
            print(dup_rows)
        
        # Remove duplicates by taking the first occurrence
        print("\nRemoving duplicates...")
        sample_df = sample_df.drop_duplicates(subset=['sample', 'pos'], keep='first')
        print(f"After removing duplicates: {len(sample_df)} rows")
    
    # Pivot data for heatmap
    heatmap_data = sample_df.pivot(index='sample', columns='pos', values='dp')
    
    # Reorder samples by site if metadata provided
    if site_data:
        # Add site information to samples
        sample_sites = []
        for sample in heatmap_data.index:
            site = site_data.get(sample, 'Unknown')
            sample_sites.append((sample, site))
        
        # Sort by site, then by sample name within each site
        sample_sites.sort(key=lambda x: (x[1], x[0]))
        ordered_samples = [s[0] for s in sample_sites]
        
        # Reorder the heatmap data
        heatmap_data = heatmap_data.reindex(ordered_samples)
        
        # Create site annotation data
        site_annotations = [s[1] for s in sample_sites]
        
        # Create figure with space for site annotations
        fig, (ax_annot, ax_heatmap) = plt.subplots(1, 2, figsize=(18, 10), 
                                                   gridspec_kw={'width_ratios': [0.05, 0.95]})
        
        # Create site annotation bar with simple colors
        unique_sites = list(dict.fromkeys(site_annotations))  # Preserve order
        simple_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
        
        # Create annotation matrix - each row represents one sample
        annot_matrix = np.zeros((len(site_annotations), 1))
        for i, site in enumerate(site_annotations):
            site_idx = unique_sites.index(site)
            annot_matrix[i, 0] = site_idx
        
        # Plot site annotations with simple colors
        im = ax_annot.imshow(annot_matrix, aspect='auto', cmap='tab10', vmin=0, vmax=len(unique_sites)-1)
        ax_annot.set_xticks([])
        ax_annot.set_yticks(range(len(site_annotations)))
        ax_annot.set_yticklabels([])
        ax_annot.set_title('Site', fontsize=10)
        
        # Add site labels on the right side of annotation bar
        for site in unique_sites:
            site_indices = [j for j, s in enumerate(site_annotations) if s == site]
            if site_indices:
                mid_point = (min(site_indices) + max(site_indices)) / 2
                ax_annot.text(-0.5, mid_point, site, ha='right', va='center', 
                             fontsize=8, rotation=0)
        
        # Plot heatmap
        sns.heatmap(heatmap_data, cmap='YlOrRd', cbar_kws={'label': 'Depth'}, 
                    vmax=100, yticklabels=False, ax=ax_heatmap)
        ax_heatmap.set_title('Gene Coverage Heatmap (Grouped by Site)')
        ax_heatmap.set_xlabel('Position')
        ax_heatmap.set_ylabel('Sample')
        ax_heatmap.set_xticklabels(ax_heatmap.get_xticklabels(), rotation=45)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Print site grouping info
        print(f"\nSite groupings:")
        for site in unique_sites:
            site_samples = [s for s in site_annotations if s == site]
            print(f"  {site}: {len(site_samples)} samples")
    
    else:
        # Original heatmap without site grouping
        plt.figure(figsize=(15, 10))
        sns.heatmap(heatmap_data, cmap='YlOrRd', cbar_kws={'label': 'Depth'}, 
                    vmax=100, yticklabels=False)
        plt.title('Gene Coverage Heatmap')
        plt.xlabel('Position')
        plt.ylabel('Sample')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

def create_variant_metrics_plots(variant_df, sample_df, output_file):
    """Create line plots for variant-level metrics"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    
    # NS across positions
    axes[0,0].plot(variant_df['pos'], variant_df['ns'])
    axes[0,0].set_title('Number of Samples with Data')
    axes[0,0].set_xlabel('Position')
    axes[0,0].set_ylabel('NS')
    
    # F_MISSING across positions
    axes[0,1].plot(variant_df['pos'], variant_df['f_missing'])
    axes[0,1].set_title('Fraction Missing')
    axes[0,1].set_xlabel('Position')
    axes[0,1].set_ylabel('F_MISSING')
    
    # MAF across positions
    axes[1,0].plot(variant_df['pos'], variant_df['maf'])
    axes[1,0].set_title('Minor Allele Frequency')
    axes[1,0].set_xlabel('Position')
    axes[1,0].set_ylabel('MAF')
    
    # Coverage statistics across positions
    coverage_stats = sample_df.groupby('pos')['dp'].agg(['mean', 'median', 'std']).reset_index()
    axes[1,1].plot(coverage_stats['pos'], coverage_stats['mean'], label='Mean')
    axes[1,1].plot(coverage_stats['pos'], coverage_stats['median'], label='Median')
    axes[1,1].set_title('Coverage Statistics')
    axes[1,1].set_xlabel('Position')
    axes[1,1].set_ylabel('Depth')
    axes[1,1].legend()
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Gene coverage analysis')
    parser.add_argument('csv_file', help='Input CSV file from bcftools query')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    parser.add_argument('--metadata', help='Sample metadata file with sample names and sites')
    
    args = parser.parse_args()
    
    # Parse site metadata if provided
    site_data = None
    if args.metadata:
        print(f"Loading site metadata from {args.metadata}")
        site_data = parse_site_metadata(args.metadata)
        print(f"Loaded site data for {len(site_data)} samples")
    
    # Parse data
    variant_df, sample_df = parse_coverage_data(args.csv_file)
    
    # Create plots
    create_coverage_heatmap(sample_df, f'{args.output_dir}/coverage_heatmap.png', site_data)
    create_variant_metrics_plots(variant_df, sample_df, f'{args.output_dir}/variant_metrics.png')
    
    print("Plots generated")

if __name__ == '__main__':
    main()
