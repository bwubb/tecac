#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import os

def parse_site_metadata(metadata_file):
    """Parse the site metadata file"""
    site_data={}
    status_data={}
    with open(metadata_file,'r') as f:
        for line in f:
            parts=line.strip().split(',')
            if len(parts)>=3 and parts[0]!='SampleID':  # Skip header
                sample=parts[0]
                site=parts[1]
                status=parts[2]
                site_data[sample]=site
                status_data[sample]=int(status)
    return site_data,status_data

def parse_coverage_data(csv_file):
    """Parse the bcftools output CSV"""
    variant_data=[]
    
    print("Parsing variant data...")
    with open(csv_file,'r') as f:
        for line in f:
            parts=line.strip().split(',')
            if len(parts)<6:
                continue
                
            chrom,pos,ref,alt,maf,ns=parts[:6]
            pos=int(pos)
            
            variant_data.append({
                'chrom':chrom,
                'pos':pos,
                'ref':ref,
                'alt':alt,
                'maf':float(maf) if maf!='.' else 0.0,
                'ns':int(ns) if ns!='.' else 0
            })
    
    variant_df=pd.DataFrame(variant_data)
    print(f"Variant data shape: {variant_df.shape}")
    
    return variant_df

def parse_sample_data_for_chunk(csv_file,chunk_variants):
    """Parse sample data for a specific chunk of variants"""
    chunk_positions=set(chunk_variants['pos'].values)
    sample_data=[]
    
    with open(csv_file,'r') as f:
        for line in f:
            parts=line.strip().split(',')
            if len(parts)<6:
                continue
                
            pos=int(parts[1])
            if pos not in chunk_positions:
                continue
                
            # Parse the remaining comma-separated sample:depth:gt triplets
            for sample_dp_gt in parts[6:]:
                if ':' in sample_dp_gt:
                    parts_sample=sample_dp_gt.split(':')
                    if len(parts_sample)>=2:
                        sample=parts_sample[0]
                        dp=int(parts_sample[1]) if parts_sample[1]!='.' else 0
                        gt=parts_sample[2] if len(parts_sample)>2 else './.'
                        
                        sample_data.append({
                            'sample':sample,
                            'pos':pos,
                            'dp':dp,
                            'gt':gt
                        })
    
    sample_df=pd.DataFrame(sample_data)
    return sample_df

def calculate_coverage_imbalance(variant_df,sample_df,site_data,status_data,site_expected_counts):
    """Calculate coverage imbalance statistics across meta groups"""
    
    # Remove duplicates
    sample_df=sample_df.drop_duplicates(subset=['sample','pos'],keep='first').copy()
    
    # Add site and status information to samples
    sample_df['site']=sample_df['sample'].map(lambda x: site_data.get(x,'Unknown'))
    sample_df['status']=sample_df['sample'].map(lambda x: status_data.get(x,-1))
    
    # Get all unique sites for column names
    all_sites=sorted(set(site_data.values()))
    
    # Calculate per-SNP statistics
    imbalance_stats=[]
    
    for _,variant in variant_df.iterrows():
        pos=variant['pos']
        pos_data=sample_df[sample_df['pos']==pos].copy()
        
        if len(pos_data)==0:
            continue
        
        # Define missing samples: no GT, GT=./., or depth=0
        pos_data['is_missing']=(pos_data['gt'].isna()) | (pos_data['gt']=='./.') | (pos_data['dp']==0)
        
        # Calculate missingness per site
        site_missingness={}
        for site in all_sites:
            expected_samples=site_expected_counts.get(site,0)
            if expected_samples>0:
                site_pos_data=pos_data[pos_data['site']==site]
                actual_samples=len(site_pos_data)
                missing_samples=len(site_pos_data[site_pos_data['is_missing']])
                site_missingness[site]=(expected_samples-actual_samples+missing_samples)/expected_samples
            else:
                site_missingness[site]=0
        
        # Case/control analysis
        case_data=pos_data[pos_data['status']==1]
        control_data=pos_data[pos_data['status']==0]
        
        # Count total cases and controls from metadata
        total_cases=len([s for s in status_data.values() if s==1])
        total_controls=len([s for s in status_data.values() if s==0])
        
        # Case missingness
        if total_cases>0:
            case_missing_samples=len(case_data[case_data['is_missing']])
            case_missingness=case_missing_samples/total_cases
        else:
            case_missingness=0
        
        # Control missingness
        if total_controls>0:
            control_missing_samples=len(control_data[control_data['is_missing']])
            control_missingness=control_missing_samples/total_controls
        else:
            control_missingness=0
        
        # Create output dictionary with basic variant info
        stats_dict={
            'chrom':variant['chrom'],
            'pos':pos,
            'ref':variant['ref'],
            'alt':variant['alt'],
            'maf':variant['maf'],
            'ns':variant['ns'],
            'case_missingness':case_missingness,
            'control_missingness':control_missingness
        }
        
        # Add per-site missingness columns
        for site in all_sites:
            stats_dict[f'{site}_missingness']=site_missingness[site]
        
        imbalance_stats.append(stats_dict)
    
    return imbalance_stats

def main():
    parser=argparse.ArgumentParser(description='SNP coverage imbalance analysis')
    parser.add_argument('csv_file',help='Input CSV file from bcftools query')
    parser.add_argument('-o','--output',default='snp_coverage_imbalance.csv',help='Output CSV file')
    parser.add_argument('--metadata',required=True,help='Sample metadata file with sample names and sites')
    parser.add_argument('--chunk-size',type=int,default=100,help='Number of SNPs to process at once')
    
    args=parser.parse_args()
    
    # Parse site metadata
    print(f"Loading site metadata from {args.metadata}")
    site_data,status_data=parse_site_metadata(args.metadata)
    print(f"Loaded site data for {len(site_data)} samples")
    print(f"Loaded status_data for {len(status_data)} samples")
    
    # Add diagnostic prints for metadata
    case_count = sum(1 for s in status_data.values() if s == 1)
    control_count = sum(1 for s in status_data.values() if s == 0)
    print(f"Metadata summary: {case_count} cases, {control_count} controls")
    
    # Show first few samples from metadata
    print("First 5 samples from metadata:")
    for i, (sample, status) in enumerate(list(status_data.items())[:5]):
        print(f"  {sample}: status={status}, site={site_data.get(sample, 'Unknown')}")
    
    # Validate site data
    unique_sites = set(site_data.values())
    print(f"Unique sites found: {unique_sites}")
    
    # Calculate expected samples per site (constant across all SNPs)
    site_expected_counts={}
    for site in unique_sites:
        site_sample_count = len([s for s in site_data.values() if s == site])
        site_expected_counts[site] = site_sample_count
        print(f"  Site {site}: {site_sample_count} samples")
    
    # Parse variant data only
    variant_df=parse_coverage_data(args.csv_file)
    
    # Process SNPs in chunks
    all_imbalance_stats=[]
    total_snps=len(variant_df)
    
    for chunk_start in range(0,total_snps,args.chunk_size):
        chunk_end=min(chunk_start+args.chunk_size,total_snps)
        chunk_variants=variant_df.iloc[chunk_start:chunk_end]
        
        # Parse sample data for this chunk only
        chunk_sample_df=parse_sample_data_for_chunk(args.csv_file,chunk_variants)
        
        # Calculate imbalance statistics for this chunk
        chunk_stats=calculate_coverage_imbalance(chunk_variants,chunk_sample_df,site_data,status_data,site_expected_counts)
        all_imbalance_stats.extend(chunk_stats)
        
        # Clear memory
        del chunk_sample_df
        del chunk_stats
    
    # Combine all results
    imbalance_df=pd.DataFrame(all_imbalance_stats)
    
    # Save results
    imbalance_df.to_csv(args.output,index=False)
    print(f"Results saved to {args.output}")
    
    # Print summary statistics
    print("\nMissingness analysis summary:")
    print(f"Total SNPs analyzed: {len(imbalance_df)}")
    
    print(f"\nCase/control analysis:")
    print(f"Mean case missingness: {imbalance_df['case_missingness'].mean():.4f}")
    print(f"Mean control missingness: {imbalance_df['control_missingness'].mean():.4f}")
    print(f"SNPs with case missingness > 0.1: {(imbalance_df['case_missingness']>0.1).sum()}")
    print(f"SNPs with control missingness > 0.1: {(imbalance_df['control_missingness']>0.1).sum()}")
    
    print(f"\nSite missingness analysis:")
    site_cols=[col for col in imbalance_df.columns if col.endswith('_missingness') and col not in ['case_missingness','control_missingness']]
    for site_col in site_cols:
        site_name=site_col.replace('_missingness','')
        mean_missing=imbalance_df[site_col].mean()
        high_missing=(imbalance_df[site_col]>0.1).sum()
        print(f"  {site_name}: mean={mean_missing:.4f}, high missingness(>0.1)={high_missing}")
    
    print(f"\nOutput columns: {list(imbalance_df.columns)}")

if __name__=='__main__':
    main()
