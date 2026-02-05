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

def get_variant_count(csv_file):
    """Count total number of variants without loading into memory"""
    count=0
    with open(csv_file,'r') as f:
        for line in f:
            parts=line.strip().split(',')
            if len(parts)>=6:
                count+=1
    return count


def parse_data(csv_file,chunk_size):
    """Parse data efficiently by reading file once and processing in chunks"""
    print("Reading CSV file and parsing data...")
    
    all_variants=[]
    all_sample_data=[]
    current_chunk=0
    current_pos=0
    chunk_variants=[]
    chunk_sample_data=[]
    
    with open(csv_file,'r') as f:
        for line_num,line in enumerate(f,1):
            if line_num%10000==0:
                print(f"Processed {line_num} lines...")
                
            parts=line.strip().split(',')
            if len(parts)<6:
                continue
                
            chrom,pos,ref,alt,maf,ns=parts[:6]
            pos=int(pos)
            
            # Add variant to current chunk
            chunk_variants.append({
                'chrom':chrom,
                'pos':pos,
                'ref':ref,
                'alt':alt,
                'maf':float(maf) if maf!='.' else 0.0,
                'ns':int(ns) if ns!='.' else 0
            })
            
            # Parse sample data for this variant
            for sample_dp_gt in parts[6:]:
                if ':' in sample_dp_gt:
                    parts_sample=sample_dp_gt.split(':')
                    if len(parts_sample)>=2:
                        sample=parts_sample[0]
                        dp=int(parts_sample[1]) if parts_sample[1]!='.' else 0
                        gt=parts_sample[2] if len(parts_sample)>2 else './.'
                        
                        chunk_sample_data.append({
                            'sample':sample,
                            'pos':pos,
                            'dp':dp,
                            'gt':gt
                        })
            
            current_pos+=1
            
            # Process chunk when it reaches chunk_size
            if len(chunk_variants)>=chunk_size:
                print(f"Processing chunk {current_chunk+1} with {len(chunk_variants)} variants...")
                
                # Convert to DataFrames
                variant_df=pd.DataFrame(chunk_variants)
                sample_df=pd.DataFrame(chunk_sample_data)
                
                # Remove duplicates
                sample_df=sample_df.drop_duplicates(subset=['sample','pos'],keep='first').copy()
                
                all_variants.append(variant_df)
                all_sample_data.append(sample_df)
                
                # Reset for next chunk
                chunk_variants=[]
                chunk_sample_data=[]
                current_chunk+=1
    
    # Process remaining data if any
    if chunk_variants:
        print(f"Processing final chunk {current_chunk+1} with {len(chunk_variants)} variants...")
        variant_df=pd.DataFrame(chunk_variants)
        sample_df=pd.DataFrame(chunk_sample_data)
        sample_df=sample_df.drop_duplicates(subset=['sample','pos'],keep='first').copy()
        all_variants.append(variant_df)
        all_sample_data.append(sample_df)
    
    print(f"Parsed {len(all_variants)} chunks total")
    return all_variants,all_sample_data

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
        
        # Define missing samples: no GT, GT=./., depth=0, or DP='.'
        pos_data['is_missing']=(pos_data['gt'].isna())|(pos_data['gt']=='./.')|(pos_data['dp']==0)|(pos_data['dp']=='.')
        
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
    
    # Parse data efficiently in one pass
    all_variants,all_sample_data=parse_data(args.csv_file,args.chunk_size)
    
    # Process each chunk and write incrementally
    header_written=False
    
    for i,(chunk_variants,chunk_sample_df) in enumerate(zip(all_variants,all_sample_data)):
        print(f"Calculating statistics for chunk {i+1}/{len(all_variants)}...")
        
        # Calculate imbalance statistics for this chunk
        chunk_stats=calculate_coverage_imbalance(chunk_variants,chunk_sample_df,site_data,status_data,site_expected_counts)
        
        # Convert to DataFrame and write immediately
        chunk_df=pd.DataFrame(chunk_stats)
        
        # Write header only on first chunk
        if not header_written:
            chunk_df.to_csv(args.output,index=False)
            header_written=True
            print(f"Started writing results to {args.output}")
        else:
            # Append without header
            chunk_df.to_csv(args.output,mode='a',header=False,index=False)
        
        print(f"  Wrote {len(chunk_stats)} variants to output file")
        
        # Clear memory
        del chunk_variants
        del chunk_sample_df
        del chunk_stats
        del chunk_df
    
    print(f"Results saved to {args.output}")
    
    # Print summary statistics
    print("\nMissingness analysis summary:")
    print(f"Results written to {args.output}")
    print("Summary statistics can be calculated from the output file if needed")

if __name__=='__main__':
    main()
