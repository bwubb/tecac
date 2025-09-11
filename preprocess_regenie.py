import pandas as pd
import argparse

def get_args():
    p=argparse.ArgumentParser(description='Preprocess data for regenie analysis')
    p.add_argument('--vep-list-file',required=True,help='File containing list of VEP CSV files (one per line)')
    p.add_argument('--eigenvec',help='PLINK eigenvec file (optional - skip if using pre-made files)')
    p.add_argument('--sites',help='Site mapping file (IID SITE format) - required if processing eigenvec')
    p.add_argument('--controls',help='Controls list file - required if processing eigenvec')
    p.add_argument('--covariates',help='Covariates file (IID and additional columns) - optional, appends additional covariates by matching IID')
    p.add_argument('--remove-PCs',action='store_true',help='Remove PCs from eigenvec file')
    p.add_argument('--annotation-only',action='store_true',help='Only generate annotation files (skip eigenvec processing)')
    p.add_argument('-O','--output-prefix',default='',help='Output prefix (optional)')
    return p.parse_args()

def main():
    args=get_args()
    
    # Set up output prefix
    prefix=f"{args.output_prefix}." if args.output_prefix else ""
    
    # Read VEP files list and validate files exist
    with open(args.vep_list_file,'r') as f:
        vep_files=[line.strip() for line in f if line.strip()]
    
    print(f"Found {len(vep_files)} VEP files in list")
    
    # Check that all files exist
    import os
    missing_files=[f for f in vep_files if not os.path.exists(f)]
    if missing_files:
        print(f"Error: The following VEP files do not exist:")
        for f in missing_files:
            print(f"  {f}")
        return
    
    # Process VEP CSV file(s) one at a time to save memory
    required_columns=['ID','Gene','Variant.LoF_level']
    
    # Initialize output files
    annotation_file=f'{prefix}regenie.annotation.txt'
    set_file=f'{prefix}regenie.set.txt'
    
    # Open output files for writing
    with open(annotation_file,'w') as ann_f, open(set_file,'w') as set_f:
        # Process each VEP file
        for vep_file in vep_files:
            print(f"Processing VEP file: {vep_file}")
            df=pd.read_csv(vep_file,usecols=required_columns)
            
            # Filter for pathogenic and VUS variants
            pathogenic_df=df[df['Variant.LoF_level']==1][['ID','Gene']]
            vus_df=df[df['Variant.LoF_level']==2][['ID','Gene']]
            
            # Write annotation file entries
            for _,row in pathogenic_df.iterrows():
                ann_f.write(f"{row['ID']} {row['Gene']} pathogenic\n")
            for _,row in vus_df.iterrows():
                ann_f.write(f"{row['ID']} {row['Gene']} vus\n")
            
            # Write set file entries (group by gene)
            all_filtered_df=df[df['Variant.LoF_level'].isin([1,2])]
            gene_groups=all_filtered_df.groupby('Gene')['ID'].agg(list)
            
            for gene,variants in gene_groups.items():
                # Get chromosome and position from first variant
                first_variant=variants[0]
                chr_pos=first_variant.split('_')[:2]
                chr_name=chr_pos[0]
                pos=int(chr_pos[1])
                
                # Format variant list
                variant_str=','.join(variants)
                set_f.write(f"{gene} {chr_name} {pos} {variant_str}\n")
    
    print(f"Processed {len(vep_files)} VEP files")
    
    # Create mask file
    with open(f'{prefix}regenie.mask.txt','w') as f:
        f.write('M1 pathogenic\n')
        f.write('M2 vus\n')
        f.write('M3 pathogenic,vus\n')
    
    # Process eigenvec only if not in annotation-only mode and eigenvec is provided
    if not args.annotation_only and args.eigenvec:
        if not args.sites or not args.controls:
            print("Error: --sites and --controls are required when processing eigenvec file")
            return
        
        # Read eigenvec file with headers
        eigenvec_df=pd.read_csv(args.eigenvec,sep='\s+')
        
        # Handle #IID column (rename to IID)
        if '#IID' in eigenvec_df.columns:
            eigenvec_df=eigenvec_df.rename(columns={'#IID':'IID'})
        
        # Add FID column if it doesn't exist (copy of IID)
        if 'FID' not in eigenvec_df.columns:
            eigenvec_df['FID']=eigenvec_df['IID']
        
        # Read key-value mapping file and create binary site indicators
        sites_df=pd.read_csv(args.sites,sep='\s+',names=['IID','SITE'])
        # Convert SITE to binary indicators (one-hot encoding) with 0/1 integers
        site_dummies=pd.get_dummies(sites_df['SITE']).astype(int)
        # Combine with IID
        sites_df=pd.concat([sites_df['IID'],site_dummies],axis=1)
        
        # Merge the dataframes on ID column
        merged_df=pd.merge(eigenvec_df,sites_df,on='IID',how='left')
        
        # Fill any NaN values with 0 (for samples not in any site)
        site_columns=[col for col in merged_df.columns if col in site_dummies.columns]
        merged_df[site_columns]=merged_df[site_columns].fillna(0).astype(int)
        
        # Add covariates if provided
        if args.covariates:
            print(f"Loading covariates file: {args.covariates}")
            covar_df=pd.read_csv(args.covariates,sep='\s+')
            # Ensure IID column exists
            if 'IID' not in covar_df.columns:
                print("Error: Covariates file must contain 'IID' column")
                return
            merged_df=pd.merge(merged_df,covar_df,on='IID',how='left')
            print(f"Covariates loaded for {len(covar_df)} samples with {len(covar_df.columns)-1} additional columns")
        
        if args.remove_PCs:
            # Remove PC columns
            pc_columns=[col for col in merged_df.columns if col.startswith('PC')]
            merged_df=merged_df.drop(columns=pc_columns)
        
        # Write output file
        merged_df.to_csv(f'{prefix}regenie.covar.txt',sep=' ',index=False)
        
        # Read control samples file
        controls=pd.read_csv(args.controls,header=None)[0].tolist()
        
        # Create phenotype dataframe from merged_df
        pheno_df=eigenvec_df[['FID','IID']].copy()
        
        # Add STATUS column - 0 for controls, 1 for cases
        pheno_df['STATUS']=pheno_df['IID'].apply(lambda x: 0 if x in controls else 1)
        
        # Write phenotype file
        pheno_df.to_csv(f'{prefix}regenie.pheno.txt',sep=' ',index=False)
    elif args.annotation_only:
        print("Annotation-only mode: Skipping eigenvec processing")
    else:
        print("No eigenvec file provided - only VEP files generated")

if __name__=='__main__':
    main()
