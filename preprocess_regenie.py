import pandas as pd
import argparse

def get_args():
    p=argparse.ArgumentParser(description='Preprocess data for regenie analysis')
    p.add_argument('--vep-csv',required=True,help='Input VEP CSV file')
    p.add_argument('--eigenvec',required=True,help='PLINK eigenvec file')
    p.add_argument('--sites',required=True,help='Site mapping file (IID SITE format)')
    p.add_argument('--controls',required=True,help='Controls list file')
    p.add_argument('--remove-PCs',action='store_true',help='Remove PCs from eigenvec file')
    p.add_argument('-O','--output-prefix',default='',help='Output prefix (optional)')
    return p.parse_args()

def main():
    args=get_args()
    
    # Set up output prefix
    prefix=f"{args.output_prefix}." if args.output_prefix else ""
    
    df=pd.read_csv(args.vep_csv)
    
    # Create annotation file with both pathogenic and vus
    pathogenic_df=df[df['Variant.LoF_level']==1][['ID','Gene']].copy()
    pathogenic_df['Status']='pathogenic'
    vus_df=df[df['Variant.LoF_level']==2][['ID','Gene']].copy()
    vus_df['Status']='vus'
    
    # Combine and write annotation file
    annotation_df=pd.concat([pathogenic_df,vus_df],ignore_index=True)
    annotation_df.to_csv(f'{prefix}regenie.annotation.txt',sep=' ',header=False,index=False)
    
    # Create set file with both categories
    # Group by gene and get all variants (both pathogenic and vus)
    all_filtered_df=df[df['Variant.LoF_level'].isin([1,2])]
    gene_groups=all_filtered_df.groupby('Gene')['ID'].agg(list).reset_index()
    
    # Extract chromosome and position from first variant ID for each gene
    def get_chr_pos(variant_list):
        if not variant_list:
            return '.',0
        first_variant=variant_list[0]
        chr_pos=first_variant.split('_')[:2]
        return chr_pos[0],int(chr_pos[1])
    
    gene_groups[['Chr','Pos']]=pd.DataFrame(gene_groups['ID'].apply(get_chr_pos).tolist())
    
    # Format variant lists as comma-separated strings
    gene_groups['Variants']=gene_groups['ID'].apply(lambda x: ','.join(x))
    
    # Create output dataframe with required columns
    output_df=gene_groups[['Gene','Chr','Pos','Variants']]
    output_df.to_csv(f'{prefix}regenie.set.txt',sep=' ',header=False,index=False)
    
    # Create mask file
    with open(f'{prefix}regenie.mask.txt','w') as f:
        f.write('M1 pathogenic\n')
        f.write('M2 vus\n')
        f.write('M3 pathogenic,vus\n')
    
    # Read eigenvec file with headers
    eigenvec_df=pd.read_csv(args.eigenvec,sep='\s+')
    eigenvec_df=eigenvec_df.rename(columns={'#FID':'FID'})
    
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

if __name__=='__main__':
    main()
