import argparse
import os

import pandas as pd


def get_args():
    p=argparse.ArgumentParser(description='Preprocess data for regenie analysis')
    p.add_argument('--vep-list-file',required=True,help='File containing list of VEP CSV files (one per line)')
    p.add_argument('--ancestry-file',help='FID/#FID, IID/#IID, then ancestry columns; renamed DM1..DMn for regenie')
    p.add_argument('--eigenvec',help=argparse.SUPPRESS)
    p.add_argument('--dm-count',type=int,default=None,help='If set, take only first dm-count ancestry columns after IID (default: keep all)')
    p.add_argument('--ancestry-columns',default=None,help='Comma-separated names to take by order; renamed DM1..DMn')
    p.add_argument('--samples',required=True,help='Samples list file (one IID per line)')
    p.add_argument('--sites',help='Site mapping file (IID SITE) optional')
    p.add_argument('--controls',help='Controls list file (only needed if covariates has no STATUS)')
    p.add_argument('--covariates',help='Covariates file with FID IID STATUS; any other columns become extra covariates (e.g. FREEZE)')
    p.add_argument('--negative-control-blacklist',default='',help='Optional file with variant IDs (one per line) to exclude from synonymous negative-control mask')
    p.add_argument('-O','--output-prefix',default='',help='Output directory prefix')
    return p.parse_args()


def is_synonymous_consequence(value):
    if pd.isna(value):
        return False
    text=str(value).strip().lower()
    if not text:
        return False
    return 'synonymous_variant' in text


def load_ancestry_table(path,dm_count,ancestry_columns):
    df=pd.read_csv(path,sep=r'\s+')
    if '#FID' in df.columns:
        df=df.rename(columns={'#FID':'FID'})
    if '#IID' in df.columns:
        df=df.rename(columns={'#IID':'IID'})
    if 'FID' not in df.columns or 'IID' not in df.columns:
        print('Error: ancestry file needs FID (or #FID) and IID (or #IID)')
        return None,0
    if df.columns[0]!='FID':
        cols=df.columns.tolist()
        cols.remove('FID')
        df=df[['FID']+cols]
    rest=[c for c in df.columns if c not in ('FID','IID')]
    if ancestry_columns:
        missing=[c for c in ancestry_columns if c not in df.columns]
        if missing:
            print(f'Error: ancestry columns not in file: {missing}')
            return None,0
        take=list(ancestry_columns)
    else:
        # Default behavior: keep all ancestry columns, unless dm_count is set.
        if dm_count is not None:
            if dm_count <= 0:
                print('Error: --dm-count must be > 0')
                return None,0
            take=rest[:dm_count]
        else:
            take=rest
    dm_n=len(take)
    rename={take[i]:f'DM{i+1}' for i in range(dm_n)}
    df=df.rename(columns=rename)
    out=df[['FID','IID']+[f'DM{i}' for i in range(1,dm_n+1)]]
    return out,dm_n


def main():
    args=get_args()
    ancestry_path=args.ancestry_file or args.eigenvec
    ancestry_col_list=None
    if args.ancestry_columns:
        ancestry_col_list=[c.strip() for c in args.ancestry_columns.split(',') if c.strip()]

    prefix=(args.output_prefix.rstrip('/')+'/') if args.output_prefix else ''

    with open(args.vep_list_file,'r') as f:
        vep_files=[line.strip() for line in f if line.strip()]

    print(f'Found {len(vep_files)} VEP files in list')

    missing_files=[f for f in vep_files if not os.path.exists(f)]
    if missing_files:
        print('Error: missing VEP files:')
        for f in missing_files:
            print(f'  {f}')
        return

    required_columns=['ID','Gene','Variant.LoF_level','Variant.Consequence']
    annotation_file=f'{prefix}regenie.annotation.txt'
    set_file=f'{prefix}regenie.set.txt'

    blacklist_ids=set()
    if args.negative_control_blacklist and os.path.exists(args.negative_control_blacklist):
        with open(args.negative_control_blacklist,'r') as f:
            blacklist_ids={line.strip() for line in f if line.strip()}
        print(f'Loaded {len(blacklist_ids)} blacklisted variant IDs for M4 negative control')
    elif args.negative_control_blacklist:
        print(f'Warning: negative-control blacklist not found: {args.negative_control_blacklist} (continuing without blacklist)')

    variant_rows=[]
    for vep_file in vep_files:
        print(f'Processing VEP file: {vep_file}')
        df=pd.read_csv(vep_file,usecols=required_columns)
        variant_rows.append(df)

    if variant_rows:
        all_variants=pd.concat(variant_rows,ignore_index=True)
    else:
        all_variants=pd.DataFrame(columns=required_columns)

    all_variants=all_variants.dropna(subset=['ID','Gene']).copy()
    all_variants['ID']=all_variants['ID'].astype(str).str.strip()
    all_variants['Gene']=all_variants['Gene'].astype(str).str.strip()
    all_variants=all_variants[(all_variants['ID']!='') & (all_variants['Gene']!='')]
    all_variants=all_variants.drop_duplicates(subset=['ID','Gene'])
    all_variants['Variant.LoF_level']=pd.to_numeric(all_variants['Variant.LoF_level'],errors='coerce')

    pathogenic_df=all_variants[all_variants['Variant.LoF_level']==1][['ID','Gene']].drop_duplicates()
    vus_df=all_variants[all_variants['Variant.LoF_level']==2][['ID','Gene']].drop_duplicates()

    syn_df=all_variants[all_variants['Variant.Consequence'].apply(is_synonymous_consequence)][['ID','Gene']].drop_duplicates()
    if len(blacklist_ids)>0:
        syn_df=syn_df[~syn_df['ID'].isin(blacklist_ids)]
    excluded_ids=set(pathogenic_df['ID'].tolist()) | set(vus_df['ID'].tolist())
    syn_df=syn_df[~syn_df['ID'].isin(excluded_ids)]
    syn_df=syn_df.drop_duplicates()

    with open(annotation_file,'w') as ann_f:
        for _,row in pathogenic_df.iterrows():
            ann_f.write(f"{row['ID']} {row['Gene']} pathogenic\n")
        for _,row in vus_df.iterrows():
            ann_f.write(f"{row['ID']} {row['Gene']} vus\n")
        for _,row in syn_df.iterrows():
            ann_f.write(f"{row['ID']} {row['Gene']} synonymous\n")

    set_variants=pd.concat([pathogenic_df,vus_df,syn_df],ignore_index=True).drop_duplicates(subset=['ID','Gene'])
    with open(set_file,'w') as set_f:
        gene_groups=set_variants.groupby('Gene')['ID'].agg(list)
        for gene,variants in gene_groups.items():
            first_variant=variants[0]
            chr_pos=first_variant.split('_')[:2]
            if len(chr_pos)<2:
                continue
            chr_name=chr_pos[0]
            try:
                pos=int(chr_pos[1])
            except ValueError:
                continue
            variant_str=','.join(variants)
            set_f.write(f'{gene} {chr_name} {pos} {variant_str}\n')

    print(f'Processed {len(vep_files)} VEP files')
    print(f"M1 pathogenic variants: {len(pathogenic_df)}")
    print(f"M2 VUS variants: {len(vus_df)}")
    print(f"M4 synonymous negative-control variants: {len(syn_df)}")

    with open(f'{prefix}regenie.mask.txt','w') as f:
        f.write('M1 pathogenic\n')
        f.write('M2 vus\n')
        f.write('M3 pathogenic,vus\n')
        f.write('M4 synonymous\n')

    samples_df=pd.read_csv(args.samples,header=None)
    if len(samples_df)>0 and samples_df.iloc[0,0] in ['IID','FID','#IID','#FID']:
        samples=samples_df.iloc[1:,0].tolist()
        print('Detected header in samples file, skipping first row')
    else:
        samples=samples_df[0].tolist()
    print(f'Loaded {len(samples)} samples from samples file')

    covar_df=None
    if ancestry_path:
        ancestry_df,dm_n=load_ancestry_table(ancestry_path,args.dm_count,ancestry_col_list)
        if ancestry_df is None:
            return
        ancestry_df=ancestry_df[ancestry_df['IID'].isin(samples)]
        print(f'Filtered ancestry table to {len(ancestry_df)} samples')
        merged_df=ancestry_df.copy()
    else:
        merged_df=pd.DataFrame({'FID':samples,'IID':samples})

    if args.sites:
        print(f'Loading sites file: {args.sites}')
        sites_df=pd.read_csv(args.sites,sep=r'\s+',names=['IID','SITE'])
        sites_df=sites_df[sites_df['IID'].isin(samples)]
        print(f'Filtered sites to {len(sites_df)} samples')
        site_dummies=pd.get_dummies(sites_df['SITE']).astype(int)
        sites_df=pd.concat([sites_df['IID'],site_dummies],axis=1)
        merged_df=pd.merge(merged_df,sites_df,on='IID',how='left')
        site_columns=[col for col in merged_df.columns if col in site_dummies.columns]
        merged_df[site_columns]=merged_df[site_columns].fillna(0).astype(int)
        print(f'Sites loaded for {len(sites_df)} samples with {len(site_dummies.columns)} site indicators')

    if args.covariates:
        print(f'Loading covariates file: {args.covariates}')
        covar_df=pd.read_csv(args.covariates,sep=r'\s+')
        if 'IID' not in covar_df.columns:
            print("Error: covariates file must contain 'IID'")
            return
        covar_df=covar_df[covar_df['IID'].isin(samples)]
        extra_cols=[c for c in covar_df.columns if c not in ('FID','IID','STATUS')]
        if extra_cols:
            covar_extra=covar_df[['IID']+extra_cols]
            merged_df=pd.merge(merged_df,covar_extra,on='IID',how='left')
            print(f'Covariates loaded: {extra_cols}')

    id_cols=[c for c in ['FID','IID'] if c in merged_df.columns]
    dm_cols=[c for c in merged_df.columns if c.startswith('DM')]
    dm_cols=sorted(dm_cols,key=lambda x:int(x[2:]))
    other_cols=[c for c in merged_df.columns if c not in id_cols and c not in dm_cols]
    merged_df=merged_df[id_cols+dm_cols+other_cols]
    merged_df.to_csv(f'{prefix}regenie.covar.txt',sep=' ',index=False)

    if args.covariates and covar_df is not None and 'STATUS' in covar_df.columns:
        pheno_df=covar_df[['FID','IID','STATUS']].copy()
        pheno_df=pheno_df[pheno_df['STATUS'].isin([1,2])]
        pheno_df['STATUS']=pheno_df['STATUS'].map({1:0,2:1})
        print(f'Pheno from covariates: n={len(pheno_df)}')
    else:
        if not args.controls:
            print('Error: need --controls when covariates has no STATUS')
            return
        if ancestry_path:
            pheno_df=ancestry_df[['FID','IID']].copy()
        else:
            pheno_df=pd.DataFrame({'FID':samples,'IID':samples})
        controls_df=pd.read_csv(args.controls,header=None)
        if len(controls_df)>0 and controls_df.iloc[0,0] in ['IID','FID','#IID','#FID']:
            controls=controls_df.iloc[1:,0].tolist()
        else:
            controls=controls_df[0].tolist()
        pheno_df['STATUS']=pheno_df['IID'].apply(lambda x:0 if x in controls else 1)
    pheno_df.to_csv(f'{prefix}regenie.pheno.txt',sep=' ',index=False)


if __name__=='__main__':
    main()
