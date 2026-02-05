import pandas as pd
from pathlib import Path
import yaml

# Read input files
df=pd.read_csv('glnexus.norm.ref.het.hom.vaf.no_sample.vep.report.csv')
eigenvec_df=pd.read_csv('TECAC_WES_PASS.all_chr.eigenvec',sep='\s+')

# Rename #FID column to FID in eigenvec_df
if '#FID' in eigenvec_df.columns:
    eigenvec_df.rename(columns={'#FID':'FID'},inplace=True)
sites_df=pd.read_csv('seq_sites.txt',sep='\s+',names=['IID','SITE'])

# Filter for pathogenic variants
filtered_df=df[df['Variant.LoF_level']==1][['ID','Gene']]

# Extract chromosome from variant ID
filtered_df['Chr']=filtered_df['ID'].str.split('_').str[0]

# Load gene locations from YAML file
with open(str(Path.home())+'/resources/refGene/hg38/refGene.hg38.locations.yml','r') as f:
    gene_locations=yaml.safe_load(f)

# Get unique list of genes from LoF level 1 variants
unique_genes=filtered_df['Gene'].unique()

# Filter to only include genes with location information
valid_genes = [gene for gene in unique_genes if gene in gene_locations]
print(f"Found {len(unique_genes)} unique genes, {len(valid_genes)} have location information")
print(f"Skipping {len(unique_genes) - len(valid_genes)} genes without location information")

# Filter the dataframe to only include variants in genes with location information
filtered_df = filtered_df[filtered_df['Gene'].isin(valid_genes)]

# Create output directory if it doesn't exist
Path('grouping_files').mkdir(exist_ok=True)

# Group by chromosome and gene
grouped=filtered_df.groupby(['Chr','Gene'])

# Get unique chromosomes
unique_chromosomes = filtered_df['Chr'].unique()

# Process each chromosome separately
for chr_name in unique_chromosomes:
    output_lines=[]
    
    # Get all groups for this chromosome
    chr_groups = [(key, group) for key, group in grouped if key[0] == chr_name]
    
    # Process each gene in the chromosome
    for (_, gene), group in chr_groups:
        variants=group['ID'].tolist()
        # Create the two lines for each gene: variants and annotations
        output_lines.append(f"{gene} var {' '.join(variants)}")
        output_lines.append(f"{gene} anno {' '.join(['group1']*len(variants))}")
    
    # Write to chromosome-specific file
    output_file=f"grouping_files/group_varID.{chr_name}.txt"
    with open(output_file,'w') as f:
        f.write('\n'.join(output_lines))

# Create gene_loc.txt file with specific headers
with open('gene_loc.txt','w') as f:
    f.write("gene_id\tchromosome\tseq_region_start\tseq_region_end\tgene_symbol\n")
    for gene in valid_genes:
        loc=gene_locations[gene]
        chromosome=loc['chromosome']
        # Add 'chr' prefix to chromosome
        chr_formatted=f"chr{chromosome}"
        start=loc['start']
        end=loc['end']
        f.write(f"{gene}\t{chr_formatted}\t{start}\t{end}\t{gene}\n")

# Process phenotype file
# Convert SITE to binary indicators (0/1)
site_dummies=pd.get_dummies(sites_df['SITE']).astype(int)

# Combine with IID
pheno_df=pd.concat([sites_df['IID'],site_dummies],axis=1)

# Read control samples file
controls=pd.read_csv('controls.txt',header=None)[0].tolist()

# Add STATUS column - 0 for controls, 1 for cases
pheno_df['STATUS']=pheno_df['IID'].apply(lambda x:0 if x in controls else 1)

# Fill any missing site values with 0
site_columns=[col for col in pheno_df.columns if col not in ['IID','STATUS']]
pheno_df[site_columns]=pheno_df[site_columns].fillna(0)

# Extract PC1-PC10 from eigenvec_df
pc_columns = [f'PC{i}' for i in range(1, 11)]
eigenvec_pc_df = eigenvec_df[['IID'] + pc_columns]

# Merge PCs into phenotype dataframe
pheno_df = pd.merge(pheno_df, eigenvec_pc_df, on='IID', how='left')

# Write phenotype file as CSV
pheno_df.to_csv('TECAC_WES.saige_pheno.csv',index=False)
