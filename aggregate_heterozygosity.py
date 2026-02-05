import pandas as pd
import argparse

p=argparse.ArgumentParser(description="Aggregate and filter sample heterozygosity across multiple input files.")
p.add_argument("--missingness-exclusions",required=True,help="Samples to exclude based on missingness (2 columns: FID IID)")
p.add_argument("--output-aggregated",required=True,help="Output file for aggregated heterozygosity per sample")
p.add_argument("--output-exclusions",required=True,help="Output file for sample exclusions (outlier Fhat)")
p.add_argument("--output-report",required=True,help="Output QC report file")
p.add_argument("input_files",nargs='+',help="Input .het files")

args=p.parse_args()

missingness_exclusions_file=args.missingness_exclusions
output_aggregated=args.output_aggregated
output_exclusions=args.output_exclusions
output_report=args.output_report
input_files=args.input_files

missingness_exclusions=pd.read_csv(missingness_exclusions_file,sep=' ',header=None,names=['IID'])

dfs=[]
for f in input_files:
    df=pd.read_csv(f,sep='\\s+')
    dfs.append(df[['#IID','OBS_CT','O(HOM)']])

combined=pd.concat(dfs)
agg=combined.groupby('#IID').sum().reset_index()
agg['N_NM']=agg['OBS_CT']
agg['O_HOM']=agg['O(HOM)']
agg['Fhat']=(agg['N_NM']-agg['O_HOM'])/agg['N_NM']

agg_filtered=agg[~agg['#IID'].isin(missingness_exclusions['IID'])]
agg_filtered.to_csv(output_aggregated,sep='\t',index=False)

mean_F=agg_filtered['Fhat'].mean()
sd_F=agg_filtered['Fhat'].std()
lower_bound=mean_F-3*sd_F
upper_bound=mean_F+3*sd_F

het_outliers=agg_filtered[(agg_filtered['Fhat']<lower_bound)|(agg_filtered['Fhat']>upper_bound)]
het_outliers['#IID'].to_csv(output_exclusions,sep=' ',index=False,header=False)
het_outliers[['#IID','Fhat']].to_csv(output_report,sep=',',index=False,header=True)
