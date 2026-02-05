import pandas as pd
import argparse

p=argparse.ArgumentParser(description="Aggregate and filter sample missingness across multiple input files.")
p.add_argument("input_files",nargs='+',help="Input missingness summary files")
p.add_argument("--mind-thr",type=float,required=True,help="Missingness threshold (F_MISS > mind_thr will be excluded)")
p.add_argument("--output-aggregated",required=True,help="Output file for aggregated missingness per sample")
p.add_argument("--output-exclusions",required=True,help="Output file for sample exclusions (exceeding missingness)")
p.add_argument("--output-report",required=True,help="Output QC report file")

args=p.parse_args()

input_files=args.input_files
mind_thr=args.mind_thr
output_aggregated=args.output_aggregated
output_exclusions=args.output_exclusions
output_report=args.output_report

dfs=[]
for f in input_files:
    df=pd.read_csv(f,sep='\\s+')
    dfs.append(df[['#IID','MISSING_CT','OBS_CT']])

combined=pd.concat(dfs)
agg=combined.groupby('#IID').sum().reset_index()
agg['F_MISS']=agg['MISSING_CT']/agg['OBS_CT']

agg.to_csv(output_aggregated,sep='\t',index=False)

high_miss=agg[agg['F_MISS']>mind_thr]
high_miss['#IID'].to_csv(output_exclusions,sep=' ',index=False,header=False)
high_miss[['#IID','F_MISS']].to_csv(output_report,sep=',',index=False,header=True)
