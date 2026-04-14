import pandas as pd
import argparse

p=argparse.ArgumentParser(description="Flag heterozygosity outliers from a single PLINK .het file (mean +/- 3*SD on F).")
p.add_argument("input_het",help="Input .het file (one per merged dataset)")
p.add_argument("--output-exclusions",required=True,help="Output file for sample exclusions (FID IID)")
p.add_argument("--output-report",required=True,help="Output QC report (n mean sd low high)")

args=p.parse_args()

input_het=args.input_het
output_exclusions=args.output_exclusions
output_report=args.output_report

df=pd.read_csv(input_het,sep=r'\s+',engine='python')
# plink2 .het: #FID IID O(HOM) E(HOM) N(NM) F
fid_col='#FID' if '#FID' in df.columns else 'FID'
if 'F' not in df.columns:
    raise SystemExit("Expected column F in .het file")
n=len(df)
if n<2:
    pd.DataFrame(columns=[fid_col,'IID']).to_csv(output_exclusions,sep=' ',index=False)
    pd.DataFrame([[0,0,0,0,0]],columns=['n','mean','sd','low','high']).to_csv(output_report,sep=' ',index=False)
else:
    mean_F=df['F'].mean()
    sd_F=df['F'].std()
    low=mean_F-3*sd_F
    high=mean_F+3*sd_F
    outliers=df[(df['F']<low)|(df['F']>high)]
    outliers[[fid_col,'IID']].to_csv(output_exclusions,sep=' ',index=False)
    pd.DataFrame([[n,mean_F,sd_F,low,high]],columns=['n','mean','sd','low','high']).to_csv(output_report,sep=' ',index=False)
