import pandas as pd
import argparse

p=argparse.ArgumentParser(description="Summarize differential missingness test results.")
p.add_argument("input_files",nargs='+',help="Input .test.missing files")
p.add_argument("--threshold",type=float,required=True,help="P-value threshold for differential missingness")
p.add_argument("--output",required=True,help="Output summary CSV file")

args=p.parse_args()

dfs=[]
for f in args.input_files:
    try:
        df=pd.read_csv(f,sep=r'\s+')
        dfs.append(df)
    except:
        continue

if len(dfs)==0:
    exit(1)

combined=pd.concat(dfs,ignore_index=True)

total=len(combined)
failed=(combined['P']<args.threshold).sum()
median_p=combined['P'].median()
mean_cases=combined['F_MISS_A'].mean()
mean_controls=combined['F_MISS_U'].mean()

summary=pd.DataFrame({
    'Metric':['Total variants tested','Variants with diff missingness P < '+str(args.threshold),'Percentage failed','Median P-value','Mean missingness in cases','Mean missingness in controls'],
    'Value':[total,failed,f"{100*failed/total:.2f}%",f"{median_p:.4e}",f"{mean_cases:.4f}",f"{mean_controls:.4f}"]
})

summary.to_csv(args.output,index=False)
print(f"Summary saved to {args.output}")
print(f"Total variants: {total}, Failed: {failed} ({100*failed/total:.2f}%)")
