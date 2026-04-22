import argparse
import csv
import math
import os
import sys
import pandas as pd

p=argparse.ArgumentParser(description="Controls-only freeze differential missingness and carrier QC per variant.")
p.add_argument("--matrix",required=True,help="Variant x sample GT matrix from bcftools query")
p.add_argument("--sample-order",required=True,help="Sample order file from bcftools query -l")
p.add_argument("--sample-freeze",required=True,help="Sample to freeze map TSV (IID TAB FREEZE)")
p.add_argument("--output",required=True,help="Output TSV path")
p.add_argument("--min-abs-missing-diff",type=float,default=0.03,help="Flag threshold for absolute missingness-rate difference")
p.add_argument("--high-missing-threshold",type=float,default=0.05,help="Flag threshold for high missingness in either freeze")
p.add_argument("--carrier-fisher-p-threshold",type=float,default=1e-4,help="Flag threshold for carrier-rate Fisher p-value")
p.add_argument("--missing-fisher-p-threshold",type=float,default=1e-6,help="Flag threshold for missingness Fisher p-value")
args=p.parse_args()
FREEZE2="2"
FREEZE3="3"

def is_missing(gt):
    g=str(gt).strip()
    return g=='' or '.' in g

def is_carrier(gt):
    if is_missing(gt):
        return False
    alleles=gt.replace('|','/').split('/')
    return any(a not in {'0',''} for a in alleles)

def logchoose(n,k):
    if k<0 or k>n:
        return float('-inf')
    return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)

def hypergeom_prob(a,r1,r2,c1,n):
    return math.exp(logchoose(r1,a)+logchoose(r2,c1-a)-logchoose(n,c1))

def fisher_exact_two_sided(a,b,c,d):
    r1=a+b
    r2=c+d
    c1=a+c
    n=r1+r2
    low=max(0,c1-r2)
    high=min(r1,c1)
    p_obs=hypergeom_prob(a,r1,r2,c1,n)
    p=0.0
    eps=1e-12
    for x in range(low,high+1):
        px=hypergeom_prob(x,r1,r2,c1,n)
        if px<=p_obs+eps:
            p+=px
    return min(max(p,0.0),1.0)

def safe_rate(num,den):
    if den<=0:
        return None
    return num/den

def should_flag(miss_f2,miss_f3,miss_p,carr_p,carr_abs,total_called_carriers):
    reasons=[]
    if miss_f2 is not None and miss_f3 is not None:
        if abs(miss_f2-miss_f3)>args.min_abs_missing_diff:
            reasons.append('missing_abs_diff')
        if max(miss_f2,miss_f3)>args.high_missing_threshold:
            reasons.append('high_missing_either_freeze')
    if miss_p is not None and miss_p<args.missing_fisher_p_threshold:
        reasons.append('missing_fisher')
    # Intentionally hardcoded for a mild first-pass screen.
    if carr_p is not None and carr_abs is not None and carr_p<args.carrier_fisher_p_threshold and carr_abs>=0.01 and total_called_carriers>=3:
        reasons.append('carrier_fisher')
    return (1 if reasons else 0,','.join(reasons))

sample_order=pd.read_csv(args.sample_order,sep='\t',header=None,dtype=str)[0].astype(str).tolist()
sample_freeze_df=pd.read_csv(args.sample_freeze,sep='\t',header=None,names=['IID','FREEZE'],dtype=str)
sample_freeze={str(r['IID']):str(r['FREEZE']) for _,r in sample_freeze_df.iterrows()}
missing_samples=[s for s in sample_order if s not in sample_freeze]
if len(missing_samples)>0:
    raise ValueError(f"sample_freeze is missing {len(missing_samples)} sample IDs from sample_order. First 10: {missing_samples[:10]}")
freeze_vec=[sample_freeze[s] for s in sample_order]
idx_f2=[i for i,f in enumerate(freeze_vec) if f==FREEZE2]
idx_f3=[i for i,f in enumerate(freeze_vec) if f==FREEZE3]
if not idx_f2 or not idx_f3:
    raise ValueError("Need controls in both freeze groups from sample files.")
print(f"Controls in sample order: {len(sample_order)}",file=sys.stderr)
print(f"Freeze 2 controls: {len(idx_f2)}",file=sys.stderr)
print(f"Freeze 3 controls: {len(idx_f3)}",file=sys.stderr)

out_cols=[
    'variant_id','variant_type','freeze2_label','freeze3_label',
    'f2_n_total','f2_n_missing','f2_n_nonmissing',
    'f3_n_total','f3_n_missing','f3_n_nonmissing',
    'total_missing','total_called',
    'f2_missing_rate','f3_missing_rate','abs_missing_rate_diff','missing_fisher_p',
    'f2_n_called','f2_n_carrier','f2_n_noncarrier',
    'f3_n_called','f3_n_carrier','f3_n_noncarrier',
    'total_called_carriers',
    'f2_carrier_rate_called','f3_carrier_rate_called','abs_carrier_rate_diff','carrier_fisher_p',
    'flag','flag_reasons'
]

os.makedirs(os.path.dirname(args.output) or '.',exist_ok=True)
with open(args.output,'w',newline='') as out:
    w=csv.writer(out,delimiter='\t')
    w.writerow(out_cols)
    n_total=0
    n_flag=0
    n_missing_fisher=0
    n_missing_abs_diff=0
    n_high_missing_either_freeze=0
    n_carrier_fisher=0
    with open(args.matrix,'r',newline='') as f:
        for line in f:
            if not line.strip():
                continue
            parts=line.rstrip('\n').split('\t')
            if len(parts)<2:
                continue
            n_total+=1
            vid=parts[0]
            vtype=parts[1]
            gts=parts[2:]
            if len(gts)!=len(sample_order):
                raise ValueError(f"Variant {vid} has {len(gts)} genotypes but expected {len(sample_order)} based on sample-order file")

            f2_missing=sum(1 for i in idx_f2 if is_missing(gts[i]))
            f3_missing=sum(1 for i in idx_f3 if is_missing(gts[i]))
            f2_total=len(idx_f2)
            f3_total=len(idx_f3)
            f2_nonmissing=f2_total-f2_missing
            f3_nonmissing=f3_total-f3_missing
            total_missing=f2_missing+f3_missing
            total_called=f2_nonmissing+f3_nonmissing

            f2_called_carrier=sum(1 for i in idx_f2 if (not is_missing(gts[i])) and is_carrier(gts[i]))
            f3_called_carrier=sum(1 for i in idx_f3 if (not is_missing(gts[i])) and is_carrier(gts[i]))
            f2_called_noncarrier=sum(1 for i in idx_f2 if (not is_missing(gts[i])) and (not is_carrier(gts[i])))
            f3_called_noncarrier=sum(1 for i in idx_f3 if (not is_missing(gts[i])) and (not is_carrier(gts[i])))
            total_called_carriers=f2_called_carrier+f3_called_carrier

            miss_rate_f2=safe_rate(f2_missing,f2_total)
            miss_rate_f3=safe_rate(f3_missing,f3_total)
            carr_rate_f2=safe_rate(f2_called_carrier,f2_nonmissing)
            carr_rate_f3=safe_rate(f3_called_carrier,f3_nonmissing)

            miss_abs=abs(miss_rate_f2-miss_rate_f3) if miss_rate_f2 is not None and miss_rate_f3 is not None else None
            carr_abs=abs(carr_rate_f2-carr_rate_f3) if carr_rate_f2 is not None and carr_rate_f3 is not None else None

            miss_p=fisher_exact_two_sided(f2_missing,f2_nonmissing,f3_missing,f3_nonmissing)
            carr_p=None
            if total_called_carriers==0:
                carr_p=None
            elif f2_nonmissing>0 and f3_nonmissing>0:
                carr_p=fisher_exact_two_sided(f2_called_carrier,f2_called_noncarrier,f3_called_carrier,f3_called_noncarrier)

            flag,flag_reasons=should_flag(miss_rate_f2,miss_rate_f3,miss_p,carr_p,carr_abs,total_called_carriers)
            reasons=set(flag_reasons.split(',')) if flag_reasons else set()
            if flag==1:
                n_flag+=1
            if 'missing_fisher' in reasons:
                n_missing_fisher+=1
            if 'missing_abs_diff' in reasons:
                n_missing_abs_diff+=1
            if 'high_missing_either_freeze' in reasons:
                n_high_missing_either_freeze+=1
            if 'carrier_fisher' in reasons:
                n_carrier_fisher+=1

            w.writerow([
                vid,vtype,FREEZE2,FREEZE3,
                f2_total,f2_missing,f2_nonmissing,
                f3_total,f3_missing,f3_nonmissing,
                total_missing,total_called,
                miss_rate_f2,miss_rate_f3,miss_abs,miss_p,
                f2_nonmissing,f2_called_carrier,f2_called_noncarrier,
                f3_nonmissing,f3_called_carrier,f3_called_noncarrier,
                total_called_carriers,
                carr_rate_f2,carr_rate_f3,carr_abs,carr_p,
                flag,flag_reasons
            ])
print(f"Total variants processed: {n_total}",file=sys.stderr)
print(f"Total variants flagged: {n_flag}",file=sys.stderr)
print(f"Total variants flagged for missing_fisher: {n_missing_fisher}",file=sys.stderr)
print(f"Total variants flagged for missing_abs_diff: {n_missing_abs_diff}",file=sys.stderr)
print(f"Total variants flagged for high_missing_either_freeze: {n_high_missing_either_freeze}",file=sys.stderr)
print(f"Total variants flagged for carrier_fisher: {n_carrier_fisher}",file=sys.stderr)
