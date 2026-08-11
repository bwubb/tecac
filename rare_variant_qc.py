import argparse
import csv
import math
import os
import sys
from itertools import combinations
import pandas as pd

# Defaults for flags and should_flag(). Edit here or pass CLI args.
MISSING_ABS_DIFF_DEFAULT=0.03
HIGH_MISSING_DEFAULT=0.05
MISSING_FISHER_P_DEFAULT=1e-6
CARRIER_FISHER_P_DEFAULT=1e-4
CARRIER_ABS_DIFF_MIN_DEFAULT=0.01

p=argparse.ArgumentParser(description="Controls-only freeze differential missingness and carrier QC per variant. "
                                      "Pairwise across all freezes present (flag if ANY freeze pair trips a threshold).")
p.add_argument("--matrix",required=True,help="Variant x sample GT matrix from bcftools query")
p.add_argument("--sample-order",required=True,help="Sample order file from bcftools query -l")
p.add_argument("--sample-freeze",required=True,help="Sample to freeze map TSV (IID TAB FREEZE)")
p.add_argument("--output",required=True,help="Output TSV path")
p.add_argument("--exclude-ids",required=True,help="Flagged variant IDs, one per line (list is built even if not applied downstream)")
p.add_argument("--min-abs-missing-diff",type=float,default=MISSING_ABS_DIFF_DEFAULT,help="Flag threshold for absolute missingness-rate difference (max over freeze pairs)")
p.add_argument("--high-missing-threshold",type=float,default=HIGH_MISSING_DEFAULT,help="Flag threshold for high missingness in any freeze")
p.add_argument("--carrier-fisher-p-threshold",type=float,default=CARRIER_FISHER_P_DEFAULT,help="Flag threshold for carrier-rate Fisher p-value (min over freeze pairs)")
p.add_argument("--missing-fisher-p-threshold",type=float,default=MISSING_FISHER_P_DEFAULT,help="Flag threshold for missingness Fisher p-value (min over freeze pairs)")
args=p.parse_args()


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
    if n==0:
        return None
    low=max(0,c1-r2)
    high=min(r1,c1)
    p_obs=hypergeom_prob(a,r1,r2,c1,n)
    pv=0.0
    eps=1e-12
    for x in range(low,high+1):
        px=hypergeom_prob(x,r1,r2,c1,n)
        if px<=p_obs+eps:
            pv+=px
    return min(max(pv,0.0),1.0)

def safe_rate(num,den):
    if den<=0:
        return None
    return num/den

def freeze_sort_key(f):
    try:
        return (0,int(f))
    except ValueError:
        return (1,f)


sample_order=pd.read_csv(args.sample_order,sep='\t',header=None,dtype=str)[0].astype(str).tolist()
sample_freeze_df=pd.read_csv(args.sample_freeze,sep='\t',header=None,names=['IID','FREEZE'],dtype=str)
sample_freeze={str(r['IID']):str(r['FREEZE']) for _,r in sample_freeze_df.iterrows()}
missing_samples=[s for s in sample_order if s not in sample_freeze]
if len(missing_samples)>0:
    raise ValueError(f"sample_freeze is missing {len(missing_samples)} sample IDs from sample_order. First 10: {missing_samples[:10]}")

freeze_vec=[sample_freeze[s] for s in sample_order]
freezes=sorted(set(freeze_vec),key=freeze_sort_key)
idx_by_freeze={f:[i for i,fv in enumerate(freeze_vec) if fv==f] for f in freezes}
if len(freezes)<2:
    raise ValueError(f"Need >=2 freeze groups in controls for pairwise QC; found: {freezes}")
freeze_pairs=list(combinations(freezes,2))

print(f"Controls in sample order: {len(sample_order)}",file=sys.stderr)
for f in freezes:
    print(f"Freeze {f} controls: {len(idx_by_freeze[f])}",file=sys.stderr)
print(f"Freeze pairs compared: {['_vs_'.join(p) for p in freeze_pairs]}",file=sys.stderr)


def should_flag(max_miss_diff,max_miss_rate,min_miss_p,min_carr_p,max_carr_diff,total_called_carriers):
    reasons=[]
    if max_miss_diff is not None and max_miss_diff>args.min_abs_missing_diff:
        reasons.append('missing_abs_diff')
    if max_miss_rate is not None and max_miss_rate>args.high_missing_threshold:
        reasons.append('high_missing_any_freeze')
    if min_miss_p is not None and min_miss_p<args.missing_fisher_p_threshold:
        reasons.append('missing_fisher')
    if (min_carr_p is not None and max_carr_diff is not None and min_carr_p<args.carrier_fisher_p_threshold
            and max_carr_diff>=CARRIER_ABS_DIFF_MIN_DEFAULT and total_called_carriers>=3):
        reasons.append('carrier_fisher')
    return (1 if reasons else 0,','.join(reasons))


# Build header: fixed + per-freeze + summary
per_freeze_cols=[]
for f in freezes:
    per_freeze_cols+= [f"f{f}_n_total",f"f{f}_n_missing",f"f{f}_missing_rate",
                       f"f{f}_n_called",f"f{f}_n_carrier",f"f{f}_carrier_rate_called"]
out_cols=(['variant_id','variant_type','freezes']
          +per_freeze_cols
          +['total_missing','total_called','total_called_carriers',
            'abs_missing_rate_diff','missing_fisher_p','missing_worst_pair','max_missing_rate',
            'abs_carrier_rate_diff','carrier_fisher_p','carrier_worst_pair',
            'flag','flag_reasons'])

os.makedirs(os.path.dirname(args.output) or '.',exist_ok=True)
os.makedirs(os.path.dirname(args.exclude_ids) or '.',exist_ok=True)

with open(args.output,'w',newline='') as out,open(args.exclude_ids,'w') as excl:
    w=csv.writer(out,delimiter='\t')
    w.writerow(out_cols)
    n_total=0
    n_flag=0
    counts={'missing_fisher':0,'missing_abs_diff':0,'high_missing_any_freeze':0,'carrier_fisher':0}
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
                raise ValueError(f"Variant {vid} has {len(gts)} genotypes but expected {len(sample_order)}")

            # per-freeze counts
            stat={}
            total_missing=0
            total_called=0
            total_called_carriers=0
            for fz in freezes:
                idx=idx_by_freeze[fz]
                n_tot=len(idx)
                n_miss=sum(1 for i in idx if is_missing(gts[i]))
                n_nonmiss=n_tot-n_miss
                n_carr=sum(1 for i in idx if (not is_missing(gts[i])) and is_carrier(gts[i]))
                n_noncarr=n_nonmiss-n_carr
                stat[fz]={'n_total':n_tot,'n_missing':n_miss,'n_nonmiss':n_nonmiss,
                          'n_carrier':n_carr,'n_noncarrier':n_noncarr,
                          'miss_rate':safe_rate(n_miss,n_tot),
                          'carr_rate':safe_rate(n_carr,n_nonmiss)}
                total_missing+=n_miss
                total_called+=n_nonmiss
                total_called_carriers+=n_carr

            # pairwise missingness + carrier
            max_miss_diff=None
            min_miss_p=None
            missing_worst_pair=''
            max_carr_diff=None
            min_carr_p=None
            carrier_worst_pair=''
            for a,b in freeze_pairs:
                sa,sb=stat[a],stat[b]
                # missingness
                if sa['miss_rate'] is not None and sb['miss_rate'] is not None:
                    d=abs(sa['miss_rate']-sb['miss_rate'])
                    if max_miss_diff is None or d>max_miss_diff:
                        max_miss_diff=d
                mp=fisher_exact_two_sided(sa['n_missing'],sa['n_nonmiss'],sb['n_missing'],sb['n_nonmiss'])
                if mp is not None and (min_miss_p is None or mp<min_miss_p):
                    min_miss_p=mp
                    missing_worst_pair=f"f{a}_vs_f{b}"
                # carrier (only if both have called samples and some carriers exist)
                if sa['n_nonmiss']>0 and sb['n_nonmiss']>0:
                    if sa['carr_rate'] is not None and sb['carr_rate'] is not None:
                        cd=abs(sa['carr_rate']-sb['carr_rate'])
                        if max_carr_diff is None or cd>max_carr_diff:
                            max_carr_diff=cd
                    if (sa['n_carrier']+sb['n_carrier'])>0:
                        cp=fisher_exact_two_sided(sa['n_carrier'],sa['n_noncarrier'],sb['n_carrier'],sb['n_noncarrier'])
                        if cp is not None and (min_carr_p is None or cp<min_carr_p):
                            min_carr_p=cp
                            carrier_worst_pair=f"f{a}_vs_f{b}"

            max_miss_rate=max((stat[fz]['miss_rate'] for fz in freezes if stat[fz]['miss_rate'] is not None),default=None)

            flag,flag_reasons=should_flag(max_miss_diff,max_miss_rate,min_miss_p,min_carr_p,max_carr_diff,total_called_carriers)
            if flag==1:
                n_flag+=1
                excl.write(vid+'\n')
            for r in (flag_reasons.split(',') if flag_reasons else []):
                if r in counts:
                    counts[r]+=1

            row=[vid,vtype,','.join(freezes)]
            for fz in freezes:
                s=stat[fz]
                row+=[s['n_total'],s['n_missing'],s['miss_rate'],
                      s['n_nonmiss'],s['n_carrier'],s['carr_rate']]
            row+=[total_missing,total_called,total_called_carriers,
                  max_miss_diff,min_miss_p,missing_worst_pair,max_miss_rate,
                  max_carr_diff,min_carr_p,carrier_worst_pair,
                  flag,flag_reasons]
            w.writerow(row)

print(f"Total variants processed: {n_total}",file=sys.stderr)
print(f"Total variants flagged: {n_flag}",file=sys.stderr)
for k,v in counts.items():
    print(f"Total variants flagged for {k}: {v}",file=sys.stderr)
