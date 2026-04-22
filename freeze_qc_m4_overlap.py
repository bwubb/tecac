import argparse
import sys

import pandas as pd

p=argparse.ArgumentParser(description="Overlap freeze-QC flagged variants with CHEK2/HEATR3 and HEATR3 synonymous proxy.")
p.add_argument("--flagged",required=True,help="rare_variant_qc flagged TSV")
p.add_argument("--pathogenic-vus",required=True,help="pathogenic_vus.csv")
p.add_argument("--output",required=True,help="Overlap summary TSV")
args=p.parse_args()

fl=pd.read_csv(args.flagged,sep='\t',engine='python')
if 'variant_id' not in fl.columns:
    raise ValueError("flagged TSV needs variant_id")
flagged=set(fl['variant_id'].astype(str))

with open(args.pathogenic_vus,'r') as f:
    hdr=f.readline()
sep=',' if hdr.count(',')>0 else '\t'
pv=pd.read_csv(args.pathogenic_vus,sep=sep,engine='python')
pv.columns=[str(c).strip() for c in pv.columns]
if 'ID' not in pv.columns or 'Gene' not in pv.columns:
    raise ValueError("pathogenic_vus needs ID and Gene")

chek2=set(pv[pv['Gene'].astype(str)=='CHEK2']['ID'].astype(str))
heatr3=set(pv[pv['Gene'].astype(str)=='HEATR3']['ID'].astype(str))
heatr3_syn=set()
cons_col=None
for c in pv.columns:
    if c.replace('.','').lower()=='variantconsequence' or c=='Variant.Consequence':
        cons_col=c
        break
if cons_col:
    h=pv[pv['Gene'].astype(str)=='HEATR3']
    m=h[h[cons_col].astype(str).str.contains('synonymous',case=False,na=False)]
    heatr3_syn=set(m['ID'].astype(str))

def summarize(name, gene_set):
    inter=sorted(flagged & gene_set)
    return {'set_name':name,'n_in_pv':len(gene_set),'n_flagged_overlap':len(inter),'overlap_ids':';'.join(inter[:200])}

rows=[
    summarize('CHEK2',chek2),
    summarize('HEATR3_all',heatr3),
    summarize('HEATR3_synonymous_proxy',heatr3_syn),
]
out=pd.DataFrame(rows)
out.to_csv(args.output,sep='\t',index=False)
print(out.to_string(index=False),file=sys.stderr)
print(f"wrote {args.output}",file=sys.stderr)
