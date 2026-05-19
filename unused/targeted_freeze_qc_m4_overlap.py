import argparse
import sys

import pandas as pd

p=argparse.ArgumentParser(description='Overlap freeze-QC flagged with HEATR3 M1/M2 (pathogenic_vus) and HEATR3 M4 synonymous (synonymous.csv). CHEK2 is not M4 and is not summarized here.')
p.add_argument('--flagged',required=True,help='rare_variant_qc flagged TSV')
p.add_argument('--pathogenic-vus',required=True,help='pathogenic_vus.csv (M1/M2 extract)')
p.add_argument('--synonymous',required=True,help='synonymous.csv from preprocess_regenie (same columns as pathogenic_vus)')
p.add_argument('--output',required=True,help='Overlap summary TSV')
args=p.parse_args()

fl=pd.read_csv(args.flagged,sep='\t',engine='python')
if 'variant_id' not in fl.columns:
    raise ValueError('flagged TSV needs variant_id')
flagged=set(fl['variant_id'].astype(str))

def read_id_gene(path,label):
    with open(path,'r',encoding='utf-8',errors='replace') as f:
        hdr=f.readline()
    sep=',' if hdr.count(',')>0 else '\t'
    df=pd.read_csv(path,sep=sep,engine='python')
    df.columns=[str(c).strip().strip('"') for c in df.columns]
    if 'ID' not in df.columns or 'Gene' not in df.columns:
        raise ValueError(f'{label} needs ID and Gene')
    return df

pv=read_id_gene(args.pathogenic_vus,'pathogenic_vus')
syn=read_id_gene(args.synonymous,'synonymous')

heatr3=set(pv[pv['Gene'].astype(str)=='HEATR3']['ID'].astype(str))
heatr3_syn=set(syn[syn['Gene'].astype(str)=='HEATR3']['ID'].astype(str))

def summarize(name,gene_set):
    inter=sorted(flagged & gene_set)
    return {'set_name':name,'n_in_pv':len(gene_set),'n_flagged_overlap':len(inter),'overlap_ids':';'.join(inter[:200])}

rows=[
    summarize('HEATR3_pathogenic_vus',heatr3),
    summarize('HEATR3_M4_synonymous',heatr3_syn),
]
out=pd.DataFrame(rows)
out.to_csv(args.output,sep='\t',index=False)
print(out.to_string(index=False),file=sys.stderr)
print(f'wrote {args.output}',file=sys.stderr)
