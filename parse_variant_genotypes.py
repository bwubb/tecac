import pandas as pd
import argparse

p=argparse.ArgumentParser(description="Summarize pathogenic/VUS genotype matrices by STATUS and FREEZE.")
p.add_argument("-c","--covariates",required=True,help="Covariates file (expects IID plus STATUS/FREEZE columns)")
p.add_argument("-p","--pathogenic-gt",required=True,help="Pathogenic genotype matrix TSV")
p.add_argument("-v","--vus-gt",required=True,help="VUS genotype matrix TSV")
p.add_argument("-P","--pathogenic-output",required=True,help="Output summary TSV for pathogenic set")
p.add_argument("-V","--vus-output",required=True,help="Output summary TSV for VUS set")
args=p.parse_args()

# PLINK-style pheno: STATUS 1 = control, 2 = case (output uses these words, not 1/2).
def status_to_case_control(s):
    t=str(s).strip()
    if t in ('1','1.0'): return 'control'
    if t in ('2','2.0'): return 'case'
    if t in ('nan','','UNKNOWN'): return 'UNKNOWN'
    return t

def load_covariates(path):
    with open(path,'r') as f:
        first=f.readline()
    sep=',' if first.count(',')>0 else r'\s+'
    cov=pd.read_csv(path,sep=sep,engine='python')
    cov.columns=[str(c).strip().upper() for c in cov.columns]
    iid_col='IID' if 'IID' in cov.columns else cov.columns[1]
    status_col='STATUS' if 'STATUS' in cov.columns else None
    freeze_col='FREEZE' if 'FREEZE' in cov.columns else None
    out=cov[[iid_col]].copy()
    out=out.rename(columns={iid_col:'sample_id'})
    out['sample_id']=out['sample_id'].astype(str)
    out['STATUS']=cov[status_col].astype(str) if status_col else 'UNKNOWN'
    out['FREEZE']=cov[freeze_col].astype(str) if freeze_col else 'UNKNOWN'
    out['STATUS']=out['STATUS'].replace({'nan':'UNKNOWN','':'UNKNOWN'})
    out['FREEZE']=out['FREEZE'].replace({'nan':'UNKNOWN','':'UNKNOWN'})
    out['STATUS']=out['STATUS'].map(status_to_case_control)
    return out

def gt_class(s):
    s=str(s).strip()
    if s=='' or '.' in s:
        return 'MISSING'
    a=s.replace('|','/').split('/')
    if any(x=='.' or x=='' for x in a):
        return 'MISSING'
    if all(x=='0' for x in a):
        return 'REF'
    if '0' in a and any(x!='0' for x in a):
        return 'HET'
    if len(set(a))==1:
        return 'HOM'
    return 'HET'

def summarize_matrix(gt_path,out_path,variant_set,cov):
    g=pd.read_csv(gt_path,sep='\t',dtype=str)
    if g.shape[1]<3:
        pd.DataFrame(columns=['variant_set','variant_id','variant_type','scope','status','freeze','n_total','n_called','n_missing','n_ref','n_het','n_hom','call_rate']).to_csv(out_path,sep='\t',index=False)
        return
    id_col=g.columns[0]
    type_col=g.columns[1]
    sample_cols=list(g.columns[2:])
    long=g.melt(id_vars=[id_col,type_col],value_vars=sample_cols,var_name='sample_id',value_name='GT')
    long=long.rename(columns={id_col:'variant_id',type_col:'variant_type'})
    long['sample_id']=long['sample_id'].astype(str)
    long=long.merge(cov,on='sample_id',how='left')
    long['STATUS']=long['STATUS'].fillna('UNKNOWN')
    long['FREEZE']=long['FREEZE'].fillna('UNKNOWN')
    long['gt_class']=long['GT'].map(gt_class)

    ct=pd.crosstab([long['variant_id'],long['variant_type']],long['gt_class']).reset_index()
    for c in ['REF','HET','HOM','MISSING']:
        if c not in ct.columns:
            ct[c]=0
    ct['scope']='overall'
    ct['status']='ALL'
    ct['freeze']='ALL'

    cs=pd.crosstab([long['variant_id'],long['variant_type'],long['STATUS']],long['gt_class']).reset_index()
    for c in ['REF','HET','HOM','MISSING']:
        if c not in cs.columns:
            cs[c]=0
    cs=cs.rename(columns={'STATUS':'status'})
    cs['scope']='status'
    cs['freeze']='ALL'

    cf=pd.crosstab([long['variant_id'],long['variant_type'],long['FREEZE']],long['gt_class']).reset_index()
    for c in ['REF','HET','HOM','MISSING']:
        if c not in cf.columns:
            cf[c]=0
    cf=cf.rename(columns={'FREEZE':'freeze'})
    cf['scope']='freeze'
    cf['status']='ALL'

    csf=pd.crosstab([long['variant_id'],long['variant_type'],long['STATUS'],long['FREEZE']],long['gt_class']).reset_index()
    for c in ['REF','HET','HOM','MISSING']:
        if c not in csf.columns:
            csf[c]=0
    csf=csf.rename(columns={'STATUS':'status','FREEZE':'freeze'})
    csf['scope']='status_freeze'

    out=pd.concat([
        ct[['variant_id','variant_type','scope','status','freeze','REF','HET','HOM','MISSING']],
        cs[['variant_id','variant_type','scope','status','freeze','REF','HET','HOM','MISSING']],
        cf[['variant_id','variant_type','scope','status','freeze','REF','HET','HOM','MISSING']],
        csf[['variant_id','variant_type','scope','status','freeze','REF','HET','HOM','MISSING']]
    ],ignore_index=True)

    out=out.rename(columns={'REF':'n_ref','HET':'n_het','HOM':'n_hom','MISSING':'n_missing'})
    out['n_total']=out['n_ref']+out['n_het']+out['n_hom']+out['n_missing']
    out['n_called']=out['n_total']-out['n_missing']
    out['call_rate']=(out['n_called']/out['n_total']).fillna(0)
    out.insert(0,'variant_set',variant_set)
    out=out[['variant_set','variant_id','variant_type','scope','status','freeze','n_total','n_called','n_missing','n_ref','n_het','n_hom','call_rate']]
    out.to_csv(out_path,sep='\t',index=False)

cov=load_covariates(args.covariates)
summarize_matrix(args.pathogenic_gt,args.pathogenic_output,'pathogenic',cov)
summarize_matrix(args.vus_gt,args.vus_output,'vus',cov)
