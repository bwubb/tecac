import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd

p=argparse.ArgumentParser(description="Per-variant case/control x freeze counts (missing,carrier) from het_miss mnp.gt BCFs.")
p.add_argument("--covariates",required=True,help="Covariates with IID,STATUS,FREEZE (1=control,2=case)")
p.add_argument("--passing-samples",required=True,help="Passing sample IDs one per line")
p.add_argument("--variant-list",default="",help="Variant IDs one per line; if empty use --pathogenic-vus + --genes + --flagged-tsv")
p.add_argument("--pathogenic-vus",default="",help="pathogenic_vus.csv for --genes pull")
p.add_argument("--genes",default="CHEK2,HEATR3",help="Comma-separated genes to pull IDs from pathogenic_vus")
p.add_argument("--flagged-tsv",default="",help="rare_variant_qc flagged TSV; top N variant_id merged when no variant-list")
p.add_argument("--top-flagged",type=int,default=50,help="How many flagged variants to add when building list")
p.add_argument("--bcf-prefix",default="data/bcftools/chr",help="Prefix before CHR")
p.add_argument("--bcf-suffix",default=".qc_filter1.het_miss.mnp.gt.bcf",help="Suffix after CHR")
p.add_argument("--output-ids",default="",help="Optional write merged variant ID list")
p.add_argument("--ids-only",action="store_true",help="Only write --output-ids then exit")
p.add_argument("--output",default="",help="Output stratified summary TSV (required unless --ids-only)")
args=p.parse_args()

FREEZE_KEEP={"2","3"}

def load_cov(path):
    with open(path,'r') as f:
        first=f.readline()
    sep=',' if first.count(',')>0 else r'\s+'
    cov=pd.read_csv(path,sep=sep,engine='python')
    cov.columns=[str(c).strip().upper() for c in cov.columns]
    iid='IID' if 'IID' in cov.columns else cov.columns[1]
    st='STATUS' if 'STATUS' in cov.columns else None
    fr='FREEZE' if 'FREEZE' in cov.columns else None
    if st is None or fr is None:
        raise ValueError("covariates need STATUS and FREEZE")
    out=pd.DataFrame({
        'sample_id':cov[iid].astype(str),
        'STATUS':cov[st].astype(str).str.strip(),
        'FREEZE':cov[fr].astype(str).str.strip()
    })
    return out

def status_cc(s):
    t=str(s).strip()
    if t in ('1','1.0'): return 'control'
    if t in ('2','2.0'): return 'case'
    return 'UNKNOWN'

def read_passing(path):
    with open(path,'r') as f:
        lines=[x.strip() for x in f if x.strip()]
    if lines and lines[0].upper() in ('IID','#IID'):
        lines=lines[1:]
    return set(lines)

def vid_chr(vid):
    p=str(vid).split('_',1)[0].replace('chr','')
    return p

def is_missing(gt):
    g=str(gt).strip()
    return g=='' or '.' in g

def is_carrier(gt):
    if is_missing(gt):
        return False
    a=gt.replace('|','/').split('/')
    return any(x not in {'0',''} for x in a)

def run_bcftools_query(bcf,sample_file,ids_file):
    cmd=['bcftools','query','-S',sample_file,'-i',f'ID=@{ids_file}','-f','%ID\t%TYPE[\t%GT]\n',bcf]
    proc=subprocess.run(cmd,capture_output=True,text=True)
    if proc.returncode!=0:
        raise RuntimeError(proc.stderr or proc.stdout or 'bcftools failed')
    return proc.stdout.splitlines()

def build_variant_ids():
    ids=[]
    if args.variant_list and os.path.isfile(args.variant_list):
        with open(args.variant_list,'r') as f:
            ids=[x.strip() for x in f if x.strip()]
    else:
        if not args.pathogenic_vus or not os.path.isfile(args.pathogenic_vus):
            raise ValueError("need --variant-list or --pathogenic-vus")
        with open(args.pathogenic_vus,'r') as f:
            hdr=f.readline()
        sep=',' if hdr.count(',')>0 else '\t'
        pv=pd.read_csv(args.pathogenic_vus,sep=sep,engine='python')
        pv.columns=[str(c).strip() for c in pv.columns]
        if 'ID' not in pv.columns or 'Gene' not in pv.columns:
            raise ValueError("pathogenic_vus needs ID and Gene columns")
        genes=[g.strip() for g in args.genes.split(',') if g.strip()]
        for g in genes:
            sub=pv[pv['Gene'].astype(str)==g]['ID'].astype(str).unique().tolist()
            ids.extend(sub)
        if args.flagged_tsv and os.path.isfile(args.flagged_tsv):
            fl=pd.read_csv(args.flagged_tsv,sep='\t')
            if 'variant_id' in fl.columns:
                top=fl.sort_values('abs_missing_rate_diff',ascending=False,na_position='last').head(args.top_flagged)
                ids.extend(top['variant_id'].astype(str).tolist())
    seen=set()
    out=[]
    for v in ids:
        v=str(v).strip()
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    if args.output_ids:
        os.makedirs(os.path.dirname(args.output_ids) or '.',exist_ok=True)
        with open(args.output_ids,'w') as f:
            for v in out:
                f.write(v+'\n')
    return out

def main():
    if args.ids_only:
        if not args.output_ids:
            raise ValueError("--ids-only needs --output-ids")
        build_variant_ids()
        print(f"wrote {args.output_ids}",file=sys.stderr)
        return
    if not args.output:
        raise ValueError("need --output unless --ids-only")
    passing=read_passing(args.passing_samples)
    cov=load_cov(args.covariates)
    cov['STATUS']=cov['STATUS'].map(status_cc)
    cov=cov[(cov['sample_id'].isin(passing))&(cov['FREEZE'].isin(FREEZE_KEEP))&(cov['STATUS'].isin({'control','case'}))]
    cov=cov.drop_duplicates(subset=['sample_id'],keep='first')
    if cov.empty:
        raise ValueError("no samples after passing+freeze2/3+case/control filter")
    variant_ids=build_variant_ids()
    if not variant_ids:
        raise ValueError("no variant IDs to process")
    by_chr={}
    for vid in variant_ids:
        c=vid_chr(vid)
        if c not in by_chr:
            by_chr[c]=[]
        by_chr[c].append(vid)
    os.makedirs(os.path.dirname(args.output) or '.',exist_ok=True)
    rows_out=[]
    for chrom,vids in sorted(by_chr.items(),key=lambda x:int(x[0]) if str(x[0]).isdigit() else 99):
        bcf=f"{args.bcf_prefix}{chrom}{args.bcf_suffix}"
        if not os.path.isfile(bcf):
            print(f"skip chr{chrom}: missing {bcf}",file=sys.stderr)
            continue
        proc=subprocess.run(['bcftools','query','-l',bcf],capture_output=True,text=True)
        if proc.returncode!=0:
            raise RuntimeError(proc.stderr)
        bcf_order=[x.strip() for x in proc.stdout.splitlines() if x.strip()]
        sample_order=[s for s in bcf_order if s in set(cov['sample_id'])]
        if not sample_order:
            raise ValueError(f"no overlapping samples for {bcf}")
        cov_idx=cov.set_index('sample_id')
        freeze_vec=[cov_idx.loc[s,'FREEZE'] for s in sample_order]
        status_vec=[cov_idx.loc[s,'STATUS'] for s in sample_order]
        idx={}
        for fr in FREEZE_KEEP:
            for st in ('control','case'):
                key=(fr,st)
                idx[key]=[i for i in range(len(sample_order)) if freeze_vec[i]==fr and status_vec[i]==st]
        with tempfile.NamedTemporaryFile('w',delete=False,suffix='.samples') as sf:
            for s in sample_order:
                sf.write(s+'\n')
            sample_file=sf.name
        with tempfile.NamedTemporaryFile('w',delete=False,suffix='.ids') as idf:
            for v in vids:
                idf.write(v+'\n')
            ids_file=idf.name
        try:
            lines=run_bcftools_query(bcf,sample_file,ids_file)
        finally:
            os.unlink(sample_file)
            os.unlink(ids_file)
        for line in lines:
            if not line.strip():
                continue
            parts=line.rstrip('\n').split('\t')
            if len(parts)<2:
                continue
            vid,vtype=parts[0],parts[1]
            gts=parts[2:]
            if len(gts)!=len(sample_order):
                raise ValueError(f"{vid}: got {len(gts)} genotypes expected {len(sample_order)}")
            for fr in FREEZE_KEEP:
                for st in ('control','case'):
                    key=(fr,st)
                    ii=idx[key]
                    n=len(ii)
                    miss=sum(1 for i in ii if is_missing(gts[i]))
                    called=n-miss
                    carr=sum(1 for i in ii if (not is_missing(gts[i])) and is_carrier(gts[i]))
                    mr=(miss/n) if n else None
                    cr=(carr/called) if called else None
                    rows_out.append({
                        'variant_id':vid,'variant_type':vtype,'freeze':fr,'status':st,
                        'n_total':n,'n_missing':miss,'n_called':called,'n_carrier':carr,
                        'missing_rate':mr,'carrier_rate_called':cr
                    })
    out_df=pd.DataFrame(rows_out)
    out_df.to_csv(args.output,sep='\t',index=False)
    print(f"wrote {len(out_df)} rows to {args.output}",file=sys.stderr)

if __name__=='__main__':
    main()
