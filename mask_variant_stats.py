#!/usr/bin/env python3
"""Full mask-site carrier counts + VEP fields for REGENIE (LoF1/2 + synonymous).

Output can be large (~hundreds of thousands of rows). The report must subset by ID
on demand — do not load this whole TSV into R.
"""
import argparse
import csv
import os
import subprocess
import tempfile

import pandas as pd

VEP_COLS=['ID','Gene','Variant.LoF_level','Variant.Class','Variant.Consequence',
          'HGVSc','HGVSp','ClinVar.SIG','AutoGVP']


def get_args():
    p=argparse.ArgumentParser(description='Full mask variant carrier counts + VEP fields')
    p.add_argument('--vep-list',required=True,help='List of VEP report CSVs (one path per line)')
    p.add_argument('--pheno',required=True,help='REGENIE pheno: FID IID STATUS (0=ctrl,1=case)')
    p.add_argument('--pfile-template',default='data/preprocess/chr{CHR}.annotation',
                   help='PLINK2 --pfile prefix with {CHR}')
    p.add_argument('--chroms',default=','.join(str(i) for i in range(1,23)))
    p.add_argument('-o','--output',required=True,help='Output TSV')
    return p.parse_args()


def is_synonymous(value):
    if pd.isna(value):
        return False
    t=str(value).strip().lower()
    return 'synonymous_variant' in t or 'stop_retained_variant' in t


def load_vep(vep_files):
    # LoF1/2 + synonymous (full mask site set); one row per ID
    rows=[]
    for fp in vep_files:
        if not os.path.isfile(fp):
            continue
        df=pd.read_csv(fp,usecols=lambda c:c in VEP_COLS,dtype=str,low_memory=False)
        if df.empty:
            continue
        df['Variant.LoF_level']=pd.to_numeric(df['Variant.LoF_level'],errors='coerce')
        keep=(df['Variant.LoF_level'].isin([1,2]))|(df['Variant.Consequence'].map(is_synonymous))
        rows.append(df.loc[keep,VEP_COLS])
    if not rows:
        return pd.DataFrame(columns=VEP_COLS)
    out=pd.concat(rows,ignore_index=True)
    out['ID']=out['ID'].astype(str).str.strip()
    out['Gene']=out['Gene'].astype(str).str.strip()
    out=out[(out['ID']!='')&(out['Gene']!='')&(out['ID']!='.')]
    out=out.sort_values(['ID','Variant.LoF_level'],na_position='last')
    return out.drop_duplicates(subset=['ID'],keep='first')


def write_keep(path,iids):
    with open(path,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['#FID','IID'])
        for iid in iids:
            w.writerow([iid,iid])


def parse_gcount(path):
    # ID -> (carriers, called); names from plink2 --geno-counts .gcount
    out={}
    if not os.path.isfile(path):
        return out
    with open(path) as f:
        header=f.readline().strip().lstrip('#').split('\t')
        idx={h:i for i,h in enumerate(header)}
        id_i=idx.get('ID',0)
        het_i=idx.get('HET_REF_ALT_CTS',idx.get('HET_REF_ALT_CT'))
        two_i=idx.get('TWO_ALT_CTS',idx.get('TWO_ALT_CT'))
        hap_i=idx.get('HAP_ALT_CTS',idx.get('HAP_ALT_CT'))
        href_i=idx.get('HOM_REF_CTS',idx.get('HOM_REF_CT'))

        def num(p,i):
            if i is None or i>=len(p) or p[i] in ('.',''):
                return 0
            try:
                return int(float(p[i]))
            except ValueError:
                return 0

        for line in f:
            p=line.rstrip('\n').split('\t')
            if len(p)<=id_i:
                continue
            carr=num(p,het_i)+num(p,two_i)+num(p,hap_i)
            called=num(p,href_i)+num(p,het_i)+num(p,two_i)+num(p,hap_i)
            out[p[id_i]]=(carr,called)
    return out


def plink_gcount(pfile,extract,keep,out_prefix):
    cmd=['plink2','--pfile',pfile,'--extract',extract,'--keep',keep,
         '--geno-counts','--out',out_prefix]
    subprocess.run(cmd,check=False,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    return parse_gcount(out_prefix+'.gcount')


def chrom_of(vid):
    c=vid.split('_',1)[0]
    return c[3:] if c.startswith('chr') else c


def main():
    args=get_args()
    with open(args.vep_list) as f:
        vep_files=[ln.strip() for ln in f if ln.strip()]
    annot=load_vep(vep_files)
    if annot.empty:
        print('Error: no mask variants in VEP files')
        return
    print(f'mask sites: {len(annot)}')

    pheno=pd.read_csv(args.pheno,sep=r'\s+')
    if 'IID' not in pheno.columns or 'STATUS' not in pheno.columns:
        print('Error: pheno needs IID and STATUS')
        return
    pheno['STATUS']=pd.to_numeric(pheno['STATUS'],errors='coerce')
    cases=pheno.loc[pheno['STATUS']==1,'IID'].astype(str).tolist()
    controls=pheno.loc[pheno['STATUS']==0,'IID'].astype(str).tolist()
    if not cases or not controls:
        print(f'Error: need cases and controls (n_case={len(cases)} n_ctrl={len(controls)})')
        return

    annot['chr']=annot['ID'].map(chrom_of)
    chroms=[c.strip() for c in args.chroms.split(',') if c.strip()]
    case_carr={}
    ctrl_carr={}
    case_called={}
    ctrl_called={}

    with tempfile.TemporaryDirectory(prefix='mask_var_') as td:
        case_keep=os.path.join(td,'cases.keep')
        ctrl_keep=os.path.join(td,'controls.keep')
        write_keep(case_keep,cases)
        write_keep(ctrl_keep,controls)

        for chrom in chroms:
            ids=annot.loc[annot['chr']==chrom,'ID'].tolist()
            if not ids:
                continue
            pfile=args.pfile_template.format(CHR=chrom)
            if not os.path.isfile(pfile+'.pgen'):
                print(f'skip chr{chrom}: missing {pfile}.pgen')
                continue
            extract=os.path.join(td,f'chr{chrom}.extract')
            with open(extract,'w') as f:
                f.write('\n'.join(ids)+'\n')
            gc=plink_gcount(pfile,extract,case_keep,os.path.join(td,f'chr{chrom}.cases'))
            gt=plink_gcount(pfile,extract,ctrl_keep,os.path.join(td,f'chr{chrom}.controls'))
            for vid,(carr,called) in gc.items():
                case_carr[vid]=carr
                case_called[vid]=called
            for vid,(carr,called) in gt.items():
                ctrl_carr[vid]=carr
                ctrl_called[vid]=called
            print(f'chr{chrom}: {len(ids)} ids case={len(gc)} ctrl={len(gt)}')

    annot['n_case_carrier']=annot['ID'].map(lambda x:case_carr.get(x,0)).astype(int)
    annot['n_control_carrier']=annot['ID'].map(lambda x:ctrl_carr.get(x,0)).astype(int)
    annot['n_case_called']=annot['ID'].map(lambda x:case_called.get(x,0)).astype(int)
    annot['n_control_called']=annot['ID'].map(lambda x:ctrl_called.get(x,0)).astype(int)
    annot['AC']=annot['n_case_carrier']+annot['n_control_carrier']
    annot['AN']=2*(annot['n_case_called']+annot['n_control_called'])
    annot['AF']=annot.apply(lambda r:(r['AC']/r['AN']) if r['AN']>0 else float('nan'),axis=1)
    annot=annot.drop(columns=['chr'])

    out_dir=os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir,exist_ok=True)
    annot.to_csv(args.output,sep='\t',index=False)
    print(f'wrote {args.output} ({len(annot)} variants)')


if __name__=='__main__':
    main()
