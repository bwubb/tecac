#!/usr/bin/env python3
# Prebuild slim REGENIE report tables (ADD top genes + contrib; SKAT/SKATO top genes;
# unique gene-level case/control carriers; CHEK2 carrier sample list for PC plots).
import argparse
import glob
import math
import os
import subprocess
import tempfile

import pandas as pd

HTML_THRESH=('0.01','0.001','singleton')
CSV_THRESH=('0.01','0.001','0.0001','singleton','All variants')
MASKS=('M1','M2','M3','M4')
SKAT_TESTS=('ADD-SKAT','ADD-SKATO')
CHEK2_MASKS=('M1',)
CHEK2_THRESH=('0.01','0.001')


def get_args():
    p=argparse.ArgumentParser(description='Prebuild REGENIE report tables')
    p.add_argument('--gene-glob',required=True,help='Glob for step2 gene-based .regenie files')
    p.add_argument('--snplist-glob',required=True,help='Glob for *_masks.snplist files')
    p.add_argument('--mask-stats',required=True,help='mask_variant_stats.tsv (full, may be large)')
    p.add_argument('--out-prefix',required=True,help='Output prefix (no extension)')
    p.add_argument('--top-n',type=int,default=20,help='Top N genes per mask×bin for contrib/CSV')
    p.add_argument('--html-n',type=int,default=10,help='Mark top N for HTML bins')
    p.add_argument('--pheno',default='',help='REGENIE pheno (FID IID STATUS); needed for unique gene carriers')
    p.add_argument('--pfile-template',default='data/preprocess/chr{CHR}.annotation',
                   help='PLINK2 --pfile prefix with {CHR} for unique carrier counts')
    return p.parse_args()


def split_gene_mask_thresh(vid):
    p=str(vid).split('.')
    if len(p)<3:
        return None,None,None
    gene=p[0]
    mask=p[1]
    raw='.'.join(p[2:])
    thresh='All variants' if raw=='all' else raw
    return gene,mask,thresh


def bh_adjust(pvals):
    p=list(pvals)
    n=len(p)
    if n==0:
        return []
    order=sorted(range(n),key=lambda i:p[i])
    adj=[0.0]*n
    prev=1.0
    for i in range(n-1,-1,-1):
        idx=order[i]
        rank=i+1
        val=p[idx]*n/rank
        if val>prev:
            val=prev
        if val>1.0:
            val=1.0
        prev=val
        adj[idx]=val
    return adj


def apply_bh(df):
    out=df.copy()
    out['PVAL_ADJ']=float('nan')
    ok=out['PVAL'].notna()&(out['PVAL']>0)&(out['PVAL']<=1)
    if ok.any():
        out.loc[ok,'PVAL_ADJ']=bh_adjust(out.loc[ok,'PVAL'].tolist())
    return out


def load_gene_based(paths):
    rows=[]
    for fp in paths:
        if not os.path.isfile(fp):
            continue
        df=pd.read_csv(fp,sep=r'\s+',comment='#')
        rows.append(df)
    if not rows:
        return pd.DataFrame()
    g=pd.concat(rows,ignore_index=True)
    for c in ('BETA','SE','LOG10P','A1FREQ','A1FREQ_CASES','A1FREQ_CONTROLS'):
        if c in g.columns:
            g[c]=pd.to_numeric(g[c],errors='coerce')
    g['PVAL']=g['LOG10P'].map(lambda x:(10**(-x)) if pd.notna(x) and math.isfinite(x) else float('nan'))
    parts=g['ID'].map(split_gene_mask_thresh)
    g['GENE']=[t[0] for t in parts]
    g['MASK']=[t[1] for t in parts]
    g['THRESH']=[t[2] for t in parts]
    return g


def add_or_ci(add):
    out=add.copy()
    out['OR']=float('nan')
    out['OR_L95']=float('nan')
    out['OR_U95']=float('nan')
    out['OR']=out['BETA'].map(lambda b:math.exp(b) if pd.notna(b) and math.isfinite(b) else float('nan'))
    out['OR_L95']=(out['BETA']-1.96*out['SE']).map(
        lambda x:math.exp(x) if pd.notna(x) and math.isfinite(x) else float('nan'))
    out['OR_U95']=(out['BETA']+1.96*out['SE']).map(
        lambda x:math.exp(x) if pd.notna(x) and math.isfinite(x) else float('nan'))
    return out


def top_genes_add(add,top_n,html_n):
    out=[]
    for mask in MASKS:
        for thresh in CSV_THRESH:
            sub=add[(add['MASK']==mask)&(add['THRESH']==thresh)].copy()
            if sub.empty:
                continue
            sub=sub.sort_values('LOG10P',ascending=False).head(top_n)
            sub['rank']=range(1,len(sub)+1)
            sub['html_bin']=thresh in HTML_THRESH
            sub['html_top']=sub['html_bin']&(sub['rank']<=html_n)
            out.append(sub)
    if not out:
        return pd.DataFrame()
    cols=['GENE','MASK','THRESH','rank','html_top','BETA','SE','OR','OR_L95','OR_U95',
          'LOG10P','PVAL','PVAL_ADJ','A1FREQ','A1FREQ_CASES','A1FREQ_CONTROLS','N','N_CASES','N_CONTROLS','ID']
    t=pd.concat(out,ignore_index=True)
    return t[[c for c in cols if c in t.columns]]


def top_genes_skat(skat,top_n,html_n):
    out=[]
    for test in SKAT_TESTS:
        for mask in MASKS:
            for thresh in CSV_THRESH:
                sub=skat[(skat['TEST']==test)&(skat['MASK']==mask)&(skat['THRESH']==thresh)].copy()
                if sub.empty:
                    continue
                sub=sub.sort_values('LOG10P',ascending=False).head(top_n)
                sub['rank']=range(1,len(sub)+1)
                sub['html_bin']=thresh in HTML_THRESH
                sub['html_top']=sub['html_bin']&(sub['rank']<=html_n)
                out.append(sub)
    if not out:
        return pd.DataFrame(columns=[
            'GENE','TEST','MASK','THRESH','rank','html_top',
            'LOG10P','PVAL','PVAL_ADJ','A1FREQ_CASES','A1FREQ_CONTROLS','N','N_CASES','N_CONTROLS','ID'
        ])
    cols=['GENE','TEST','MASK','THRESH','rank','html_top',
          'LOG10P','PVAL','PVAL_ADJ','A1FREQ','A1FREQ_CASES','A1FREQ_CONTROLS','N','N_CASES','N_CONTROLS','ID']
    t=pd.concat(out,ignore_index=True)
    return t[[c for c in cols if c in t.columns]]


def needed_keys(top):
    return set(zip(top['GENE'].tolist(),top['MASK'].tolist(),top['THRESH'].tolist()))


def load_snplist_for_keys(paths,keys):
    rows=[]
    for fp in paths:
        if not os.path.isfile(fp):
            continue
        with open(fp) as f:
            for line in f:
                line=line.strip()
                if not line or line.startswith('#'):
                    continue
                parts=line.split('\t')
                if len(parts)<2:
                    continue
                mid=parts[0]
                gene,mask,thresh=split_gene_mask_thresh(mid)
                if gene is None:
                    continue
                if (gene,mask,thresh) not in keys:
                    continue
                for vid in parts[1].split(','):
                    vid=vid.strip()
                    if vid:
                        rows.append({'GENE':gene,'MASK':mask,'THRESH':thresh,'mask_id':mid,'ID':vid})
    if not rows:
        return pd.DataFrame(columns=['GENE','MASK','THRESH','mask_id','ID'])
    return pd.DataFrame(rows)


def load_snplist_for_gene(paths,gene,masks,threshs):
    rows=[]
    prefix=gene+'.'
    mask_set=set(masks)
    thresh_set=set(threshs)
    for fp in paths:
        if not os.path.isfile(fp):
            continue
        with open(fp) as f:
            for line in f:
                line=line.strip()
                if not line or line.startswith('#') or not line.startswith(prefix):
                    continue
                parts=line.split('\t')
                if len(parts)<2:
                    continue
                mid=parts[0]
                g,mask,thresh=split_gene_mask_thresh(mid)
                if g!=gene or mask not in mask_set or thresh not in thresh_set:
                    continue
                for vid in parts[1].split(','):
                    vid=vid.strip()
                    if vid:
                        rows.append({'GENE':g,'MASK':mask,'THRESH':thresh,'ID':vid})
    if not rows:
        return pd.DataFrame(columns=['GENE','MASK','THRESH','ID'])
    return pd.DataFrame(rows)


def id_aliases(vid):
    s=str(vid)
    if s.startswith('chr'):
        return [s,s[3:]]
    return [s,'chr'+s]


def unique_keep_order(ids):
    seen=set()
    out=[]
    for x in ids:
        if x in seen:
            continue
        seen.add(x)
        out.append(x)
    return out


def load_stats_for_ids(path,ids):
    want=set()
    for vid in ids:
        want.update(id_aliases(vid))
    if not want:
        return pd.DataFrame()
    chunks=[]
    for chunk in pd.read_csv(path,sep='\t',chunksize=100000,dtype=str,low_memory=False):
        chunk=chunk[chunk['ID'].isin(want)]
        if len(chunk):
            chunks.append(chunk)
    if not chunks:
        return pd.DataFrame()
    out=pd.concat(chunks,ignore_index=True).drop_duplicates(subset=['ID'],keep='first')
    for c in ('n_case_carrier','n_control_carrier','n_case_called','n_control_called','AC','AN','AF'):
        if c in out.columns:
            out[c]=pd.to_numeric(out[c],errors='coerce')
    return out


def reindex_stats_to_ids(stats,ids):
    # Map mask_stats rows onto snplist/query IDs (chr prefix may differ)
    if stats is None or stats.empty or not len(ids):
        return pd.DataFrame()
    by=stats.drop_duplicates(subset=['ID']).set_index('ID',drop=False)
    rows=[]
    for vid in unique_keep_order(ids):
        hit=None
        for a in id_aliases(vid):
            if a in by.index:
                hit=by.loc[a]
                break
        if hit is None:
            continue
        r=hit.copy()
        r['ID']=vid
        rows.append(r)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).drop_duplicates(subset=['ID'],keep='first')


def consequence_matrix(contrib):
    if contrib.empty or 'Variant.Consequence' not in contrib.columns:
        return pd.DataFrame(columns=['MASK','THRESH','GENE','Variant.Consequence','n'])
    c=contrib[['MASK','THRESH','GENE','Variant.Consequence']].copy()
    c=c[c['Variant.Consequence'].notna()&(c['Variant.Consequence']!='')&(c['Variant.Consequence']!='.')]
    if c.empty:
        return pd.DataFrame(columns=['MASK','THRESH','GENE','Variant.Consequence','n'])
    c['Variant.Consequence']=c['Variant.Consequence'].astype(str).str.split('&')
    c=c.explode('Variant.Consequence')
    c['Variant.Consequence']=c['Variant.Consequence'].str.strip()
    c=c[c['Variant.Consequence']!='']
    if c.empty:
        return pd.DataFrame(columns=['MASK','THRESH','GENE','Variant.Consequence','n'])
    g=c.groupby(['MASK','THRESH','GENE','Variant.Consequence'],as_index=False).size()
    return g.rename(columns={'size':'n'})


def chrom_of(vid):
    c=str(vid).split('_',1)[0]
    return c[3:] if c.startswith('chr') else c


def alt_of(vid):
    # IDs are CHR_POS_REF_ALT
    return str(vid).rsplit('_',1)[-1]


def load_pheno_status(path):
    ph=pd.read_csv(path,sep=r'\s+')
    if 'IID' not in ph.columns or 'STATUS' not in ph.columns:
        raise ValueError('pheno needs IID and STATUS')
    ph['IID']=ph['IID'].astype(str)
    ph['STATUS']=pd.to_numeric(ph['STATUS'],errors='coerce')
    ph=ph[ph['STATUS'].isin([0,1])]
    return dict(zip(ph['IID'],ph['STATUS'].astype(int)))


def plink_export_a(pfile,extract,out_prefix,ids):
    # Force counted allele = ALT from ID (CHR_POS_REF_ALT). Default --export A can
    # count REF and mark nearly every sample as a "carrier".
    allele_fp=out_prefix+'.export_allele.txt'
    with open(allele_fp,'w') as f:
        for vid in ids:
            f.write(f'{vid}\t{alt_of(vid)}\n')
    cmd=['plink2','--pfile',pfile,'--extract',extract,
         '--export','A','--export-allele',allele_fp,'--out',out_prefix]
    subprocess.run(cmd,check=False,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    raw=out_prefix+'.raw'
    if not os.path.isfile(raw):
        return None
    return pd.read_csv(raw,sep=r'\s+',dtype=str,low_memory=False)


def match_raw_cols(raw_cols,want_ids):
    """Map variant ID -> .raw dosage column; only ALT-counted columns (ID_ALT)."""
    want=set(want_ids)
    alt_by={vid:alt_of(vid) for vid in want}
    out={}
    for col in raw_cols:
        if col in ('FID','IID','PAT','MAT','SEX','PHENOTYPE'):
            continue
        if '_' not in col:
            continue
        base,allele=col.rsplit('_',1)
        if base not in want:
            continue
        if allele!=alt_by[base]:
            continue
        if base not in out:
            out[base]=col
    return out


def compute_gene_carriers(mem,status_by_iid,pfile_template,extra_sets=None):
    """
    Unique case/control carriers per (GENE,MASK,THRESH) using plink2 --export A (ALT).
    mem: snplist rows for top genes. extra_sets: optional list of dicts with
    GENE,MASK,THRESH,ID rows (e.g. CHEK2 M1) merged into the dosage export.
    Returns (counts_df, chek2_carrier_iids_df).
    """
    empty_counts=pd.DataFrame(columns=['GENE','MASK','THRESH','n_case_carriers','n_control_carriers'])
    empty_chek2=pd.DataFrame(columns=['IID','STATUS'])
    if mem is None or mem.empty:
        mem=pd.DataFrame(columns=['GENE','MASK','THRESH','ID'])
    work=mem[['GENE','MASK','THRESH','ID']].copy()
    if extra_sets is not None and len(extra_sets):
        work=pd.concat([work,extra_sets[['GENE','MASK','THRESH','ID']]],ignore_index=True)
    work=work.drop_duplicates()
    if work.empty:
        return empty_counts,empty_chek2

    work['chr']=work['ID'].map(chrom_of)
    sets={}
    for (gene,mask,thresh),sub in work.groupby(['GENE','MASK','THRESH']):
        sets[(gene,mask,thresh)]=sub['ID'].tolist()

    iid_to_status=status_by_iid
    carriers={k:set() for k in sets}

    with tempfile.TemporaryDirectory(prefix='gene_carr_') as td:
        for chrom,sub in work.groupby('chr'):
            ids=sorted(set(sub['ID'].tolist()))
            if not ids:
                continue
            pfile=pfile_template.format(CHR=chrom)
            if not os.path.isfile(pfile+'.pgen'):
                print(f'skip chr{chrom}: missing {pfile}.pgen')
                continue
            extract=os.path.join(td,f'chr{chrom}.extract')
            with open(extract,'w') as f:
                f.write('\n'.join(ids)+'\n')
            raw=plink_export_a(pfile,extract,os.path.join(td,f'chr{chrom}'),ids)
            if raw is None:
                print(f'skip chr{chrom}: plink --export A failed')
                continue
            colmap=match_raw_cols(list(raw.columns),ids)
            print(f'chr{chrom}: {len(ids)} ids, {len(colmap)} ALT raw cols matched')
            for key,vids in sets.items():
                cols=[colmap[v] for v in vids if v in colmap]
                if not cols:
                    continue
                dos=raw[['IID']+cols].copy()
                for c in cols:
                    dos[c]=pd.to_numeric(dos[c],errors='coerce').fillna(0)
                hit=dos[cols].max(axis=1)>0
                carriers[key].update(dos.loc[hit,'IID'].astype(str).tolist())

    rows=[]
    for (gene,mask,thresh),iids in carriers.items():
        n_case=sum(1 for i in iids if iid_to_status.get(i)==1)
        n_ctrl=sum(1 for i in iids if iid_to_status.get(i)==0)
        rows.append({
            'GENE':gene,'MASK':mask,'THRESH':thresh,
            'n_case_carriers':n_case,'n_control_carriers':n_ctrl,
        })
    counts=pd.DataFrame(rows) if rows else empty_counts

    # CHEK2: M1 only (union across CHEK2 AAF bins; not singleton)
    chek2_iids=set()
    for (gene,mask,thresh),iids in carriers.items():
        if gene=='CHEK2' and mask=='M1' and thresh in CHEK2_THRESH:
            chek2_iids.update(iids)
    chek2_rows=[{'IID':i,'STATUS':iid_to_status.get(i)} for i in sorted(chek2_iids) if i in iid_to_status]
    chek2_df=pd.DataFrame(chek2_rows) if chek2_rows else empty_chek2
    return counts,chek2_df


def main():
    args=get_args()
    gene_files=sorted(glob.glob(args.gene_glob))
    snp_files=sorted(glob.glob(args.snplist_glob))
    print(f'gene files: {len(gene_files)}')
    print(f'snplist files: {len(snp_files)}')

    g=load_gene_based(gene_files)
    if g.empty:
        print('Error: no gene-based rows')
        return

    add=apply_bh(g[g['TEST']=='ADD'].copy())
    add=add_or_ci(add)
    skat=apply_bh(g[g['TEST'].isin(SKAT_TESTS)].copy())
    print(f'gene rows: {len(g)} ADD: {len(add)} SKAT/SKATO: {len(skat)}')

    top=top_genes_add(add,args.top_n,args.html_n)
    top_skat=top_genes_skat(skat,args.top_n,args.html_n)
    print(f'top ADD gene rows: {len(top)}')
    print(f'top SKAT/SKATO gene rows: {len(top_skat)}')
    keys=needed_keys(top) if len(top) else set()
    print(f'snplist keys: {len(keys)}')

    mem=load_snplist_for_keys(snp_files,keys)
    print(f'snplist variant rows: {len(mem)}')

    stats_raw=load_stats_for_ids(args.mask_stats,mem['ID'].tolist()) if len(mem) else pd.DataFrame()
    stats=reindex_stats_to_ids(stats_raw,mem['ID'].tolist()) if len(mem) else pd.DataFrame()
    print(f'stats matched: {len(stats)} (raw hits {len(stats_raw)})')

    gene_cols=['GENE','MASK','THRESH','rank','html_top','BETA','SE','OR','OR_L95','OR_U95',
               'LOG10P','PVAL','PVAL_ADJ','A1FREQ','A1FREQ_CASES','A1FREQ_CONTROLS']
    gene_cols=[c for c in gene_cols if c in top.columns]
    if len(top) and len(mem):
        contrib=mem.merge(top[gene_cols],on=['GENE','MASK','THRESH'],how='left')
        if len(stats):
            contrib=contrib.merge(stats,on='ID',how='left')
        contrib=contrib.rename(columns={
            'BETA':'gene_BETA','SE':'gene_SE','OR':'gene_OR','OR_L95':'gene_OR_L95','OR_U95':'gene_OR_U95',
            'LOG10P':'gene_LOG10P','PVAL':'gene_PVAL','PVAL_ADJ':'gene_PVAL_ADJ',
            'A1FREQ':'gene_A1FREQ','A1FREQ_CASES':'gene_A1FREQ_CASES','A1FREQ_CONTROLS':'gene_A1FREQ_CONTROLS',
            'rank':'gene_rank'
        })
        contrib=contrib.sort_values(['MASK','THRESH','gene_LOG10P','GENE','ID'],ascending=[True,True,False,True,True])
    else:
        contrib=pd.DataFrame()

    cons=consequence_matrix(contrib)

    # Unique gene-level carriers for every gene×mask×bin we show in the report:
    # all top ADD (CSV bins, incl. singleton), top SKAT, suggestive ADD, CHEK2 M1
    gene_carr=pd.DataFrame(columns=['GENE','MASK','THRESH','n_case_carriers','n_control_carriers'])
    chek2_carr=pd.DataFrame(columns=['IID','STATUS'])
    if args.pheno and os.path.isfile(args.pheno):
        status_by_iid=load_pheno_status(args.pheno)
        carr_keys=set()
        if len(top):
            carr_keys|=needed_keys(top)
        if len(top_skat):
            carr_keys|=needed_keys(top_skat)
        suggestive=add[(add['LOG10P']>3)&(add['PVAL_ADJ']>=0.05)].copy() if len(add) else pd.DataFrame()
        if len(suggestive):
            # top-ish suggestive per mask (same spirit as report)
            sug_bits=[]
            for mask in MASKS:
                sub=suggestive[suggestive['MASK']==mask].sort_values('LOG10P',ascending=False).head(15)
                if len(sub):
                    sug_bits.append(sub)
            if sug_bits:
                sug=pd.concat(sug_bits,ignore_index=True)
                carr_keys|=needed_keys(sug)
        # FDR-significant ADD (risk/protective tables)
        if len(add):
            sig_add=add[add['PVAL_ADJ']<0.05]
            if len(sig_add):
                carr_keys|=needed_keys(sig_add)

        carr_mem=load_snplist_for_keys(snp_files,carr_keys) if carr_keys else pd.DataFrame()
        chek2_extra=load_snplist_for_gene(snp_files,'CHEK2',CHEK2_MASKS,CHEK2_THRESH)
        print(f'gene-carrier snplist keys: {len(carr_keys)}; rows: {len(carr_mem)}')
        print(f'CHEK2 M1 snplist rows (bins {CHEK2_THRESH}): {len(chek2_extra)}')
        gene_carr,chek2_carr=compute_gene_carriers(
            carr_mem,status_by_iid,args.pfile_template,
            extra_sets=chek2_extra if len(chek2_extra) else None,
        )
        print(f'gene carrier count rows: {len(gene_carr)}; CHEK2 carrier samples: {len(chek2_carr)}')
        if len(gene_carr):
            by_th=gene_carr.groupby('THRESH').size().to_dict()
            print(f'  carrier rows by THRESH: {by_th}')
    else:
        print('Skipping unique gene carriers (--pheno not set or missing)')

    out_dir=os.path.dirname(args.out_prefix)
    if out_dir:
        os.makedirs(out_dir,exist_ok=True)
    top_fp=args.out_prefix+'.top_genes_ADD.tsv'
    skat_fp=args.out_prefix+'.top_genes_SKAT.tsv'
    contrib_fp=args.out_prefix+'.variant_contrib.tsv'
    cons_fp=args.out_prefix+'.consequence_matrix.tsv'
    carr_fp=args.out_prefix+'.gene_carriers.tsv'
    chek2_fp=args.out_prefix+'.chek2_carriers.tsv'
    top.to_csv(top_fp,sep='\t',index=False)
    top_skat.to_csv(skat_fp,sep='\t',index=False)
    contrib.to_csv(contrib_fp,sep='\t',index=False)
    cons.to_csv(cons_fp,sep='\t',index=False)
    gene_carr.to_csv(carr_fp,sep='\t',index=False)
    chek2_carr.to_csv(chek2_fp,sep='\t',index=False)
    print(f'wrote {top_fp}')
    print(f'wrote {skat_fp}')
    print(f'wrote {contrib_fp}')
    print(f'wrote {cons_fp}')
    print(f'wrote {carr_fp}')
    print(f'wrote {chek2_fp}')


if __name__=='__main__':
    main()
