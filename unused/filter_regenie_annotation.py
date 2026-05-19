import argparse
import os

import pandas as pd


def parse_args():
    p=argparse.ArgumentParser(description="Filter regenie annotation/set files by flagged variant IDs.")
    p.add_argument("--flagged-tsv",required=True,help="rare_variant_qc flagged TSV with variant_id column")
    p.add_argument("--annotation",required=True,help="Input regenie.annotation.txt")
    p.add_argument("--annotation-out",required=True,help="Output filtered regenie.annotation.txt")
    p.add_argument("--set-file",required=True,help="Input regenie.set.txt")
    p.add_argument("--set-out",required=True,help="Output filtered regenie.set.txt")
    return p.parse_args()


def load_flagged_ids(flagged_tsv):
    fl=pd.read_csv(flagged_tsv,sep='\t',engine='python')
    if 'variant_id' not in fl.columns:
        raise ValueError("--flagged-tsv needs a variant_id column")
    return set(fl['variant_id'].astype(str).str.strip())


def filter_annotation(in_path,out_path,exclude):
    os.makedirs(os.path.dirname(out_path) or '.',exist_ok=True)
    kept=0
    dropped=0
    with open(in_path,'r',encoding='utf-8',errors='replace') as inf,open(out_path,'w',encoding='utf-8') as outf:
        for line in inf:
            line=line.rstrip('\n')
            if not line.strip():
                continue
            vid=line.split()[0]
            if vid in exclude:
                dropped+=1
                continue
            outf.write(line+'\n')
            kept+=1
    return kept,dropped


def filter_set_file(in_path,out_path,exclude):
    os.makedirs(os.path.dirname(out_path) or '.',exist_ok=True)
    kept_lines=0
    dropped_lines=0
    with open(in_path,'r',encoding='utf-8',errors='replace') as inf,open(out_path,'w',encoding='utf-8') as outf:
        for raw in inf:
            line=raw.strip()
            if not line:
                continue
            parts=line.split()
            if len(parts)<4:
                continue
            gene,chrom,pos=parts[0],parts[1],parts[2]
            ids=[x for x in parts[3].split(',') if x]
            ids=[x for x in ids if x not in exclude]
            if not ids:
                dropped_lines+=1
                continue
            outf.write(f"{gene} {chrom} {pos} {','.join(ids)}\n")
            kept_lines+=1
    return kept_lines,dropped_lines


def main():
    args=parse_args()
    exclude=load_flagged_ids(args.flagged_tsv)
    if not exclude:
        raise ValueError("no variant IDs found in --flagged-tsv")
    ann_kept,ann_dropped=filter_annotation(args.annotation,args.annotation_out,exclude)
    set_kept,set_dropped=filter_set_file(args.set_file,args.set_out,exclude)
    print(f"annotation kept={ann_kept} dropped={ann_dropped} wrote {args.annotation_out}")
    print(f"set kept_lines={set_kept} dropped_lines={set_dropped} wrote {args.set_out}")


if __name__=='__main__':
    main()
