#!/usr/bin/env python3
"""Filter PASS TSV to read-check supported MNPs.

Two modes:
  pair-level (default): keep every carrier row of a pair if ANY read-check row
      for that pair is Strong_Cis. (legacy set_cover behavior)
  --per-sample: keep a carrier row only if THAT (id1,id2,sample) was read-checked
      and is Strong_Cis. Use this when you read-checked every case/control carrier
      and only want to rewrite GT for samples you actually validated.
"""
import argparse
import csv
import sys


def load_keep_pairs(read_check_csv,keep_col="Strong_Cis_Default",keep_val="1"):
    """Return set of (id1,id2) pairs (both orders) with at least one Strong_Cis row."""
    keep=set()
    with open(read_check_csv,newline="") as f:
        r=csv.DictReader(f)
        missing=[c for c in ("id1","id2",keep_col) if c not in (r.fieldnames or [])]
        if missing:
            raise RuntimeError(f"Missing columns in {read_check_csv}: {missing}")
        for row in r:
            if (row.get(keep_col) or "").strip()!=keep_val:
                continue
            id1=(row.get("id1") or "").strip()
            id2=(row.get("id2") or "").strip()
            if not id1 or not id2:
                continue
            keep.add((id1,id2))
            keep.add((id2,id1))
    return keep


def load_target_pairs(path):
    """Read a targets TSV with id1,id2 columns; return set of (id1,id2) both orders.
    Restricts MNP creation to these pairs (e.g. cohort AAF <= cutoff)."""
    keep=set()
    with open(path,newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        missing=[c for c in ("id1","id2") if c not in (r.fieldnames or [])]
        if missing:
            raise RuntimeError(f"Missing columns in {path}: {missing}")
        for row in r:
            id1=(row.get("id1") or "").strip()
            id2=(row.get("id2") or "").strip()
            if not id1 or not id2:
                continue
            keep.add((id1,id2))
            keep.add((id2,id1))
    return keep


def load_keep_samples(read_check_csv,keep_col="Strong_Cis_Default",keep_val="1"):
    """Return set of (id1,id2,sample) (both id orders) that are Strong_Cis."""
    keep=set()
    with open(read_check_csv,newline="") as f:
        r=csv.DictReader(f)
        missing=[c for c in ("id1","id2","sample",keep_col) if c not in (r.fieldnames or [])]
        if missing:
            raise RuntimeError(f"Missing columns in {read_check_csv}: {missing}")
        for row in r:
            if (row.get(keep_col) or "").strip()!=keep_val:
                continue
            id1=(row.get("id1") or "").strip()
            id2=(row.get("id2") or "").strip()
            sample=(row.get("sample") or "").strip()
            if not id1 or not id2 or not sample:
                continue
            keep.add((id1,id2,sample))
            keep.add((id2,id1,sample))
    return keep


def filter_pass(pass_in,pass_out,keep,ids_out,per_sample=False,target_pairs=None):
    kept=0
    total=0
    all_ids=set()
    with open(pass_in,newline="") as fin:
        r=csv.DictReader(fin,delimiter="\t")
        fieldnames=r.fieldnames or []
        need=("id1","id2","sample") if per_sample else ("id1","id2")
        missing=[c for c in need if c not in fieldnames]
        if missing:
            raise RuntimeError(f"Missing columns in {pass_in}: {missing}")
        with open(pass_out,"w",newline="") as fout:
            w=csv.DictWriter(fout,fieldnames=fieldnames,delimiter="\t",extrasaction="ignore")
            w.writeheader()
            for row in r:
                total+=1
                id1=(row.get("id1") or "").strip()
                id2=(row.get("id2") or "").strip()
                if target_pairs is not None and (id1,id2) not in target_pairs:
                    continue
                if per_sample:
                    sample=(row.get("sample") or "").strip()
                    hit=(id1,id2,sample) in keep
                else:
                    hit=(id1,id2) in keep
                if hit:
                    w.writerow(row)
                    kept+=1
                    all_ids.add(id1)
                    all_ids.add(id2)
    with open(ids_out,"w") as f:
        for vid in sorted(all_ids):
            f.write(vid+"\n")
    print(f"IDs (one per line for bcftools ID=@file): {len(all_ids)} -> {ids_out}",file=sys.stderr)
    return kept,total


def main():
    ap=argparse.ArgumentParser(description="Filter PASS TSV to read-check supported MNP pairs.")
    ap.add_argument("-r","--read-check",required=True,help="check_mnp_reads results CSV")
    ap.add_argument("-i","--pass-in",required=True,help="Input PASS TSV")
    ap.add_argument("-o","--pass-out",required=True,help="Filtered PASS TSV output")
    ap.add_argument("--ids-out",required=True,help="Every distinct id1 and id2 from kept rows; one ID per line; sorted (bcftools -i ID=@file)")
    ap.add_argument("--keep-col",default="Strong_Cis_Default",help="Column in read-check CSV used for keep decision")
    ap.add_argument("--keep-val",default="1",help="Value in --keep-col to keep")
    ap.add_argument("--per-sample",action="store_true",help="Keep a carrier row only if that (id1,id2,sample) is Strong_Cis (only rewrite GT for validated samples)")
    ap.add_argument("--target-mnps",default=None,help="TSV with id1,id2 columns; restrict MNP creation to these pairs (e.g. cohort AAF <= cutoff)")
    args=ap.parse_args()

    target_pairs=load_target_pairs(args.target_mnps) if args.target_mnps else None
    if target_pairs is not None:
        print(f"Target MNPs (restrict to): {len(target_pairs)//2}",file=sys.stderr)

    if args.per_sample:
        keep=load_keep_samples(args.read_check,args.keep_col,args.keep_val)
        kept,total=filter_pass(args.pass_in,args.pass_out,keep,args.ids_out,per_sample=True,target_pairs=target_pairs)
        print(f"Keep (pair,sample) Strong_Cis: {len(keep)//2}",file=sys.stderr)
    else:
        keep=load_keep_pairs(args.read_check,args.keep_col,args.keep_val)
        kept,total=filter_pass(args.pass_in,args.pass_out,keep,args.ids_out,per_sample=False,target_pairs=target_pairs)
        print(f"Keep-pairs: {len(keep)//2}",file=sys.stderr)
    print(f"Filtered {args.pass_in}: kept {kept}/{total} rows -> {args.pass_out}",file=sys.stderr)


if __name__=="__main__":
    main()
