#!/usr/bin/env python3
"""Filter PASS TSV to pairs with strong read-check cis evidence."""
import argparse
import csv
import sys


def load_keep_pairs(read_check_csv,keep_col="Strong_Cis_Default",keep_val="1",pair_rule="any"):
    """Return set of (id1,id2) pairs to keep. pair_rule any: one strong row; all: every row per pair strong."""
    rows_by_pair={}
    with open(read_check_csv,newline="") as f:
        r=csv.DictReader(f)
        missing=[c for c in ("id1","id2",keep_col) if c not in (r.fieldnames or [])]
        if missing:
            raise RuntimeError(f"Missing columns in {read_check_csv}: {missing}")
        for row in r:
            id1=(row.get("id1") or "").strip()
            id2=(row.get("id2") or "").strip()
            if not id1 or not id2:
                continue
            key=(id1,id2) if id1<=id2 else (id2,id1)
            strong=(row.get(keep_col) or "").strip()==keep_val
            rows_by_pair.setdefault(key,[]).append(strong)
    keep=set()
    for key,flags in rows_by_pair.items():
        if pair_rule=="all":
            if flags and all(flags):
                keep.add(key)
                keep.add((key[1],key[0]))
        else:
            if any(flags):
                keep.add(key)
                keep.add((key[1],key[0]))
    return keep


def filter_pass(pass_in,pass_out,keep_pairs,ids_out):
    kept=0
    total=0
    all_ids=set()
    with open(pass_in,newline="") as fin:
        r=csv.DictReader(fin,delimiter="\t")
        fieldnames=r.fieldnames or []
        missing=[c for c in ("id1","id2") if c not in fieldnames]
        if missing:
            raise RuntimeError(f"Missing columns in {pass_in}: {missing}")
        with open(pass_out,"w",newline="") as fout:
            w=csv.DictWriter(fout,fieldnames=fieldnames,delimiter="\t",extrasaction="ignore")
            w.writeheader()
            for row in r:
                total+=1
                id1=(row.get("id1") or "").strip()
                id2=(row.get("id2") or "").strip()
                if (id1,id2) in keep_pairs:
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
    ap.add_argument("--pair-rule",choices=("any","all"),default="any",
        help="any: keep pair if any read-check row is strong (default); all: every case read-check must be strong")
    args=ap.parse_args()

    keep_pairs=load_keep_pairs(args.read_check,args.keep_col,args.keep_val,args.pair_rule)
    kept,total=filter_pass(args.pass_in,args.pass_out,keep_pairs,args.ids_out)
    print(f"Keep-pairs: {len(keep_pairs)//2}",file=sys.stderr)
    print(f"Filtered {args.pass_in}: kept {kept}/{total} rows -> {args.pass_out}",file=sys.stderr)


if __name__=="__main__":
    main()

