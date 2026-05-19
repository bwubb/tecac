#!/usr/bin/env python3
"""Summarize MNP validation metrics across candidate,VAF-pass,and read-check stages."""
import argparse
import csv
import sys


def pair_key(id1,id2):
    a=(id1 or "").strip()
    b=(id2 or "").strip()
    if not a or not b:
        return None
    return (a,b) if a<=b else (b,a)


def load_candidate_pairs(pair_files):
    pairs=set()
    for path in pair_files:
        with open(path,"r") as f:
            for line in f:
                line=line.strip()
                if not line or line.startswith("#"):
                    continue
                toks=line.split("\t")
                if len(toks)<2:
                    continue
                key=pair_key(toks[0],toks[1])
                if key:
                    pairs.add(key)
    return pairs


def load_pass_metrics(pass_files):
    pairs=set()
    rows=0
    for path in pass_files:
        with open(path,newline="") as f:
            r=csv.DictReader(f,delimiter="\t")
            for row in r:
                id1=(row.get("id1") or "").strip()
                id2=(row.get("id2") or "").strip()
                key=pair_key(id1,id2)
                if not key:
                    continue
                rows+=1
                pairs.add(key)
    return {"rows":rows,"pairs":pairs}


def count_tsv_rows(path):
    n=0
    with open(path,newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        for _ in r:
            n+=1
    return n


def load_read_check_metrics(path,strong_col="Strong_Cis_Default"):
    pairs=set()
    strong_pairs=set()
    n_rows=0
    n_strong=0
    n_not_strong=0
    n_cis_only=0
    n_mixed=0
    n_trans_only=0
    n_no_phaseable=0
    n_zero_total_reads=0
    with open(path,newline="") as f:
        r=csv.DictReader(f)
        needed={"id1","id2","Reads_with_Both","Trans_Evidence","Phaseable_Reads","Total_Reads",strong_col}
        missing=needed-set(r.fieldnames or [])
        if missing:
            raise RuntimeError(f"Missing columns in {path}: {sorted(missing)}")
        for row in r:
            n_rows+=1
            id1=(row.get("id1") or "").strip()
            id2=(row.get("id2") or "").strip()
            key=pair_key(id1,id2)
            if key:
                pairs.add(key)
            strong=(row.get(strong_col) or "").strip()=="1"
            if strong:
                n_strong+=1
                if key:
                    strong_pairs.add(key)
            else:
                n_not_strong+=1
            both=int(float((row.get("Reads_with_Both") or "0")))
            trans=int(float((row.get("Trans_Evidence") or "0")))
            phaseable=int(float((row.get("Phaseable_Reads") or "0")))
            total_reads=int(float((row.get("Total_Reads") or "0")))
            if both>0 and trans==0:
                n_cis_only+=1
            elif both>0 and trans>0:
                n_mixed+=1
            elif both==0 and trans>0:
                n_trans_only+=1
            if phaseable==0:
                n_no_phaseable+=1
            if total_reads==0:
                n_zero_total_reads+=1
    return {
        "rows":n_rows,
        "pairs":pairs,
        "strong_pairs":strong_pairs,
        "n_strong":n_strong,
        "n_not_strong":n_not_strong,
        "n_cis_only":n_cis_only,
        "n_mixed":n_mixed,
        "n_trans_only":n_trans_only,
        "n_no_phaseable":n_no_phaseable,
        "n_zero_total_reads":n_zero_total_reads,
    }


def write_metrics(path,metrics):
    with open(path,"w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["metric","value"])
        for k,v in metrics:
            w.writerow([k,v])


def main():
    ap=argparse.ArgumentParser(description="Build MNP validation summary metrics table.")
    ap.add_argument("--candidate-pairs",nargs="+",required=True,help="chr*.annotation.mnp.pairs.txt")
    ap.add_argument("--pass-files",nargs="+",required=True,help="chr*.mnp_sample_info.PASS.txt")
    ap.add_argument("--assignments",required=True,help="mnp_read_check.assignments.tsv")
    ap.add_argument("--read-check",required=True,help="mnp_read_check.results.mapq20.csv")
    ap.add_argument("-o","--output",required=True,help="Output metrics TSV")
    args=ap.parse_args()

    candidate_pairs=load_candidate_pairs(args.candidate_pairs)
    pass_metrics=load_pass_metrics(args.pass_files)
    n_assignments=count_tsv_rows(args.assignments)
    rc=load_read_check_metrics(args.read_check)

    metrics=[
        ("n_mnp_candidates_identified_pairs",len(candidate_pairs)),
        ("n_vaf_pass_pairs",len(pass_metrics["pairs"])),
        ("n_vaf_pass_rows",pass_metrics["rows"]),
        ("n_read_validation_targets_assigned",n_assignments),
        ("n_read_validations_performed_rows",rc["rows"]),
        ("n_unique_pairs_validated",len(rc["pairs"])),
        ("n_strong_cis_rows",rc["n_strong"]),
        ("n_not_strong_cis_rows",rc["n_not_strong"]),
        ("n_strong_cis_unique_pairs",len(rc["strong_pairs"])),
        ("n_cis_only_rows",rc["n_cis_only"]),
        ("n_mixed_cis_trans_rows",rc["n_mixed"]),
        ("n_trans_only_rows",rc["n_trans_only"]),
        ("n_no_phaseable_rows",rc["n_no_phaseable"]),
        ("n_zero_total_reads_rows",rc["n_zero_total_reads"]),
    ]
    write_metrics(args.output,metrics)
    print(f"Wrote metrics -> {args.output}",file=sys.stderr)


if __name__=="__main__":
    main()

