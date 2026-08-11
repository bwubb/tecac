#!/usr/bin/env python3
"""Every (MNP, case) PASS row -> read-check targets for check_mnp_reads.py."""
import argparse
import csv
import sys


def load_sample_list(path):
    keep=set()
    with open(path) as f:
        for line in f:
            sample=line.strip()
            if not sample or sample.startswith('#'):
                continue
            keep.add(sample)
    return keep


def load_pass_rows(pass_paths):
    rows=[]
    for path in pass_paths:
        with open(path,newline='') as f:
            r=csv.DictReader(f,delimiter='\t')
            for row in r:
                id1=(row.get('id1') or '').strip()
                id2=(row.get('id2') or '').strip()
                sample=(row.get('sample') or '').strip()
                if not id1 or not id2 or not sample:
                    continue
                rows.append(row)
    return rows


def main():
    ap=argparse.ArgumentParser(description='Read-check targets: every case carrier per MNP.')
    ap.add_argument('pass_files',nargs='+',help='chr*.mnp_sample_info.PASS.txt')
    ap.add_argument('-o','--targets',required=True,dest='targets')
    ap.add_argument('-s','--samples-out',required=True)
    ap.add_argument('--case-list',required=True,help='One case ID per line')
    args=ap.parse_args()

    cases=load_sample_list(args.case_list)
    target_fields=['mnp_id','sample','id1','id2','VAF1','VAF2','VAF_diff','Cis_likely']
    out_rows=[]
    samples=set()

    for row in load_pass_rows(args.pass_files):
        sample=(row.get('sample') or '').strip()
        if sample not in cases:
            continue
        id1=(row.get('id1') or '').strip()
        id2=(row.get('id2') or '').strip()
        mnp_id=id1+'__'+id2
        out_rows.append({
            'mnp_id':mnp_id,
            'sample':sample,
            'id1':id1,
            'id2':id2,
            'VAF1':row.get('VAF1','.'),
            'VAF2':row.get('VAF2','.'),
            'VAF_diff':row.get('VAF_diff','.'),
            'Cis_likely':row.get('Cis_likely','.'),
        })
        samples.add(sample)

    with open(args.targets,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=target_fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(out_rows)

    with open(args.samples_out,'w') as f:
        for s in sorted(samples):
            f.write(s+'\n')

    n_mnp=len({(r['id1'],r['id2']) for r in out_rows})
    print(
        f"MNPs: {n_mnp} | read-check rows: {len(out_rows)} | case BAMs: {len(samples)}",
        file=sys.stderr,
    )
    print(f"Wrote targets -> {args.targets}",file=sys.stderr)
    print(f"Wrote sample list -> {args.samples_out}",file=sys.stderr)


if __name__=='__main__':
    main()
