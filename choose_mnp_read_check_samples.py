#!/usr/bin/env python3
"""Greedy set cover on mnp_sample_info.PASS -> minimal sample list + targets TSV for check_mnp_reads.py."""
import argparse
import csv
import sys


def load_pass_files(pass_paths):
    """
    Read test_mnp_sample_info PASS TSVs (pair_id id1 id2 sample ...).
    One MNP = unique (id1,id2). Carriers = distinct samples with a PASS row for that pair.
    Returns (mnp_keys, carrier_sets) aligned by index.
    """
    mnp_keys=[]
    carriers_by_key={}
    seen_key=set()
    for path in pass_paths:
        with open(path,newline='') as f:
            r=csv.DictReader(f,delimiter='\t')
            if r.fieldnames:
                names={fn.lstrip('\ufeff') for fn in r.fieldnames}
                missing={'id1','id2','sample'}-names
                if missing:
                    print(f"Warning: {path} missing {missing}; got {r.fieldnames}",file=sys.stderr)
            for row in r:
                id1=(row.get('id1') or '').strip()
                id2=(row.get('id2') or '').strip()
                sample=(row.get('sample') or '').strip()
                if not id1 or not id2 or not sample:
                    continue
                key=(id1,id2)
                if key not in seen_key:
                    seen_key.add(key)
                    mnp_keys.append(key)
                    carriers_by_key[key]=set()
                carriers_by_key[key].add(sample)
    carrier_sets=[carriers_by_key[k] for k in mnp_keys]
    return mnp_keys,carrier_sets


def greedy_set_cover(n_mnp,sample_to_mnps):
    """Ordered samples, each once, union of sample_to_mnps[s] covers 0..n_mnp-1."""
    uncovered=set(range(n_mnp))
    remaining=set(sample_to_mnps)
    chosen=[]
    while uncovered:
        best_s=None
        best_hit=0
        for s in remaining:
            hit=len(sample_to_mnps[s]&uncovered)
            if hit>best_hit:
                best_hit=hit
                best_s=s
        if best_s is None or best_hit==0:
            bad=next(iter(uncovered))
            raise RuntimeError(
                f"Greedy set cover stalled, {len(uncovered)} uncovered (e.g. index {bad})."
            )
        chosen.append(best_s)
        uncovered-=sample_to_mnps[best_s]
        remaining.discard(best_s)
    return chosen


def main():
    ap=argparse.ArgumentParser(
        description='Greedy minimal sample set + per-MNP assignment for check_mnp_reads.py'
    )
    ap.add_argument('pass_files',nargs='+',help='chr*.mnp_sample_info.PASS.txt (tab from test_mnp_sample_info.py)')
    ap.add_argument('-o','--targets',required=True,dest='targets',help='Minimal TSV: mnp_id, sample, id1, id2 (for check_mnp_reads.py --targets-tsv)')
    ap.add_argument('-s','--samples-out',required=True,help='One sample per line (minimal BAM set)')
    ap.add_argument('--also-list-redundant',default=None,help='Write all distinct PASS samples for comparison')
    args=ap.parse_args()

    mnp_keys,carrier_sets=load_pass_files(args.pass_files)
    n=len(mnp_keys)
    if n==0:
        print('No PASS rows / MNPs loaded.',file=sys.stderr)
        with open(args.targets,'w',newline='') as f:
            f.write('mnp_id\tsample\tid1\tid2\n')
        open(args.samples_out,'w').close()
        return

    sample_to_mnps={}
    all_carriers=set()
    for i,carriers in enumerate(carrier_sets):
        all_carriers.update(carriers)
        for s in carriers:
            sample_to_mnps.setdefault(s,set()).add(i)

    chosen_order=greedy_set_cover(n,sample_to_mnps)
    chosen_set=set(chosen_order)

    target_fields=['mnp_id','sample','id1','id2']
    out_rows=[]
    for i,key in enumerate(mnp_keys):
        id1,id2=key
        carriers=carrier_sets[i]
        pick=sorted(carriers&chosen_set)
        if not pick:
            raise RuntimeError(f"MNP {i} id1={id1!r} id2={id2!r} no carrier in chosen set (bug).")
        sample=pick[0]
        mnp_id=id1+'__'+id2
        out_rows.append({'mnp_id':mnp_id,'sample':sample,'id1':id1,'id2':id2})

    with open(args.targets,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=target_fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(out_rows)

    with open(args.samples_out,'w') as f:
        for s in sorted(chosen_set):
            f.write(s+'\n')

    if args.also_list_redundant:
        with open(args.also_list_redundant,'w') as out:
            for s in sorted(all_carriers):
                out.write(s+'\n')

    print(
        f"MNPs: {n} | distinct carriers: {len(all_carriers)} | greedy BAM set: {len(chosen_set)}",
        file=sys.stderr,
    )
    print(f"Wrote targets -> {args.targets}",file=sys.stderr)
    print(f"Wrote sample list -> {args.samples_out}",file=sys.stderr)


if __name__=='__main__':
    main()
