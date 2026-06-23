#!/usr/bin/env python3
"""Build read-check targets from mnp_sample_info.PASS for check_mnp_reads.py."""
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


def build_set_cover_targets(mnp_keys,carrier_sets):
    n=len(mnp_keys)
    sample_to_mnps={}
    all_carriers=set()
    for i,carriers in enumerate(carrier_sets):
        all_carriers.update(carriers)
        for s in carriers:
            sample_to_mnps.setdefault(s,set()).add(i)
    chosen_order=greedy_set_cover(n,sample_to_mnps)
    chosen_set=set(chosen_order)
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
    return out_rows,sorted(chosen_set),all_carriers


def build_all_case_targets(pass_rows,case_list_path):
  cases=load_sample_list(case_list_path)
  out_rows=[]
  samples=set()
  for row in pass_rows:
      sample=(row.get('sample') or '').strip()
      if sample not in cases:
          continue
      id1=(row.get('id1') or '').strip()
      id2=(row.get('id2') or '').strip()
      if not id1 or not id2:
          continue
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
  return out_rows,sorted(samples)


def main():
    ap=argparse.ArgumentParser(
        description='Build read-check targets for check_mnp_reads.py'
    )
    ap.add_argument('pass_files',nargs='+',help='chr*.mnp_sample_info.PASS.txt (tab from test_mnp_sample_info.py)')
    ap.add_argument('-o','--targets',required=True,dest='targets',help='TSV: mnp_id, sample, id1, id2 (for check_mnp_reads.py --targets-tsv)')
    ap.add_argument('-s','--samples-out',required=True,help='One sample per line (BAMs to open)')
    ap.add_argument('--mode',choices=('set_cover','all_cases'),default='all_cases',
        help='set_cover: one carrier per MNP, minimal BAM set; all_cases: every case carrier per MNP (default)')
    ap.add_argument('--case-list',default=None,help='Required for all_cases: one case ID per line')
    ap.add_argument('--also-list-redundant',default=None,help='Write all distinct PASS carriers for comparison')
    args=ap.parse_args()

    target_fields=['mnp_id','sample','id1','id2','VAF1','VAF2','VAF_diff','Cis_likely']

    if args.mode=='set_cover':
        mnp_keys,carrier_sets=load_pass_files(args.pass_files)
        n=len(mnp_keys)
        if n==0:
            print('No PASS rows / MNPs loaded.',file=sys.stderr)
            with open(args.targets,'w',newline='') as f:
                f.write('mnp_id\tsample\tid1\tid2\n')
            open(args.samples_out,'w').close()
            return
        out_rows,chosen_samples,all_carriers=build_set_cover_targets(mnp_keys,carrier_sets)
        print(
            f"MODE set_cover | MNPs: {n} | distinct carriers: {len(all_carriers)} | greedy BAM set: {len(chosen_samples)}",
            file=sys.stderr,
        )
    else:
        if not args.case_list:
            ap.error('--case-list is required for --mode all_cases')
        pass_rows=load_pass_rows(args.pass_files)
        if not pass_rows:
            print('No PASS rows loaded.',file=sys.stderr)
            with open(args.targets,'w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=target_fields,delimiter='\t',extrasaction='ignore')
                w.writeheader()
            open(args.samples_out,'w').close()
            return
        out_rows,chosen_samples=build_all_case_targets(pass_rows,args.case_list)
        all_carriers={row['sample'] for row in out_rows}
        n_mnp=len({(r['id1'],r['id2']) for r in out_rows})
        print(
            f"MODE all_cases | MNPs with case carriers: {n_mnp} | read-check rows: {len(out_rows)} | case BAMs: {len(chosen_samples)}",
            file=sys.stderr,
        )

    with open(args.targets,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=target_fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(out_rows)

    with open(args.samples_out,'w') as f:
        for s in chosen_samples:
            f.write(s+'\n')

    if args.also_list_redundant:
        with open(args.also_list_redundant,'w') as out:
            for s in sorted(all_carriers):
                out.write(s+'\n')

    print(f"Wrote targets -> {args.targets}",file=sys.stderr)
    print(f"Wrote sample list -> {args.samples_out}",file=sys.stderr)


if __name__=='__main__':
    main()
