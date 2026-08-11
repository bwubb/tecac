#!/usr/bin/env python3
"""Per-MNP case/control carrier stats + validation targeting, for threshold choice.

"Carrier" = VAF-cis carrier (a row in the PASS file), i.e. "possibly has the MNP"
before CRAM read-check. Case PASS files (test_mnp_sample_info --sample-list case.list)
contain all carrier rows (case + control) for every pair with >=1 case carrier, so
they alone give both counts per case MNP.

Cohort AAF (allele freq, het-approx) = total carriers / (2 * N_total), which is what
decides the REGENIE aaf-bin a variant lands in. N_total defaults to cases+controls;
override with --n-total or --passing-samples. Control carriers are the *download cost*;
cohort AAF is the *include/exclude gate*.

Outputs:
  -o TSV: one row per unique MNP with counts, freqs, cohort AAF, bin, optional Fisher,
          and (if --read-check) per-MNP case/control cis-fraction from read-check.
  stderr: distribution + a TARGETING table showing, per candidate AAF cutoff, the
          UNIQUE control CRAMs to download (deduped; one sample carrying several MNPs
          counts once), and NEW ones after subtracting --already-samples.
  --controls-out / --targets-out (with --aaf-max): write the actual download/target
          lists for the chosen cutoff.
"""
import argparse
import csv
import os
import sys
from statistics import quantiles

try:
    from scipy.stats import fisher_exact
    HAVE_SCIPY=True
except Exception:
    HAVE_SCIPY=False

AAF_BINS=(0.0001,0.001,0.01)


def load_id_list(path):
    keep=set()
    with open(path) as f:
        for line in f:
            s=line.strip().strip('"')
            if not s or s.startswith('#'):
                continue
            parts=s.split()
            keep.add(parts[-1] if len(parts)>=2 else parts[0])
    return keep


def load_id_lists(paths):
    keep=set()
    for p in paths or []:
        if p and os.path.isfile(p):
            keep|=load_id_list(p)
    return keep


def load_pairs(pass_paths):
    """(id1,id2) -> set(samples). De-dupes carrier rows per pair."""
    pairs={}
    for path in pass_paths:
        with open(path,newline='') as f:
            r=csv.DictReader(f,delimiter='\t')
            for row in r:
                id1=(row.get('id1') or '').strip().strip('"')
                id2=(row.get('id2') or '').strip().strip('"')
                sample=(row.get('sample') or '').strip().strip('"')
                if not id1 or not id2 or not sample:
                    continue
                pairs.setdefault((id1,id2),set()).add(sample)
    return pairs


def load_read_check(path,cases,controls):
    """frozenset({id1,id2}) -> dict of case/control checked & strong counts."""
    if not path or not os.path.isfile(path):
        return {}
    agg={}
    with open(path,newline='') as f:
        r=csv.DictReader(f)
        for row in r:
            id1=(row.get('id1') or '').strip().strip('"')
            id2=(row.get('id2') or '').strip().strip('"')
            sample=(row.get('sample') or '').strip().strip('"')
            if not id1 or not id2 or not sample:
                continue
            strong=(row.get('Strong_Cis_Default') or '').strip().strip('"')=='1'
            key=frozenset((id1,id2))
            d=agg.setdefault(key,{'case_checked':0,'case_strong':0,'ctrl_checked':0,'ctrl_strong':0})
            if sample in cases:
                d['case_checked']+=1
                d['case_strong']+=1 if strong else 0
            elif sample in controls:
                d['ctrl_checked']+=1
                d['ctrl_strong']+=1 if strong else 0
    return agg


def pct(sorted_vals,p):
    if not sorted_vals:
        return 0
    if len(sorted_vals)==1:
        return sorted_vals[0]
    qs=quantiles(sorted_vals,n=100,method='inclusive')
    return qs[min(max(int(p)-1,0),len(qs)-1)]


def aaf_bin_label(aaf):
    for b in AAF_BINS:
        if aaf<=b:
            return f"<={b:g}"
    return f">{AAF_BINS[-1]:g}"


def main():
    ap=argparse.ArgumentParser(description='Per-MNP case/control carrier stats + validation targeting.')
    ap.add_argument('pass_files',nargs='+',help='chr*.mnp_sample_info.PASS.txt (case PASS files)')
    ap.add_argument('--case-list',required=True)
    ap.add_argument('--control-list',required=True)
    ap.add_argument('-o','--output',required=True)
    ap.add_argument('--n-total',type=int,default=None,help='Cohort N for AAF denominator (default cases+controls)')
    ap.add_argument('--passing-samples',default=None,help='File of passing IIDs; its line count sets N_total')
    ap.add_argument('--read-check',default=None,help='check_mnp_reads results CSV to join per-MNP cis-fraction')
    ap.add_argument('--already-samples',nargs='*',default=None,help='Files of control IIDs already downloaded/validated (subtracted from NEW)')
    ap.add_argument('--aaf-max',type=float,default=None,help='If set with the *-out options, write lists for this cohort-AAF cutoff')
    ap.add_argument('--controls-out',default=None,help='Write unique NEW control IIDs to download for --aaf-max (minus --already-samples)')
    ap.add_argument('--targets-out',default=None,help='Write target MNP rows (TSV) for --aaf-max')
    ap.add_argument('--assignments-out',default=None,help='Write control read-check assignments (mnp_id sample id1 id2) for ALL control carriers of target MNPs (feeds check_mnp_reads)')
    args=ap.parse_args()

    cases=load_id_list(args.case_list)
    controls=load_id_list(args.control_list)
    n_cases=len(cases)
    n_controls=len(controls)
    if n_cases==0 or n_controls==0:
        sys.exit(f"Empty list: cases={n_cases} controls={n_controls}")

    if args.passing_samples and os.path.isfile(args.passing_samples):
        n_total=len(load_id_list(args.passing_samples))
    elif args.n_total:
        n_total=args.n_total
    else:
        n_total=n_cases+n_controls

    already=load_id_lists(args.already_samples)
    pairs=load_pairs(args.pass_files)
    rc=load_read_check(args.read_check,cases,controls)

    rows=[]
    for (id1,id2),samples in pairs.items():
        case_set={s for s in samples if s in cases}
        ctrl_set={s for s in samples if s in controls}
        n_case=len(case_set)
        n_ctrl=len(ctrl_set)
        n_other=len(samples)-n_case-n_ctrl
        n_carrier=len(samples)
        cohort_aaf=n_carrier/(2*n_total)
        case_freq=n_case/n_cases
        ctrl_freq=n_ctrl/n_controls
        ratio=(case_freq/ctrl_freq) if ctrl_freq>0 else float('inf')
        fp='.'
        orat='.'
        if HAVE_SCIPY:
            try:
                orr,p=fisher_exact([[n_case,n_cases-n_case],[n_ctrl,n_controls-n_ctrl]],alternative='greater')
                fp=f"{p:.3e}"
                orat=f"{orr:.3f}" if orr!=float('inf') else 'inf'
            except Exception:
                pass
        row={
            'mnp_id':f"{id1}__{id2}",'id1':id1,'id2':id2,
            'n_case':n_case,'n_control':n_ctrl,'n_other':n_other,'carrier_total':n_carrier,
            'case_freq':f"{case_freq:.6f}",'control_freq':f"{ctrl_freq:.6f}",
            'cohort_aaf':f"{cohort_aaf:.6g}",'aaf_bin':aaf_bin_label(cohort_aaf),
            'case_or_ctrl_ratio':('inf' if ratio==float('inf') else f"{ratio:.3f}"),
            'fisher_p':fp,'odds_ratio':orat,
        }
        if args.read_check:
            d=rc.get(frozenset((id1,id2)),{})
            cc=d.get('case_checked',0); csg=d.get('case_strong',0)
            tc=d.get('ctrl_checked',0); tsg=d.get('ctrl_strong',0)
            row.update({
                'case_checked':cc,'case_strong':csg,
                'case_cis_frac':(f"{csg/cc:.3f}" if cc else '.'),
                'ctrl_checked':tc,'ctrl_strong':tsg,
                'ctrl_cis_frac':(f"{tsg/tc:.3f}" if tc else '.'),
            })
        row['_ctrl_set']=ctrl_set
        row['_aaf']=cohort_aaf
        rows.append(row)

    rows.sort(key=lambda r:(r['_aaf'],r['n_control']),reverse=True)

    fields=['mnp_id','id1','id2','n_case','n_control','n_other','carrier_total',
            'case_freq','control_freq','cohort_aaf','aaf_bin',
            'case_or_ctrl_ratio','fisher_p','odds_ratio']
    if args.read_check:
        fields+=['case_checked','case_strong','case_cis_frac','ctrl_checked','ctrl_strong','ctrl_cis_frac']
    with open(args.output,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

    # ---------------- summary ----------------
    n_mnp=len(rows)
    ctrl_counts=sorted(r['n_control'] for r in rows)
    print(f"Cohort: {n_cases} cases, {n_controls} controls, N_total={n_total} (AAF denom 2N={2*n_total})",file=sys.stderr)
    print(f"Unique case MNPs: {n_mnp}    Wrote -> {args.output}",file=sys.stderr)
    if already:
        print(f"Already-validated controls supplied: {len(already)}",file=sys.stderr)
    print("",file=sys.stderr)
    print("n_control carriers per MNP: "
          f"min={ctrl_counts[0]} median={ctrl_counts[n_mnp//2]} "
          f"p90={pct(ctrl_counts,90)} p95={pct(ctrl_counts,95)} max={ctrl_counts[-1]}",file=sys.stderr)
    print("",file=sys.stderr)

    # MNPs per aaf bin
    print("MNPs by cohort AAF bin:",file=sys.stderr)
    for b in AAF_BINS:
        nb=sum(1 for r in rows if r['_aaf']<=b)
        print(f"  AAF <= {b:<7g}: {nb:>5} MNPs",file=sys.stderr)
    nb=sum(1 for r in rows if r['_aaf']>AAF_BINS[-1])
    print(f"  AAF >  {AAF_BINS[-1]:<7g}: {nb:>5} MNPs (left as separate SNVs)",file=sys.stderr)
    print("",file=sys.stderr)

    # TARGETING: unique control downloads per cutoff (deduped + minus already)
    print("Validation targeting (full case+control, MNP record only below cutoff):",file=sys.stderr)
    print(f"  {'AAF cutoff':>11} | {'MNPs':>5} | {'uniq ctrl':>9} | {'NEW ctrl':>8}",file=sys.stderr)
    for b in AAF_BINS:
        tgt=[r for r in rows if r['_aaf']<=b]
        uniq=set()
        for r in tgt:
            uniq|=r['_ctrl_set']
        new=uniq-already
        print(f"  {('<= %g'%b):>11} | {len(tgt):>5} | {len(uniq):>9} | {len(new):>8}",file=sys.stderr)

    # optional: write lists for chosen cutoff
    if args.aaf_max is not None and (args.controls_out or args.targets_out or args.assignments_out):
        tgt=[r for r in rows if r['_aaf']<=args.aaf_max]
        uniq=set()
        for r in tgt:
            uniq|=r['_ctrl_set']
        new=sorted(uniq-already)
        if args.controls_out:
            with open(args.controls_out,'w') as f:
                for s in new:
                    f.write(s+'\n')
            print(f"\nWrote {len(new)} NEW control IIDs (AAF<={args.aaf_max:g}) -> {args.controls_out}",file=sys.stderr)
        if args.targets_out:
            with open(args.targets_out,'w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
                w.writeheader()
                w.writerows(tgt)
            print(f"Wrote {len(tgt)} target MNPs (AAF<={args.aaf_max:g}) -> {args.targets_out}",file=sys.stderr)
        if args.assignments_out:
            n_assign=0
            with open(args.assignments_out,'w') as f:
                f.write("mnp_id\tsample\tid1\tid2\n")
                for r in tgt:
                    for s in sorted(r['_ctrl_set']):
                        f.write(f"{r['mnp_id']}\t{s}\t{r['id1']}\t{r['id2']}\n")
                        n_assign+=1
            print(f"Wrote {n_assign} control read-check assignments (AAF<={args.aaf_max:g}) -> {args.assignments_out}",file=sys.stderr)


if __name__=='__main__':
    main()
