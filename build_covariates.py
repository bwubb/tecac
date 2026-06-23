#!/usr/bin/env python3
"""Build covariates.txt (FID IID FREEZE STATUS) for ExWAS / REGENIE pipelines.

FREEZE: integer parsed from metadata strings like "Freeze 4".
STATUS: 1 = control (in --controls), 2 = case (everyone else in --samples).
FID and IID are both set to the PMBB sample ID.
"""
import argparse
import csv
import re
import sys

FREEZE_RE=re.compile(r'Freeze\s*(\d+)',re.IGNORECASE)


def load_id_list(path):
    ids=[]
    with open(path) as f:
        for line in f:
            s=line.strip()
            if not s or s.startswith('#'):
                continue
            parts=s.split()
            ids.append(parts[-1] if len(parts)>=2 else parts[0])
    return ids


def parse_freeze_from_row(row):
    for cell in row:
        m=FREEZE_RE.search(str(cell))
        if m:
            return int(m.group(1))
    return None


def load_metadata(metadata_paths):
    """Map PMBB IID -> freeze int from one or more PARCC CSV files."""
    freeze_by_iid={}
    for path in metadata_paths:
        with open(path,newline='') as f:
            for row in csv.reader(f):
                if not row:
                    continue
                iid=(row[0] or '').strip()
                if not iid.startswith('PMBB'):
                    continue
                fz=parse_freeze_from_row(row)
                if fz is None:
                    continue
                if iid in freeze_by_iid and freeze_by_iid[iid]!=fz:
                    print(f'Warning: {iid} freeze {freeze_by_iid[iid]} -> {fz} ({path})',file=sys.stderr)
                freeze_by_iid[iid]=fz
    return freeze_by_iid


def load_legacy_covariates(path):
    """FID IID FREEZE STATUS from prior ExWAS covariates.txt."""
    out={}
    with open(path) as f:
        for i,line in enumerate(f):
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            parts=line.split()
            if i==0 and parts[0].upper() in ('FID','#FID'):
                continue
            if len(parts)<4:
                continue
            fid,iid,freeze,status=parts[0],parts[1],parts[2],parts[3]
            out[iid]={'freeze':int(freeze),'status':int(status),'fid':fid}
    return out


def main():
    ap=argparse.ArgumentParser(description='Build covariates.txt for ExWAS/REGENIE.')
    ap.add_argument('--samples',required=True,help='Sample list (one IID per line, or FID IID)')
    ap.add_argument('--controls',required=True,help='Control list (IIDs present -> STATUS 1)')
    ap.add_argument('--metadata',nargs='+',required=True,help='PARCC metadata CSV(s); freeze parsed from "Freeze N"')
    ap.add_argument('--legacy-covariates',default=None,help='Optional prior covariates.txt for FREEZE fallback')
    ap.add_argument('-o','--output',default='covariates.txt',help='Output path')
    ap.add_argument('--default-freeze',type=int,default=None,help='Freeze if not in metadata or legacy (else skip sample)')
    args=ap.parse_args()

    samples=load_id_list(args.samples)
    controls=set(load_id_list(args.controls))
    freeze_by_iid=load_metadata(args.metadata)
    legacy=load_legacy_covariates(args.legacy_covariates) if args.legacy_covariates else {}

    rows=[]
    missing_freeze=[]
    for iid in samples:
        fid=iid
        if iid in freeze_by_iid:
            freeze=freeze_by_iid[iid]
        elif iid in legacy:
            freeze=legacy[iid]['freeze']
        elif args.default_freeze is not None:
            freeze=args.default_freeze
        else:
            missing_freeze.append(iid)
            continue
        status=1 if iid in controls else 2
        rows.append((fid,iid,freeze,status))

    rows.sort(key=lambda r:r[1])
    with open(args.output,'w',newline='') as f:
        f.write('FID\tIID\tFREEZE\tSTATUS\n')
        for fid,iid,freeze,status in rows:
            f.write(f'{fid}\t{iid}\t{freeze}\t{status}\n')

    n_ctrl=sum(1 for *_,s in rows if s==1)
    n_case=sum(1 for *_,s in rows if s==2)
    by_freeze={}
    for *_,fz,s in rows:
        by_freeze.setdefault(fz,{'ctrl':0,'case':0})
        by_freeze[fz]['ctrl' if s==1 else 'case']+=1

    print(f'Wrote {len(rows)} samples -> {args.output}',file=sys.stderr)
    print(f'  controls (STATUS=1): {n_ctrl}',file=sys.stderr)
    print(f'  cases    (STATUS=2): {n_case}',file=sys.stderr)
    for fz in sorted(by_freeze):
        c=by_freeze[fz]
        print(f'  FREEZE {fz}: {c["ctrl"]} controls, {c["case"]} cases',file=sys.stderr)
    if missing_freeze:
        print(f'Missing FREEZE for {len(missing_freeze)} sample(s); not written:',file=sys.stderr)
        for iid in missing_freeze[:20]:
            print(f'  {iid}',file=sys.stderr)
        if len(missing_freeze)>20:
            print(f'  ... and {len(missing_freeze)-20} more',file=sys.stderr)
        sys.exit(1)


if __name__=='__main__':
    main()
