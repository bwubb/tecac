import csv
import argparse
import sys

def _parse_id_pos(vid):
    """From ID chr22_15528747_G_A return (chr,pos) or None."""
    parts=vid.strip().split('_')
    if len(parts)<2:
        return None
    try:
        return (parts[0],int(parts[1]))
    except ValueError:
        return None

def load_pairs(pairs_path,chrom_filter=None):
    """Load pairs from .pairs.txt (id1\\tid2 per line) or find_mnp_variants CSV. Returns list of (pair_id,id1,id2,chr,pos1,pos2)."""
    pairs=[]
    with open(pairs_path,'r') as f:
        first=f.readline()
    toks=first.strip().split('\t')
    if len(toks)>=2 and toks[0].strip()!='Chr' and 'Chr' not in first:
        with open(pairs_path) as f:
            pair_id=0
            for line in f:
                line=line.strip()
                if not line or line.startswith('#'):
                    continue
                parts=line.split('\t')
                if len(parts)<2:
                    continue
                id1,id2=parts[0].strip(),parts[1].strip()
                if not id1 or not id2:
                    continue
                p1=_parse_id_pos(id1)
                p2=_parse_id_pos(id2)
                if not p1 or not p2:
                    continue
                chr1,pos1=p1
                chr2,pos2=p2
                if chr1!=chr2:
                    continue
                if chrom_filter and chr1!=chrom_filter:
                    continue
                pair_id+=1
                pairs.append((pair_id,id1,id2,chr1,pos1,pos2))
        return pairs
    with open(pairs_path,'r') as f:
        r=csv.DictReader(f)
        rows=list(r)
    pair_id=0
    i=0
    while i+1<len(rows):
        r1,r2=rows[i],rows[i+1]
        chr1=(r1.get('Chr') or '').strip()
        chr2=(r2.get('Chr') or '').strip()
        id1=(r1.get('ID') or '').strip()
        id2=(r2.get('ID') or '').strip()
        pos1_str=(r1.get('Start') or '').strip()
        pos2_str=(r2.get('Start') or '').strip()
        if not chr1 or not id1 or not chr2 or not id2:
            i+=1
            continue
        if chr1!=chr2:
            i+=1
            continue
        try:
            pos1,pos2=int(pos1_str),int(pos2_str)
        except ValueError:
            i+=1
            continue
        if chrom_filter and chr1!=chrom_filter:
            i+=2
            continue
        pair_id+=1
        pairs.append((pair_id,id1,id2,chr1,pos1,pos2))
        i+=2
    return pairs

def load_sample_info(sample_info_path):
    """Parse space-separated ID SAMPLE GT AD DP VAF. Returns (calls_by_(id,sample), chrom_seen)."""
    calls={}
    chrom_seen=None
    with open(sample_info_path) as f:
        for line in f:
            toks=line.strip().split()
            if len(toks)<6:
                continue
            vid,sample,gt,ad,dp,vaf_str=toks[0],toks[1],toks[2],toks[3],toks[4],toks[5]
            if chrom_seen is None and '_' in vid:
                chrom_seen=vid.split('_')[0]
            try:
                vaf=float(vaf_str) if vaf_str and vaf_str!='.' else None
            except ValueError:
                vaf=None
            key=(vid,sample)
            calls[key]={'gt':gt,'ad':ad,'dp':dp,'vaf':vaf}
    return calls,chrom_seen

def main():
    ap=argparse.ArgumentParser(description='Test MNP pairs against cohort sample info: which samples have both variants and cis-like VAF.')
    ap.add_argument('-p','--pairs',required=True,help='MNP pairs CSV from find_mnp_variants (consecutive rows = pair)')
    ap.add_argument('-i','--input',required=True,help='Sample info file: space-separated ID SAMPLE GT AD DP VAF')
    ap.add_argument('-o','--pass-out',dest='pass_out',required=True,help='Output: PASS samples (both variants, VAF close)')
    ap.add_argument('-f','--fail-out',dest='fail_out',required=True,help='Output: FAIL samples (missing variant or VAF not close)')
    ap.add_argument('--vaf-diff',type=float,default=0.10,help='Max |VAF1-VAF2| to call cis (default 0.10)')
    args=ap.parse_args()

    calls,chrom_seen=load_sample_info(args.input)
    if not calls:
        print('No sample info rows.',file=sys.stderr)
        with open(args.pass_out,'w') as f:
            f.write('')
        with open(args.fail_out,'w') as f:
            f.write('')
        return

    pairs=load_pairs(args.pairs,chrom_filter=chrom_seen)
    pass_rows=[]
    fail_rows=[]
    pass_fields=['pair_id','id1','id2','sample','GT1','AD1','DP1','VAF1','GT2','AD2','DP2','VAF2','VAF_diff','Cis_likely']
    fail_fields=['pair_id','id1','id2','sample','reason','GT1','AD1','DP1','VAF1','GT2','AD2','DP2','VAF2','VAF_diff','Cis_likely']

    for pair_id,id1,id2,chr_,pos1,pos2 in pairs:
        samples_id1={s for (vid,s) in calls if vid==id1}
        samples_id2={s for (vid,s) in calls if vid==id2}
        both=samples_id1&samples_id2
        only1=samples_id1-samples_id2
        only2=samples_id2-samples_id1

        for sample in both:
            c1=calls[(id1,sample)]
            c2=calls[(id2,sample)]
            vaf1,vaf2=c1.get('vaf'),c2.get('vaf')
            vaf_diff=None
            if vaf1 is not None and vaf2 is not None:
                vaf_diff=abs(vaf1-vaf2)
            cis='No'
            if vaf_diff is not None and vaf_diff<=args.vaf_diff:
                cis='Yes'
            row={
                'pair_id':pair_id,'id1':id1,'id2':id2,'sample':sample,
                'GT1':c1.get('gt','.'),'AD1':c1.get('ad','.'),'DP1':c1.get('dp','.'),
                'VAF1':f"{vaf1:.4f}" if vaf1 is not None else '.',
                'GT2':c2.get('gt','.'),'AD2':c2.get('ad','.'),'DP2':c2.get('dp','.'),
                'VAF2':f"{vaf2:.4f}" if vaf2 is not None else '.',
                'VAF_diff':f"{vaf_diff:.4f}" if vaf_diff is not None else '.',
                'Cis_likely':cis,
            }
            if cis=='Yes':
                pass_rows.append(row)
            else:
                row['reason']='vaf_not_close'
                fail_rows.append(row)

        for sample in only1:
            c1=calls[(id1,sample)]
            fail_rows.append({
                'pair_id':pair_id,'id1':id1,'id2':id2,'sample':sample,'reason':'missing_id2',
                'GT1':c1.get('gt','.'),'AD1':c1.get('ad','.'),'DP1':c1.get('dp','.'),
                'VAF1':f"{c1.get('vaf'):.4f}" if c1.get('vaf') is not None else '.',
                'GT2':'.','AD2':'.','DP2':'.','VAF2':'.','VAF_diff':'.','Cis_likely':'.',
            })

        for sample in only2:
            c2=calls[(id2,sample)]
            fail_rows.append({
                'pair_id':pair_id,'id1':id1,'id2':id2,'sample':sample,'reason':'missing_id1',
                'GT1':'.','AD1':'.','DP1':'.','VAF1':'.',
                'GT2':c2.get('gt','.'),'AD2':c2.get('ad','.'),'DP2':c2.get('dp','.'),
                'VAF2':f"{c2.get('vaf'):.4f}" if c2.get('vaf') is not None else '.',
                'VAF_diff':'.','Cis_likely':'.',
            })

    with open(args.pass_out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=pass_fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(pass_rows)

    with open(args.fail_out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fail_fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(fail_rows)

    print(f"PASS: {len(pass_rows)} -> {args.pass_out}",file=sys.stderr)
    print(f"FAIL: {len(fail_rows)} -> {args.fail_out}",file=sys.stderr)

if __name__=='__main__':
    main()
