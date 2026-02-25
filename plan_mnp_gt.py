#!/usr/bin/env python3
"""Build MNP GT plan from tested-pairs table (same input/REF logic as build_mnp). Output for manage_gt."""
import argparse
import subprocess
import sys

def parse_id(vid):
    """Parse ID like chr22_16591614_C_T into (chr,pos,ref,alt)."""
    parts=vid.strip().split('_')
    if len(parts)<4:
        return None
    try:
        pos=int(parts[1])
        return (parts[0],pos,parts[2],parts[3])
    except (ValueError,IndexError):
        return None

def chrom_for_ref(chr_,strip_chr=False):
    """Return contig name for ref/plan: strip leading 'chr' if strip_chr else as-is."""
    if strip_chr and (chr_ or '').lower().startswith('chr'):
        return chr_[3:]
    return chr_ or ''

def fetch_ref_region(ref_fasta,chr_,start,end,strip_chr=False):
    """Return reference sequence for chr:start-end (1-based inclusive)."""
    contig=chrom_for_ref(chr_,strip_chr)
    region=f"{contig}:{start}-{end}"
    result=subprocess.run(
        ['samtools','faidx',ref_fasta,region],
        capture_output=True,text=True,check=True
    )
    seq=result.stdout.strip().split('\n')[-1]
    return seq

def load_pairs_and_carriers(input_path):
    """Read tab-separated pairs table (header: pair_id id1 id2 sample ...). Return dict (id1,id2) -> set(sample)."""
    carriers={}
    with open(input_path) as f:
        for i,line in enumerate(f):
            toks=line.strip().split('\t')
            if len(toks)<4:
                continue
            if i==0 and toks[0]=='pair_id':
                continue
            id1,id2,sample=toks[1],toks[2],toks[3]
            if id1 and id2 and sample:
                key=(id1,id2)
                carriers.setdefault(key,set()).add(sample)
    return carriers

def main():
    ap=argparse.ArgumentParser(description='Build MNP GT plan from tested-pairs table (for manage_gt).')
    ap.add_argument('-r','--ref',required=True,help='Reference FASTA (indexed with samtools faidx)')
    ap.add_argument('-i','--input',required=True,help='Tab-separated pairs table: pair_id id1 id2 sample ...')
    ap.add_argument('-o','--output',required=True,help='Output GT plan path (tab: mnp_id chrom pos ref alt id1 id2 carriers)')
    ap.add_argument('--no-chr',action='store_true',help='Strip leading chr for ref contig and plan chrom (e.g. GRCh38 no-chr fasta)')
    args=ap.parse_args()

    pair_carriers=load_pairs_and_carriers(args.input)
    n_adjacent=0
    n_2apart=0
    rows=[]

    for (id1,id2),samples in pair_carriers.items():
        p1=parse_id(id1)
        p2=parse_id(id2)
        if not p1 or not p2:
            continue
        chr1,pos1,ref1,alt1=p1
        chr2,pos2,ref2,alt2=p2
        if chr1!=chr2:
            continue
        if pos1>pos2:
            pos1,pos2=pos2,pos1
            ref1,ref2=ref2,ref1
            alt1,alt2=alt2,alt1

        if pos2==pos1+1:
            ref=ref1+ref2
            alt=alt1+alt2
            pos=pos1
            n_adjacent+=1
        elif pos2==pos1+2:
            seq=fetch_ref_region(args.ref,chr1,pos1,pos2,strip_chr=args.no_chr)
            if len(seq)!=3:
                continue
            if seq[0]!=ref1 or seq[2]!=ref2:
                continue
            ref=seq
            alt=alt1+seq[1]+alt2
            pos=pos1
            n_2apart+=1
        else:
            continue

        mnp_id=f"{chr1}_{pos1}_{ref}_{alt}"
        out_chr=chrom_for_ref(chr1,args.no_chr)
        carriers_str=','.join(sorted(samples))
        rows.append((mnp_id,out_chr,pos,ref,alt,id1,id2,carriers_str))

    with open(args.output,'w') as f:
        f.write('mnp_id\tchrom\tpos\tref\talt\tid1\tid2\tcarriers\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r)+'\n')

    print(f"Wrote {len(rows)} MNPs ({n_adjacent} adjacent + {n_2apart} 2-apart) -> {args.output}",file=sys.stderr)

if __name__=='__main__':
    main()
