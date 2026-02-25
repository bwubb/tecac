#!/usr/bin/env python3
"""Apply MNP GT plan to cohort VCF: set id1/id2 to 0/0 in carriers, add MNP record with 0/1 carriers, 0/0 others, ./. when either SNV missing."""
import argparse
import sys
import vcfpy

PROGRESS_INTERVAL=1000
READ_PROGRESS_INTERVAL=10000

def norm_chrom(c):
    if (c or '').lower().startswith('chr'):
        return c[3:]
    return c or ''

def load_plan(plan_path):
    """Return (id_to_pair, pair_to_info). id_to_pair: id -> (id1,id2). pair_to_info: (id1,id2) -> {mnp_id,chrom,pos,ref,alt,carriers}."""
    id_to_pair={}
    pair_to_info={}
    with open(plan_path) as f:
        next(f)
        for line in f:
            toks=line.strip().split('\t')
            if len(toks)<8:
                continue
            mnp_id,chrom,pos_s,ref,alt,id1,id2,carriers_s=toks[0],toks[1],toks[2],toks[3],toks[4],toks[5],toks[6],toks[7]
            key=(id1,id2)
            id_to_pair[id1]=key
            id_to_pair[id2]=key
            pair_to_info[key]={
                'mnp_id':mnp_id,'chrom':chrom,'pos':int(pos_s),'ref':ref,'alt':alt,
                'carriers':set(c for c in carriers_s.split(',') if c)
            }
    return id_to_pair,pair_to_info

def is_missing_gt(gt):
    return not gt or gt in ('.','./.','.|.')

def main():
    ap=argparse.ArgumentParser(description='Apply MNP GT plan: set pair GTs to 0/0 in carriers, add MNP record.')
    ap.add_argument('-p','--plan',required=True,help='GT plan file (from plan_mnp_gt.py)')
    ap.add_argument('-i','--input',required=True,help='Input VCF')
    ap.add_argument('-o','--output',required=True,help='Output VCF')
    args=ap.parse_args()

    id_to_pair,pair_to_info=load_plan(args.plan)
    n_mnps=len(pair_to_info)
    print(f"Plan: {n_mnps} MNPs",file=sys.stderr)
    reader=vcfpy.Reader.from_path(args.input)
    writer=vcfpy.Writer.from_path(args.output,reader.header)
    samples=reader.header.samples.names
    buffer=[]
    id1_gt_store={}
    insert_after={}
    cur_chrom=None
    n_written=0
    n_read=0

    def write_rec(rec):
        writer.write_record(rec)
        nonlocal n_written
        n_written+=1
        if n_written==1 or n_written%PROGRESS_INTERVAL==0:
            print(f"  ... {n_written} written",file=sys.stderr)

    for record in reader:
        n_read+=1
        if n_read%READ_PROGRESS_INTERVAL==0:
            print(f"  ... {n_read} read",file=sys.stderr)
        chrom_norm=norm_chrom(record.CHROM)
        if chrom_norm!=cur_chrom:
            if buffer:
                buffer.sort(key=lambda x:(x[0],len(x[1].REF)))
                for _pos,rec in buffer:
                    write_rec(rec)
                buffer.clear()
                id1_gt_store.clear()
                insert_after.clear()
            cur_chrom=chrom_norm
            print(f"Processing {record.CHROM}...",file=sys.stderr)

        ids=list(record.ID) if record.ID else []
        vid=ids[0] if ids else None
        pair=id_to_pair.get(vid) if vid else None
        if not pair:
            if not buffer:
                write_rec(record)
            else:
                buffer.append((record.POS,record))
            continue

        id1,id2=pair
        info=pair_to_info[pair]
        carriers=info['carriers']

        if vid==id1:
            for call in record.calls:
                if call.sample in carriers:
                    call.set_genotype('0/0')
            id1_gt_store[pair]={c.sample:(c.data.get('GT') or './.') for c in record.calls}
            buffer.append((record.POS,record))
            insert_after[pair]=len(buffer)-1
            continue

        if vid==id2:
            for call in record.calls:
                if call.sample in carriers:
                    call.set_genotype('0/0')
            idx=insert_after.get(pair)
            if idx is None:
                if not buffer:
                    write_rec(record)
                else:
                    buffer.append((record.POS,record))
                id1_gt_store.pop(pair,None)
                continue

            id1_gts=id1_gt_store.get(pair,{})
            mnp_gt_by_sample={}
            for call in record.calls:
                s=call.sample
                gt2=call.data.get('GT') or './.'
                gt1=id1_gts.get(s,'./.')
                if is_missing_gt(gt1) or is_missing_gt(gt2):
                    mnp_gt_by_sample[s]='./.'
                elif s in carriers:
                    mnp_gt_by_sample[s]='0/1'
                else:
                    mnp_gt_by_sample[s]='0/0'
            for s in samples:
                if s not in mnp_gt_by_sample:
                    mnp_gt_by_sample[s]='./.'

            mnp_calls=[vcfpy.Call(s,{'GT':mnp_gt_by_sample[s]},site=None) for s in samples]
            alt_rec=vcfpy.Substitution(type_='MNV',value=info['alt'])
            mnp_record=vcfpy.Record(
                CHROM=record.CHROM,
                POS=info['pos'],
                ID=[info['mnp_id']],
                REF=info['ref'],
                ALT=[alt_rec],
                QUAL=None,
                FILTER=['PASS'],
                INFO={},
                FORMAT=['GT'],
                calls=mnp_calls
            )
            buffer.insert(idx+1,(info['pos'],mnp_record))
            buffer.append((record.POS,record))

            def is_id1(rec):
                ids=list(rec.ID) if rec.ID else []
                vid=ids[0] if ids else None
                p=id_to_pair.get(vid) if vid else None
                return p and vid==p[0]

            block=[buffer[idx],buffer[idx+1]]
            for i in range(idx+2,len(buffer)-1):
                if not is_id1(buffer[i][1]):
                    block.append(buffer[i])
            block.append(buffer[-1])
            block.sort(key=lambda x:(x[0],len(x[1].REF)))
            for _pos,rec in block:
                write_rec(rec)
            emitted_idx={idx,idx+1,len(buffer)-1}
            for i in range(idx+2,len(buffer)-1):
                if not is_id1(buffer[i][1]):
                    emitted_idx.add(i)
            new_buffer=[x for i,x in enumerate(buffer) if i not in emitted_idx]
            buffer[:]=new_buffer
            for k in list(insert_after):
                if k==pair:
                    del insert_after[k]
                else:
                    old=insert_after[k]
                    insert_after[k]=old-sum(1 for j in emitted_idx if j<old)
            id1_gt_store.pop(pair,None)
            continue

        buffer.append((record.POS,record))

    if buffer:
        buffer.sort(key=lambda x:(x[0],len(x[1].REF)))
        for _pos,rec in buffer:
            write_rec(rec)
    print(f"Done. {n_written} records written.",file=sys.stderr)
    reader.close()
    writer.stream.flush()
    writer.close()

if __name__=='__main__':
    main()
