#!/usr/bin/env python3
"""
Check whether two SNVs fall on the same reads using samtools view + SAM parsing only
(no pysam — avoids OpenSSL/FIPS issues on some clusters).

Phasing summary for downstream MNP review/reporting:
- Cis_Evidence / Reads_with_Both: count of reads carrying both SNV alts (same-read support).
- Trans_Evidence: reads supporting one alt but not both (SNV1+SNV2-2*Both).
- Phaseable_Reads: Cis_Evidence+Trans_Evidence (reads informative for cis vs trans).
- Cis_Fraction_Phaseable: Cis_Evidence/Phaseable_Reads (0 if no phaseable reads).
- Strong_Cis_Default: "1" when Reads_with_Both>=3 and Phaseable_Reads>=5 and
  Cis_Fraction_Phaseable>=0.90, else "0".

Default intent: keep all evidence in output, then threshold by Strong_Cis_Default
or custom cutoffs in downstream analysis/reporting.
"""
import argparse
import csv
import os
import re
import subprocess
import sys

# SAM FLAG bits
FLAG_UNMAP=0x4
FLAG_SECONDARY=0x100
FLAG_SUPPLEMENTARY=0x800
PROGRESS_EVERY_N=50


def parse_variant_id(vid):
    """VCF-style ID: chr22_15528747_G_A -> (chrom, pos, ref, alt)."""
    if not vid or not str(vid).strip():
        return None
    parts=str(vid).strip().split("_")
    if len(parts)<4:
        return None
    try:
        chrom=parts[0]
        pos=int(parts[1])
        if len(parts)==4:
            ref,alt=parts[2],parts[3]
        else:
            ref="_".join(parts[2:-1])
            alt=parts[-1]
        return chrom,pos,ref,alt
    except (ValueError,IndexError):
        return None


def load_bam_table(path):
    """sample<TAB>bam_or_cram_path per line."""
    bams={}
    with open(path,"r") as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"):
                continue
            toks=line.split("\t")
            if len(toks)>=2:
                bams[toks[0].strip()]=toks[1].strip()
    return bams


def parse_cigar_string(cigar_str):
    """Return list of (length, op_char) for CIGAR; empty if missing or '*'."""
    if not cigar_str or cigar_str=="*":
        return []
    return [(int(m.group(1)),m.group(2)) for m in re.finditer(r"(\d+)([MIDNSHP=X])",cigar_str)]


def parse_cigar_get_query_pos(cigar_ops,ref_pos_1based,read_start_0based):
    """
    Map 1-based reference position to 0-based index into query SEQ.
    cigar_ops: list of (length, op_char) from parse_cigar_string.
    """
    if not cigar_ops:
        return (None,None)
    ref_offset=ref_pos_1based-1-read_start_0based
    if ref_offset<0:
        # Target position lies before this read's aligned start.
        return (None,None)
    query_pos=0
    for length,op in cigar_ops:
        if op in "M=X":
            if ref_offset<length:
                return (query_pos+ref_offset,None)
            query_pos+=length
            ref_offset-=length
        elif op=="I":
            query_pos+=length
        elif op in "DN":
            if ref_offset<length:
                return (None,None)
            ref_offset-=length
        elif op=="S":
            query_pos+=length
        elif op in "HP":
            pass
        if ref_offset<0:
            break
    return (None,None)


def _read_skip_flag(flag_int):
    return bool(flag_int & (FLAG_UNMAP | FLAG_SECONDARY | FLAG_SUPPLEMENTARY))


def check_read_for_snvs(seq,cigar_ops,read_start_0based,pos1,alt1,pos2,alt2):
    """SNV-oriented: compare single query base to expected alt (first char if multi)."""
    if not seq or seq=="*":
        return (False,False,None,None)
    a1=(alt1 or "").strip().upper()[:1]
    a2=(alt2 or "").strip().upper()[:1]
    if not a1 or not a2:
        return (False,False,None,None)
    q1,_=parse_cigar_get_query_pos(cigar_ops,pos1,read_start_0based)
    q2,_=parse_cigar_get_query_pos(cigar_ops,pos2,read_start_0based)
    has1=has2=False
    b1=b2=None
    if q1 is not None and 0<=q1<len(seq):
        b1=seq[q1].upper()
        has1=b1==a1
    if q2 is not None and 0<=q2<len(seq):
        b2=seq[q2].upper()
        has2=b2==a2
    return (has1,has2,b1,b2)


def iter_samtools_view(aln_path,region,min_mapq=20,reference_fasta=None):
    """
    Yield SAM alignment lines (not headers) from: samtools view [-q] [-T ref] aln region
    """
    cmd=["samtools","view"]
    if min_mapq>0:
        cmd.extend(["-q",str(min_mapq)])
    if reference_fasta:
        cmd.extend(["-T",reference_fasta])
    cmd.extend([aln_path,region])
    proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    try:
        for line in proc.stdout:
            line=line.rstrip("\n")
            if not line or line.startswith("@"):
                continue
            yield line
    finally:
        proc.stdout.close()
        err=proc.stderr.read()
        code=proc.wait()
        if code!=0:
            raise RuntimeError(f"samtools view failed (exit {code}): {err.strip() or cmd}")


def check_mnp_in_alignment(aln_path,chrom,pos1,ref1,alt1,pos2,ref2,alt2,min_mapq=20,reference_fasta=None):
    """
    Stream overlapping reads with samtools view; count same-read support for both alts.
    ref1/ref2 are not used for base checks (alt alleles only); kept for API compatibility.
    """
    lo=min(pos1,pos2)
    hi=max(pos1,pos2)
    region=f"{chrom}:{max(1,lo-1)}-{hi+1}"

    reads_with_snv1=reads_with_snv2=reads_with_both=reads_with_neither=total_reads=0
    read_details=[]

    for line in iter_samtools_view(aln_path,region,min_mapq=min_mapq,reference_fasta=reference_fasta):
        fields=line.split("\t")
        if len(fields)<11:
            continue
        qname,flag_s,rname,pos_s,mapq_s,cigar_s=(fields[0],fields[1],fields[2],fields[3],fields[4],fields[5])
        try:
            flag=int(flag_s)
            mapq=int(mapq_s)
            pos_1based=int(pos_s)
        except ValueError:
            continue
        if _read_skip_flag(flag):
            continue
        if rname!=chrom:
            continue
        # MAPQ filter: samtools -q when min_mapq > 0

        read_start_0=pos_1based-1
        cigar_ops=parse_cigar_string(cigar_s)
        seq=fields[9]

        total_reads+=1
        has_snv1,has_snv2,base1,base2=check_read_for_snvs(seq,cigar_ops,read_start_0,pos1,alt1,pos2,alt2)
        if has_snv1:
            reads_with_snv1+=1
        if has_snv2:
            reads_with_snv2+=1
        if has_snv1 and has_snv2:
            reads_with_both+=1
        if not has_snv1 and not has_snv2:
            reads_with_neither+=1
        read_details.append({"read_name":qname,"has_snv1":has_snv1,"has_snv2":has_snv2,"base1":base1,"base2":base2,"mapq":mapq})

    return {
        "total_reads":total_reads,
        "reads_with_snv1":reads_with_snv1,
        "reads_with_snv2":reads_with_snv2,
        "reads_with_both":reads_with_both,
        "reads_with_neither":reads_with_neither,
        "cis_evidence":reads_with_both,
        "trans_evidence":reads_with_snv1+reads_with_snv2-2*reads_with_both,
        "read_details":read_details,
    }


def process_mnp_pairs(mnp_csv,aln_path,output_csv,min_mapq=20,reference_fasta=None):
    mnp_pairs=[]
    with open(mnp_csv,"r") as infile:
        reader=csv.DictReader(infile)
        var1=None
        for row in reader:
            if var1 is None:
                var1=row
            else:
                var2=row
                if var1.get("Chr")==var2.get("Chr") and (
                    var1.get("Sample.ID")==var2.get("Sample.ID")
                    or var1.get("Tumor.ID")==var2.get("Tumor.ID")
                ):
                    mnp_pairs.append((var1,var2))
                var1=var2
                var2=None

    results=[]
    for i,(var1,var2) in enumerate(mnp_pairs):
        c=var1.get("Chr","")
        p1=int(var1.get("Start",0))
        ref1=var1.get("REF","")
        alt1=var1.get("ALT","")
        p2=int(var2.get("Start",0))
        ref2=var2.get("REF","")
        alt2=var2.get("ALT","")
        print(
            f"Processing pair {i+1}/{len(mnp_pairs)}: {c}:{p1} {ref1}>{alt1} and {c}:{p2} {ref2}>{alt2}",
            file=sys.stderr,
        )
        stats=check_mnp_in_alignment(aln_path,c,p1,ref1,alt1,p2,ref2,alt2,min_mapq,reference_fasta)
        result={
            "Chr":c,
            "Pos1":p1,
            "Ref1":ref1,
            "Alt1":alt1,
            "Pos2":p2,
            "Ref2":ref2,
            "Alt2":alt2,
            "Total_Reads":stats["total_reads"],
            "Reads_with_SNV1":stats["reads_with_snv1"],
            "Reads_with_SNV2":stats["reads_with_snv2"],
            "Reads_with_Both":stats["reads_with_both"],
            "Reads_with_Neither":stats["reads_with_neither"],
            "Cis_Evidence":stats["cis_evidence"],
            "Trans_Evidence":stats["trans_evidence"],
            "Cis_Ratio": f"{stats['reads_with_both']/max(stats['total_reads'],1):.3f}"
            if stats["total_reads"]>0
            else "0.000",
        }
        for key in ["Gene","Sample.ID","Tumor.ID","HGVSc","HGVSp"]:
            if key in var1:
                result[f"Var1_{key}"]=var1.get(key,".")
            if key in var2:
                result[f"Var2_{key}"]=var2.get(key,".")
        results.append(result)

    _write_results_csv(results,output_csv)


def _write_results_csv(results,output_csv):
    if output_csv:
        with open(output_csv,"w",newline="") as outfile:
            if results:
                fieldnames=list(results[0].keys())
                writer=csv.DictWriter(
                    outfile,
                    fieldnames=fieldnames,
                    delimiter=",",
                    restval=".",
                    extrasaction="ignore",
                    quoting=csv.QUOTE_NONNUMERIC,
                    dialect="excel",
                )
                writer.writeheader()
                for result in results:
                    writer.writerow(result)
        print(f"Processed {len(results)} MNP pairs -> {output_csv}",file=sys.stderr)
        return

    outfile=sys.stdout
    if results:
        fieldnames=list(results[0].keys())
        writer=csv.DictWriter(
            outfile,
            fieldnames=fieldnames,
            delimiter=",",
            restval=".",
            extrasaction="ignore",
            quoting=csv.QUOTE_NONNUMERIC,
            dialect="excel",
        )
        writer.writeheader()
        for result in results:
            writer.writerow(result)
    print(f"Processed {len(results)} MNP pairs",file=sys.stderr)


def _write_read_details_tsv(read_detail_rows,out_tsv):
    with open(out_tsv,"w",newline="") as outfile:
        fieldnames=[
            "pair_id",
            "sample",
            "id1",
            "id2",
            "Chr",
            "Pos1",
            "Ref1",
            "Alt1",
            "Pos2",
            "Ref2",
            "Alt2",
            "read_name",
            "mapq",
            "base1",
            "base2",
            "has_snv1",
            "has_snv2",
            "read_class",
        ]
        writer=csv.DictWriter(outfile,fieldnames=fieldnames,delimiter="\t",extrasaction="ignore")
        writer.writeheader()
        for row in read_detail_rows:
            writer.writerow(row)


def _classify_read(has_snv1,has_snv2):
    if has_snv1 and has_snv2:
        return "both_alt"
    if has_snv1:
        return "snv1_only"
    if has_snv2:
        return "snv2_only"
    return "neither_alt"


def process_targets_tsv(
    targets_tsv,bam_table,output_csv,min_mapq=20,reference_fasta=None,read_details_tsv=None
):
    """
    Tabular per-sample targets: requires columns id1, id2, sample (VCF-style IDs for the two SNVs).
    Optional mnp_id or pair_id (label in output). Optional VAF/GT columns echoed if present.
    Input is typically from choose_mnp_read_check_samples.py (mnp_id, sample, id1, id2) or legacy
    full PASS from test_mnp_sample_info.py.
    """
    bams=load_bam_table(bam_table)
    results=[]
    read_detail_rows=[]
    skipped=0
    processed=0
    seen=0
    with open(targets_tsv,"r",newline="") as count_in:
        total_targets=sum(1 for _ in csv.DictReader(count_in,delimiter="\t"))
    print(f"Targets to evaluate: {total_targets}",file=sys.stderr)
    with open(targets_tsv,"r",newline="") as infile:
        reader=csv.DictReader(infile,delimiter="\t")
        for row in reader:
            seen+=1
            id1=(row.get("id1") or "").strip()
            id2=(row.get("id2") or "").strip()
            sample=(row.get("sample") or "").strip()
            if not id1 or not id2 or not sample:
                skipped+=1
                continue
            p1,p2=parse_variant_id(id1),parse_variant_id(id2)
            if not p1 or not p2:
                print(f"Skip unparseable id: {id1!r} {id2!r}",file=sys.stderr)
                skipped+=1
                continue
            chr1,pos1,ref1,alt1=p1
            chr2,pos2,ref2,alt2=p2
            if chr1!=chr2:
                skipped+=1
                continue
            aln_path=bams.get(sample)
            if not aln_path or not os.path.isfile(aln_path):
                print(f"Skip sample {sample}: missing alignment path",file=sys.stderr)
                skipped+=1
                continue
            stats=check_mnp_in_alignment(aln_path,chr1,pos1,ref1,alt1,pos2,ref2,alt2,min_mapq,reference_fasta)
            processed+=1
            row_label=(row.get("mnp_id") or row.get("pair_id") or ".").strip() or "."
            phaseable_reads=stats["reads_with_both"]+stats["trans_evidence"]
            cis_fraction_phaseable=stats["reads_with_both"]/phaseable_reads if phaseable_reads>0 else 0.0
            strong_cis_default=stats["reads_with_both"]>=3 and phaseable_reads>=5 and cis_fraction_phaseable>=0.90
            results.append(
                {
                    "pair_id": row_label,
                    "sample": sample,
                    "id1": id1,
                    "id2": id2,
                    "Chr": chr1,
                    "Pos1": pos1,
                    "Ref1": ref1,
                    "Alt1": alt1,
                    "Pos2": pos2,
                    "Ref2": ref2,
                    "Alt2": alt2,
                    "VAF1":row.get("VAF1","."),
                    "VAF2":row.get("VAF2","."),
                    "VAF_diff":row.get("VAF_diff","."),
                    "Cis_likely":row.get("Cis_likely","."),
                    "Total_Reads": stats["total_reads"],
                    "Reads_with_SNV1": stats["reads_with_snv1"],
                    "Reads_with_SNV2": stats["reads_with_snv2"],
                    "Reads_with_Both": stats["reads_with_both"],
                    "Reads_with_Neither": stats["reads_with_neither"],
                    "Cis_Evidence": stats["cis_evidence"],
                    "Trans_Evidence": stats["trans_evidence"],
                    "Phaseable_Reads": phaseable_reads,
                    "Cis_Fraction_Phaseable": f"{cis_fraction_phaseable:.3f}",
                    "Cis_Ratio": f"{stats['reads_with_both']/max(stats['total_reads'],1):.3f}"
                    if stats["total_reads"]>0
                    else "0.000",
                    "Strong_Cis_Default": "1" if strong_cis_default else "0",
                }
            )
            if read_details_tsv:
                for detail in stats["read_details"]:
                    has_snv1=bool(detail["has_snv1"])
                    has_snv2=bool(detail["has_snv2"])
                    read_class=_classify_read(has_snv1,has_snv2)
                    read_detail_rows.append(
                        {
                            "pair_id": row_label,
                            "sample": sample,
                            "id1": id1,
                            "id2": id2,
                            "Chr": chr1,
                            "Pos1": pos1,
                            "Ref1": ref1,
                            "Alt1": alt1,
                            "Pos2": pos2,
                            "Ref2": ref2,
                            "Alt2": alt2,
                            "read_name":detail.get("read_name","."),
                            "mapq":detail.get("mapq","."),
                            "base1":detail.get("base1","."),
                            "base2":detail.get("base2","."),
                            "has_snv1": "1" if has_snv1 else "0",
                            "has_snv2": "1" if has_snv2 else "0",
                            "read_class": read_class,
                        }
                    )
            if seen%PROGRESS_EVERY_N==0:
                pct=(seen/total_targets*100.0) if total_targets>0 else 100.0
                print(
                    f"Processed {seen}/{total_targets} targets ({pct:.1f}%); "
                    f"evaluated={processed}, skipped={skipped}",
                    file=sys.stderr,
                )
    _write_results_csv(results,output_csv)
    if read_details_tsv:
        _write_read_details_tsv(read_detail_rows,read_details_tsv)
        print(f"Read-level details: {len(read_detail_rows)} rows -> {read_details_tsv}",file=sys.stderr)
    print(f"Targets read-check: {len(results)} rows, skipped {skipped}",file=sys.stderr)


def main():
    parser=argparse.ArgumentParser(
        description="MNP same-read check via samtools view (no pysam). "
        "Requires `samtools` on PATH."
    )
    parser.add_argument("-o","--output_csv",default=None,help="Output CSV (default: stdout)")
    parser.add_argument("-q","--min_mapq",type=int,default=20,help="Minimum MAPQ (default: 20)")
    parser.add_argument("-T","--reference",default=None,help="Reference FASTA for CRAM decoding (-T passed to samtools view)")
    src=parser.add_mutually_exclusive_group(required=True)
    src.add_argument("-t","--targets-tsv",default=None,help="Per-sample targets TSV: mnp_id, sample, id1, id2 (choose_mnp_read_check_samples.py); or legacy full PASS from test_mnp_sample_info.py (extra columns ignored).")
    src.add_argument("-i","--mnp-csv",default=None,help="MNP pair CSV from find_mnp_variants.py")
    parser.add_argument("-b","--bam",default=None,help="BAM/CRAM (with --mnp-csv only)")
    parser.add_argument("-B","--bam-table",default=None,help="sample<TAB>path (with --targets-tsv / --pass-tsv only)")
    parser.add_argument("-d","--read-details-tsv",default=None,help="Optional TSV of per-read evidence for each target pair")
    args=parser.parse_args()

    if args.targets_tsv:
        if not args.bam_table:
            parser.error("--targets-tsv requires --bam-table")
        process_targets_tsv(args.targets_tsv,args.bam_table,args.output_csv,args.min_mapq,args.reference,args.read_details_tsv)
    else:
        if not args.bam:
            parser.error("--mnp-csv requires -b/--bam")
        process_mnp_pairs(args.mnp_csv,args.bam,args.output_csv,args.min_mapq,args.reference)


if __name__=="__main__":
    main()
