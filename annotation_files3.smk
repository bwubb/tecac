#needs vep_vcf_parser2.py

# =============================================================================
# ANNOTATION FILES (annotation_files3.smk)
# =============================================================================
# - Passing-sample subset → het_miss BCF + no-sample VCF; VEP that track for QC CSVs.
# - Post-MNP: PASS / rv-qc exclusions → bcftools -e ID=@file on MNP GT BCF + merged no-sample MNP VEP. include rare_variant_qc.smk.
# =============================================================================

include: "mnp.smk"
include: "rare_variant_qc.smk"

CHROMOSOMES_AUTOSOMAL=list(range(1,23))

wildcard_constraints:
    CHR='[0-9]+'

rule annotation_files:
    input:
        expand("data/preprocess/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL),
        postmnp_stats=expand("data/qc/reports/chr{CHR}.postmnp.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        rare_qc_exclusions="data/qc/exclusions/rare_variant_exclusions.txt",

# -----------------------------------------------------------------------------
# 1. Subset to passing samples only (from build). Output: BCF + no-sample VCF.
# -----------------------------------------------------------------------------
rule bcftools_subset_passing_het_miss:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        samples="data/qc/passing_samples.txt"
    output:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.bcf",
        vcf="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vcf"
    shell:
        """
        bcftools view -S {input.samples} -a -Ob -W=csi -o {output.bcf} {input.bcf}
        bcftools view -G -Ov -o {output.vcf} {output.bcf}
        """

# -----------------------------------------------------------------------------
# 2. First VEP annotation (used directly as final annotation VEP in NO-MNP mode).
# -----------------------------------------------------------------------------
rule first_variant_annotation:
    input:
        "data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vcf"
    output:
        "data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.vcf"
    resources:
        mem_mb=32000
    shell:
        """
        export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache
        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i $PWD/{input} \
        -o $PWD/{output} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical --mane \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta /opt/vep/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.fasta \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
        --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
        --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
        """

rule parse_first_variant_annotation:
    input:
        "data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.vcf"
    output:
        "data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.report.csv"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m no_sample"


#MNP pipeline: Happens elswhere, leave it alone.

rule adjust_mnp_vep:
    input:
        ids="data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.id",
        vcf="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.vcf",
        mnp_bcf="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.bcf"
    output:
        not_mnp="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.not_mnp.bcf",
        mnp_bcf="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.bcf"
    shell:
        """
        bcftools view -e 'ID=@{input.ids}' -Ob -W=csi -o {output.not_mnp} {input.vcf}
        bcftools concat -a {input.mnp_bcf} {output.not_mnp} | bcftools sort -Ob -o {output.mnp_bcf}
        """

rule parse_adjust_mnp_vep:
    input:
        "data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.bcf"
    output:
        "data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.report.csv"
    params:
        vcf="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.vcf"
    shell:
        """
        bcftools view -Ov -o {params.vcf} {input}
        python vep_vcf_parser2.py -i {params.vcf} -o {output} -m no_sample
        """

rule preprocess_pathogenic_vus:
    input:
        expand("data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.report.csv",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/preprocess/pathogenic_vus.csv"
    shell:
        """
        echo '"ID","Gene","Variant.LoF_level","HGVSc","HGVSp"' > {output}
        for vep_file in {input}; do
            if [ -f "$vep_file" ]; then
                awk -F',' '($9 == "\\"1\\"" || $9 == "\\"2\\"") {{print $6 "," $7 "," $9 "," $13 "," $14}}' "$vep_file" >> {output}
            fi
        done
        """


rule adjust_mnp_gt:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.bcf",
        ids="data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.id",
        mnp="data/mnp/chr{CHR}.mnp.pass.gt.sorted.bcf"
    output:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.mnp.gt.bcf",
        not_mnp="data/bcftools/chr{CHR}.site-qc.het_miss.not_mnp.bcf"
    shell:
        """
        bcftools view -e 'ID=@{input.ids}' -W=csi -Ob -o {output.not_mnp} {input.bcf}
        bcftools concat -Ob -a {input.mnp} {output.not_mnp} | bcftools sort -Ob -W=csi -o {output.bcf}
        """

rule count_variant_types_postmnp:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.mnp.gt.bcf"
    output:
        "data/qc/reports/chr{CHR}.postmnp.variant_types.txt"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# Drop variants flagged by rare_variant_qc (per-chr TSV) before annotation PGEN / counts.
rule bcftools_exclude_rare_variant_qc_flags:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.mnp.gt.bcf",
        vep="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.bcf",
        ids="data/qc/exclusions/rare_variant_exclusions.txt"
    output:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.mnp.gt.rv-qc.bcf",
        vep="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.rv-qc.vcf",
    params:
        # The exclusion list is always built (input dependency above). Whether it is
        # actually applied is gated by config rare_variant_qc.apply_exclusions (default False),
        # so you can flip it on later without re-running the QC. When False, variants pass
        # through unchanged.
        excl=lambda wildcards, input: ("-e 'ID=@%s'" % input.ids
              if config.get("rare_variant_qc",{}).get("apply_exclusions",False) else "")
    shell:
        """
        bcftools view {params.excl} -Ob -W=csi -o {output.bcf} {input.bcf}
        bcftools view {params.excl} -Ov -o {output.vep} {input.vep}
        """

rule plink2_annotation_pgen:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.mnp.gt.rv-qc.bcf",
        sex_file="data/plink/chrX.sex_update.txt"
    output:
        pgen="data/preprocess/chr{CHR}.annotation.pgen",
        pvar="data/preprocess/chr{CHR}.annotation.pvar",
        psam="data/preprocess/chr{CHR}.annotation.psam"
    params:
        output_prefix="data/preprocess/chr{CHR}.annotation" 
    shell:
        """
        plink2 --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}
        """

# -----------------------------------------------------------------------------
# 5. Copy post-MNP VEP to preprocess path (stable filename for downstream).
# -----------------------------------------------------------------------------
rule stage_preprocess_annotation_no_sample_vep:
    input:
        "data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.mnp.rv-qc.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    shell:
        "cp -f {input} {output}"

# -----------------------------------------------------------------------------
# 6. Parse final annotation VEP to report CSV.
# -----------------------------------------------------------------------------
rule parse_annotation_no_sample_vep:
    input:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv"
    resources:
        mem_mb=32000,
        lsf_err="logs/lsf/parse_annotation_no_sample_vep.chr{CHR}.e",
        lsf_out="logs/lsf/parse_annotation_no_sample_vep.chr{CHR}.o"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m no_sample"
