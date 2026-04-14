# =============================================================================
# ANNOTATION FILES: subset to build passing samples → VEP → exclude MNPs from annotation → MNP genotype path → two concat+sorts
# =============================================================================
# Pipeline order:
#   1. Subset qc_filter1.bcf to passing samples only (from build: data/qc/passing_samples.txt)
#   2. First VEP on full no-sample VCF (for MNP detection and non-MNP annotations)
#   3. Parse first VEP report
#   4. Find MNP variants (pairs, regions, id_file)
#   5a. Exclude MNP from first no_sample.vep.vcf → non-MNP annotation VCF
#   5b–5i. Genotype path: exclude/filter MNP pairs, manage_mnp_gt, concat, sort → sorted BCF
#   5j–5l. Annotation path: mnp_gt → no_sample VCF → VEP (small); then concat+sort non-MNP.vep.vcf + MNP.vep.vcf → annotation.no_sample.vep.vcf
#   6. Un-annotated no_sample.vcf from sorted BCF; PLINK annotation pgen from sorted BCF
#   7. parse_vep on final annotation.no_sample.vep.vcf → report CSV
# =============================================================================
# Expects: build_files (or main) produces data/qc/passing_samples.txt (IID per line).
# =============================================================================
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
        postmnp_stats=expand("data/qc/reports/chr{CHR}.postmnp.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),

# -----------------------------------------------------------------------------
# 1. Subset to passing samples only (from build). Output: BCF + no-sample VCF for VEP.
# -----------------------------------------------------------------------------
rule bcftools_annotation_subset:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        samples="data/qc/passing_samples.txt"
    output:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.bcf",
        vcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vcf"
    resources:
        lsf_err="logs/lsf/bcftools_annotation_subset.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_annotation_subset.chr{CHR}.o"
    shell:
        """
        bcftools view -S {input.samples} -a -Ob -W=csi -o {output.bcf} {input.bcf}
        bcftools view -G -Ov -o {output.vcf} {output.bcf}
        """

# -----------------------------------------------------------------------------
# 2. First VEP annotation (no-sample VCF) for MNP detection and variant typing.
# -----------------------------------------------------------------------------
rule first_variants_annotation:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
       "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vcf"
    output:
        "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.vcf"
    resources:
        mem_mb=32000,
        lsf_err="logs/lsf/first_variants_annotation.chr{CHR}.e",
        lsf_out="logs/lsf/first_variants_annotation.chr{CHR}.o"
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
        --vcf --everything --canonical \
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

# -----------------------------------------------------------------------------
# 3. Parse first VEP report → CSV.
# -----------------------------------------------------------------------------
rule parse_first_vep:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.vcf"
    output:
        "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.report.csv"
    resources:
        lsf_err="logs/lsf/parse_first_vep.chr{CHR}.e",
        lsf_out="logs/lsf/parse_first_vep.chr{CHR}.o"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m no_sample"

# -----------------------------------------------------------------------------
# 4. Find MNP variants; split into non-MNP vs MNP-pairs-only (MNP split may be revised).
# -----------------------------------------------------------------------------
rule find_mnp_variants:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        csv="data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.report.csv"
    output:
        csv="data/mnp/chr{CHR}.annotation.mnp.csv",
        pairs="data/mnp/chr{CHR}.annotation.mnp.pairs.txt",
        regions="data/mnp/chr{CHR}.regions.txt",
        id_file="data/mnp/chr{CHR}.id.txt"
    resources:
        lsf_err="logs/lsf/find_mnp_variants.chr{CHR}.e",
        lsf_out="logs/lsf/find_mnp_variants.chr{CHR}.o"
    shell:
        "python find_mnp_variants.py -i {input.csv} -o {output.csv} -m no_sample -r {output.regions} -I {output.id_file}"

# 5a. Exclude MNP variants from first no_sample.vep.vcf → non-MNP annotation VCF (for annotation concat later)
rule bcftools_exclude_mnp_from_vep_vcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.vcf",
        id_file="data/mnp/chr{CHR}.id.txt"
    output:
        vcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.non_mnp.vep.vcf"
    resources:
        lsf_err="logs/lsf/bcftools_exclude_mnp_from_vep_vcf.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_exclude_mnp_from_vep_vcf.chr{CHR}.o"
    shell:
        """
        bcftools view -e 'ID=@{input.id_file}' -Ov -o {output.vcf} {input.vcf}
        """

# 5b. Extract variants that are NOT MNP pairs (genotype track; for concat after manage_mnp_gt)
rule bcftools_mnp_exclude_pairs:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.bcf",
        id_file="data/mnp/chr{CHR}.id.txt"
    output:
        vcf="data/mnp/chr{CHR}.het_miss.excl_mnp_pairs.vcf"
    resources:
        lsf_err="logs/lsf/bcftools_mnp_exclude_pairs.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_mnp_exclude_pairs.chr{CHR}.o"
    shell:
        """
        bcftools view -e 'ID=@{input.id_file}' -Ov -o {output.vcf} {input.bcf}
        """

# 5c. Extract ONLY MNP pair variants (id1, id2) - small VCF for manage_mnp_gt
rule bcftools_mnp_filter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.bcf",
        regions="data/mnp/chr{CHR}.regions.txt",
        id_file="data/mnp/chr{CHR}.id.txt"
    output:
        vcf="data/mnp/chr{CHR}.mnp_pairs_only.vcf"
    resources:
        lsf_err="logs/lsf/bcftools_mnp_filter.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_mnp_filter.chr{CHR}.o"
    shell:
        """
        bcftools view -R {input.regions} -i 'ID=@{input.id_file}' -Ov -o {output.vcf} {input.vcf}
        """

# 5d. Sample info for MNP pairs (for plan_mnp_gt / manage_mnp_gt)
rule bcftools_mnp_sample_info:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/mnp/chr{CHR}.mnp_pairs_only.vcf"
    output:
        "data/mnp/chr{CHR}.mnp_sample_info.txt"
    resources:
        lsf_err="logs/lsf/bcftools_mnp_sample_info.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_mnp_sample_info.chr{CHR}.o"
    shell:
        """
        bcftools query -f '[%ID %SAMPLE %GT %AD %DP %VAF\n]' {input} |
        awk -F' ' '$3=="0/1" || $3=="1/1"' > {output}
        """

# 5e. Test MNP sample info (PASS/FAIL)
rule test_mnp_sample_info:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pairs="data/mnp/chr{CHR}.annotation.mnp.pairs.txt",
        sample_info="data/mnp/chr{CHR}.mnp_sample_info.txt"
    output:
        pass_out="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",
        fail_out="data/mnp/chr{CHR}.mnp_sample_info.FAIL.txt"
    resources:
        lsf_err="logs/lsf/test_mnp_sample_info.chr{CHR}.e",
        lsf_out="logs/lsf/test_mnp_sample_info.chr{CHR}.o"
    shell:
        "python test_mnp_sample_info.py -p {input.pairs} -i {input.sample_info} -o {output.pass_out} -f {output.fail_out}"

# 5f. Plan MNP GT corrections
rule plan_mnp_gt:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        sample_info="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt"
    output:
        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt"
    params:
        ref="/home/bwubb/resources/Genomes/Human/hg38/Homo_sapiens_assembly38.fasta.gz"
    resources:
        lsf_err="logs/lsf/plan_mnp_gt.chr{CHR}.e",
        lsf_out="logs/lsf/plan_mnp_gt.chr{CHR}.o"
    shell:
        "python plan_mnp_gt.py -i {input.sample_info} -o {output.plan} -r {params.ref}"

# 5g. Apply GT plan to MNP pairs only (small input = low memory)
rule manage_mnp_gt:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt",
        vcf="data/mnp/chr{CHR}.mnp_pairs_only.vcf"
    output:
        vcf="data/mnp/chr{CHR}.mnp_gt.vcf"
    resources:
        mem_mb=32000,
        lsf_err="logs/lsf/manage_mnp_gt.chr{CHR}.e",
        lsf_out="logs/lsf/manage_mnp_gt.chr{CHR}.o"
    shell:
        "python manage_mnp_gt.py -p {input.plan} -i {input.vcf} -o {output.vcf}"

# 5h. Concat genotype: non-MNP + mnp_gt → combined VCF
rule bcftools_mnp_concat:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        excl="data/mnp/chr{CHR}.het_miss.excl_mnp_pairs.vcf",
        mnp_gt="data/mnp/chr{CHR}.mnp_gt.vcf"
    output:
        vcf="data/mnp/chr{CHR}.het_miss.mnp_gt.combined.vcf"
    resources:
        lsf_err="logs/lsf/bcftools_mnp_concat.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_mnp_concat.chr{CHR}.o"
    shell:
        """
        bcftools concat {input.excl} {input.mnp_gt} -Ov -o {output.vcf}
        """

# 5i. Sort combined genotype VCF → BCF
rule bcftools_sort_mnp_gt_bcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/mnp/chr{CHR}.het_miss.mnp_gt.combined.vcf"
    output:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.mnp_gt.sorted.bcf"
    resources:
        lsf_err="logs/lsf/bcftools_sort_mnp_gt_bcf.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_sort_mnp_gt_bcf.chr{CHR}.o"
    shell:
        """
        bcftools sort -W=csi -Ob -o {output.bcf} {input.vcf}
        """

# 5j. No-sample VCF from mnp_gt.vcf (for small VEP run)
rule bcftools_mnp_gt_no_sample:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/mnp/chr{CHR}.mnp_gt.vcf"
    output:
        vcf="data/mnp/chr{CHR}.mnp_gt.no_sample.vcf"
    resources:
        lsf_err="logs/lsf/bcftools_mnp_gt_no_sample.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_mnp_gt_no_sample.chr{CHR}.o"
    shell:
        "bcftools view -G -Ov -o {output.vcf} {input.vcf}"

# 5k. VEP on small MNP-only no_sample VCF
rule vep_mnp_only:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/mnp/chr{CHR}.mnp_gt.no_sample.vcf"
    output:
        "data/mnp/chr{CHR}.mnp_gt.no_sample.vep.vcf"
    resources:
        mem_mb=32000,
        lsf_err="logs/lsf/vep_mnp_only.chr{CHR}.e",
        lsf_out="logs/lsf/vep_mnp_only.chr{CHR}.o"
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
        --vcf --everything --canonical \
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

# 5l. Concat + sort annotation VCFs: non-MNP.vep.vcf + MNP.vep.vcf → final annotation.no_sample.vep.vcf
rule bcftools_concat_sort_annotation_vep:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        non_mnp="data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.non_mnp.vep.vcf",
        mnp="data/mnp/chr{CHR}.mnp_gt.no_sample.vep.vcf"
    output:
        vcf="data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    resources:
        lsf_err="logs/lsf/bcftools_concat_sort_annotation_vep.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_concat_sort_annotation_vep.chr{CHR}.o"
    params:
        tmp="data/preprocess/chr{CHR}.annotation.no_sample.vep.concat_tmp.vcf"
    shell:
        """
        bcftools concat {input.non_mnp} {input.mnp} -Ov -o {params.tmp}
        bcftools sort -Ov -o {output.vcf} {params.tmp}
        rm -f {params.tmp}
        """

# -----------------------------------------------------------------------------
# 6. Un-annotated no_sample VCF from sorted BCF (optional); PLINK annotation pgen/pvar/psam.
# -----------------------------------------------------------------------------
rule bcftools_annotation2_vcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.mnp_gt.sorted.bcf"
    output:
        vcf="data/preprocess/chr{CHR}.annotation.no_sample.vcf"  # un-annotated no_sample from combined genotype
    resources:
        lsf_err="logs/lsf/bcftools_annotation2_vcf.chr{CHR}.e",
        lsf_out="logs/lsf/bcftools_annotation2_vcf.chr{CHR}.o"
    shell:
        """
        bcftools view -G -Ov -o {output.vcf} {input.bcf}
        """

# Count variant types after MNP correction (same format as postfilter: TYPE COUNT)
rule count_variant_types_postmnp:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.mnp_gt.sorted.bcf"
    output:
        "data/qc/reports/chr{CHR}.postmnp.variant_types.txt"
    resources:
        lsf_err="logs/lsf/count_variant_types_postmnp.chr{CHR}.e",
        lsf_out="logs/lsf/count_variant_types_postmnp.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

rule plink2_annotation_pgen:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.mnp_gt.sorted.bcf",
        sex_file="data/plink/chrX.sex_update.txt"
    output:
        pgen="data/preprocess/chr{CHR}.annotation.pgen",
        pvar="data/preprocess/chr{CHR}.annotation.pvar",
        psam="data/preprocess/chr{CHR}.annotation.psam"
    params:
        output_prefix="data/preprocess/chr{CHR}.annotation"
    resources:
        lsf_err="logs/lsf/plink2_annotation_pgen.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_annotation_pgen.chr{CHR}.o"
    shell:
        """
        plink2 --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}
        """


# -----------------------------------------------------------------------------
# 7. Parse final annotation VCF (from concat+sort of non-MNP + MNP VEP) → report CSV.
# -----------------------------------------------------------------------------
rule parse_vep:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv"
    resources:
        mem_mb=32000,
        lsf_err="logs/lsf/parse_vep.chr{CHR}.e",
        lsf_out="logs/lsf/parse_vep.chr{CHR}.o"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m no_sample"