
#needs find_mnp_variants.py
#needs test_mnp_sample_info.py
#needs filter_mnp_pass_by_read_check.py
#needs plan_mnp_gt.py
#needs manage_mnp_gt.py
#needs check_mnp_reads.py
#needs choose_mnp_read_check_samples.py
#needs mnp_validation_metrics.py
#needs vep_vcf_parser2.py

import datetime

DATE=datetime.date.today().strftime('%Y%m%d')

# Read-check CSV for filter_mnp_pass_by_read_check. Optional sidestep: set input.mnp in config
# to an existing results file; otherwise use output from check_mnp_reads.
MNP_READ_CHECK=config.get("input",{}).get("mnp") or "data/mnp/mnp_read_check.results.mapq20.csv"
# minimal: one read-check per MNP (greedy BAM set); pair-level filter.
# all_cases: read-check every (MNP, case) PASS row; per-sample filter.
MNP_READ_CHECK_MODE=config.get("mnp",{}).get("read_check_mode","minimal")
MNP_ALL_CARRIERS="--all-carriers" if MNP_READ_CHECK_MODE=="all_cases" else ""
MNP_FILTER_PER_SAMPLE="--per-sample" if MNP_READ_CHECK_MODE=="all_cases" else ""
MNP_CASE_LIST=config.get("input",{}).get("cases","case.list")
# set_cover: one BAM per MNP (minimal); all_cases: every case carrier per MNP (default)
MNP_READ_CHECK_MODE=config.get("mnp",{}).get("read_check_mode","all_cases")
MNP_PAIR_RULE=config.get("mnp",{}).get("pair_rule","any")


# Target: same outputs as mnp for GT/VEP handoff, but the DAG does not require assignments/samples.
# Put read-check outputs in place yourself (same paths as rule check_mnp_reads) so that rule is skipped.
rule mnp_gt_for_downstream:
    input:
        ids=expand("data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.id",CHR=range(1,23)),
        bcf=expand("data/mnp/chr{CHR}.mnp.pass.gt.sorted.bcf",CHR=range(1,23)),
        vep_vcf=expand("data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.vcf",CHR=range(1,23)),
        vep_bcf=expand("data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.bcf",CHR=range(1,23))


# Stop here before samtools read validation (check_mnp_reads).
# Run: snakemake -s qc_pipeline.smk mnp_before_read_check --configfile config.yml
rule mnp_before_read_check:
    input:
        assignments="data/mnp/mnp_read_check.assignments.tsv",
        samples="data/mnp/mnp_read_check.samples.txt",
        pass_files=expand("data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",CHR=range(1,23)),


# Target: full path including choose_mnp_read_check_samples, check_mnp_reads, validation metrics.
rule mnp_check_and_blacklist:
    input:
        "data/mnp/mnp_read_check.assignments.tsv",
        "data/mnp/mnp_read_check.samples.txt",
        "data/mnp/mnp_read_check.results.mapq20.csv",
        expand("data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.txt",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.id",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp_gt.plan.txt",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp.pass.vcf",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp.pass.gt.vcf",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp.pass.gt.sorted.bcf",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vcf",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.vcf",CHR=range(1,23)),
        expand("data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.report.csv",CHR=range(1,23)),
        "data/mnp/mnp.validation.metrics.tsv",
        "data/mnp/mnp_blacklist.txt",


rule find_mnp_variants:
    input:
        csv="data/bcftools/chr{CHR}.site-qc.het_miss.no_sample.vep.report.csv"
    output:
        csv="data/mnp/chr{CHR}.annotation.mnp.csv",
        pairs="data/mnp/chr{CHR}.annotation.mnp.pairs.txt",
        regions="data/mnp/chr{CHR}.regions.txt",
        id_file="data/mnp/chr{CHR}.id.txt"
    shell:
        "python find_mnp_variants.py -i {input.csv} -o {output.csv} -m no_sample -r {output.regions} -I {output.id_file}"

rule bcftools_mnp_filter:
    input:
        vcf="data/bcftools/chr{CHR}.site-qc.het_miss.bcf",
        regions="data/mnp/chr{CHR}.regions.txt",
        id_file="data/mnp/chr{CHR}.id.txt"
    output:
        vcf="data/mnp/chr{CHR}.mnp.vcf"
    shell:
        """
        bcftools view -R {input.regions} -i 'ID=@{input.id_file}' -Ov -o {output.vcf} {input.vcf}
        """

rule bcftools_mnp_sample_info:
    input:
        vcf="data/mnp/chr{CHR}.mnp.vcf"
    output:
        "data/mnp/chr{CHR}.mnp_sample_info.txt"
    shell:
        """
        bcftools query -f '[%ID %SAMPLE %GT %AD %DP %VAF\n]' {input} |
        awk -F' ' '$3=="0/1" || $3=="1/1"' > {output}
        """

rule test_mnp_sample_info:
    input:
        pairs="data/mnp/chr{CHR}.annotation.mnp.pairs.txt",
        sample_info="data/mnp/chr{CHR}.mnp_sample_info.txt",
    output:
        pass_out="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",
        fail_out="data/mnp/chr{CHR}.mnp_sample_info.FAIL.txt"
    params:
        sample_list=MNP_CASE_LIST
    shell:
        "python test_mnp_sample_info.py -p {input.pairs} -i {input.sample_info} -o {output.pass_out} -f {output.fail_out} --sample-list {params.sample_list} --min-samples 1"


rule choose_mnp_read_check_samples:
    input:
        pass_files=expand("data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",CHR=range(1,23)),
    output:
        assignments="data/mnp/mnp_read_check.assignments.tsv",
        samples="data/mnp/mnp_read_check.samples.txt",
    params:
        case_list=MNP_CASE_LIST,
        mode=MNP_READ_CHECK_MODE,
        # set_cover (minimal BAMs): --mode set_cover  (drop --case-list)
    shell:
        "python choose_mnp_read_check_samples.py {input.pass_files} "
        "-o {output.assignments} -s {output.samples} "
        "--mode {params.mode} --case-list {params.case_list}"

rule check_mnp_reads:
    input:
        assignments="data/mnp/mnp_read_check.assignments.tsv",
        samples="data/mnp/mnp_read_check.samples.txt",
    output:
        csv="data/mnp/mnp_read_check.results.mapq20.csv",
        read_details="data/mnp/mnp_read_check.read_details.mapq20.tsv"
    params:
        bam_table=config.get("mnp",{}).get("bam_table","mnp_bam.table"),
        reference="resources/GRCh38_full_analysis_set_plus_decoy_hla.fa",
        min_mapq=20
    shell:
        "python check_mnp_reads.py "
        "--targets-tsv {input.assignments} "
        "--bam-table {params.bam_table} "
        "-q {params.min_mapq} "
        "-T {params.reference} "
        "-o {output.csv} "
        "--read-details-tsv {output.read_details}"

rule filter_mnp_pass_by_read_check:
    input:
        pass_in="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",
        read_check=MNP_READ_CHECK,
    params:
        pair_rule=MNP_PAIR_RULE,
    output:
        pass_out="data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.txt",
        ids_out="data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.id"
    shell:
        "python filter_mnp_pass_by_read_check.py "
        "-r {input.read_check} "
        "-i {input.pass_in} "
        "-o {output.pass_out} "
        "--ids-out {output.ids_out} "
        "--pair-rule {params.pair_rule}"

rule plan_mnp_gt:
    input:
        sample_info="data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.txt"
    output:
        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt"
    params:
        ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    shell:
        "python plan_mnp_gt.py -i {input.sample_info} -o {output.plan} -r {params.ref}"

rule subset_mnp_vcf:
    input:
        vcf="data/mnp/chr{CHR}.mnp.vcf",
        ids="data/mnp/chr{CHR}.mnp_sample_info.read_check_filtered.PASS.id"
    output:
        vcf="data/mnp/chr{CHR}.mnp.pass.vcf"
    shell:
        "bcftools view -i 'ID=@{input.ids}' -Ov -o {output.vcf} {input.vcf}"

rule manage_mnp_gt:
    input:
        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt",
        vcf="data/mnp/chr{CHR}.mnp.pass.vcf"
    output:
        vcf="data/mnp/chr{CHR}.mnp.pass.gt.vcf"
    shell:
        "python manage_mnp_gt.py -p {input.plan} -i {input.vcf} -o {output.vcf}"

rule bcftools_mnp_pass_gt_to_bcf_and_no_sample:
    input:
        vcf="data/mnp/chr{CHR}.mnp.pass.gt.vcf"
    output:
        vcf="data/mnp/chr{CHR}.mnp.pass.gt.sorted.vcf",
        bcf="data/mnp/chr{CHR}.mnp.pass.gt.sorted.bcf",
        no_sample_vcf="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vcf"
    params:
        tmp="data/mnp/"
    shell:
        """
        bcftools sort -Ov -o {output.vcf} -T {params.tmp} {input.vcf}
        bcftools view -W=csi -Ob -o {output.bcf} {output.vcf}
        bcftools view -G -Ov -o {output.no_sample_vcf} {output.vcf}
        """

##FIGURE OUT A WAY TO CRUSH THESE SCRIPTS DOWN INTO ONE SCRIPT.

rule annotate_mnp_gt:
    input:
        vcf="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vcf"
    output:
        vcf="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.vcf"
    shell:
        """
        export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache
        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i $PWD/{input.vcf} \
        -o $PWD/{output.vcf} \
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

rule bcftools_mnp_bcf:
    input:
        vcf="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.vcf"
    output:
        bcf="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.bcf"
    shell:
        "bcftools view -W=csi -Ob -o {output.bcf} {input.vcf}"

rule parse_mnp_gt_vep:
    input:
        vcf="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.vcf"
    output:
        csv="data/mnp/chr{CHR}.mnp.pass.no_sample.sorted.vep.report.csv"
    shell:
        "python vep_vcf_parser2.py -i {input.vcf} -o {output.csv} -m cohort"

rule mnp_validation_metrics:
    input:
        candidate_pairs=expand("data/mnp/chr{CHR}.annotation.mnp.pairs.txt",CHR=range(1,23)),
        pass_files=expand("data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",CHR=range(1,23)),
        assignments="data/mnp/mnp_read_check.assignments.tsv",
        read_check="data/mnp/mnp_read_check.results.mapq20.csv"
    output:
        tsv="data/mnp/mnp.validation.metrics.tsv"
    shell:
        "python mnp_validation_metrics.py "
        "--candidate-pairs {input.candidate_pairs} "
        "--pass-files {input.pass_files} "
        "--assignments {input.assignments} "
        "--read-check {input.read_check} "
        "-o {output.tsv}"

rule mnp_blacklist_from_pairs:
    input:
        pairs=expand("data/mnp/chr{CHR}.annotation.mnp.pairs.txt",CHR=range(1,23))
    output:
        txt="data/mnp/mnp_blacklist.txt"
    shell:
        "awk -F'\\t' 'NF>=2{{print $1\"\\n\"$2}}' {input.pairs} | awk 'NF>0{{print}}' | sort -u > {output.txt}"
