# =============================================================================
# MISSINGNESS BIAS: differential missingness (case vs control, freeze 2 vs 3)
# =============================================================================
# Generates and parses PLINK --test-missing results (case/control and freeze 2 vs 3).
# Does NOT exclude samples; keeps thresholds and builds SNP exclusion lists for
# optional use (e.g. reporting, or future filtering).
#
# Outputs:
#   - diff_miss_fail_variants.txt       SNPs failing case/control (P < thr)
#   - diff_miss_freeze_fail_variants.txt SNPs failing freeze 2 vs 3 (P < thr)
#   - diff_miss_snp_exclusion_list.txt  Union (either test)
#   - diff_miss_both_fail_variants.txt  Intersection (both tests) — what the two lists have in common
#   - diff_miss_bias_snps_for_plots.tsv Table for QC report: ID, CHR, POS, fail_cc, fail_freeze,
#       P_cc, P_freeze, F_MISS_* and neglog10_P_* for genome-position and loading plots.
#
# QC report (future): e.g. bias SNP locations along the genome; loadings (-log10(P) or F_MISS) along the genome.
#
# Expects when included:
#   - config['qc']['diff_miss_thr'] (default 1e-5)
#   - CHROMOSOMES_AUTOSOMAL
#   - Inputs: data/plink/chr{CHR}.qc_filter2.pgen/.pvar/.psam (from build_files),
#             config input covariates file (FID IID FREEZE STATUS; STATUS 1/2)
#   - Internal intermediates: chr*.bed_diffmiss_temp.*, chr*.bed_freeze_temp.*,
#             chr*.pheno.txt, chr*.pheno_freeze.txt, chr*.freeze_sample_list.txt
#   - Plot prep pvar source: data/plink/chr{CHR}.qc_filter2.pvar
# =============================================================================

CHROMOSOMES_AUTOSOMAL=list(range(1,23))
DIFF_MISS_THR=config.get('qc', {}).get('diff_miss_thr', 1e-5)
COVARIATES_FILE=config.get('input', {}).get('covariates', '')

wildcard_constraints:
    CHR='[0-9]+'

rule missingness_bias:
    input:
        expand("data/plink/chr{CHR}.test.missing", CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.test_freeze.missing", CHR=CHROMOSOMES_AUTOSOMAL),
        "data/qc/reports/diff_miss_freeze_summary.csv",
        "data/qc/reports/diff_miss_summary.csv",

# -----------------------------------------------------------------------------
# Generate case/control differential missingness (STATUS: 1=control, 2=case)
# -----------------------------------------------------------------------------
rule plink_test_missing:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.qc_filter2.pgen",
        pvar="data/plink/chr{CHR}.qc_filter2.pvar",
        psam="data/plink/chr{CHR}.qc_filter2.psam",
        covariates=COVARIATES_FILE
    output:
        "data/plink/chr{CHR}.test.missing"
    params:
        pfile_prefix="data/plink/chr{CHR}.qc_filter2",
        bed_prefix="data/plink/chr{CHR}.bed_diffmiss_temp",
        pheno_file="data/plink/chr{CHR}.pheno.txt",
        outputname="data/plink/chr{CHR}.test"
    resources:
        lsf_err="logs/lsf/plink_test_missing.chr{CHR}.e",
        lsf_out="logs/lsf/plink_test_missing.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.pfile_prefix} --make-bed --out {params.bed_prefix}
        awk 'BEGIN {{OFS="\\t"}} NR==1 {{print "FID","IID","STATUS"; next}} $4==1 || $4==2 {{print $1,$2,$4}}' {input.covariates} > {params.pheno_file}
        plink --bfile {params.bed_prefix} --test-missing --pheno {params.pheno_file} --allow-no-sex --out {params.outputname}
        rm -f {params.bed_prefix}.bed {params.bed_prefix}.bim {params.bed_prefix}.fam {params.bed_prefix}.log {params.pheno_file}
        """

# -----------------------------------------------------------------------------
# Generate freeze 2 vs 3 differential missingness (FREEZE 2->1, 3->2)
# -----------------------------------------------------------------------------
rule plink_test_missing_freeze:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.qc_filter2.pgen",
        pvar="data/plink/chr{CHR}.qc_filter2.pvar",
        psam="data/plink/chr{CHR}.qc_filter2.psam",
        covariates=COVARIATES_FILE
    output:
        "data/plink/chr{CHR}.test_freeze.missing"
    params:
        pfile_prefix="data/plink/chr{CHR}.qc_filter2",
        bed_prefix="data/plink/chr{CHR}.bed_freeze_temp",
        sample_list="data/plink/chr{CHR}.freeze_sample_list.txt",
        pheno_file="data/plink/chr{CHR}.pheno_freeze.txt",
        outputname="data/plink/chr{CHR}.test_freeze"
    resources:
        lsf_err="logs/lsf/plink_test_missing_freeze.chr{CHR}.e",
        lsf_out="logs/lsf/plink_test_missing_freeze.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.pfile_prefix} --make-bed --out {params.bed_prefix}
        awk 'NR>1 && ($3==2 || $3==3) {{print $1, $2}}' {input.covariates} > {params.sample_list}
        awk 'BEGIN {{OFS="\\t"}} NR==1 {{print "FID","IID","FREEZE_GROUP"; next}} $3==2 || $3==3 {{g=($3==2?1:2); print $1,$2,g}}' {input.covariates} > {params.pheno_file}
        plink --bfile {params.bed_prefix} --keep {params.sample_list} --test-missing --pheno {params.pheno_file} --allow-no-sex --out {params.outputname}
        rm -f {params.bed_prefix}.bed {params.bed_prefix}.bim {params.bed_prefix}.fam {params.bed_prefix}.log {params.sample_list} {params.pheno_file}
        """

# -----------------------------------------------------------------------------
# Case/control: parse test.missing → SNPs with P < threshold (exclusion list)
# -----------------------------------------------------------------------------
rule parse_diff_miss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/plink/chr{CHR}.test.missing"
    output:
        "data/qc/exclusions/chr{CHR}.diff_miss_fail.txt"
    params:
        threshold=DIFF_MISS_THR
    resources:
        lsf_err="logs/lsf/parse_diff_miss.chr{CHR}.e",
        lsf_out="logs/lsf/parse_diff_miss.chr{CHR}.o"
    shell:
        """
        awk 'NR>1 && $5 < {params.threshold} {{print $2}}' {input} > {output}
        touch {output}
        """

# -----------------------------------------------------------------------------
# Freeze 2 vs 3: parse test_freeze.missing → SNPs with P < threshold
# -----------------------------------------------------------------------------
rule parse_diff_miss_freeze:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/plink/chr{CHR}.test_freeze.missing"
    output:
        "data/qc/exclusions/chr{CHR}.diff_miss_freeze_fail.txt"
    params:
        threshold=DIFF_MISS_THR
    resources:
        lsf_err="logs/lsf/parse_diff_miss_freeze.chr{CHR}.e",
        lsf_out="logs/lsf/parse_diff_miss_freeze.chr{CHR}.o"
    shell:
        """
        awk 'NR>1 && $5 < {params.threshold} {{print $2}}' {input} > {output}
        touch {output}
        """

# -----------------------------------------------------------------------------
# Aggregate per-chr exclusion lists → one SNP list per test (and combined)
# -----------------------------------------------------------------------------
rule aggregate_diff_miss:
    input:
        expand("data/qc/exclusions/chr{CHR}.diff_miss_fail.txt", CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/exclusions/diff_miss_fail_variants.txt"
    resources:
        lsf_err="logs/lsf/aggregate_diff_miss.e",
        lsf_out="logs/lsf/aggregate_diff_miss.o"
    shell:
        "cat {input} | sort -u > {output}"

rule aggregate_diff_miss_freeze:
    input:
        expand("data/qc/exclusions/chr{CHR}.diff_miss_freeze_fail.txt", CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/exclusions/diff_miss_freeze_fail_variants.txt"
    resources:
        lsf_err="logs/lsf/aggregate_diff_miss_freeze.e",
        lsf_out="logs/lsf/aggregate_diff_miss_freeze.o"
    shell:
        "cat {input} | sort -u > {output}"

# Combined SNP exclusion list: SNPs that fail either case/control or freeze (for optional use)
rule aggregate_diff_miss_combined:
    input:
        case_control="data/qc/exclusions/diff_miss_fail_variants.txt",
        freeze="data/qc/exclusions/diff_miss_freeze_fail_variants.txt"
    output:
        "data/qc/exclusions/diff_miss_snp_exclusion_list.txt"
    resources:
        lsf_err="logs/lsf/aggregate_diff_miss_combined.e",
        lsf_out="logs/lsf/aggregate_diff_miss_combined.o"
    shell:
        "cat {input.case_control} {input.freeze} | sort -u > {output}"

# SNPs that fail BOTH case/control and freeze (overlap of the two lists)
rule aggregate_diff_miss_both_fail:
    input:
        case_control="data/qc/exclusions/diff_miss_fail_variants.txt",
        freeze="data/qc/exclusions/diff_miss_freeze_fail_variants.txt"
    output:
        "data/qc/exclusions/diff_miss_both_fail_variants.txt"
    resources:
        lsf_err="logs/lsf/aggregate_diff_miss_both_fail.e",
        lsf_out="logs/lsf/aggregate_diff_miss_both_fail.o"
    shell:
        "comm -12 <(sort {input.case_control}) <(sort {input.freeze}) > {output}"

# -----------------------------------------------------------------------------
# Table for QC report plots: bias SNP locations and "loadings" along the genome
# (ID, CHR, POS, fail_cc, fail_freeze, P_cc, P_freeze, F_MISS cols → genome + loading plots)
# -----------------------------------------------------------------------------
rule prepare_diff_miss_plot_data:
    input:
        cc_fail="data/qc/exclusions/diff_miss_fail_variants.txt",
        freeze_fail="data/qc/exclusions/diff_miss_freeze_fail_variants.txt",
        cc_tests=expand("data/plink/chr{CHR}.test.missing",CHR=CHROMOSOMES_AUTOSOMAL),
        freeze_tests=expand("data/plink/chr{CHR}.test_freeze.missing",CHR=CHROMOSOMES_AUTOSOMAL),
        pvars=expand("data/plink/chr{CHR}.qc_filter2.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
    output:
        "data/qc/reports/diff_miss_bias_snps_for_plots.tsv"
    params:
        threshold=DIFF_MISS_THR
    resources:
        lsf_err="logs/lsf/prepare_diff_miss_plot_data.e",
        lsf_out="logs/lsf/prepare_diff_miss_plot_data.o"
    shell:
        """
        python prepare_diff_miss_plot_data.py \
          --cc-fail {input.cc_fail} --freeze-fail {input.freeze_fail} \
          --cc-missing {input.cc_tests} --freeze-missing {input.freeze_tests} \
          --pvar {input.pvars} --threshold {params.threshold} --output {output}
        """

# -----------------------------------------------------------------------------
# Summaries (for reporting; same threshold)
# -----------------------------------------------------------------------------
rule summarize_diff_miss:
    input:
        expand("data/plink/chr{CHR}.test.missing",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/reports/diff_miss_summary.csv"
    params:
        threshold=DIFF_MISS_THR
    resources:
        lsf_err="logs/lsf/summarize_diff_miss.e",
        lsf_out="logs/lsf/summarize_diff_miss.o"
    shell:
        "python summarize_diff_miss.py {input} --threshold {params.threshold} --output {output}"

rule summarize_diff_miss_freeze:
    input:
        expand("data/plink/chr{CHR}.test_freeze.missing",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/reports/diff_miss_freeze_summary.csv"
    params:
        threshold=DIFF_MISS_THR
    resources:
        lsf_err="logs/lsf/summarize_diff_miss_freeze.e",
        lsf_out="logs/lsf/summarize_diff_miss_freeze.o"
    shell:
        "python summarize_diff_miss.py {input} --threshold {params.threshold} --output {output}"
