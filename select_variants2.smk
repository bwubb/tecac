# select_variants2.smk
# - Applies differential missingness as a VARIANT filter (passing_variants excludes diff_miss failures)
# - Reports Path A vs Path B sample exclusion comparison:
#   Path A: MIND+HET computed on full variants_filtered (current order) -> N samples excluded
#   Path B: diff_miss filter first, then MIND+HET on variants_filtered_nodiffmiss -> N samples excluded
# Does NOT change which samples we exclude (still use Path A for now); report is for comparison only.

include: "select_variants.smk"

# Override rule all to require the comparison report
rule all:
    input:
        expand("data/preprocess/{PROJECT}.build.pgen",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.pvar",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.psam",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.snplist",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.eigenvec",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.eigenval",PROJECT=PROJECT_NAME),
        expand("data/qc/reports/{PROJECT}_exwas_qc_report.html",PROJECT=PROJECT_NAME),
        expand("data/qc/reports/{PROJECT}.sample_exclusion_comparison.txt",PROJECT=PROJECT_NAME),

# --- Override: Apply diff_miss as variant filter (passing_variants = pvar minus diff_miss_fail) ---
rule extract_variant_ids:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar",
        diff_miss_fail="data/qc/exclusions/{PROJECT}.chr{CHR}.diff_miss_fail.txt"
    output:
        "data/qc/{PROJECT}.chr{CHR}.passing_variants.txt"
    log:
        stderr="logs/lsf/extract_variant_ids.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/extract_variant_ids.{PROJECT}.chr{CHR}.o"
    shell:
        """
        awk 'NR>1 {{print $3}}' {input.pvar} | sort -u > tmp_pvar.ids
        sort -u {input.diff_miss_fail} 2>/dev/null > tmp_diffmiss.ids || touch tmp_diffmiss.ids
        comm -23 tmp_pvar.ids tmp_diffmiss.ids > {output}
        rm -f tmp_pvar.ids tmp_diffmiss.ids
        """

# --- Path B: variants_filtered with diff_miss removed, then MIND+HET on that set ---
rule plink2_filter_diffmiss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered.psam",
        diff_miss_fail="data/qc/exclusions/{PROJECT}.chr{CHR}.diff_miss_fail.txt"
    output:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.psam"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered",
        outputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss"
    log:
        stderr="logs/lsf/plink2_filter_diffmiss.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink2_filter_diffmiss.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --exclude {input.diff_miss_fail} --make-pgen --out {params.outputname}
        """

rule plink2_sample_missing_nodiffmiss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.psam"
    output:
        smiss="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.smiss"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss",
        outputname="data/plink/{PROJECT}.chr{CHR}.nodiffmiss_missing"
    log:
        stderr="logs/lsf/plink2_sample_missing_nodiffmiss.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink2_sample_missing_nodiffmiss.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --missing --out {params.outputname}
        mv {params.outputname}.smiss {output.smiss}
        """

rule aggregate_sample_missingness_nodiffmiss:
    input:
        expand("data/plink/{{PROJECT}}.chr{CHR}.variants_filtered_nodiffmiss.smiss",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        aggregated="data/plink/{PROJECT}.aggregated_missingness_nodiffmiss.txt",
        exclusions="data/qc/exclusions/{PROJECT}.high_missingness_nodiffmiss.txt",
        report="data/qc/reports/{PROJECT}.missingness_report_nodiffmiss.txt"
    params:
        mind_thr=MIND_THR
    log:
        stderr="logs/lsf/aggregate_sample_missingness_nodiffmiss.{PROJECT}.e",
        stdout="logs/lsf/aggregate_sample_missingness_nodiffmiss.{PROJECT}.o"
    shell:
        """
        python aggregate_missingness.py --mind-thr {params.mind_thr} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input}
        """

rule calculate_heterozygosity_nodiffmiss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss.psam"
    output:
        het="data/plink/{PROJECT}.chr{CHR}.nodiffmiss.het"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered_nodiffmiss",
        outputname="data/plink/{PROJECT}.chr{CHR}.nodiffmiss"
    log:
        stderr="logs/lsf/calculate_heterozygosity_nodiffmiss.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/calculate_heterozygosity_nodiffmiss.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --het --out {params.outputname}
        """

rule aggregate_heterozygosity_nodiffmiss:
    input:
        het_files=expand("data/plink/{{PROJECT}}.chr{CHR}.nodiffmiss.het",CHR=CHROMOSOMES_AUTOSOMAL),
        missingness_exclusions="data/qc/exclusions/{PROJECT}.high_missingness_nodiffmiss.txt"
    output:
        aggregated="data/plink/{PROJECT}.aggregated_het_nodiffmiss.txt",
        exclusions="data/qc/exclusions/{PROJECT}.het_outliers_nodiffmiss.txt",
        report="data/qc/reports/{PROJECT}.heterozygosity_report_nodiffmiss.txt"
    log:
        stderr="logs/lsf/aggregate_heterozygosity_nodiffmiss.{PROJECT}.e",
        stdout="logs/lsf/aggregate_heterozygosity_nodiffmiss.{PROJECT}.o"
    shell:
        """
        python aggregate_heterozygosity.py --missingness-exclusions {input.missingness_exclusions} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input.het_files}
        """

rule het_missingness_exclusions_nodiffmiss:
    input:
        missingness_exclusions="data/qc/exclusions/{PROJECT}.high_missingness_nodiffmiss.txt",
        het_exclusions="data/qc/exclusions/{PROJECT}.het_outliers_nodiffmiss.txt"
    output:
        exclusions="data/qc/exclusions/{PROJECT}.het_missingness_failures_nodiffmiss.txt"
    log:
        stderr="logs/lsf/het_missingness_exclusions_nodiffmiss.{PROJECT}.e",
        stdout="logs/lsf/het_missingness_exclusions_nodiffmiss.{PROJECT}.o"
    shell:
        """
        cat {input.missingness_exclusions} {input.het_exclusions} | sort -u > {output.exclusions}
        """

rule sample_exclusion_comparison:
    input:
        path_a="data/qc/exclusions/{PROJECT}.het_missingness_failures.txt",
        path_b="data/qc/exclusions/{PROJECT}.het_missingness_failures_nodiffmiss.txt"
    output:
        "data/qc/reports/{PROJECT}.sample_exclusion_comparison.txt"
    log:
        stderr="logs/lsf/sample_exclusion_comparison.{PROJECT}.e",
        stdout="logs/lsf/sample_exclusion_comparison.{PROJECT}.o"
    shell:
        """
        echo "Sample exclusion comparison: MIND+HET thresholds" > {output}
        echo "================================================" >> {output}
        echo "" >> {output}
        echo "Path A (current): MIND+HET computed on full variants_filtered (before diff_miss filter)" >> {output}
        echo "  Samples excluded:" $(wc -l < {input.path_a}) >> {output}
        echo "" >> {output}
        echo "Path B (diff_miss first): diff_miss variant filter applied, then MIND+HET on remaining variants" >> {output}
        echo "  Samples excluded:" $(wc -l < {input.path_b}) >> {output}
        """

# Override generate_qc_report to require sample exclusion comparison (so report shows it)
rule generate_qc_report:
    input:
        sex_report="data/qc/reports/{PROJECT}.sex_inference.txt",
        miss_report="data/qc/reports/{PROJECT}.missingness_report.txt",
        het_report="data/qc/reports/{PROJECT}.heterozygosity_report.txt",
        diff_miss_variants="data/qc/exclusions/{PROJECT}.diff_miss_fail_variants.txt",
        diff_miss_summary="data/qc/reports/{PROJECT}.diff_miss_summary.csv",
        pca_eigenvec="data/preprocess/{PROJECT}.build.eigenvec",
        pca_eigenval="data/preprocess/{PROJECT}.build.eigenval",
        sample_exclusion_comparison="data/qc/reports/{PROJECT}.sample_exclusion_comparison.txt",
        prefilter_stats=expand("data/qc/reports/{{PROJECT}}.chr{CHR}.prefilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postfilter_stats=expand("data/qc/reports/{{PROJECT}}.chr{CHR}.postfilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        plink_filter_metrics=expand("data/qc/reports/{{PROJECT}}.chr{CHR}.plink_filter_metrics.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postplink_counts=expand("data/qc/reports/{{PROJECT}}.chr{CHR}.postplink.variant_count.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postdiffmiss_stats=expand("data/qc/reports/{{PROJECT}}.chr{CHR}.postdiffmiss.variant_count.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_afreq=expand("data/plink/{{PROJECT}}.chr{CHR}.afreq",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_vmiss=expand("data/plink/{{PROJECT}}.chr{CHR}.vmiss",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_hardy=expand("data/plink/{{PROJECT}}.chr{CHR}.hardy",CHR=CHROMOSOMES_AUTOSOMAL),
        diff_miss_tests=expand("data/plink/{{PROJECT}}.chr{CHR}.test.missing",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/reports/{PROJECT}_exwas_qc_report.html"
    params:
        project="{PROJECT}",
        controls=config.get('input',{}).get('controls','controls.EUR.txt'),
        data_dir="data"
    log:
        stderr="logs/lsf/generate_qc_report.{PROJECT}.e",
        stdout="logs/lsf/generate_qc_report.{PROJECT}.o"
    shell:
        """
        Rscript -e "rmarkdown::render('exwas_qc_report.Rmd', output_file='{output}', params=list(project='{params.project}', controls_file='{params.controls}', data_dir='{params.data_dir}'))"
        DATE=$(date +%Y%m%d)
        cp {output} {wildcards.PROJECT}_exwas_qc_report.${{DATE}}.html
        """
