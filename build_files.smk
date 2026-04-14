# =============================================================================
# BUILD FILES: pipeline to produce final analysis dataset (build) and PCA
# =============================================================================
# Pipeline order:
#   1. Variant QC: geno, maf, hwe, SNPs only per chr → chr*.qc_filter2 (no LD on chr)
#   2. Sample missingness on qc_filter2 per chr → .smiss
#   3. Aggregate missingness + apply mind → high_missingness.txt (exclude list)
#   4. Mind-passing keep list → mind_passing_keep.txt (samples to keep)
#   5. Merge chr*.qc_filter2, keep only mind-passing samples → merge_mind (not build)
#   6. LD prune merge_mind (chr8 + long-LD excluded) → merge_mind.prune.in/.out (SNP list only)
#   7. Heterozygosity on merge_mind.pruned (pruned SNP set) → merge_mind_pruned.het
#   8. Identify het outliers from merge_mind_pruned.het → het_exclusions.txt
#   9. Build = merge_mind minus het outliers → data/preprocess/build.* (final dataset for REGENIE; full variant set, chr8 retained)
#  10. PCA on build using merge_mind.prune.in, chr 1-7,9-22, long-LD excluded (defensive) → build_pca.eigenvec/.eigenval
# =============================================================================
# Include from main Snakefile: include: "build_files.smk"
# Expects: config, CHROMOSOMES_AUTOSOMAL, GENO_THR, MAF_THR, HWE_THR, MIND_THR.
# Main workflow provides: data/plink/chr{CHR}.pgen, data/plink/chrX.sex_update.txt.
# Aux scripts: aggregate_missingness.py, het_outliers.py, identify_pca_outliers.R
# =============================================================================

CHROMOSOMES_AUTOSOMAL=list(range(1,23))
GENO_THR=config.get('qc',{}).get('geno_thr',0.01)
MAF_THR=config.get('qc',{}).get('maf_thr',0.05)
HWE_THR=config.get('qc',{}).get('hwe_thr',1e-6)
MIND_THR=config.get('qc',{}).get('mind_thr',0.05)#0.01

wildcard_constraints:
    CHR='[0-9]+'


rule build_files:
    input:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam",
        samples="data/qc/passing_samples.txt",
        snplist="data/preprocess/build.snplist",
        eigenvec="data/preprocess/build_pca.eigenvec",
        eigenval="data/preprocess/build_pca.eigenval",
        eigenvec_var="data/preprocess/build_pca.eigenvec.allele",
        pca_outliers="data/qc/exclusions/pca_outliers.txt",
        pca_outlier_plot="data/qc/reports/pca_outliers_pc1_pc2.png",
        pca_outlier_report="data/qc/reports/pca_outlier_report.txt",
        clean_eigenvec="data/preprocess/build_pca_clean.eigenvec",
        clean_eigenval="data/preprocess/build_pca_clean.eigenval",
        clean_eigenvec_var="data/preprocess/build_pca_clean.eigenvec.allele",
        prefilter_stats=expand("data/qc/reports/chr{CHR}.prefilter.variant_types.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        postfilter_stats=expand("data/qc/reports/chr{CHR}.postfilter.variant_types.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        postplink_counts=expand("data/qc/reports/chr{CHR}.postplink.variant_count.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        plink_filter_metrics=expand("data/qc/reports/chr{CHR}.plink_filter_metrics.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        variant_afreq=expand("data/plink/chr{CHR}.afreq", CHR=CHROMOSOMES_AUTOSOMAL),
        variant_vmiss=expand("data/plink/chr{CHR}.vmiss", CHR=CHROMOSOMES_AUTOSOMAL),
        variant_hardy=expand("data/plink/chr{CHR}.hardy", CHR=CHROMOSOMES_AUTOSOMAL),

# -----------------------------------------------------------------------------
# 1. Variant QC: geno, maf, hwe, SNPs only per chr. No LD pruning on chr files.
#    Output: chr{CHR}.qc_filter2.pgen/pvar/psam
# -----------------------------------------------------------------------------
rule plink2_filter_variants:
    input:
        pgen="data/plink/chr{CHR}.pgen",
        pvar="data/plink/chr{CHR}.pvar",
        psam="data/plink/chr{CHR}.psam"
    output:
        pgen="data/plink/chr{CHR}.qc_filter2.pgen",
        pvar="data/plink/chr{CHR}.qc_filter2.pvar",
        psam="data/plink/chr{CHR}.qc_filter2.psam",
        log="data/plink/chr{CHR}.qc_filter2.log"
    params:
        inputname="data/plink/chr{CHR}",
        outputname="data/plink/chr{CHR}.qc_filter2",
        geno=GENO_THR,
        maf=MAF_THR,
        hwe=HWE_THR
    resources:
        lsf_err="logs/lsf/plink2_filter_variants.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_filter_variants.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --snps-only \
        --geno {params.geno} \
        --maf {params.maf} \
        --hwe {params.hwe} midp \
        --make-pgen \
        --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# QC report inputs (prefilter, postfilter, postplink, plink_filter_metrics, afreq/vmiss/hardy).
# Same paths as report expects; built from build's pvar/log only (no main pipeline).
# -----------------------------------------------------------------------------
rule count_variant_types_prefilter_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.pvar"
    output:
        "data/qc/reports/chr{CHR}.prefilter.variant_types.txt"
    resources:
        lsf_err="logs/lsf/count_variant_types_prefilter_build.chr{CHR}.e",
        lsf_out="logs/lsf/count_variant_types_prefilter_build.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        awk -v out={output} '$0 !~ /^#/ && NF>=5 {{tot++; ref=length($4); alt=length($5); if(ref==1 && alt==1) snp++; else indel++}} END {{print "TOTAL", tot+0 >> out; print "SNP", snp+0 >> out; print "INDEL", indel+0 >> out; print "OTHER", 0 >> out}}' {input.pvar}
        """

rule count_variant_types_postfilter_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.qc_filter2.pvar"
    output:
        "data/qc/reports/chr{CHR}.postfilter.variant_types.txt"
    resources:
        lsf_err="logs/lsf/count_variant_types_postfilter_build.chr{CHR}.e",
        lsf_out="logs/lsf/count_variant_types_postfilter_build.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        awk -v out={output} '$0 !~ /^#/ && NF>=5 {{tot++; ref=length($4); alt=length($5); if(ref==1 && alt==1) snp++; else indel++}} END {{print "TOTAL", tot+0 >> out; print "SNP", snp+0 >> out; print "INDEL", indel+0 >> out; print "OTHER", 0 >> out}}' {input.pvar}
        """

rule count_variants_postplink_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.qc_filter2.pvar"
    output:
        "data/qc/reports/chr{CHR}.postplink.variant_count.txt"
    resources:
        lsf_err="logs/lsf/count_variants_postplink_build.chr{CHR}.e",
        lsf_out="logs/lsf/count_variants_postplink_build.chr{CHR}.o"
    shell:
        """
        TOTAL=$(awk '$0 !~ /^#/' {input.pvar} | wc -l)
        echo "CHR VARIANTS" > {output}
        echo "{wildcards.CHR} $TOTAL" >> {output}
        """

rule parse_plink_filter_metrics_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        log="data/plink/chr{CHR}.qc_filter2.log"
    output:
        "data/qc/reports/chr{CHR}.plink_filter_metrics.txt"
    resources:
        lsf_err="logs/lsf/parse_plink_filter_metrics_build.chr{CHR}.e",
        lsf_out="logs/lsf/parse_plink_filter_metrics_build.chr{CHR}.o"
    shell:
        """
        echo "CHR FILTER VARIANTS_REMOVED" > {output}
        GENO=$(sed -n 's/.*\\([0-9][0-9]*\\) variant(s) removed due to missing genotype data.*/\\1/p' {input.log} 2>/dev/null | head -1); GENO=${{GENO:-0}}
        MAF=$(sed -n 's/.*\\([0-9][0-9]*\\) variant(s) removed due to allele frequency threshold.*/\\1/p' {input.log} 2>/dev/null | head -1); MAF=${{MAF:-0}}
        HWE=$(sed -n 's/.*\\([0-9][0-9]*\\) variant(s) removed due to Hardy-Weinberg.*/\\1/p' {input.log} 2>/dev/null | head -1); HWE=${{HWE:-0}}
        echo "{wildcards.CHR} GENO $GENO" >> {output}
        echo "{wildcards.CHR} MAF $MAF" >> {output}
        echo "{wildcards.CHR} HWE $HWE" >> {output}
        echo "{wildcards.CHR} PRUNE 0" >> {output}
        """

rule plink2_freq_missing_hardy_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.qc_filter2.pgen",
        pvar="data/plink/chr{CHR}.qc_filter2.pvar",
        psam="data/plink/chr{CHR}.qc_filter2.psam"
    output:
        afreq="data/plink/chr{CHR}.afreq",
        vmiss="data/plink/chr{CHR}.vmiss",
        hardy="data/plink/chr{CHR}.hardy"
    params:
        inputname="data/plink/chr{CHR}.qc_filter2",
        outputname="data/plink/chr{CHR}.qc_filter2_report"
    resources:
        lsf_err="logs/lsf/plink2_freq_missing_hardy_build.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_freq_missing_hardy_build.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --freq --missing --hardy midp --out {params.outputname}
        mv {params.outputname}.afreq {output.afreq}
        mv {params.outputname}.vmiss {output.vmiss}
        mv {params.outputname}.hardy {output.hardy}
        """

# -----------------------------------------------------------------------------
# 2. Sample missingness on qc_filter2 per chr → chr{CHR}.qc_filter2.smiss
# -----------------------------------------------------------------------------
rule plink2_sample_missing_postfilter:
    input:
        pgen="data/plink/chr{CHR}.qc_filter2.pgen",
        pvar="data/plink/chr{CHR}.qc_filter2.pvar",
        psam="data/plink/chr{CHR}.qc_filter2.psam"
    output:
        smiss="data/plink/chr{CHR}.qc_filter2.smiss"
    params:
        inputname="data/plink/chr{CHR}.qc_filter2",
        outputname="data/plink/chr{CHR}.qc_filter2_missing"
    resources:
        lsf_err="logs/lsf/plink2_sample_missing_postfilter.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_sample_missing_postfilter.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --missing --out {params.outputname}
        mv {params.outputname}.smiss {output.smiss}
        """

# -----------------------------------------------------------------------------
# 3. Aggregate missingness across chr1-22; apply mind threshold.
#    Output: high_missingness.txt (samples to exclude), missingness_report.txt
# -----------------------------------------------------------------------------
rule aggregate_sample_missingness:
    input:
        expand("data/plink/chr{CHR}.qc_filter2.smiss", CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        aggregated="data/plink/aggregated_missingness.txt",
        exclusions="data/qc/exclusions/high_missingness.txt",
        report="data/qc/reports/missingness_report.txt"
    params:
        mind_thr=MIND_THR
    resources:
        lsf_err="logs/lsf/aggregate_sample_missingness.e",
        lsf_out="logs/lsf/aggregate_sample_missingness.o"
    shell:
        """
        python aggregate_missingness.py --mind-thr {params.mind_thr} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input}
        """

# -----------------------------------------------------------------------------
# 4. Mind-passing keep list: samples not in high_missingness.txt (FID IID for --keep).
#    Output: mind_passing_keep.txt
# -----------------------------------------------------------------------------
rule mind_passing_keep_samples:
    input:
        psam="data/plink/chr1.qc_filter2.psam",
        exclusions="data/qc/exclusions/high_missingness.txt"
    output:
        "data/plink/mind_passing_keep.txt"
    resources:
        lsf_err="logs/lsf/mind_passing_keep_samples.e",
        lsf_out="logs/lsf/mind_passing_keep_samples.o"
    shell:
        """
        echo "#FID IID" > {output}
        join -1 2 -2 1 -v 1 <(awk 'NR>1 {{print $1,$2}}' {input.psam} | sort -k2,2) <(sort {input.exclusions}) | sort -k2,2 >> {output}
        """

# -----------------------------------------------------------------------------
# 5. Merge chr*.qc_filter2 and keep only mind-passing samples.
#    Output: merge_mind.pgen/pvar/psam (this is NOT build; name deliberately not "build")
# -----------------------------------------------------------------------------
rule plink2_merge_mind:
    input:
        expand("data/plink/chr{CHR}.qc_filter2.pgen", CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.qc_filter2.pvar", CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.qc_filter2.psam", CHR=CHROMOSOMES_AUTOSOMAL),
        sex_file="data/plink/chrX.sex_update.txt",
        keep="data/plink/mind_passing_keep.txt"
    output:
        pgen="data/plink/merge_mind.pgen",
        pvar="data/plink/merge_mind.pvar",
        psam="data/plink/merge_mind.psam",
        snplist="data/plink/merge_mind.snplist",
        samples="data/plink/merge_mind.samples"
    params:
        merge_prefix="data/plink/merge_mind_tmp",
        output_prefix="data/plink/merge_mind"
    resources:
        lsf_err="logs/lsf/plink2_merge_mind.e",
        lsf_out="logs/lsf/plink2_merge_mind.o"
    shell:
        """
        rm -f file_list.txt
        for chrom in {{2..22}}; do echo "data/plink/chr${{chrom}}.qc_filter2.pgen data/plink/chr${{chrom}}.qc_filter2.pvar data/plink/chr${{chrom}}.qc_filter2.psam" >> file_list.txt; done
        plink2 --pfile data/plink/chr1.qc_filter2 --pmerge-list file_list.txt --make-pgen --double-id --write-snplist --out {params.merge_prefix}
        plink2 --pfile {params.merge_prefix} --update-sex {input.sex_file} --make-pgen --out {params.merge_prefix}
        plink2 --pfile {params.merge_prefix} --keep {input.keep} --make-pgen --write-snplist --out {params.output_prefix}
        awk 'NR>1 {{print $2}}' {output.psam} > {output.samples}
        rm -f {params.merge_prefix}.pgen {params.merge_prefix}.pvar {params.merge_prefix}.psam {params.merge_prefix}.log {params.merge_prefix}.snplist
        """

# -----------------------------------------------------------------------------
# 6. LD pruning on merge_mind. Produces SNP list only; does not change dataset.
#    Exclude chr8 and long-LD (HLA, 17q21.31). Output: merge_mind.prune.in, merge_mind.prune.out
# -----------------------------------------------------------------------------
rule plink2_ld_prune_merge_mind:
    input:
        pgen="data/plink/merge_mind.pgen",
        pvar="data/plink/merge_mind.pvar",
        psam="data/plink/merge_mind.psam",
        long_ld=config.get("input", {}).get("long_ld_bed", "long_ld_regions.bed")
    output:
        prune_in="data/plink/merge_mind.prune.in",
        prune_out="data/plink/merge_mind.prune.out"
    params:
        inputname="data/plink/merge_mind",
        outputname="data/plink/merge_mind"
    resources:
        lsf_err="logs/lsf/plink2_ld_prune_merge_mind.e",
        lsf_out="logs/lsf/plink2_ld_prune_merge_mind.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --exclude range {input.long_ld} \
        --chr 1-7,9-22 \
        --indep-pairwise 200 20 0.2 \
        --out {params.outputname}
        """

rule plink2_apply_ld_prune:
    input:
        pgen="data/plink/merge_mind.pgen",
        pvar="data/plink/merge_mind.pvar",
        psam="data/plink/merge_mind.psam",
        prune_in="data/plink/merge_mind.prune.in"
    output:
        pgen="data/plink/merge_mind.pruned.pgen",
        pvar="data/plink/merge_mind.pruned.pvar",
        psam="data/plink/merge_mind.pruned.psam"
    params:
        inputname="data/plink/merge_mind",
        outputname="data/plink/merge_mind.pruned"
    resources:
        lsf_err="logs/lsf/plink2_apply_ld_prune.e",
        lsf_out="logs/lsf/plink2_apply_ld_prune.o"
    shell:
        """
        plink2 --pfile {params.inputname} --extract {input.prune_in} --make-pgen --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# 7. Heterozygosity on merge_mind.pruned (pruned SNP set only; not full merge_mind).
#    Output: merge_mind_pruned.het
# -----------------------------------------------------------------------------
rule plink2_het:
    input:
        pgen="data/plink/merge_mind.pruned.pgen",
        pvar="data/plink/merge_mind.pruned.pvar",
        psam="data/plink/merge_mind.pruned.psam"
    output:
        "data/plink/merge_mind_pruned.het"
    params:
        inputname="data/plink/merge_mind.pruned",
        outputname="data/plink/merge_mind_pruned"
    resources:
        lsf_err="logs/lsf/plink2_het.e",
        lsf_out="logs/lsf/plink2_het.o"
    shell:
        """
        plink2 --pfile {params.inputname} --het --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# 8. Identify het outliers from merge_mind_pruned.het (F = inbreeding; flag beyond mean ± 3*SD).
#    Output: het_exclusions.txt (samples to exclude), het_report.txt
# -----------------------------------------------------------------------------
rule plink2_het_outliers:
    input:
        "data/plink/merge_mind_pruned.het"
    output:
        exclusions="data/qc/exclusions/het_exclusions.txt",
        report="data/qc/reports/het_report.txt"
    resources:
        lsf_err="logs/lsf/plink2_het_outliers.e",
        lsf_out="logs/lsf/plink2_het_outliers.o"
    shell:
        """
        python het_outliers.py {input} --output-exclusions {output.exclusions} --output-report {output.report}
        """

# -----------------------------------------------------------------------------
# 9. Build = merge_mind minus het outliers. Final analysis dataset for REGENIE.
#    Output: data/preprocess/build.pgen/pvar/psam, data/qc/passing_samples.txt
# -----------------------------------------------------------------------------
rule plink2_create_build:
    input:
        pgen="data/plink/merge_mind.pgen",
        pvar="data/plink/merge_mind.pvar",
        psam="data/plink/merge_mind.psam",
        het_exclusions="data/qc/exclusions/het_exclusions.txt"
    output:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam",
        samples="data/qc/passing_samples.txt",
        snplist="data/preprocess/build.snplist"
    params:
        inputname="data/plink/merge_mind",
        output_prefix="data/preprocess/build"
    resources:
        lsf_err="logs/lsf/plink2_create_build.e",
        lsf_out="logs/lsf/plink2_create_build.o"
    shell:
        """
        plink2 --pfile {params.inputname} --remove {input.het_exclusions} --make-pgen --out {params.output_prefix} --write-snplist
        awk 'NR>1 {{print $2}}' {output.psam} > {output.samples}
        """

# -----------------------------------------------------------------------------
# 10. PCA on build (already het-filtered) using merge_mind.prune.in, chr 1-7,9-22, long-LD excluded (defensive).
#     Output: build_pca.eigenvec, build_pca.eigenval (build unchanged; chr8 retained for REGENIE)
# -----------------------------------------------------------------------------
rule calculate_pca:
    input:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam",
        pruned_snps="data/plink/merge_mind.prune.in",
        long_ld=config.get("input", {}).get("long_ld_bed", "long_ld_regions.bed")
    output:
        eigenvec="data/preprocess/build_pca.eigenvec",
        eigenvec_var="data/preprocess/build_pca.eigenvec.allele",
        eigenval="data/preprocess/build_pca.eigenval"
    params:
        inputname="data/preprocess/build",
        outputname="data/preprocess/build_pca"
    resources:
        lsf_err="logs/lsf/calculate_pca.e",
        lsf_out="logs/lsf/calculate_pca.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --extract {input.pruned_snps} \
        --exclude range {input.long_ld} \
        --chr 1-7,9-22 \
        --pca 20 allele-wts \
        --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# 11. PCA outliers: PC4 > 0.3 on build_pca.eigenvec (see identify_pca_outliers.R).
#     Writes pca_outliers.txt; pca_outlier_report.txt; pca_outliers_pc1_pc2*.png
# -----------------------------------------------------------------------------
rule identify_pca_outliers:
    input:
        eigenvec="data/preprocess/build_pca.eigenvec"
    output:
        outliers="data/qc/exclusions/pca_outliers.txt",
        plot="data/qc/reports/pca_outliers_pc1_pc2.png",
        report="data/qc/reports/pca_outlier_report.txt"
    resources:
        lsf_err="logs/lsf/identify_pca_outliers.e",
        lsf_out="logs/lsf/identify_pca_outliers.o"
    shell:
        """
        Rscript identify_pca_outliers.R \
        --eigenvec {input.eigenvec} \
        --outliers {output.outliers} \
        --plot {output.plot} \
        --report {output.report}
        """

# -----------------------------------------------------------------------------
# 12. Re-run PCA on build after removing pca_outliers.txt.
#     Output prefix: build_pca_clean
# -----------------------------------------------------------------------------
rule calculate_pca_clean:
    input:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam",
        pruned_snps="data/plink/merge_mind.prune.in",
        pca_outliers="data/qc/exclusions/pca_outliers.txt",
        long_ld=config.get("input", {}).get("long_ld_bed", "long_ld_regions.bed")
    output:
        eigenvec="data/preprocess/build_pca_clean.eigenvec",
        eigenvec_var="data/preprocess/build_pca_clean.eigenvec.allele",
        eigenval="data/preprocess/build_pca_clean.eigenval"
    params:
        inputname="data/preprocess/build",
        outputname="data/preprocess/build_pca_clean"
    resources:
        lsf_err="logs/lsf/calculate_pca_clean.e",
        lsf_out="logs/lsf/calculate_pca_clean.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --extract {input.pruned_snps} \
        --exclude range {input.long_ld} \
        --chr 1-7,9-22 \
        --remove {input.pca_outliers} \
        --pca 20 allele-wts \
        --out {params.outputname}
        """
