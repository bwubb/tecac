import os

os.makedirs('logs/lsf',exist_ok=True)


# Extract configuration values
PROJECT_NAME=config['project']['name']
CHROMOSOMES_AUTOSOMAL=list(range(1,23))
CHROMOSOMES_ALL=list(range(1,23))+['X']

# Individual chromosome files (REQUIRED)
BCF_INPUT=config['input']['chromosomes']

# QC thresholds
GENO_THR=config.get('qc',{}).get('geno_thr',0.01)
MAF_THR=config.get('qc',{}).get('maf_thr',0.05)
HWE_THR=config.get('qc',{}).get('hwe_thr',1e-6)
MIND_THR=config.get('qc',{}).get('mind_thr',0.05)#0.01
DIFF_MISS_THR=config.get('qc',{}).get('diff_miss_thr',0.01)

wildcard_constraints:
    CHR='[0-9XY]+'


rule all:
    input:
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        "data/preprocess/build.pgen",
        "data/preprocess/build.pvar",
        "data/preprocess/build.psam",
        "data/preprocess/build.snplist",
        "data/preprocess/build.eigenvec",
        "data/preprocess/build.eigenval",
        "data/qc/reports/exwas_qc_report.html",
       
rule finalize:
    input:
        "data/qc/reports/exwas_qc_report.html",
        expand("data/final/chr{CHR}.annotation.no_sample.vep.bcf",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        "data/final/build.pgen",
        "data/final/build.eigenvec",
        "data/final/samples.txt",
        "data/final/controls.txt",
        "data/final/cases.txt"


# Convert chrX BCF to PLINK format for sex check
# Note: We need to provide unknown sex (0) for all samples to import chrX
rule bcf_to_plink_chrX:
    input:
        lambda wildcards: BCF_INPUT['chrX']
    output:
        pgen="data/plink/chrX.pgen",
        pvar="data/plink/chrX.pvar",
        psam="data/plink/chrX.psam"
    params:
        outputname="data/plink/chrX",
        temp_psam="data/plink/chrX.temp.psam"
        #add temp() for cleanup
    log:
        stderr="logs/lsf/bcf_to_plink_chrX.e",
        stdout="logs/lsf/bcf_to_plink_chrX.o"
    shell:
        """
        bcftools query -l {input} | awk '{{OFS="\\t"}} {{print $1,$1,"0"}}' | cat <(echo -e "#FID\\tIID\\tSEX") - > {params.temp_psam}
        plink2 --bcf {input} --psam {params.temp_psam} --split-par hg38 --vcf-half-call haploid --make-pgen --out {params.outputname}
        """

# Calculate X chromosome F-statistic for sex imputation
rule plink2_chrX_het:
    input:
        pgen="data/plink/chrX.pgen",
        pvar="data/plink/chrX.pvar",
        psam="data/plink/chrX.psam"
    output:
        het="data/plink/chrX.het"
    params:
        inputname="data/plink/chrX",
        outputname="data/plink/chrX"
    log:
        stderr="logs/lsf/plink2_chrX_het.e",
        stdout="logs/lsf/plink2_chrX_het.o"
    shell:
        """
        plink2 --pfile {params.inputname} --het --out {params.outputname}
        """

# Infer sex from X chromosome F-statistic and create exclusion list
# F > 0.8: likely male (mostly homozygous/hemizygous on X)
# F < 0.2: likely female (heterozygous on X)
# 0.2 <= F <= 0.8: ambiguous/problem
rule infer_sex_from_het:
    input:
        "data/plink/chrX.het"
    output:
        exclusions="data/qc/exclusions/sexcheck_fail.txt",
        report="data/qc/reports/sex_inference.txt",
        sex_file="data/plink/chrX.sex_update.txt"
    log:
        stderr="logs/lsf/infer_sex_from_het.e",
        stdout="logs/lsf/infer_sex_from_het.o"
    shell:
        """
        awk 'NR>1 {{
            f = $3 / $5;
            if (f >= 0.2 && f <= 0.8) {{
                print $1, $2 > "{output.exclusions}";
            }}
        }}' {input}
        
        # Create report
        echo "Sample_ID O_HOM E_HOM OBS_CT F_inbreed F_sex Status" > {output.report}
        awk 'NR>1 {{
            f_sex = $3 / $5;
            status = (f_sex > 0.8) ? "MALE" : (f_sex < 0.2) ? "FEMALE" : "AMBIGUOUS";
            print $1, $3, $4, $5, $6, f_sex, status
        }}' {input} >> {output.report}
        
        # Touch exclusions file if empty (in case no problems found)
        touch {output.exclusions}

        # Create sex update file (only for non-ambiguous samples)
        # Format: #FID IID SEX (plink2 requires FID for --update-sex)
        echo "#FID IID SEX" > {output.sex_file}
        awk 'NR>1 && $NF!="AMBIGUOUS" {{print $1, $1, ($NF=="MALE"?1:2)}}' {output.report} >> {output.sex_file}
        """

# Generate stats before bcftools filtering
rule bcftools_stats_prefilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf=lambda wildcards: BCF_INPUT[f'chr{wildcards.CHR}'],
        exclusions="data/qc/exclusions/sexcheck_fail.txt"
    output:
        bcf="data/bcftools/chr{CHR}.prefilter.bcf",
        csi="data/bcftools/chr{CHR}.prefilter.bcf.csi",
        stats="data/qc/reports/chr{CHR}.prefilter.stats"
    log:
        stderr="logs/lsf/bcftools_stats_prefilter.chr{CHR}.e",
        stdout="logs/lsf/bcftools_stats_prefilter.chr{CHR}.o"
    shell:
        """
        bcftools view -S ^{input.exclusions} -W=csi -Ob -o {output.bcf} {input.bcf}
        bcftools stats -s - {output.bcf} > {output.stats}
        """

# Apply bcftools filters and remove sex check failures
rule bcftools_filter_per_chr:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/bcftools/chr{CHR}.prefilter.bcf"
    output:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/chr{CHR}.qc_filter1.bcf.csi"
    log:
        stderr="logs/lsf/bcftools_filter_per_chr.chr{CHR}.e",
        stdout="logs/lsf/bcftools_filter_per_chr.chr{CHR}.o"
    shell:
        """
        bcftools filter -s NoHQHet -e 'COUNT(FORMAT/GT="0/1" && FORMAT/DP>=10 && FORMAT/GQ>=20 && FORMAT/VAF>0.2)=0' -m + {input} |
        bcftools filter -s LowCallRate -e 'INFO/F_MISSING > 0.05' -m + |
        bcftools view -f PASS -W=csi -Ob -o {output.bcf}
        """

# Generate stats after bcftools filtering
rule bcftools_stats_postfilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/chr{CHR}.qc_filter1.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.postfilter.stats"
    log:
        stderr="logs/lsf/bcftools_stats_postfilter.chr{CHR}.e",
        stdout="logs/lsf/bcftools_stats_postfilter.chr{CHR}.o"
    shell:
        "bcftools stats -s - {input.bcf} > {output}"

# Count variant types before filtering
rule count_variant_types_prefilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.prefilter.bcf",
        csi="data/bcftools/chr{CHR}.prefilter.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.prefilter.variant_types.txt"
    log:
        stderr="logs/lsf/count_variant_types_prefilter.chr{CHR}.e",
        stdout="logs/lsf/count_variant_types_prefilter.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# Count variant types after bcftools filtering
rule count_variant_types_postfilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/chr{CHR}.qc_filter1.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.postfilter.variant_types.txt"
    log:
        stderr="logs/lsf/count_variant_types_postfilter.chr{CHR}.e",
        stdout="logs/lsf/count_variant_types_postfilter.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# ---------------------------------------------------------------------------
# Pipeline flow (high level):
#   chrX sex check -> prefilter (drop sex fails) -> qc_filter1 (HQHet, call rate)
#   -> bcf_to_plink_per_chr -> plink2_missing_freq
#   -> plink2_filter_variants (geno/maf/hwe/LD) -> [plink_test_missing*, sample missing, het, ...]
#   -> aggregate missingness/het -> exclusions -> two branches:
#       (A) qc_filter2 + build merge + PCA  (B) het_miss VCF -> MNP pipeline -> annotation pgen + VEP
#   -> report, final_annotation, final_build, final_lists
# Requires: config input.covariates = path to file (cols: FID, IID, Status, Freeze) for diff-miss and report.
# ---------------------------------------------------------------------------

# Convert filtered BCF to PLINK format
rule bcf_to_plink_per_chr:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/chr{CHR}.qc_filter1.bcf.csi",
        sex_file="data/plink/chrX.sex_update.txt"
    output:
        pgen="data/plink/chr{CHR}.pgen",
        pvar="data/plink/chr{CHR}.pvar",
        psam="data/plink/chr{CHR}.psam"
    params:
        tmp="data/plink/chr{CHR}.temp",
        outputname="data/plink/chr{CHR}"
    log:
        stderr="logs/lsf/bcf_to_plink_per_chr.chr{CHR}.e",
        stdout="logs/lsf/bcf_to_plink_per_chr.chr{CHR}.o"
    shell:
        """
        # Convert BCF to PLINK
        plink2 --bcf {input.bcf} --vcf-half-call reference --make-pgen --out {params.tmp}
        plink2 --pfile {params.tmp} --update-sex {input.sex_file} --make-pgen --out {params.outputname}
        
        #rm -f {params.tmp}.pgen {params.tmp}.pvar {params.tmp}.psam
        """

# Calculate missingness and frequency per chromosome (chr1-22)
rule plink2_missing_freq:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.pgen",
        pvar="data/plink/chr{CHR}.pvar",
        psam="data/plink/chr{CHR}.psam"
    output:
        vmiss="data/plink/chr{CHR}.vmiss",
        smiss="data/plink/chr{CHR}.smiss",
        afreq="data/plink/chr{CHR}.afreq",
        hardy="data/plink/chr{CHR}.hardy"
    params:
        inputname="data/plink/chr{CHR}",
        outputname="data/plink/chr{CHR}"
    log:
        stderr="logs/lsf/plink2_missing_freq.chr{CHR}.e",
        stdout="logs/lsf/plink2_missing_freq.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --missing --out {params.outputname}
        plink2 --pfile {params.inputname} --freq --out {params.outputname}
        plink2 --pfile {params.inputname} --hardy --out {params.outputname}
        """

# Apply variant-level filters per chromosome (geno, maf, hwe, LD prune)
rule plink2_filter_variants:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.pgen",
        pvar="data/plink/chr{CHR}.pvar",
        psam="data/plink/chr{CHR}.psam"
    output:
        pgen="data/plink/chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/chr{CHR}.variants_filtered.pvar",
        psam="data/plink/chr{CHR}.variants_filtered.psam",
        log="data/plink/chr{CHR}.variants_filtered.log",
        prune_out="data/plink/chr{CHR}.variants_filtered.prune.out"
    params:
        inputname="data/plink/chr{CHR}",
        outputname="data/plink/chr{CHR}.variants_filtered",
        geno=GENO_THR,
        maf=MAF_THR,
        hwe=HWE_THR
    log:
        stderr="logs/lsf/plink2_filter_variants.chr{CHR}.e",
        stdout="logs/lsf/plink2_filter_variants.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --snps-only \
        --geno {params.geno} \
        --maf {params.maf} \
        --hwe {params.hwe} midp \
        --indep-pairwise 200 20 0.2 \
        --make-pgen \
        --out {params.outputname}
        """

# Count variants after plink filtering
rule count_variants_postplink:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.variants_filtered.pvar"
    output:
        "data/qc/reports/chr{CHR}.postplink.variant_count.txt"
    log:
        stderr="logs/lsf/count_variants_postplink.chr{CHR}.e",
        stdout="logs/lsf/count_variants_postplink.chr{CHR}.o"
    shell:
        """
        TOTAL=$(awk 'NR>1' {input.pvar} | wc -l)
        echo "CHR VARIANTS" > {output}
        echo "{wildcards.CHR} $TOTAL" >> {output}
        """

# Parse plink2 log to extract filter statistics
rule parse_plink_filter_metrics:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        log="data/plink/chr{CHR}.variants_filtered.log",
        prune_out="data/plink/chr{CHR}.variants_filtered.prune.out"
    output:
        "data/qc/reports/chr{CHR}.plink_filter_metrics.txt"
    log:
        stderr="logs/lsf/parse_plink_filter_metrics.chr{CHR}.e",
        stdout="logs/lsf/parse_plink_filter_metrics.chr{CHR}.o"
    shell:
        """
        echo "CHR FILTER VARIANTS_REMOVED" > {output}
        
        # Extract variant removal counts from log (plink2 wording may vary)
        GENO=$(grep -oP '\\K[0-9]+(?= variants removed due to missing genotype data)' {input.log} 2>/dev/null || echo 0)
        MAF=$(grep -oP '\\K[0-9]+(?= variants removed due to allele frequency threshold)' {input.log} 2>/dev/null || echo 0)
        HWE=$(grep -oP '\\K[0-9]+(?= variants removed due to Hardy-Weinberg exact test)' {input.log} 2>/dev/null || echo 0)
        # LD pruning: plink2 writes excluded IDs to .prune.out; use line count (no header)
        PRUNE=$(wc -l < {input.prune_out} 2>/dev/null || echo 0)
        
        echo "{wildcards.CHR} GENO $GENO" >> {output}
        echo "{wildcards.CHR} MAF $MAF" >> {output}
        echo "{wildcards.CHR} HWE $HWE" >> {output}
        echo "{wildcards.CHR} PRUNE $PRUNE" >> {output}
        """

# Test for differential missingness case vs control (plink 1.9). Needs covariates file: cols 1=FID, 2=IID, 3=Status, 4=Freeze.
rule plink_test_missing:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/chr{CHR}.variants_filtered.pvar",
        psam="data/plink/chr{CHR}.variants_filtered.psam",
        covariates=config.get('input',{}).get('covariates','')
    output:
        "data/plink/chr{CHR}.test.missing"
    params:
        inputname="data/plink/chr{CHR}.variants_filtered",
        bed_prefix="data/plink/chr{CHR}.bed_temp",
        pheno_file="data/plink/chr{CHR}.pheno.txt",
        outputname="data/plink/chr{CHR}.test"
    log:
        stderr="logs/lsf/plink_test_missing.chr{CHR}.e",
        stdout="logs/lsf/plink_test_missing.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --make-bed --out {params.bed_prefix}
        awk 'NR==1 {{next}} {{g=($3==1||$3=="Control"||$3=="control"?1:($3==2||$3=="Case"||$3=="case"?2:0)); if(g>0) print $1,$2,g}}' {input.covariates} | cat <(echo -e "FID\\tIID\\tCASE") - > {params.pheno_file}
        plink --bfile {params.bed_prefix} --test-missing --pheno {params.pheno_file} --allow-no-sex --out {params.outputname}
        rm -f {params.bed_prefix}.bed {params.bed_prefix}.bim {params.bed_prefix}.fam {params.bed_prefix}.log {params.pheno_file}
        """

# Differential missingness freeze 2 vs 3 (covariate only; same covariates file, col 4=Freeze).
rule plink_test_missing_freeze:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/chr{CHR}.variants_filtered.pvar",
        psam="data/plink/chr{CHR}.variants_filtered.psam",
        covariates=config.get('input',{}).get('covariates','')
    output:
        "data/plink/chr{CHR}.test_freeze.missing"
    params:
        inputname="data/plink/chr{CHR}.variants_filtered",
        bed_prefix="data/plink/chr{CHR}.bed_freeze_temp",
        pheno_file="data/plink/chr{CHR}.pheno_freeze.txt",
        sample_list="data/plink/chr{CHR}.freeze_sample_list.txt",
        outputname="data/plink/chr{CHR}.test_freeze"
    log:
        stderr="logs/lsf/plink_test_missing_freeze.chr{CHR}.e",
        stdout="logs/lsf/plink_test_missing_freeze.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --make-bed --out {params.bed_prefix}
        awk 'NR==1 {{next}} $4==2||$4==3 {{print $1,$2}}' {input.covariates} > {params.sample_list}
        awk 'NR==1 {{next}} $4==2||$4==3 {{g=($4==2?1:2); print $1,$2,g}}' {input.covariates} | cat <(echo -e "FID\\tIID\\tFREEZE_GROUP") - > {params.pheno_file}
        plink --bfile {params.bed_prefix} --keep {params.sample_list} --test-missing --pheno {params.pheno_file} --allow-no-sex --out {params.outputname}
        rm -f {params.bed_prefix}.bed {params.bed_prefix}.bim {params.bed_prefix}.fam {params.bed_prefix}.log {params.pheno_file} {params.sample_list}
        """

# Recompute sample missingness on filtered variants (chr1-22 only)
rule plink2_sample_missing_postfilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/chr{CHR}.variants_filtered.pvar",
        psam="data/plink/chr{CHR}.variants_filtered.psam"
    output:
        smiss="data/plink/chr{CHR}.variants_filtered.smiss"
    params:
        inputname="data/plink/chr{CHR}.variants_filtered",
        outputname="data/plink/chr{CHR}.variants_filtered_missing"
    log:
        stderr="logs/lsf/plink2_sample_missing_postfilter.chr{CHR}.e",
        stdout="logs/lsf/plink2_sample_missing_postfilter.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --missing --out {params.outputname}
        mv {params.outputname}.smiss data/plink/chr{wildcards.CHR}.variants_filtered.smiss
        """

# Aggregate sample missingness across chr1-22
rule aggregate_sample_missingness:
    input:
        expand("data/plink/chr{CHR}.variants_filtered.smiss",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        aggregated="data/plink/aggregated_missingness.txt",
        exclusions="data/qc/exclusions/high_missingness.txt",
        report="data/qc/reports/missingness_report.txt"
    params:
        mind_thr=MIND_THR
    log:
        stderr="logs/lsf/aggregate_sample_missingness.e",
        stdout="logs/lsf/aggregate_sample_missingness.o"
    shell:
        """
        python aggregate_missingness.py --mind-thr {params.mind_thr} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input}
        """

# Calculate heterozygosity on missingness-filtered samples
rule calculate_heterozygosity:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/chr{CHR}.variants_filtered.pvar",
        psam="data/plink/chr{CHR}.variants_filtered.psam"
    output:
        het="data/plink/chr{CHR}.het"
    params:
        inputname="data/plink/chr{CHR}.variants_filtered",
        outputname="data/plink/chr{CHR}"
    log:
        stderr="logs/lsf/calculate_heterozygosity.chr{CHR}.e",
        stdout="logs/lsf/calculate_heterozygosity.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --het --out {params.outputname}
        """

# Aggregate heterozygosity across chr1-22 and identify outliers
rule aggregate_heterozygosity:
    input:
        het_files=expand("data/plink/chr{CHR}.het",CHR=CHROMOSOMES_AUTOSOMAL),
        missingness_exclusions="data/qc/exclusions/high_missingness.txt"
    output:
        aggregated="data/plink/aggregated_het.txt",
        exclusions="data/qc/exclusions/het_outliers.txt",
        report="data/qc/reports/heterozygosity_report.txt"
    log:
        stderr="logs/lsf/aggregate_heterozygosity.e",
        stdout="logs/lsf/aggregate_heterozygosity.o"
    shell:
        """
        python aggregate_heterozygosity.py --missingness-exclusions {input.missingness_exclusions} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input.het_files}
        """

# Parse differential missingness results
rule parse_diff_miss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/plink/chr{CHR}.test.missing"
    output:
        "data/qc/exclusions/chr{CHR}.diff_miss_fail.txt"
    params:
        threshold=DIFF_MISS_THR
    log:
        stderr="logs/lsf/parse_diff_miss.chr{CHR}.e",
        stdout="logs/lsf/parse_diff_miss.chr{CHR}.o"
    shell:
        """
        # Extract variants with differential missingness P < threshold
        # Column 7 is the P-value in plink .test.missing output
        awk 'NR>1 && $7 < {params.threshold} {{print $2}}' {input} > {output}
        
        # Touch file if empty (no failures)
        touch {output}
        """

# Generate variant ID list from filtered plink files
# NOTE: Differential missingness is calculated for QC reporting but NOT used for filtering
rule extract_variant_ids:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.variants_filtered.pvar"
    output:
        "data/qc/chr{CHR}.passing_variants.txt"
    log:
        stderr="logs/lsf/extract_variant_ids.chr{CHR}.e",
        stdout="logs/lsf/extract_variant_ids.chr{CHR}.o"
    shell:
        """
        awk 'NR>1 {{print $3}}' {input.pvar} > {output}
        """

# Count variants after differential missingness filtering
rule count_variants_postdiffmiss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        passing="data/qc/chr{CHR}.passing_variants.txt",
        diff_miss_fail="data/qc/exclusions/chr{CHR}.diff_miss_fail.txt"
    output:
        "data/qc/reports/chr{CHR}.postdiffmiss.variant_count.txt"
    log:
        stderr="logs/lsf/count_variants_postdiffmiss.chr{CHR}.e",
        stdout="logs/lsf/count_variants_postdiffmiss.chr{CHR}.o"
    shell:
        """
        PASSING=$(wc -l < {input.passing})
        DIFF_MISS=$(wc -l < {input.diff_miss_fail})
        echo "CHR PASSING DIFF_MISS_REMOVED" > {output}
        echo "{wildcards.CHR} $PASSING $DIFF_MISS" >> {output}
        """

# Aggregate differential missingness exclusions across chromosomes
rule aggregate_diff_miss:
    input:
        expand("data/qc/exclusions/chr{CHR}.diff_miss_fail.txt",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/exclusions/diff_miss_fail_variants.txt"
    log:
        stderr="logs/lsf/aggregate_diff_miss.e",
        stdout="logs/lsf/aggregate_diff_miss.o"
    shell:
        """
        cat {input} | sort -u > {output}
        """

# Summarize differential missingness statistics
rule summarize_diff_miss:
    input:
        expand("data/plink/chr{CHR}.test.missing",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/reports/diff_miss_summary.csv"
    params:
        threshold=DIFF_MISS_THR
    log:
        stderr="logs/lsf/summarize_diff_miss.e",
        stdout="logs/lsf/summarize_diff_miss.o"
    shell:
        """
        python summarize_diff_miss.py {input} --threshold {params.threshold} --output {output}
        """

rule het_missingness_exclusions:
    input:
        missingness_exclusions="data/qc/exclusions/high_missingness.txt",
        het_exclusions="data/qc/exclusions/het_outliers.txt"
    output:
        exclusions="data/qc/exclusions/het_missingness_failures.txt"
    log:
        stderr="logs/lsf/het_missingness_exclusions.e",
        stdout="logs/lsf/het_missingness_exclusions.o"
    shell:
        """
        cat {input.missingness_exclusions} {input.het_exclusions} | sort -u > {output.exclusions}
        """

rule bcftools_filter2_bcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/chr{CHR}.qc_filter1.bcf.csi",
        passing_variants="data/qc/chr{CHR}.passing_variants.txt",
        exclusions="data/qc/exclusions/het_missingness_failures.txt"
    output:
        bcf="data/bcftools/chr{CHR}.qc_filter2.bcf",
        csi="data/bcftools/chr{CHR}.qc_filter2.bcf.csi"
    log:
        stderr="logs/lsf/bcftools_filter2_bcf.chr{CHR}.e",
        stdout="logs/lsf/bcftools_filter2_bcf.chr{CHR}.o"
    shell:
        """
        bcftools view -S ^{input.exclusions} -i 'ID=@{input.passing_variants}' -W=csi -Ob -o {output.bcf}  {input.bcf}
        """

rule bcftools_het_miss_vcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        exclusions="data/qc/exclusions/het_missingness_failures.txt"
    output:
        vcf1="data/bcftools/chr{CHR}.qc_filter1.het_miss.vcf",
        vcf2="data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vcf"
    log:
        stderr="logs/lsf/bcftools_het_miss_vcf.chr{CHR}.e",
        stdout="logs/lsf/bcftools_het_miss_vcf.chr{CHR}.o"
    shell:
        """
        bcftools view -S ^{input.exclusions} -a -Ov -o {output.vcf1}  {input.bcf}
        bcftools view -G -Ov -o {output.vcf2} {output.vcf1}
        """

rule first_variants_annotation:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
       "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vcf"
    output:
        "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.vcf"
    log:
        stderr="logs/lsf/first_variants_annotation.chr{CHR}.e",
        stdout="logs/lsf/first_variants_annotation.chr{CHR}.o"
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

rule parse_first_vep:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.vcf"
    output:
        "data/bcftools/chr{CHR}.qc_filter1.het_miss.no_sample.vep.report.csv"
    log:
        stderr="logs/lsf/parse_first_vep.chr{CHR}.e",
        stdout="logs/lsf/parse_first_vep.chr{CHR}.o"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m no_sample"

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
    log:
        stderr="logs/lsf/find_mnp_variants.chr{CHR}.e",
        stdout="logs/lsf/find_mnp_variants.chr{CHR}.o"
    shell:
        "python find_mnp_variants.py -i {input.csv} -o {output.csv} -m no_sample -r {output.regions} -I {output.id_file}"

rule bcftools_mnp_filter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.vcf",
        regions="data/mnp/chr{CHR}.regions.txt",
        id_file="data/mnp/chr{CHR}.id.txt"
    output:
        vcf="data/mnp/chr{CHR}.mnp.vcf"
    log:
        stderr="logs/lsf/bcftools_mnp_filter.chr{CHR}.e",
        stdout="logs/lsf/bcftools_mnp_filter.chr{CHR}.o"
    shell:
        """
        bcftools view -R {input.regions} -i 'ID=@{input.id_file}' -Ov -o {output.vcf} {input.vcf}
        """

rule bcftools_mnp_sample_info:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/mnp/chr{CHR}.mnp.vcf"
    output:
        "data/mnp/chr{CHR}.mnp_sample_info.txt"
    log:
        stderr="logs/lsf/bcftools_mnp_sample_info.chr{CHR}.e",
        stdout="logs/lsf/bcftools_mnp_sample_info.chr{CHR}.o"
    shell:
        """
        bcftools query -f '[%ID %SAMPLE %GT %AD %DP %VAF\n]' {input} |
        awk -F' ' '$3=="0/1" || $3=="1/1"' > {output}
        """

rule test_mnp_sample_info:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pairs="data/mnp/chr{CHR}.annotation.mnp.pairs.txt",
        sample_info="data/mnp/chr{CHR}.mnp_sample_info.txt"
    output:
        pass_out="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt",
        fail_out="data/mnp/chr{CHR}.mnp_sample_info.FAIL.txt"
    log:
        stderr="logs/lsf/test_mnp_sample_info.chr{CHR}.e",
        stdout="logs/lsf/test_mnp_sample_info.chr{CHR}.o"
    shell:
        "python test_mnp_sample_info.py -p {input.pairs} -i {input.sample_info} -o {output.pass_out} -f {output.fail_out}"

rule plan_mnp_gt:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        sample_info="data/mnp/chr{CHR}.mnp_sample_info.PASS.txt"
    output:
        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt"
    params:
        ref="/home/bwubb/resources/Genomes/Human/hg38/Homo_sapiens_assembly38.fasta.gz"
    log:
        stderr="logs/lsf/plan_mnp_gt.chr{CHR}.e",
        stdout="logs/lsf/plan_mnp_gt.chr{CHR}.o"
    shell:
        "python plan_mnp_gt.py -i {input.sample_info} -o {output.plan} -r {params.ref}"

rule manage_mnp_gt:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        plan="data/mnp/chr{CHR}.mnp_gt.plan.txt",
        vcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.vcf"
    output:
        vcf="data/mnp/chr{CHR}.qc_filter1.het_miss.mnp_gt.vcf"
    log:
        stderr="logs/lsf/manage_mnp_gt.chr{CHR}.e",
        stdout="logs/lsf/manage_mnp_gt.chr{CHR}.o"
    shell:
        "python manage_mnp_gt.py -p {input.plan} -i {input.vcf} -o {output.vcf}"

rule bcftools_sort_mnp_gt_bcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/mnp/chr{CHR}.qc_filter1.het_miss.mnp_gt.vcf"
    output:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.mnp_gt.sorted.bcf"
    log:
        stderr="logs/lsf/bcftools_sort_mnp_gt_bcf.chr{CHR}.e",
        stdout="logs/lsf/bcftools_sort_mnp_gt_bcf.chr{CHR}.o"
    shell:
        """
        bcftools sort -W=csi -Ob -o {output.bcf}  {input.vcf}
        """

rule bcftools_annotation2_vcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.mnp_gt.sorted.bcf"
    output:
        vcf="data/preprocess/chr{CHR}.annotation.no_sample.vcf"
    log:
        stderr="logs/lsf/bcftools_annotation2_vcf.chr{CHR}.e",
        stdout="logs/lsf/bcftools_annotation2_vcf.chr{CHR}.o"
    shell:
        """
        bcftools view -G -Ov -o {output.vcf} {input.bcf}
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
    log:
        stderr="logs/lsf/plink2_annotation_pgen.chr{CHR}.e",
        stdout="logs/lsf/plink2_annotation_pgen.chr{CHR}.o"
    shell:
        """
        plink2 --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}
        """

rule plink2_build_pgen:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter2.bcf",
        sex_file="data/plink/chrX.sex_update.txt"
    output:
        pgen="data/plink/chr{CHR}.build.pgen",
        pvar="data/plink/chr{CHR}.build.pvar",
        psam="data/plink/chr{CHR}.build.psam",
    params:
        output_prefix="data/plink/chr{CHR}.build"
    log:
        stderr="logs/lsf/plink2_build_pgen.chr{CHR}.e",
        stdout="logs/lsf/plink2_build_pgen.chr{CHR}.o"
    shell:
        """
        plink2 --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}
        """

rule plink2_build_merge_snplist_samples:
    input:
        expand("data/plink/chr{CHR}.build.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        sex_file="data/plink/chrX.sex_update.txt"
    output:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam",
        snplist="data/preprocess/build.snplist",
        samples="data/qc/passing_samples.txt"
    params:
        output_prefix="data/preprocess/build"
    log:
        stderr="logs/lsf/plink2_build_merge_snplist_samples.e",
        stdout="logs/lsf/plink2_build_merge_snplist_samples.o"
    shell:
        """
        rm -f file_list.txt
        for chrom in {{2..22}}; do echo \"data/plink/chr${{chrom}}.build.pgen data/plink/chr${{chrom}}.build.pvar data/plink/chr${{chrom}}.build.psam\" >> file_list.txt; done
        plink2 --pfile data/plink/chr1.build --pmerge-list file_list.txt --make-pgen --double-id --write-snplist --out {params.output_prefix}
        
        # Explicitly update sex after merge to ensure it's preserved
        plink2 --pfile {params.output_prefix} --update-sex {input.sex_file} --make-pgen --out {params.output_prefix}
        
        awk 'NR>1 {{print $2}}' {output.psam} > {output.samples}
        """

rule calculate_pca:
    input:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam"
    output:
        eigenvec="data/preprocess/build.eigenvec",
        eigenval="data/preprocess/build.eigenval"
    params:
        inputname="data/preprocess/build",
        outputname="data/preprocess/build"
    log:
        stderr="logs/lsf/calculate_pca.e",
        stdout="logs/lsf/calculate_pca.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --chr 1-7,9-22 \
        --pca 20 approx \
        --out {params.outputname}
        """

rule annotate_variants:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
       "data/preprocess/chr{CHR}.annotation.no_sample.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    log:
        stderr="logs/lsf/annotate_variants.chr{CHR}.e",
        stdout="logs/lsf/annotate_variants.chr{CHR}.o"
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

rule parse_vep:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv"
    log:
        stderr="logs/lsf/parse_vep.chr{CHR}.e",
        stdout="logs/lsf/parse_vep.chr{CHR}.o"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m no_sample"

#OK WE HAVE A MAJOR REVISION TO MAKE HERE. I HAVE CODE AND STUFF FOR MNP
#we need to go back to qc filter 1.annoation.bcf and we will rename that. we need to
#take that , annotate it, and then look for mnps and do that mnp pipeline, which I have code for.
#Then once we write new vcf files that have adjusted genotypes for the mnp and variants that make it up
#We will make the annotation pgen, etc from that.
#and annotate that one.




# Generate QC report
rule generate_qc_report:
    input:
        sex_report="data/qc/reports/sex_inference.txt",
        miss_report="data/qc/reports/missingness_report.txt",
        het_report="data/qc/reports/heterozygosity_report.txt",
        diff_miss_variants="data/qc/exclusions/diff_miss_fail_variants.txt",
        diff_miss_summary="data/qc/reports/diff_miss_summary.csv",
        pca_eigenvec="data/preprocess/build.eigenvec",
        pca_eigenval="data/preprocess/build.eigenval",
        prefilter_stats=expand("data/qc/reports/chr{CHR}.prefilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postfilter_stats=expand("data/qc/reports/chr{CHR}.postfilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        plink_filter_metrics=expand("data/qc/reports/chr{CHR}.plink_filter_metrics.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postplink_counts=expand("data/qc/reports/chr{CHR}.postplink.variant_count.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postdiffmiss_stats=expand("data/qc/reports/chr{CHR}.postdiffmiss.variant_count.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_afreq=expand("data/plink/chr{CHR}.afreq",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_vmiss=expand("data/plink/chr{CHR}.vmiss",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_hardy=expand("data/plink/chr{CHR}.hardy",CHR=CHROMOSOMES_AUTOSOMAL),
        diff_miss_tests=expand("data/plink/chr{CHR}.test.missing",CHR=CHROMOSOMES_AUTOSOMAL),
        diff_miss_freeze_tests=expand("data/plink/chr{CHR}.test_freeze.missing",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/reports/exwas_qc_report.html"
    params:
        project=PROJECT_NAME,
        data_dir="data",
        covariate_file=config.get('input',{}).get('covariates','')
    log:
        stderr="logs/lsf/generate_qc_report.e",
        stdout="logs/lsf/generate_qc_report.o"
    shell:
        """
        Rscript -e "rmarkdown::render('exwas_qc_report.Rmd', output_file='{output}', params=list(project='{params.project}', data_dir='{params.data_dir}', covariate_file='{params.covariate_file}'))"
        DATE=$(date +%Y%m%d)
        cp {output} exwas_qc_report.${{DATE}}.html
        """

rule final_annotation:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf",
        csv="data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",
        pgen="data/preprocess/chr{CHR}.annotation.pgen",
        pvar="data/preprocess/chr{CHR}.annotation.pvar",
        psam="data/preprocess/chr{CHR}.annotation.psam"
    output:
        bcf="data/final/chr{CHR}.annotation.no_sample.vep.bcf",
        csv="data/final/chr{CHR}.annotation.no_sample.vep.report.csv",
        pgen="data/final/chr{CHR}.annotation.pgen",
        pvar="data/final/chr{CHR}.annotation.pvar",
        psam="data/final/chr{CHR}.annotation.psam"
    log:
        stderr="logs/lsf/final_annotation.chr{CHR}.e",
        stdout="logs/lsf/final_annotation.chr{CHR}.o"
    shell:
        """
        bcftools view -Ob -o {output.bcf} {input.vcf}
        rsync -av {input.csv} {output.csv}
        rsync -av {input.pgen} {output.pgen}
        rsync -av {input.pvar} {output.pvar}
        rsync -av {input.psam} {output.psam}
        """
rule final_build:
    input:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam",
        snplist="data/preprocess/build.snplist",
        eigenvec="data/preprocess/build.eigenvec",
        eigenval="data/preprocess/build.eigenval"
    output:
        pgen="data/final/build.pgen",
        pvar="data/final/build.pvar",
        psam="data/final/build.psam",
        snplist="data/final/build.snplist",
        eigenvec="data/final/build.eigenvec",
        eigenval="data/final/build.eigenval"
    log:
        stderr="logs/lsf/final_build.e",
        stdout="logs/lsf/final_build.o"
    shell:
        """
        rsync -av {input.pgen} {output.pgen}
        rsync -av {input.pvar} {output.pvar}
        rsync -av {input.psam} {output.psam}
        rsync -av {input.snplist} {output.snplist}
        rsync -av {input.eigenvec} {output.eigenvec}
        rsync -av {input.eigenval} {output.eigenval}
        """
    
# Splits passing samples into controls vs cases. Prefers covariate file (col3=Status); else separate controls file.
rule final_lists:
    input:
        samples="data/qc/passing_samples.txt",
        exclusions1="data/qc/exclusions/sexcheck_fail.txt",
        exclusions2="data/qc/exclusions/het_missingness_failures.txt"
    output:
        samples="data/final/samples.txt",
        controls="data/final/controls.txt",
        cases="data/final/cases.txt"
    params:
        covariates=config.get('input',{}).get('covariates',''),
        controls_file=config.get('input',{}).get('controls','')
    log:
        stderr="logs/lsf/final_lists.e",
        stdout="logs/lsf/final_lists.o"
    run:
        with open(input.exclusions1,'r') as e1:
            exclusions = e1.read().splitlines()
        with open(input.exclusions2,'r') as e2:
            exclusions.extend(e2.read().splitlines())

        # Controls: from covariate file (col2=IID, col3=Status) if present, else separate controls file
        controls = set()
        if params.covariates and os.path.isfile(params.covariates):
            with open(params.covariates,'r') as f:
                lines = f.readlines()
            if len(lines) > 1:
                for line in lines[1:]:
                    parts = line.strip().split()
                    if len(parts) >= 3 and parts[2] in ('1', 'Control', 'control', 'CONTROL'):
                        controls.add(parts[1])  # IID
        elif params.controls_file and os.path.isfile(params.controls_file):
            with open(params.controls_file,'r') as c:
                controls = set(s.strip() for s in c.read().splitlines())

        with open(input.samples,'r') as s, open(output.samples,'w') as f1, open(output.controls,'w') as f2, open(output.cases,'w') as f3:
            for sample in s.read().splitlines():
                if sample not in exclusions:
                    f1.write(sample + '\n')
                    if sample in controls:
                        f2.write(sample + '\n')
                    else:
                        f3.write(sample + '\n')



            
