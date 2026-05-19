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
    CHR='[0-9XY]+',
    PROJECT='[A-Za-z0-9_-]+'


rule all:
    input:
        expand("data/preprocess/{PROJECT}.build.pgen",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.pvar",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.psam",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.snplist",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.eigenvec",PROJECT=PROJECT_NAME),
        expand("data/preprocess/{PROJECT}.build.eigenval",PROJECT=PROJECT_NAME),
        expand("data/qc/reports/{PROJECT}_exwas_qc_report.html",PROJECT=PROJECT_NAME),
        #expand("data/preprocess/{PROJECT}.chr{CHR}.annotation.pgen",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        #expand("data/preprocess/{PROJECT}.chr{CHR}.annotation.pvar",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
       # expand("data/preprocess/{PROJECT}.chr{CHR}.annotation.psam",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
       # expand("data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vep.report.csv",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
       
rule finalize:
    input:
        expand("data/qc/reports/{PROJECT}_exwas_qc_report.html",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.chr{CHR}.annotation.no_sample.vep.bcf",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/{PROJECT}.chr{CHR}.annotation.no_sample.vep.report.csv",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/{PROJECT}.chr{CHR}.annotation.pgen",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/{PROJECT}.build.pgen",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.build.eigenvec",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.samples.txt",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.controls.txt",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.cases.txt",PROJECT=PROJECT_NAME)


# Convert chrX BCF to PLINK format for sex check
# Note: We need to provide unknown sex (0) for all samples to import chrX
rule bcf_to_plink_chrX:
    input:
        lambda wildcards: BCF_INPUT['chrX']
    output:
        pgen="data/plink/{PROJECT}.chrX.pgen",
        pvar="data/plink/{PROJECT}.chrX.pvar",
        psam="data/plink/{PROJECT}.chrX.psam"
    params:
        outputname="data/plink/{PROJECT}.chrX",
        temp_psam="data/plink/{PROJECT}.chrX.temp.psam"
        #add temp() for cleanup
    log:
        stderr="logs/lsf/bcf_to_plink_chrX.{PROJECT}.e",
        stdout="logs/lsf/bcf_to_plink_chrX.{PROJECT}.o"
    shell:
        """
        bcftools query -l {input} | awk '{{OFS="\\t"}} {{print $1,$1,"0"}}' | cat <(echo -e "#FID\\tIID\\tSEX") - > {params.temp_psam}
        plink2 --bcf {input} --psam {params.temp_psam} --split-par hg38 --vcf-half-call haploid --make-pgen --out {params.outputname}
        """

# Calculate X chromosome F-statistic for sex imputation
rule plink2_chrX_het:
    input:
        pgen="data/plink/{PROJECT}.chrX.pgen",
        pvar="data/plink/{PROJECT}.chrX.pvar",
        psam="data/plink/{PROJECT}.chrX.psam"
    output:
        het="data/plink/{PROJECT}.chrX.het"
    params:
        inputname="data/plink/{PROJECT}.chrX",
        outputname="data/plink/{PROJECT}.chrX"
    log:
        stderr="logs/lsf/plink2_chrX_het.{PROJECT}.e",
        stdout="logs/lsf/plink2_chrX_het.{PROJECT}.o"
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
        "data/plink/{PROJECT}.chrX.het"
    output:
        exclusions="data/qc/exclusions/{PROJECT}.sexcheck_fail.txt",
        report="data/qc/reports/{PROJECT}.sex_inference.txt",
        sex_file="data/plink/{PROJECT}.chrX.sex_update.txt"
    log:
        stderr="logs/lsf/infer_sex_from_het.{PROJECT}.e",
        stdout="logs/lsf/infer_sex_from_het.{PROJECT}.o"
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
        exclusions="data/qc/exclusions/{PROJECT}.sexcheck_fail.txt"
    output:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.prefilter.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.prefilter.bcf.csi",
        stats="data/qc/reports/{PROJECT}.chr{CHR}.prefilter.stats"
    log:
        stderr="logs/lsf/bcftools_stats_prefilter.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/bcftools_stats_prefilter.{PROJECT}.chr{CHR}.o"
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
        "data/bcftools/{PROJECT}.chr{CHR}.prefilter.bcf"
    output:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf.csi"
    log:
        stderr="logs/lsf/bcftools_filter_per_chr.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/bcftools_filter_per_chr.{PROJECT}.chr{CHR}.o"
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
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf.csi"
    output:
        "data/qc/reports/{PROJECT}.chr{CHR}.postfilter.stats"
    log:
        stderr="logs/lsf/bcftools_stats_postfilter.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/bcftools_stats_postfilter.{PROJECT}.chr{CHR}.o"
    shell:
        "bcftools stats -s - {input.bcf} > {output}"

# Count variant types before filtering
rule count_variant_types_prefilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.prefilter.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.prefilter.bcf.csi"
    output:
        "data/qc/reports/{PROJECT}.chr{CHR}.prefilter.variant_types.txt"
    log:
        stderr="logs/lsf/count_variant_types_prefilter.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/count_variant_types_prefilter.{PROJECT}.chr{CHR}.o"
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
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf.csi"
    output:
        "data/qc/reports/{PROJECT}.chr{CHR}.postfilter.variant_types.txt"
    log:
        stderr="logs/lsf/count_variant_types_postfilter.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/count_variant_types_postfilter.{PROJECT}.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# Convert filtered BCF to PLINK format
rule bcf_to_plink_per_chr:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf.csi",
        sex_file="data/plink/{PROJECT}.chrX.sex_update.txt"
    output:
        pgen="data/plink/{PROJECT}.chr{CHR}.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.psam"
    params:
        tmp="data/plink/{PROJECT}.chr{CHR}.temp",
        outputname="data/plink/{PROJECT}.chr{CHR}"
    log:
        stderr="logs/lsf/bcf_to_plink_per_chr.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/bcf_to_plink_per_chr.{PROJECT}.chr{CHR}.o"
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
        pgen="data/plink/{PROJECT}.chr{CHR}.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.psam"
    output:
        vmiss="data/plink/{PROJECT}.chr{CHR}.vmiss",
        smiss="data/plink/{PROJECT}.chr{CHR}.smiss",
        afreq="data/plink/{PROJECT}.chr{CHR}.afreq",
        hardy="data/plink/{PROJECT}.chr{CHR}.hardy"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}",
        outputname="data/plink/{PROJECT}.chr{CHR}"
    log:
        stderr="logs/lsf/plink2_missing_freq.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink2_missing_freq.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --missing --out {params.outputname}
        plink2 --pfile {params.inputname} --freq --out {params.outputname}
        plink2 --pfile {params.inputname} --hardy --out {params.outputname}
        """

# Test for differential missingness between cases and controls (requires plink 1.9)
rule plink_test_missing:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered.psam",
        controls=config.get('input',{}).get('controls','controls.EUR.txt')
    output:
        "data/plink/{PROJECT}.chr{CHR}.test.missing"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered",
        bed_prefix="data/plink/{PROJECT}.chr{CHR}.bed_temp",
        pheno_file="data/plink/{PROJECT}.chr{CHR}.pheno.txt",
        outputname="data/plink/{PROJECT}.chr{CHR}.test"
    log:
        stderr="logs/lsf/plink_test_missing.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink_test_missing.{PROJECT}.chr{CHR}.o"
    shell:
        """
        # Convert pgen to bed for plink 1.9
        plink2 --pfile {params.inputname} --make-bed --out {params.bed_prefix}
        
        # Generate pheno file from .fam (matches FID/IID format automatically)
        awk 'NR==FNR {{controls[$1]=1; next}} {{print $1, $2, ($2 in controls ? 1 : 2)}}' {input.controls} {params.bed_prefix}.fam | \
        cat <(echo -e "FID\\tIID\\tCASE") - > {params.pheno_file}
        
        # Run test-missing with plink 1.9
        plink --bfile {params.bed_prefix} --test-missing --pheno {params.pheno_file} --allow-no-sex --out {params.outputname}
        
        # Clean up temp files
        rm -f {params.bed_prefix}.bed {params.bed_prefix}.bim {params.bed_prefix}.fam {params.bed_prefix}.log {params.pheno_file}
        """

# Apply variant-level filters per chromosome
rule plink2_filter_variants:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/{PROJECT}.chr{CHR}.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.psam"
    output:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered.psam",
        log="data/plink/{PROJECT}.chr{CHR}.variants_filtered.log"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}",
        outputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered",
        geno=GENO_THR,
        maf=MAF_THR,
        hwe=HWE_THR
    log:
        stderr="logs/lsf/plink2_filter_variants.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink2_filter_variants.{PROJECT}.chr{CHR}.o"
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
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar"
    output:
        "data/qc/reports/{PROJECT}.chr{CHR}.postplink.variant_count.txt"
    log:
        stderr="logs/lsf/count_variants_postplink.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/count_variants_postplink.{PROJECT}.chr{CHR}.o"
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
        log="data/plink/{PROJECT}.chr{CHR}.variants_filtered.log"
    output:
        "data/qc/reports/{PROJECT}.chr{CHR}.plink_filter_metrics.txt"
    log:
        stderr="logs/lsf/parse_plink_filter_metrics.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/parse_plink_filter_metrics.{PROJECT}.chr{CHR}.o"
    shell:
        """
        echo "CHR FILTER VARIANTS_REMOVED" > {output}
        
        # Extract variant removal counts from log
        GENO=$(grep -oP '\\K[0-9]+(?= variants removed due to missing genotype data)' {input.log} || echo 0)
        MAF=$(grep -oP '\\K[0-9]+(?= variants removed due to allele frequency threshold)' {input.log} || echo 0)
        HWE=$(grep -oP '\\K[0-9]+(?= variants removed due to Hardy-Weinberg exact test)' {input.log} || echo 0)
        PRUNE=$(grep -oP '\\K[0-9]+(?= variants removed due to linkage disequilibrium pruning)' {input.log} || echo 0)
        
        echo "{wildcards.CHR} GENO $GENO" >> {output}
        echo "{wildcards.CHR} MAF $MAF" >> {output}
        echo "{wildcards.CHR} HWE $HWE" >> {output}
        echo "{wildcards.CHR} PRUNE $PRUNE" >> {output}
        """

# Recompute sample missingness on filtered variants (chr1-22 only)
rule plink2_sample_missing_postfilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered.psam"
    output:
        smiss="data/plink/{PROJECT}.chr{CHR}.variants_filtered.smiss"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered",
        outputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered_missing"
    log:
        stderr="logs/lsf/plink2_sample_missing_postfilter.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink2_sample_missing_postfilter.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --missing --out {params.outputname}
        mv {params.outputname}.smiss data/plink/{wildcards.PROJECT}.chr{wildcards.CHR}.variants_filtered.smiss
        """

# Aggregate sample missingness across chr1-22
rule aggregate_sample_missingness:
    input:
        expand("data/plink/{{PROJECT}}.chr{CHR}.variants_filtered.smiss",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        aggregated="data/plink/{PROJECT}.aggregated_missingness.txt",
        exclusions="data/qc/exclusions/{PROJECT}.high_missingness.txt",
        report="data/qc/reports/{PROJECT}.missingness_report.txt"
    params:
        mind_thr=MIND_THR
    log:
        stderr="logs/lsf/aggregate_sample_missingness.{PROJECT}.e",
        stdout="logs/lsf/aggregate_sample_missingness.{PROJECT}.o"
    shell:
        """
        python aggregate_missingness.py --mind-thr {params.mind_thr} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input}
        """

# Calculate heterozygosity on missingness-filtered samples
rule calculate_heterozygosity:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.variants_filtered.psam"
    output:
        het="data/plink/{PROJECT}.chr{CHR}.het"
    params:
        inputname="data/plink/{PROJECT}.chr{CHR}.variants_filtered",
        outputname="data/plink/{PROJECT}.chr{CHR}"
    log:
        stderr="logs/lsf/calculate_heterozygosity.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/calculate_heterozygosity.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --het --out {params.outputname}
        """

# Aggregate heterozygosity across chr1-22 and identify outliers
rule aggregate_heterozygosity:
    input:
        het_files=expand("data/plink/{{PROJECT}}.chr{CHR}.het",CHR=CHROMOSOMES_AUTOSOMAL),
        missingness_exclusions="data/qc/exclusions/{PROJECT}.high_missingness.txt"
    output:
        aggregated="data/plink/{PROJECT}.aggregated_het.txt",
        exclusions="data/qc/exclusions/{PROJECT}.het_outliers.txt",
        report="data/qc/reports/{PROJECT}.heterozygosity_report.txt"
    log:
        stderr="logs/lsf/aggregate_heterozygosity.{PROJECT}.e",
        stdout="logs/lsf/aggregate_heterozygosity.{PROJECT}.o"
    shell:
        """
        python aggregate_heterozygosity.py --missingness-exclusions {input.missingness_exclusions} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input.het_files}
        """

# Parse differential missingness results
rule parse_diff_miss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        "data/plink/{PROJECT}.chr{CHR}.test.missing"
    output:
        "data/qc/exclusions/{PROJECT}.chr{CHR}.diff_miss_fail.txt"
    params:
        threshold=DIFF_MISS_THR
    log:
        stderr="logs/lsf/parse_diff_miss.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/parse_diff_miss.{PROJECT}.chr{CHR}.o"
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
        pvar="data/plink/{PROJECT}.chr{CHR}.variants_filtered.pvar"
    output:
        "data/qc/{PROJECT}.chr{CHR}.passing_variants.txt"
    log:
        stderr="logs/lsf/extract_variant_ids.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/extract_variant_ids.{PROJECT}.chr{CHR}.o"
    shell:
        """
        awk 'NR>1 {{print $3}}' {input.pvar} > {output}
        """

# Count variants after differential missingness filtering
rule count_variants_postdiffmiss:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        passing="data/qc/{PROJECT}.chr{CHR}.passing_variants.txt",
        diff_miss_fail="data/qc/exclusions/{PROJECT}.chr{CHR}.diff_miss_fail.txt"
    output:
        "data/qc/reports/{PROJECT}.chr{CHR}.postdiffmiss.variant_count.txt"
    log:
        stderr="logs/lsf/count_variants_postdiffmiss.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/count_variants_postdiffmiss.{PROJECT}.chr{CHR}.o"
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
        expand("data/qc/exclusions/{{PROJECT}}.chr{CHR}.diff_miss_fail.txt",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/exclusions/{PROJECT}.diff_miss_fail_variants.txt"
    log:
        stderr="logs/lsf/aggregate_diff_miss.{PROJECT}.e",
        stdout="logs/lsf/aggregate_diff_miss.{PROJECT}.o"
    shell:
        """
        cat {input} | sort -u > {output}
        """

# Summarize differential missingness statistics
rule summarize_diff_miss:
    input:
        expand("data/plink/{{PROJECT}}.chr{CHR}.test.missing",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/reports/{PROJECT}.diff_miss_summary.csv"
    params:
        threshold=DIFF_MISS_THR
    log:
        stderr="logs/lsf/summarize_diff_miss.{PROJECT}.e",
        stdout="logs/lsf/summarize_diff_miss.{PROJECT}.o"
    shell:
        """
        python summarize_diff_miss.py {input} --threshold {params.threshold} --output {output}
        """

rule het_missingness_exclusions:
    input:
        missingness_exclusions="data/qc/exclusions/{PROJECT}.high_missingness.txt",
        het_exclusions="data/qc/exclusions/{PROJECT}.het_outliers.txt"
    output:
        exclusions="data/qc/exclusions/{PROJECT}.het_missingness_failures.txt"
    log:
        stderr="logs/lsf/het_missingness_exclusions.{PROJECT}.e",
        stdout="logs/lsf/het_missingness_exclusions.{PROJECT}.o"
    shell:
        """
        cat {input.missingness_exclusions} {input.het_exclusions} | sort -u > {output.exclusions}
        """


rule bcftools_filter2_bcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf.csi",
        passing_variants="data/qc/{PROJECT}.chr{CHR}.passing_variants.txt",
        exclusions="data/qc/exclusions/{PROJECT}.het_missingness_failures.txt"
    output:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter2.bcf",
        csi="data/bcftools/{PROJECT}.chr{CHR}.qc_filter2.bcf.csi"
    log:
        stderr="logs/lsf/bcftools_filter2_bcf.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/bcftools_filter2_bcf.{PROJECT}.chr{CHR}.o"
    shell:
        """
        bcftools view -S ^{input.exclusions} -i 'ID=@{input.passing_variants}' -W=csi -Ob -o {output.bcf}  {input.bcf}
        """

rule bcftools_annotation_set:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.bcf",
        exclusions="data/qc/exclusions/{PROJECT}.het_missingness_failures.txt"
    output:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.annotation.bcf",
    log:
        stderr="logs/lsf/bcftools_annotation_set.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/bcftools_annotation_set.{PROJECT}.chr{CHR}.o"
    shell:
        """
        bcftools view -S ^{input.exclusions} -W=csi -Ob -o {output.bcf}  {input.bcf}
        """
#Note: My previous method filtered the rare variants to make sure there was less missing genotypes
#      I may have to do that here as well.

rule bcftools_annotation_vcf:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.annotation.bcf"
    output:
        vcf="data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vcf"
    log:
        stderr="logs/lsf/bcftools_annotation_vcf.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/bcftools_annotation_vcf.{PROJECT}.chr{CHR}.o"
    shell:
        """
        bcftools view -G -Ov -o {output.vcf} {input.bcf}
        """

rule plink2_annotation_pgen:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter1.annotation.bcf",
        sex_file="data/plink/{PROJECT}.chrX.sex_update.txt"
    output:
        pgen="data/preprocess/{PROJECT}.chr{CHR}.annotation.pgen",
        pvar="data/preprocess/{PROJECT}.chr{CHR}.annotation.pvar",
        psam="data/preprocess/{PROJECT}.chr{CHR}.annotation.psam"
    params:
        output_prefix="data/preprocess/{PROJECT}.chr{CHR}.annotation"
    log:
        stderr="logs/lsf/plink2_annotation_pgen.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink2_annotation_pgen.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}
        """

rule plink2_build_pgen:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/{PROJECT}.chr{CHR}.qc_filter2.bcf",
        sex_file="data/plink/{PROJECT}.chrX.sex_update.txt"
    output:
        pgen="data/plink/{PROJECT}.chr{CHR}.build.pgen",
        pvar="data/plink/{PROJECT}.chr{CHR}.build.pvar",
        psam="data/plink/{PROJECT}.chr{CHR}.build.psam",
    params:
        output_prefix="data/plink/{PROJECT}.chr{CHR}.build"
    log:
        stderr="logs/lsf/plink2_build_pgen.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/plink2_build_pgen.{PROJECT}.chr{CHR}.o"
    shell:
        """
        plink2 --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}
        """

rule plink2_build_merge_snplist_samples:
    input:
        expand("data/plink/{{PROJECT}}.chr{CHR}.build.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        sex_file="data/plink/{PROJECT}.chrX.sex_update.txt"
    output:
        pgen="data/preprocess/{PROJECT}.build.pgen",
        pvar="data/preprocess/{PROJECT}.build.pvar",
        psam="data/preprocess/{PROJECT}.build.psam",
        snplist="data/preprocess/{PROJECT}.build.snplist",
        samples="data/qc/{PROJECT}.passing_samples.txt"
    params:
        output_prefix="data/preprocess/{PROJECT}.build"
    log:
        stderr="logs/lsf/plink2_build_merge_snplist_samples.{PROJECT}.e",
        stdout="logs/lsf/plink2_build_merge_snplist_samples.{PROJECT}.o"
    shell:
        """
        rm -f file_list.txt
        for chrom in {{2..22}}; do echo \"data/plink/{wildcards.PROJECT}.chr${{chrom}}.build.pgen data/plink/{wildcards.PROJECT}.chr${{chrom}}.build.pvar data/plink/{wildcards.PROJECT}.chr${{chrom}}.build.psam\" >> file_list.txt; done
        plink2 --pfile data/plink/{wildcards.PROJECT}.chr1.build --pmerge-list file_list.txt --make-pgen --double-id --write-snplist --out {params.output_prefix}
        
        # Explicitly update sex after merge to ensure it's preserved
        plink2 --pfile {params.output_prefix} --update-sex {input.sex_file} --make-pgen --out {params.output_prefix}
        
        awk 'NR>1 {{print $2}}' {output.psam} > {output.samples}
        """

rule calculate_pca:
    input:
        pgen="data/preprocess/{PROJECT}.build.pgen",
        pvar="data/preprocess/{PROJECT}.build.pvar",
        psam="data/preprocess/{PROJECT}.build.psam"
    output:
        eigenvec="data/preprocess/{PROJECT}.build.eigenvec",
        eigenval="data/preprocess/{PROJECT}.build.eigenval"
    params:
        inputname="data/preprocess/{PROJECT}.build",
        outputname="data/preprocess/{PROJECT}.build"
    log:
        stderr="logs/lsf/calculate_pca.{PROJECT}.e",
        stdout="logs/lsf/calculate_pca.{PROJECT}.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --pca 20 approx \
        --out {params.outputname}
        """

rule annotate_variants:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
       "data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vcf"
    output:
        "data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vep.vcf"
    log:
        stderr="logs/lsf/annotate_variants.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/annotate_variants.{PROJECT}.chr{CHR}.o"
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
        "data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vep.vcf"
    output:
        "data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vep.report.csv"
    log:
        stderr="logs/lsf/parse_vep.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/parse_vep.{PROJECT}.chr{CHR}.o"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m no_sample"

# Generate QC report
rule generate_qc_report:
    input:
        sex_report="data/qc/reports/{PROJECT}.sex_inference.txt",
        miss_report="data/qc/reports/{PROJECT}.missingness_report.txt",
        het_report="data/qc/reports/{PROJECT}.heterozygosity_report.txt",
        diff_miss_variants="data/qc/exclusions/{PROJECT}.diff_miss_fail_variants.txt",
        diff_miss_summary="data/qc/reports/{PROJECT}.diff_miss_summary.csv",
        pca_eigenvec="data/preprocess/{PROJECT}.build.eigenvec",
        pca_eigenval="data/preprocess/{PROJECT}.build.eigenval",
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

rule final_annotation:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        vcf="data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vep.vcf",
        csv="data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vep.report.csv",
        pgen="data/preprocess/{PROJECT}.chr{CHR}.annotation.pgen",
        pvar="data/preprocess/{PROJECT}.chr{CHR}.annotation.pvar",
        psam="data/preprocess/{PROJECT}.chr{CHR}.annotation.psam"
    output:
        bcf="data/final/{PROJECT}.chr{CHR}.annotation.no_sample.vep.bcf",
        csv="data/final/{PROJECT}.chr{CHR}.annotation.no_sample.vep.report.csv",
        pgen="data/final/{PROJECT}.chr{CHR}.annotation.pgen",
        pvar="data/final/{PROJECT}.chr{CHR}.annotation.pvar",
        psam="data/final/{PROJECT}.chr{CHR}.annotation.psam"
    log:
        stderr="logs/lsf/final_annotation.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/final_annotation.{PROJECT}.chr{CHR}.o"
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
        pgen="data/preprocess/{PROJECT}.build.pgen",
        pvar="data/preprocess/{PROJECT}.build.pvar",
        psam="data/preprocess/{PROJECT}.build.psam",
        snplist="data/preprocess/{PROJECT}.build.snplist",
        eigenvec="data/preprocess/{PROJECT}.build.eigenvec",
        eigenval="data/preprocess/{PROJECT}.build.eigenval"
    output:
        pgen="data/final/{PROJECT}.build.pgen",
        pvar="data/final/{PROJECT}.build.pvar",
        psam="data/final/{PROJECT}.build.psam",
        snplist="data/final/{PROJECT}.build.snplist",
        eigenvec="data/final/{PROJECT}.build.eigenvec",
        eigenval="data/final/{PROJECT}.build.eigenval"
    log:
        stderr="logs/lsf/final_build.{PROJECT}.e",
        stdout="logs/lsf/final_build.{PROJECT}.o"
    shell:
        """
        rsync -av {input.pgen} {output.pgen}
        rsync -av {input.pvar} {output.pvar}
        rsync -av {input.psam} {output.psam}
        rsync -av {input.snplist} {output.snplist}
        rsync -av {input.eigenvec} {output.eigenvec}
        rsync -av {input.eigenval} {output.eigenval}
        """
    
rule final_lists:
    input:
        samples="data/qc/{PROJECT}.passing_samples.txt",
        exclusions1="data/qc/exclusions/{PROJECT}.sexcheck_fail.txt",
        exclusions2="data/qc/exclusions/{PROJECT}.het_missingness_failures.txt"
    output:
        samples="data/final/{PROJECT}.samples.txt",
        controls="data/final/{PROJECT}.controls.txt",
        cases="data/final/{PROJECT}.cases.txt"
    params:
        controls=config.get('input',{}).get('controls','controls.txt')
    log:
        stderr="logs/lsf/final_lists.{PROJECT}.e",
        stdout="logs/lsf/final_lists.{PROJECT}.o"
    run:
        with open(input.exclusions1,'r') as e1:
            exclusions=e1.read().splitlines()
        
        with open(input.exclusions2,'r') as e2:
            exclusions.extend(e2.read().splitlines())
        
        with open(params.controls,'r') as c:
            controls=c.read().splitlines()
        
        with open(input.samples,'r') as s, open(output.samples,'w') as f1, open(output.controls,'w') as f2, open(output.cases,'w') as f3:
            for sample in s.read().splitlines():
                if sample not in exclusions:
                    f1.write(sample+'\n')
                    if sample in controls:
                        f2.write(sample+'\n')
                    else:
                        f3.write(sample+'\n')



            
