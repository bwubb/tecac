import os

os.makedirs('logs/lsf',exist_ok=True)
os.makedirs('data/bcftools',exist_ok=True)
os.makedirs('data/plink',exist_ok=True)
os.makedirs('data/qc/reports',exist_ok=True)
os.makedirs('data/qc/exclusions',exist_ok=True)


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
DIFF_MISS_THR=config.get('qc',{}).get('diff_miss_thr',1e-5)

wildcard_constraints:
    CHR='[0-9XY]+'


rule select_variants:
    input:
        "data/bcftools/chrX.input.bcf",
        "data/plink/chrX.sex_update.txt",
        expand("data/qc/reports/chr{CHR}.prefilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/qc/reports/chr{CHR}.postfilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/bcftools/chr{CHR}.qc_filter1.bcf",CHR=CHROMOSOMES_AUTOSOMAL)

rule bcftools_samples_region_tags:
    input:
        bcf=lambda wildcards: BCF_INPUT[f'chr{wildcards.CHR}']
    output:
        bcf="data/bcftools/chr{CHR}.input.bcf"
    params:
        samples=config['input']['samples'],
        regions=config['input']['targets']
    shell:
        """
        bcftools view -S {params.samples} -a -R {params.regions} {input.bcf} |
        bcftools norm -m-both |
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' |
        bcftools +fill-tags - -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE |
        bcftools view -Ob -W=csi -o {output.bcf}
        """


######### INITIAL QC #########

# Convert chrX BCF to PLINK format for sex check
# Note: We need to provide unknown sex (0) for all samples to import chrX
rule bcf_to_plink_chrX:
    input:
        bcf="data/bcftools/chrX.input.bcf"
    output:
        pgen="data/plink/chrX.pgen",
        pvar="data/plink/chrX.pvar",
        psam="data/plink/chrX.psam"
    params:
        outputname="data/plink/chrX",
        temp_psam="data/plink/chrX.temp.psam"
        #add temp() for cleanup
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
        bcf="data/bcftools/chr{CHR}.input.bcf",
        exclusions="data/qc/exclusions/sexcheck_fail.txt"
    output:
        bcf="data/bcftools/chr{CHR}.prefilter.bcf",
        csi="data/bcftools/chr{CHR}.prefilter.bcf.csi",
        stats="data/qc/reports/chr{CHR}.prefilter.stats"
    shell:
        """
        bcftools view -S ^{input.exclusions} -W=csi -Ob -o {output.bcf} {input.bcf}
        bcftools stats -s - {output.bcf} > {output.stats}
        """

# Count variant types before filtering
rule count_variant_types_prefilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.prefilter.bcf",
        csi="data/bcftools/chr{CHR}.prefilter.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.prefilter.variant_types.txt"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
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
    shell:
        "bcftools stats -s - {input.bcf} > {output}"

# Count variant types after bcftools filtering
rule count_variant_types_postfilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.bcf",
        csi="data/bcftools/chr{CHR}.qc_filter1.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.postfilter.variant_types.txt"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# Convert filtered BCF to PLINK format. Psam has FID=IID (double ID) and SEX; sex updated from chrX.
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
    shell:
        """
        # psam with FID=IID and SEX=0 (all plink files use double ID; sex updated below from chrX)
        bcftools query -l {input.bcf} | awk 'BEGIN {{OFS="\\t"}} {{print $1,$1,"0"}}' | cat <(echo -e "#FID\\tIID\\tSEX") - > {params.tmp}.psam
        plink2 --bcf {input.bcf} --psam {params.tmp}.psam --vcf-half-call reference --make-pgen --out {params.tmp}
        plink2 --pfile {params.tmp} --update-sex {input.sex_file} --make-pgen --out {params.outputname}
        """