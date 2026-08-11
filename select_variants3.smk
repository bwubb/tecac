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
        expand("data/qc/reports/chr{CHR}.prefilter.variant_types.select_variants.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/qc/reports/chr{CHR}.postfilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.site-qc.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.site-qc.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.site-qc.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/bcftools/chr{CHR}.site-qc.bcf",CHR=CHROMOSOMES_AUTOSOMAL)

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
        mkdir -p data/plink
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
        mkdir -p data/qc/exclusions data/qc/reports data/plink
        : > {output.exclusions}

        # PLINK2 default .het cols (with FID): #FID IID O(HOM) E(HOM) OBS_CT F
        # F_sex = O(HOM)/OBS_CT for chrX sex inference (report Status = MALE/FEMALE/AMBIGUOUS).
        # Females are NOT auto-removed; only AMBIGUOUS are written to sexcheck_fail.
        echo "Sample_ID O_HOM E_HOM OBS_CT F_inbreed F_sex Status" > {output.report}
        awk 'NR>1 {{
            oh=$3; eh=$4; obs=$5; f_inbreed=$6;
            if (obs+0 == 0) {{
                status="AMBIGUOUS"; f_sex="NA";
                print $1, $2 >> "{output.exclusions}";
            }} else {{
                f_sex = oh / obs;
                if (f_sex > 0.8) status="MALE";
                else if (f_sex < 0.2) status="FEMALE";
                else {{
                    status="AMBIGUOUS";
                    print $1, $2 >> "{output.exclusions}";
                }}
            }}
            print $1, oh, eh, obs, f_inbreed, f_sex, status
        }}' {input} >> {output.report}

        n=$(awk 'NR>1' {output.report} | wc -l)
        if [ "$n" -eq 0 ]; then
            echo "ERROR: sex_inference.txt has header only (0 samples). Check {input}." >&2
            exit 1
        fi

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
        "data/qc/reports/chr{CHR}.prefilter.variant_types.select_variants.txt"
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
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
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
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.postfilter.stats"
    shell:
        "bcftools stats -s - {input.bcf} > {output}"

# Count variant types after bcftools filtering
rule count_variant_types_postfilter:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
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
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi",
        sex_file="data/plink/chrX.sex_update.txt"
    output:
        pgen="data/plink/chr{CHR}.site-qc.pgen",
        pvar="data/plink/chr{CHR}.site-qc.pvar",
        psam="data/plink/chr{CHR}.site-qc.psam"
    params:
        tmp="data/plink/chr{CHR}.site-qc.temp",
        outputname="data/plink/chr{CHR}.site-qc"
    shell:
        """
        # psam with FID=IID and SEX=0 (all plink files use double ID; sex updated below from chrX)
        bcftools query -l {input.bcf} | awk 'BEGIN {{OFS="\\t"}} {{print $1,$1,"0"}}' | cat <(echo -e "#FID\\tIID\\tSEX") - > {params.tmp}.psam
        plink2 --bcf {input.bcf} --psam {params.tmp}.psam --vcf-half-call reference --make-pgen --out {params.tmp}
        plink2 --pfile {params.tmp} --update-sex {input.sex_file} --make-pgen --out {params.outputname}
        """