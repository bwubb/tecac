from datetime import datetime
import os

os.makedirs("logs/cluster/regenie2",exist_ok=True)

# Extract configuration values
PROJECT_NAME=config['project']['name']
CHROMOSOMES=list(range(1,23)) # Chromosomes 1-22

# Determine VCF input method
if 'all_chr' in config['input'] and config['input']['all_chr']:
    # Single all-chromosome file
    VCF_INPUT=config['input']['all_chr']
    USE_INDIVIDUAL_CHR=False
elif 'chromosomes' in config['input'] and config['input']['chromosomes']:
    # Individual chromosome files
    VCF_INPUT=config['input']['chromosomes']
    USE_INDIVIDUAL_CHR=True
else:
    raise ValueError("Either 'all_chr' or 'chromosomes' must be specified in config")

date=datetime.now().strftime("%Y%m%d")

localrules:merge_chr_snplist,create_vep_list

rule run_regenie:
    input:
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.bcf",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT_NAME),
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.filtered.pgen",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT_NAME),
        f"data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.pgen",
        f"{PROJECT_NAME}.regenie.annotation.txt",
        f"{PROJECT_NAME}.regenie.set.txt",
        f"{PROJECT_NAME}.regenie.mask.txt",
        f"{PROJECT_NAME}.regenie.covar.txt",
        f"{PROJECT_NAME}.regenie.pheno.txt",
        f"{PROJECT_NAME}.step1_1.loco.gz",
        f"{PROJECT_NAME}.step2_single_variant_STATUS.regenie",
        f"{PROJECT_NAME}.step2_gene_based_STATUS.regenie"

rule run_report:
    input:
        expand("TECAC_{date}_regenie_report.html",date=datetime.now().strftime("%Y%m%d"))


rule preprocess_vcf:
    input:
        lambda wildcards: VCF_INPUT[f"chr{wildcards.CHR}"] if USE_INDIVIDUAL_CHR else VCF_INPUT
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.bcf"
    params:
        samples=config['input']['samples']
    shell:
        "bcftools view -f 'PASS' -S {params.samples} {input} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ob -o {output}"

#--vcf-half-call <mode>

#The current VCF standard does not specify how '0/.' and similar GT values should be interpreted.
#By default (mode 'error'/'e'), PLINK 1.9 errors out and reports the line number of the anomaly.
#Should the half-call be intentional, though (this can be the case with Complete Genomics data),
#you can request the following other modes:
#
#    'haploid'/'h': Treat half-calls as haploid/homozygous (the PLINK 1 file format does not distinguish between the two). This maximizes similarity between the VCF and BCF2 parsers.
#    'missing'/'m': Treat half-calls as missing.
#    'reference'/'r': Treat the missing part as reference.

rule chr_variant_list:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.bcf"
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.filtered.pgen",
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.filtered.snplist"
    params:
        output="data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.filtered"
    shell:
        "plink2 --bcf {input} --vcf-half-call r --maf 0.05 --snps-only --geno 0.1 --hwe 1e-15 --indep-pairwise 500 50 0.4 --write-snplist --make-pgen --out {params.output}"

rule merge_chr_pgen:
    input:
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.filtered.pgen",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT_NAME)
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.pgen",
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.eigenvec"
    params:
        output="data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered"
    shell:
        """
        rm -f file_list.txt
        for chrom in {{2..22}}; do echo \"data/work/regenie/preprocess/{PROJECT_NAME}.chr${{chrom}}.filtered.pgen data/work/regenie/preprocess/{PROJECT_NAME}.chr${{chrom}}.filtered.pvar data/work/regenie/preprocess/{PROJECT_NAME}.chr${{chrom}}.filtered.psam\" >> file_list.txt; done
        plink2 --pfile data/work/regenie/preprocess/{PROJECT_NAME}.chr1.filtered --pmerge-list file_list.txt --make-pgen --pca --out {params.output}
        """

rule merge_chr_snplist:
    input:
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.filtered.snplist",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT_NAME)
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.snplist"
    run:
        with open(output[0],'w') as f:
            for i in range(1,23):
                with open(f"data/work/regenie/preprocess/{wildcards.PROJECT_NAME}.chr{i}.filtered.snplist",'r') as f2:
                    f.write(f2.read())

rule create_vep_list:
    output:
        "vep_files.list"
    run:
        with open(output[0],'w') as f:
            for i in config['input']['vep_csv']:
                f.write(i+'\n')

rule preprocess_regenie:
    input:
        vep_list="vep_files.list",
        eigenvec="data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.eigenvec",
        sites=config['input']['sites'],
        controls=config['input']['controls']
    output:
        annotation="{PROJECT_NAME}.regenie.annotation.txt",
        set_file="{PROJECT_NAME}.regenie.set.txt",
        mask="{PROJECT_NAME}.regenie.mask.txt",
        covar="{PROJECT_NAME}.regenie.covar.txt",
        pheno="{PROJECT_NAME}.regenie.pheno.txt"
    params:
        output_prefix=PROJECT_NAME
    shell:
        """
        COVAR_ARG=""
        if [ -n "{config[input][covariates]}" ]; then
            COVAR_ARG="--covariates {config[input][covariates]}"
        fi

        python preprocess_regenie.py \
        --vep-list-file {input.vep_list} \
        --eigenvec {input.eigenvec} \
        --sites {input.sites} \
        --controls {input.controls} \
        $COVAR_ARG \
        -O {params.output_prefix}
        """

rule run_step1_regenie:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.snplist"
    output:
        "{PROJECT_NAME}.step1_1.loco.gz",
        "{PROJECT_NAME}.step1_pred.list"
    params:
        input_basename="data/work/regenie/preprocess/{PROJECT_NAME}.all_chr",
        output_basename="{PROJECT_NAME}.step1",
        lowmem_prefix="data/work/regenie/tmp_rg_"
    shell:
        "regenie --step 1 --pgen {params.input_basename} --covarFile {input[1]} --phenoFile {input[2]} --extract {input[3]} "
        "--bsize 1000 --gz --bt --lowmem --lowmem-prefix {params.lowmem_prefix} --out {params.output_basename}"
        #--strict --loocv

rule run_step2_single_variant:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "{PROJECT_NAME}.step1_pred.list"
    output:
        "{PROJECT_NAME}.step2_single_variant_STATUS.regenie"
    params:
        input_basename="data/work/regenie/preprocess/{PROJECT_NAME}.all_chr",
        output_basename="{PROJECT_NAME}.step2_single_variant"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc "
        "--out {params.output_basename}"
        #--strict --minMAC 10 --test additive

rule run_step2_gene_based:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.filtered.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "{PROJECT_NAME}.regenie.annotation.txt",
        "{PROJECT_NAME}.regenie.set.txt",
        "{PROJECT_NAME}.regenie.mask.txt",
        "{PROJECT_NAME}.step1_pred.list"
    output:
        "{PROJECT_NAME}.step2_gene_based_STATUS.regenie"
    params:
        input_basename="data/work/regenie/preprocess/{PROJECT_NAME}.all_chr",
        output_basename="{PROJECT_NAME}.step2_gene_based"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato --out {params.output_basename}"
        #--aaf-bins 0.01,0.001,0.0001 --strict-check-burden --test additive --minMAC 10

rule report_regenie:
    input:
        "{PROJECT_NAME}.regenie.pheno.txt",
        "{PROJECT_NAME}.step2_single_variant_STATUS.regenie",
        "{PROJECT_NAME}.step2_gene_based_STATUS.regenie"
    output:
        "{PROJECT_NAME}_{date}_regenie_report.html"
    shell:
        """
        Rscript -e "rmarkdown::render('report_regenie.Rmd',output_file='{output}')"
        """
