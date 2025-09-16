from datetime import datetime
import os

os.makedirs("logs/cluster/regenie2",exist_ok=True)

# Extract configuration values
PROJECT=config['project']['name']
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
wildcard_constraints:
    CHR='[1-22]'

rule run_regenie:
    input:
        expand("data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_single_variant_STATUS.regenie",PROJECT_NAME=PROJECT,CHR=CHROMOSOMES),
        expand("data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_gene_based_STATUS.regenie",PROJECT_NAME=PROJECT,CHR=CHROMOSOMES)
        

rule preprocess_vcf:
    input:
        lambda wildcards: VCF_INPUT[f"chr{wildcards.CHR}"] if USE_INDIVIDUAL_CHR else VCF_INPUT
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.bcf"
    params:
        samples=config['input']['samples']
    shell:
        "bcftools view -f 'PASS' -r chr{wildcards.CHR} -S {params.samples} {input} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ob -o {output}"

#--vcf-half-call <mode>

#The current VCF standard does not specify how '0/.' and similar GT values should be interpreted.
#By default (mode 'error'/'e'), PLINK 1.9 errors out and reports the line number of the anomaly.
#Should the half-call be intentional, though (this can be the case with Complete Genomics data),
#you can request the following other modes:
#
#    'haploid'/'h': Treat half-calls as haploid/homozygous (the PLINK 1 file format does not distinguish between the two). This maximizes similarity between the VCF and BCF2 parsers.
#    'missing'/'m': Treat half-calls as missing.
#    'reference'/'r': Treat the missing part as reference.

rule chr_pgen:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.bcf"
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.pgen"
    params:
        sex_info=config['input']['sex_info'],
        output_prefix="data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}"
    shell:
        "plink2 --update-sex {params.sex_info} --bcf {input} --double-id --vcf-half-call r --make-pgen --out {params.output_prefix}"

rule chr_variant_list:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.pgen"
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.pgen",
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.pvar",
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.psam",
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.snplist"
    params:
        pfile="data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}",
        out="data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1",
    shell:
        "plink2 --pfile {params.pfile} --double-id  --maf 0.05 --snps-only --geno 0.1 --hwe 1e-6 --indep-pairwise 500 50 0.4 --write-snplist --make-pgen --out {params.out}"

rule merge_chr_pgen:
    input:
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.pgen",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT),
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.pvar",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT),
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.psam",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT)
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.pgen",
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.pvar",
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.psam"
        #"data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.eigenvec"
    params:
        out="data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1"
    shell:
        """
        rm -f file_list.txt
        for chrom in {{2..22}}; do echo \"data/work/regenie/preprocess/{wildcards.PROJECT_NAME}.chr${{chrom}}.step1.pgen data/work/regenie/preprocess/{wildcards.PROJECT_NAME}.chr${{chrom}}.step1.pvar data/work/regenie/preprocess/{wildcards.PROJECT_NAME}.chr${{chrom}}.step1.psam\" >> file_list.txt; done
        plink2 --pfile data/work/regenie/preprocess/{wildcards.PROJECT_NAME}.chr1.step1 --double-id --pmerge-list file_list.txt --make-pgen --out {params.out}
        """
        #--pca

rule merge_chr_snplist:
    input:
        expand("data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.step1.snplist",CHR=CHROMOSOMES,PROJECT_NAME=PROJECT)
    output:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.snplist"
    run:
        with open(output[0],'w') as f:
            for i in range(1,23):
                with open(f"data/work/regenie/preprocess/{PROJECT}.chr{i}.step1.snplist",'r') as f2:
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
        #eigenvec=config['eigenvec'],
        sites=config['input']['sites'],
        controls=config['input']['controls']
    output:
        annotation="{PROJECT_NAME}.regenie.annotation.txt",
        set_file="{PROJECT_NAME}.regenie.set.txt",
        mask="{PROJECT_NAME}.regenie.mask.txt",
        covar="{PROJECT_NAME}.regenie.covar.txt",
        pheno="{PROJECT_NAME}.regenie.pheno.txt"
    params:
        output_prefix=PROJECT,
        eigenvec="subset_withFID.eigenvec" #Hard code for now
    shell:
        """
        COVAR_ARG=""
        if [ -n "{config[input][covariates]}" ] && [ -f "{config[input][covariates]}" ]; then
            COVAR_ARG="--covariates {config[input][covariates]}"
        fi

        python preprocess_regenie.py \
        --vep-list-file {input.vep_list} \
        --eigenvec {params.eigenvec} \
        --sites {input.sites} \
        --controls {input.controls} \
        $COVAR_ARG \
        -O {params.output_prefix}
        """

rule run_step1_regenie:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.snplist"
    output:
        "{PROJECT_NAME}.step1_1.loco.gz",
        "{PROJECT_NAME}.step1_pred.list"
    params:
        input_basename="data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1",
        output_basename="{PROJECT_NAME}.step1",
        lowmem_prefix="data/work/regenie/tmp_rg_"
    shell:
        "regenie --step 1 --pgen {params.input_basename} --covarFile {input[1]} --phenoFile {input[2]} --extract {input[3]} "
        "--bsize 1000 --gz --bt --lowmem --lowmem-prefix {params.lowmem_prefix} --out {params.output_basename}"
        #--strict --loocv

rule run_step2_single_variant:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "{PROJECT_NAME}.step1_pred.list"
    output:
        "data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_single_variant_STATUS.regenie"
    params:
        input_basename="data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}",
        output_basename="data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_single_variant"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc "
        "--out {params.output_basename}"
        #--strict --minMAC 10 --test additive

rule run_step2_gene_based:
    input:
        "data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "{PROJECT_NAME}.regenie.annotation.txt",
        "{PROJECT_NAME}.regenie.set.txt",
        "{PROJECT_NAME}.regenie.mask.txt",
        "{PROJECT_NAME}.step1_pred.list"
    output:
        "data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_gene_based_STATUS.regenie"
    params:
        input_basename="data/work/regenie/preprocess/{PROJECT_NAME}.chr{CHR}",
        output_basename="data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_gene_based"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato --out {params.output_basename}"
        #--aaf-bins 0.01,0.001,0.0001 --strict-check-burden --test additive --minMAC 10
