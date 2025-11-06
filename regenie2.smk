from datetime import datetime
import os

os.makedirs("logs/cluster/regenie2",exist_ok=True)
os.makedirs("data/work/regenie",exist_ok=True)

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

#date=datetime.now().strftime("%Y%m%d")

localrules:merge_chr_snplist,create_vep_list
wildcard_constraints:
    CHR='\d+'

def get_covariates(config):
    cmd=""
    if config['input'].get('covariates',False):
        cmd=f"--covariates {config['input']['covariates']} "
    if config['input'].get('sites',False):
        cmd+=f"--sites {config['input']['sites']} "
    if config['input'].get('eigenvec',False):
        cmd+=f"--eigenvec {config['input']['eigenvec']} "
    return cmd

rule run_regenie:
    input:
        expand("data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_single_variant_STATUS.regenie",PROJECT_NAME=PROJECT,CHR=CHROMOSOMES),
        expand("data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_gene_based_STATUS.regenie",PROJECT_NAME=PROJECT,CHR=CHROMOSOMES),

rule report_regenie:
    input:
        expand("{PROJECT_NAME}.regenie_report.{DATE}.html",PROJECT_NAME=PROJECT,DATE=datetime.now().strftime("%Y%m%d"))

rule preprocess_vcf:
    input:
        lambda wildcards: VCF_INPUT[f"chr{wildcards.CHR}"] if USE_INDIVIDUAL_CHR else VCF_INPUT
    output:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.bcf"
    params:
        samples=config['input']['samples']
    shell:
        """
        bcftools view -r chr{wildcards.CHR} -S {params.samples} {input} |
        bcftools norm -m-both |
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT'|
        bcftools +fill-tags -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE |
        bcftools filter -s ExcessHet -e 'INFO/ExcHet=0' -m + |
        bcftools filter -s LowCallRate -e 'INFO/F_MISSING > 0.05' -m + |
        bcftools filter -s LowSampleCount -e 'INFO/NS < 1000' -m + |
        bcftools filter -s NoHQHet -e 'COUNT(FORMAT/GT="0/1" && FORMAT/DP>=10 && FORMAT/GQ>=20 && FORMAT/VAF>0.2)=0' -m + -W=csi -Ob -o {output}
        """

rule vcf_filter:
    input:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.bcf"
    output:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.bcf"
    shell:
        "bcftools view -f PASS {input} -Ob -o {output}" #simple

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
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.bcf"
    output:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.pgen"
    params:
        sex_info=config['input']['sex_info'],
        output_prefix="data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc"
    shell:
        "plink2 --update-sex {params.sex_info} --bcf {input} --double-id --vcf-half-call r --make-pgen --out {params.output_prefix}"

rule chr_variant_list:
    input:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.pgen"
    output:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.regenie_step1.pgen",
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.regenie_step1.snplist"
    params:
        pfile="data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc",
        out="data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.regenie_step1",
    shell:
        "plink2 --pfile {params.pfile} --double-id  --maf 0.05 --snps-only --geno 0.1 --hwe 1e-6 --indep-pairwise 500 50 0.4 --write-snplist --make-pgen --out {params.out}"

rule merge_chr_pgen:
    input:
        expand("data/work/preprocess/{{PROJECT_NAME}}.chr{CHR}.qc.regenie_step1.pgen",CHR=CHROMOSOMES)
    output:
        "data/work/preprocess/{PROJECT_NAME}.all_chr.qc.regenie_step1.pgen"
        #"data/work/regenie/preprocess/{PROJECT_NAME}.all_chr.step1.eigenvec"
    params:
        out="data/work/preprocess/{PROJECT_NAME}.all_chr.qc.regenie_step1"
    shell:
        """
        rm -f file_list.txt
        for chrom in {{2..22}}; do echo \"data/work/preprocess/{wildcards.PROJECT_NAME}.chr${{chrom}}.regenie_step1.pgen data/work/preprocess/{wildcards.PROJECT_NAME}.chr${{chrom}}.regenie_step1.pvar data/work/preprocess/{wildcards.PROJECT_NAME}.chr${{chrom}}.regenie_step1.psam\" >> file_list.txt; done
        plink2 --pfile data/work/preprocess/{wildcards.PROJECT_NAME}.chr1.regenie_step1 --double-id --pmerge-list file_list.txt --make-pgen --out {params.out}
        """
        #--pca

rule merge_chr_snplist:
    input:
        expand("data/work/preprocess/{{PROJECT_NAME}}.chr{CHR}.qc.regenie_step1.snplist",CHR=CHROMOSOMES)
    output:
        "data/work/preprocess/{PROJECT_NAME}.all_chr.qc.regenie_step1.snplist"
    run:
        with open(output[0],'w') as f:
            for i in range(1,23):
                with open(f"data/work/preprocess/{PROJECT}.chr{i}.qc.regenie_step1.snplist",'r') as f2:
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
        samples=config['input']['samples'],
        controls=config['input']['controls']
    output:
        annotation="{PROJECT_NAME}.regenie.annotation.txt", #Im getting duplicates annotations. Need to fix this.
        set_file="{PROJECT_NAME}.regenie.set.txt",
        mask="{PROJECT_NAME}.regenie.mask.txt",
        covar="{PROJECT_NAME}.regenie.covar.txt",
        pheno="{PROJECT_NAME}.regenie.pheno.txt"
    params:
        optional_cmd=get_covariates(config), #Handles sites, covariates, eigenvec
        output_prefix=PROJECT,
    shell:
        """
        python preprocess_regenie.py \
        --vep-list-file {input.vep_list} \
        --samples {input.samples} \
        --controls {input.controls} \
        -O {params.output_prefix} \
        {params.optional_cmd}
        """

rule run_step1_regenie:
    input:
        "data/work/preprocess/{PROJECT_NAME}.all_chr.qc.regenie_step1.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "data/work/preprocess/{PROJECT_NAME}.all_chr.qc.regenie_step1.snplist"
    output:
        "{PROJECT_NAME}.step1_1.loco.gz",
        "{PROJECT_NAME}.step1_pred.list"
    params:
        input_basename="data/work/preprocess/{PROJECT_NAME}.all_chr.qc.regenie_step1",
        output_basename="{PROJECT_NAME}.step1",
        lowmem_prefix="data/work/regenie/tmp_rg_"
    shell:
        "regenie --step 1 --pgen {params.input_basename} --covarFile {input[1]} --phenoFile {input[2]} --extract {input[3]} "
        "--bsize 1000 --gz --bt --lowmem --lowmem-prefix {params.lowmem_prefix} --out {params.output_basename}"

rule run_step2_single_variant:
    input:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "{PROJECT_NAME}.step1_pred.list"
    output:
        "data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_single_variant_STATUS.regenie"
    params:
        input_basename="data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc",
        output_basename="data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_single_variant"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --covarFile {input[1]} --phenoFile {input[2]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc --minMAC 10 "
        "--out {params.output_basename}"

rule step2_single_variant_aux:
    input:
        "vep_files.list"
    output:
        "{PROJECT_NAME}.pathogenic_vus.csv"
    shell:
        """
        echo '"ID","Gene","Variant.LoF_level","HGVSc","HGVSp"' > {output}

        while IFS= read -r vep_file; do
            if [ -f "$vep_file" ]; then
                awk -F',' '($9 == "\"1\"" || $9 == "\"2\"") {print $6","$7","$9","$13","$14}' "$vep_file" >> {output}
            fi
        done < {input}
        """

rule run_step2_gene_based:
    input:
        "data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc.pgen",
        "{PROJECT_NAME}.regenie.covar.txt",
        "{PROJECT_NAME}.regenie.pheno.txt",
        "{PROJECT_NAME}.regenie.annotation.txt",
        "{PROJECT_NAME}.regenie.set.txt",
        "{PROJECT_NAME}.regenie.mask.txt",
        "{PROJECT_NAME}.step1_pred.list"
    output:
        "data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_gene_based_STATUS.regenie"
    params:
        input_basename="data/work/preprocess/{PROJECT_NAME}.chr{CHR}.qc",
        output_basename="data/work/regenie/{PROJECT_NAME}.chr{CHR}.step2_gene_based"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--aaf-bins 0.01,0.001,0.0001 --strict-check-burden "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato "
        "--out {params.output_basename}"

rule run_regenie_report:
    input:
        expand("data/work/regenie/{{PROJECT_NAME}}.chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES),
        expand("data/work/regenie/{{PROJECT_NAME}}.chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES),
        "{PROJECT_NAME}.pathogenic_vus.csv"
    output:
        "{PROJECT_NAME}.regenie_report.{DATE}.html"
    shell:
        """
        R -e "rmarkdown::render('report_regenie.Rmd',params=list(project_name='{wildcards.PROJECT_NAME}'),output_file='{wildcards.PROJECT_NAME}.regenie_report.{wildcards.DATE}.html')"
        """