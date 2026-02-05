from datetime import datetime
import os

os.makedirs("logs/cluster/regenie3",exist_ok=True)
os.makedirs("logs/lsf",exist_ok=True)
os.makedirs("data/regenie",exist_ok=True)

# Extract configuration values
PROJECT_NAME=config['project']['name']
CHROMOSOMES_AUTOSOMAL=list(range(1,23)) # Chromosomes 1-22

localrules:create_vep_list
wildcard_constraints:
    CHR=r'\d+',
    PROJECT='[A-Za-z0-9_-]+'

def get_covariates(config):
    cmd=""
    if config['input'].get('covariates',False):
        cmd=f"--covariates {config['input']['covariates']} "
    if config['input'].get('sites',False):
        cmd+=f"--sites {config['input']['sites']} "
    # Default to QC pipeline PCA output from select_variants
    eigenvec=config['input'].get('eigenvec',f"data/preprocess/{PROJECT_NAME}.build.eigenvec")
    cmd+=f"--eigenvec {eigenvec} "
    return cmd

rule run_regenie_report:
    input:
        expand("data/final/{PROJECT}.regenie_report.html",PROJECT=PROJECT_NAME)

rule run_regenie:
    input:
        expand("data/final/{PROJECT}.chr{CHR}.step2_single_variant_STATUS.regenie",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/{PROJECT}.chr{CHR}.step2_gene_based_STATUS.regenie",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/{PROJECT}.regenie_report.html",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.regenie.covar.txt",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.regenie.pheno.txt",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.pathogenic_vus.csv",PROJECT=PROJECT_NAME)

rule create_vep_list:
    input:
        expand("data/preprocess/{PROJECT}.chr{CHR}.annotation.no_sample.vep.report.csv",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/regenie/vep_files.list"
    run:
        with open(output[0],'w') as f:
            for chr in CHROMOSOMES_AUTOSOMAL:
                vep_file=f"data/preprocess/{PROJECT_NAME}.chr{chr}.annotation.no_sample.vep.report.csv"
                if os.path.exists(vep_file):
                    f.write(vep_file+'\n')

rule preprocess_regenie:
    input:
        vep_list="data/regenie/vep_files.list",
        samples="data/qc/{PROJECT}.passing_samples.txt",
        controls=config['input'].get('controls','controls.txt')
    output:
        annotation="data/regenie/{PROJECT}.regenie.annotation.txt",
        set_file="data/regenie/{PROJECT}.regenie.set.txt",
        mask="data/regenie/{PROJECT}.regenie.mask.txt",
        covar="data/regenie/{PROJECT}.regenie.covar.txt",
        pheno="data/regenie/{PROJECT}.regenie.pheno.txt"
    params:
        optional_cmd=get_covariates(config),
        output_prefix="data/regenie/{PROJECT}",
    log:
        stderr="logs/lsf/preprocess_regenie.{PROJECT}.e",
        stdout="logs/lsf/preprocess_regenie.{PROJECT}.o"
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
        "data/preprocess/{PROJECT}.build.pgen",
        "data/regenie/{PROJECT}.regenie.covar.txt",
        "data/regenie/{PROJECT}.regenie.pheno.txt",
        "data/preprocess/{PROJECT}.build.snplist"
    output:
        "data/regenie/{PROJECT}.step1_1.loco.gz",
        "data/regenie/{PROJECT}.step1_pred.list"
    params:
        input_basename="data/preprocess/{PROJECT}.build",
        output_basename="data/regenie/{PROJECT}.step1",
        lowmem_prefix="data/regenie/tmp_rg_"
    log:
        stderr="logs/lsf/run_step1_regenie.{PROJECT}.e",
        stdout="logs/lsf/run_step1_regenie.{PROJECT}.o"
    shell:
        "regenie --step 1 --pgen {params.input_basename} --covarFile {input[1]} --phenoFile {input[2]} --extract {input[3]} "
        "--bsize 1000 --gz --bt --lowmem --lowmem-prefix {params.lowmem_prefix} --out {params.output_basename}"

rule run_step2_single_variant:
    input:
        "data/preprocess/{PROJECT}.chr{CHR}.annotation.pgen",
        "data/regenie/{PROJECT}.regenie.covar.txt",
        "data/regenie/{PROJECT}.regenie.pheno.txt",
        "data/regenie/{PROJECT}.step1_pred.list"
    output:
        "data/regenie/{PROJECT}.chr{CHR}.step2_single_variant_STATUS.regenie"
    params:
        input_basename="data/preprocess/{PROJECT}.chr{CHR}.annotation",
        output_basename="data/regenie/{PROJECT}.chr{CHR}.step2_single_variant"
    log:
        stderr="logs/lsf/run_step2_single_variant.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/run_step2_single_variant.{PROJECT}.chr{CHR}.o"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --covarFile {input[1]} --phenoFile {input[2]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc --minMAC 10 "
        "--out {params.output_basename}"

rule step2_single_variant_aux:
    input:
        "data/regenie/vep_files.list"
    output:
        "data/regenie/{PROJECT}.pathogenic_vus.csv"
    log:
        stderr="logs/lsf/step2_single_variant_aux.{PROJECT}.e",
        stdout="logs/lsf/step2_single_variant_aux.{PROJECT}.o"
    shell:
        """
        echo '"ID","Gene","Variant.LoF_level","HGVSc","HGVSp"' > {output}

        while IFS= read -r vep_file; do
            if [ -f "$vep_file" ]; then
                awk -F',' '($9 == "\\"1\\"" || $9 == "\\"2\\"") {{print $6 "," $7 "," $9 "," $13 "," $14}}' "$vep_file" >> {output}
            fi
        done < {input}
        """

rule run_step2_gene_based:
    input:
        "data/preprocess/{PROJECT}.chr{CHR}.annotation.pgen",
        "data/regenie/{PROJECT}.regenie.covar.txt",
        "data/regenie/{PROJECT}.regenie.pheno.txt",
        "data/regenie/{PROJECT}.regenie.annotation.txt",
        "data/regenie/{PROJECT}.regenie.set.txt",
        "data/regenie/{PROJECT}.regenie.mask.txt",
        "data/regenie/{PROJECT}.step1_pred.list"
    output:
        "data/regenie/{PROJECT}.chr{CHR}.step2_gene_based_STATUS.regenie"
    params:
        input_basename="data/preprocess/{PROJECT}.chr{CHR}.annotation",
        output_basename="data/regenie/{PROJECT}.chr{CHR}.step2_gene_based"
    log:
        stderr="logs/lsf/run_step2_gene_based.{PROJECT}.chr{CHR}.e",
        stdout="logs/lsf/run_step2_gene_based.{PROJECT}.chr{CHR}.o"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--aaf-bins 0.01,0.001,0.0001 --strict-check-burden "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato "
        "--out {params.output_basename}"

rule final_regenie_results:
    input:
        single_variant=expand("data/regenie/{{PROJECT}}.chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/regenie/{{PROJECT}}.chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        covar="data/regenie/{PROJECT}.regenie.covar.txt",
        pheno="data/regenie/{PROJECT}.regenie.pheno.txt",
        pathogenic="data/regenie/{PROJECT}.pathogenic_vus.csv"
    output:
        single_variant=expand("data/final/{{PROJECT}}.chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/final/{{PROJECT}}.chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        covar="data/final/{PROJECT}.regenie.covar.txt",
        pheno="data/final/{PROJECT}.regenie.pheno.txt",
        pathogenic="data/final/{PROJECT}.pathogenic_vus.csv"
    log:
        stderr="logs/lsf/final_regenie_results.{PROJECT}.e",
        stdout="logs/lsf/final_regenie_results.{PROJECT}.o"
    shell:
        """
        rsync -av data/regenie/{wildcards.PROJECT}.chr*.step2_single_variant_STATUS.regenie data/final/
        rsync -av data/regenie/{wildcards.PROJECT}.chr*.step2_gene_based_STATUS.regenie data/final/
        rsync -av {input.covar} {output.covar}
        rsync -av {input.pheno} {output.pheno}
        rsync -av {input.pathogenic} {output.pathogenic}
        """

rule regenie_report:
    input:
        single_variant=expand("data/final/{{PROJECT}}.chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/final/{{PROJECT}}.chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        pathogenic="data/final/{PROJECT}.pathogenic_vus.csv"
    output:
        "data/final/{PROJECT}.regenie_report.html"
    log:
        stderr="logs/lsf/regenie_report.{PROJECT}.e",
        stdout="logs/lsf/regenie_report.{PROJECT}.o"
    shell:
        """
        R -e "rmarkdown::render('report_regenie.Rmd',params=list(project_name='{wildcards.PROJECT}'),output_file='{output}')"
        """