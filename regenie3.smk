
from datetime import datetime
import os

os.makedirs("logs/cluster/regenie3",exist_ok=True)
os.makedirs("logs/lsf",exist_ok=True)
os.makedirs("data/regenie",exist_ok=True)

# Extract configuration values
PROJECT_NAME=config['project']['name']
CHROMOSOMES_AUTOSOMAL=list(range(1,23)) # Chromosomes 1-22

with open(config['input'].get('ancestry_file','data/preprocess/build_pca.eigenvec'),'r') as f:
    header=f.readline().strip()
    DM_COUNT=len(header.split())-2
    print(f"{DM_COUNT} Ancestry dimensions found.")

localrules:create_vep_list
wildcard_constraints:
    CHR=r'\d+'


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
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/regenie/vep_files.list"
    run:
        with open(output[0],'w') as f:
            for chr in CHROMOSOMES_AUTOSOMAL:
                vep_file=f"data/preprocess/chr{chr}.annotation.no_sample.vep.report.csv"
                if os.path.exists(vep_file):
                    f.write(vep_file+'\n')

rule preprocess_regenie:
    input:
        vep_list="data/regenie/vep_files.list",
        samples="data/qc/passing_samples.txt",
    output:
        annotation="data/regenie/regenie.annotation.txt",
        set_file="data/regenie/regenie.set.txt",
        mask="data/regenie/regenie.mask.txt",
        covar="data/regenie/regenie.covar.txt",
        pheno="data/regenie/regenie.pheno.txt"
    params:
        covariates=config['input'].get('covariates','covariates.txt'),
        ancestry=config['input'].get('ancestry_file','data/preprocess/build_pca.eigenvec'),
        output_prefix="data/regenie",
    shell:
        """
        python preprocess_regenie.py \
        --vep-list-file {input.vep_list} \
        --samples {input.samples} \
        -O {params.output_prefix} \
        --covariates {params.covariates} \
        --ancestry-file {params.ancestry}
        """

rule run_step1_regenie:
    input:
        "data/preprocess/build.pgen",
        "data/regenie/regenie.covar.txt",
        "data/regenie/regenie.pheno.txt",
        "data/preprocess/build.snplist"
    output:
        "data/regenie/step1_1.loco.gz",
        "data/regenie/step1_pred.list"
    params:
        input_basename="data/preprocess/build",
        output_basename="data/regenie/step1",
        lowmem_prefix="data/regenie/tmp_rg_",
        dm_covar_flags=f'DM{{1:{DM_COUNT}}}',
    shell:
        "regenie --step 1 --pgen {params.input_basename} --covarFile {input[1]} --covarCol FREEZE --covarCol {params.dm_covar_flags} "
        "--phenoFile {input[2]} --phenoCol STATUS --extract {input[3]} "
        "--bsize 1000 --gz --bt --lowmem --lowmem-prefix {params.lowmem_prefix} --out {params.output_basename}"

rule run_step2_single_variant:
    input:
        "data/preprocess/chr{CHR}.annotation.pgen",
        "data/regenie/regenie.covar.txt",
        "data/regenie/regenie.pheno.txt",
        "data/regenie/step1_pred.list"
    output:
        "data/regenie/chr{CHR}.step2_single_variant_STATUS.regenie"
    params:
        input_basename="data/preprocess/chr{CHR}.annotation",
        output_basename="data/regenie/chr{CHR}.step2_single_variant",
        dm_covar_flags=f'DM{{1:{DM_COUNT}}}'
    shell:
        "regenie --step 2 --pgen {params.input_basename} --covarFile {input[1]} --covarCol FREEZE --covarCol {params.dm_covar_flags} "
        "--phenoFile {input[2]} --phenoCol STATUS --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc --minMAC 10 "
        "--out {params.output_basename}"

rule step2_single_variant_aux:
    input:
        "data/regenie/vep_files.list"
    output:
        "data/regenie/pathogenic_vus.csv"
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
        "data/preprocess/chr{CHR}.annotation.pgen",
        "data/regenie/regenie.covar.txt",
        "data/regenie/regenie.pheno.txt",
        "data/regenie/regenie.annotation.txt",
        "data/regenie/regenie.set.txt",
        "data/regenie/regenie.mask.txt",
        "data/regenie/step1_pred.list"
    output:
        "data/regenie/chr{CHR}.step2_gene_based_STATUS.regenie"
    params:
        input_basename="data/preprocess/chr{CHR}.annotation",
        output_basename="data/regenie/chr{CHR}.step2_gene_based",
        dm_covar_flags=f'DM{{1:{DM_COUNT}}}'
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --phenoCol STATUS "
        "--covarFile {input[1]} --covarCol FREEZE --covarCol {params.dm_covar_flags} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--aaf-bins 0.01,0.001,0.0001 --strict-check-burden "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato "
        "--out {params.output_basename}"

rule final_regenie_results:
    input:
        single_variant=expand("data/regenie/chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/regenie/chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        covar="data/regenie/regenie.covar.txt",
        pheno="data/regenie/regenie.pheno.txt",
        pathogenic="data/regenie/pathogenic_vus.csv"
    output:
        single_variant=expand("data/final/{{PROJECT}}.chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/final/{{PROJECT}}.chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        covar="data/final/{PROJECT}.regenie.covar.txt",
        pheno="data/final/{PROJECT}.regenie.pheno.txt",
        pathogenic="data/final/{PROJECT}.pathogenic_vus.csv"
    shell:
        """
        for f in data/regenie/chr*.step2_single_variant_STATUS.regenie; do cp "${{f}}" "data/final/{wildcards.PROJECT}.$(basename "${{f}}")"; done
        for f in data/regenie/chr*.step2_gene_based_STATUS.regenie; do cp "${{f}}" "data/final/{wildcards.PROJECT}.$(basename "${{f}}")"; done
        cp {input.covar} {output.covar}
        cp {input.pheno} {output.pheno}
        cp {input.pathogenic} {output.pathogenic}
        """

rule regenie_report:
    input:
        single_variant=expand("data/final/{{PROJECT}}.chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/final/{{PROJECT}}.chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        pathogenic="data/final/{PROJECT}.pathogenic_vus.csv"
    output:
        "data/final/{PROJECT}.regenie_report.html"
    shell:
        """
        R -e "rmarkdown::render('report_regenie.Rmd',params=list(project_name='{wildcards.PROJECT}'),output_file='{output}')"
        """