
#needs preprocess_regenie.py
#needs report_regenie.Rmd

from datetime import datetime
import os

os.makedirs("logs/cluster/regenie3",exist_ok=True)
os.makedirs("logs/lsf",exist_ok=True)
os.makedirs("data/regenie",exist_ok=True)

_raw_af_delta=config.get("qc",{}).get("af_delta_thr", None)
try:
    AF_DELTA_THR=float(_raw_af_delta) if _raw_af_delta is not None and str(_raw_af_delta).strip()!="" else 0.0
except (TypeError, ValueError):
    AF_DELTA_THR=0.0
AF_DELTA_EXCLUDE="data/qc/exclusions/af_delta_exclude.snplist"

def af_delta_exclude_input(wildcards=None):
    return [AF_DELTA_EXCLUDE] if AF_DELTA_THR>0 else []

def snp_id_exclude_flag(path):
    if not path:
        return ""
    if not isinstance(path, str):
        path=str(path[0]) if len(path) else ""
    if not path or not os.path.isfile(path):
        return ""
    with open(path) as f:
        for line in f:
            s=line.strip()
            if s and not s.startswith("#"):
                return f"--exclude {path}"
    return ""

# Extract configuration values
PROJECT_NAME=config['project']['name']
CHROMOSOMES_AUTOSOMAL=list(range(1,23)) # Chromosomes 1-22

DM_TARGET=int(config.get('regenie',{}).get('dm_count',10))
if DM_TARGET<=0:
    DM_FOUND=0
    DM_COUNT=0
    print("regenie.dm_count<=0: FREEZE-only run, no ancestry/PC covariates.")
else:
    with open(config['input'].get('ancestry_file','data/preprocess/build_pca_clean.eigenvec'),'r') as f:
        header=f.readline().strip()
        DM_FOUND=len(header.split())-2
    DM_COUNT=min(DM_FOUND,DM_TARGET)
    print(f"{DM_FOUND} ancestry dimensions found; using {DM_COUNT} for REGENIE covariates.")

localrules:create_vep_list,build_regenie_report_tables,regenie_report
wildcard_constraints:
    CHR=r'\d+'

rule run_regenie:
    input:
        expand("data/final/{PROJECT}.chr{CHR}.step2_single_variant_STATUS.regenie",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/{PROJECT}.chr{CHR}.step2_gene_based_STATUS.regenie",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/final/{PROJECT}.regenie.covar.txt",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.regenie.pheno.txt",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.pathogenic_vus.csv",PROJECT=PROJECT_NAME),
        expand("data/final/{PROJECT}.synonymous.csv",PROJECT=PROJECT_NAME)

rule run_regenie_report:
    input:
        expand("data/final/{PROJECT}.regenie_report.html",PROJECT=PROJECT_NAME)

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
        gene_consequence_exclude=config.get('regenie',{}).get(
            'gene_consequence_exclude','mask_gene_consequence_exclude.txt'),
    output:
        annotation="data/regenie/regenie.annotation.txt",
        set_file="data/regenie/regenie.set.txt",
        mask="data/regenie/regenie.mask.txt",
        covar="data/regenie/regenie.covar.txt",
        pheno="data/regenie/regenie.pheno.txt",
        synonymous="data/regenie/synonymous.csv"
    params:
        covariates=config['input'].get('covariates','covariates.txt'),
        # FREEZE-only run (DM_COUNT==0): omit ancestry so covar = FID IID FREEZE.
        ancestry_args=(f"--ancestry-file {config['input'].get('ancestry_file','data/preprocess/build_pca_clean.eigenvec')} --dm-count {DM_COUNT}" if DM_COUNT>0 else ""),
        negative_control_blacklist="data/mnp/mnp_blacklist.txt",
        output_prefix="data/regenie",
    shell:
        """
        python preprocess_regenie.py \
        --vep-list-file {input.vep_list} \
        --samples {input.samples} \
        -O {params.output_prefix} \
        --synonymous-out {output.synonymous} \
        --covariates {params.covariates} \
        {params.ancestry_args} \
        --negative-control-blacklist {params.negative_control_blacklist} \
        --gene-consequence-exclude {input.gene_consequence_exclude}
        """

rule make_step1_snplist:
    input:
        pgen="data/preprocess/build.pgen",
        pvar="data/preprocess/build.pvar",
        psam="data/preprocess/build.psam",
        long_ld=config.get("input", {}).get("long_ld_bed", "long_ld_regions.bed"),
        id_excl=af_delta_exclude_input,
    output:
        snplist="data/regenie/step1.snplist",
    params:
        input_basename="data/preprocess/build",
        output_basename="data/regenie/step1_extract",
        id_excl_flag=lambda wildcards, input: snp_id_exclude_flag(input.id_excl),
    shell:
        """
        plink2 --pfile {params.input_basename} \
          --exclude range {input.long_ld} \
          --write-snplist --out {params.output_basename}_ld
        if [ -n "{params.id_excl_flag}" ]; then
          plink2 --pfile {params.input_basename} \
            --extract {params.output_basename}_ld.snplist \
            {params.id_excl_flag} \
            --write-snplist --out {params.output_basename}
          mv {params.output_basename}.snplist {output.snplist}
          rm -f {params.output_basename}.log
        else
          mv {params.output_basename}_ld.snplist {output.snplist}
        fi
        rm -f {params.output_basename}_ld.snplist {params.output_basename}_ld.log
        """

rule run_step1_regenie:
    input:
        "data/preprocess/build.pgen",
        "data/regenie/regenie.covar.txt",
        "data/regenie/regenie.pheno.txt",
        "data/regenie/step1.snplist",
    output:
        "data/regenie/step1_1.loco.gz",
        "data/regenie/step1_pred.list"
    params:
        input_basename="data/preprocess/build",
        output_basename="data/regenie/step1",
        lowmem_prefix="data/regenie/tmp_rg_",
        covar_args=(f'--covarCol FREEZE --covarCol DM{{1:{DM_COUNT}}}' if DM_COUNT>0 else '--covarCol FREEZE'),
    shell:
        "regenie --step 1 --pgen {params.input_basename} --covarFile {input[1]} {params.covar_args} "
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
        covar_args=(f'--covarCol FREEZE --covarCol DM{{1:{DM_COUNT}}}' if DM_COUNT>0 else '--covarCol FREEZE')
    shell:
        "regenie --step 2 --pgen {params.input_basename} --covarFile {input[1]} {params.covar_args} "
        "--phenoFile {input[2]} --phenoCol STATUS --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc --minMAC 1 "
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
        "data/regenie/chr{CHR}.step2_gene_based_STATUS.regenie",
        "data/regenie/chr{CHR}.step2_gene_based_masks.snplist",
    params:
        input_basename="data/preprocess/chr{CHR}.annotation",
        output_basename="data/regenie/chr{CHR}.step2_gene_based",
        covar_args=(f'--covarCol FREEZE --covarCol DM{{1:{DM_COUNT}}}' if DM_COUNT>0 else '--covarCol FREEZE')
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --phenoCol STATUS "
        "--covarFile {input[1]} {params.covar_args} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--aaf-bins 0.01,0.001,0.0001 --strict-check-burden --minMAC 1 "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato "
        "--out {params.output_basename}"

rule mask_variant_stats:
    # Full mask site set (large). Report subsets by ID; do not read whole TSV in R.
    input:
        vep_list="data/regenie/vep_files.list",
        pheno="data/regenie/regenie.pheno.txt",
        pgen=expand("data/preprocess/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/regenie/mask_variant_stats.tsv"
    shell:
        """
        python mask_variant_stats.py \
          --vep-list {input.vep_list} \
          --pheno {input.pheno} \
          --pfile-template data/preprocess/chr{{CHR}}.annotation \
          -o {output}
        """

rule build_regenie_report_tables:
    # Small tables for the HTML/CSV report (top genes + snplist variants + consequence counts
    # + unique gene-level case/control carriers + CHEK2 carrier IIDs for PC plots)
    input:
        gene=expand("data/final/{PROJECT}.chr{CHR}.step2_gene_based_STATUS.regenie",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        stats="data/regenie/mask_variant_stats.tsv",
        snplist=expand("data/regenie/chr{CHR}.step2_gene_based_masks.snplist",CHR=CHROMOSOMES_AUTOSOMAL),
        pheno="data/regenie/regenie.pheno.txt",
        pgen=expand("data/preprocess/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
    output:
        top=f"data/final/{PROJECT_NAME}.report.top_genes_ADD.tsv",
        top_skat=f"data/final/{PROJECT_NAME}.report.top_genes_SKAT.tsv",
        contrib=f"data/final/{PROJECT_NAME}.report.variant_contrib.tsv",
        consequence=f"data/final/{PROJECT_NAME}.report.consequence_matrix.tsv",
        gene_carriers=f"data/final/{PROJECT_NAME}.report.gene_carriers.tsv",
        chek2_carriers=f"data/final/{PROJECT_NAME}.report.chek2_carriers.tsv",
    params:
        gene_glob=f"data/final/{PROJECT_NAME}.chr*.step2_gene_based_STATUS.regenie",
        snplist_glob="data/regenie/chr*.step2_gene_based_masks.snplist",
        out_prefix=f"data/final/{PROJECT_NAME}.report"
    shell:
        """
        python build_regenie_report_tables.py \
          --gene-glob '{params.gene_glob}' \
          --snplist-glob '{params.snplist_glob}' \
          --mask-stats {input.stats} \
          --pheno {input.pheno} \
          --pfile-template data/preprocess/chr{{CHR}}.annotation \
          --out-prefix {params.out_prefix}
        """

rule final_regenie_results:
    input:
        single_variant=expand("data/regenie/chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/regenie/chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        covar="data/regenie/regenie.covar.txt",
        pheno="data/regenie/regenie.pheno.txt",
        pathogenic="data/regenie/pathogenic_vus.csv",
        synonymous="data/regenie/synonymous.csv"
    output:
        single_variant=expand("data/final/{{PROJECT}}.chr{CHR}.step2_single_variant_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/final/{{PROJECT}}.chr{CHR}.step2_gene_based_STATUS.regenie",CHR=CHROMOSOMES_AUTOSOMAL),
        covar="data/final/{PROJECT}.regenie.covar.txt",
        pheno="data/final/{PROJECT}.regenie.pheno.txt",
        pathogenic="data/final/{PROJECT}.pathogenic_vus.csv",
        synonymous="data/final/{PROJECT}.synonymous.csv"
    shell:
        """
        for f in data/regenie/chr*.step2_single_variant_STATUS.regenie; do cp "${{f}}" "data/final/{wildcards.PROJECT}.$(basename "${{f}}")"; done
        for f in data/regenie/chr*.step2_gene_based_STATUS.regenie; do cp "${{f}}" "data/final/{wildcards.PROJECT}.$(basename "${{f}}")"; done
        cp {input.covar} {output.covar}
        cp {input.pheno} {output.pheno}
        cp {input.pathogenic} {output.pathogenic}
        cp {input.synonymous} {output.synonymous}
        """

rule regenie_report:
    input:
        single_variant=expand("data/final/{PROJECT}.chr{CHR}.step2_single_variant_STATUS.regenie",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        gene_based=expand("data/final/{PROJECT}.chr{CHR}.step2_gene_based_STATUS.regenie",PROJECT=PROJECT_NAME,CHR=CHROMOSOMES_AUTOSOMAL),
        pathogenic=f"data/final/{PROJECT_NAME}.pathogenic_vus.csv",
        blacklist="data/mnp/mnp_blacklist.txt",
        report_top=f"data/final/{PROJECT_NAME}.report.top_genes_ADD.tsv",
        report_skat=f"data/final/{PROJECT_NAME}.report.top_genes_SKAT.tsv",
        report_contrib=f"data/final/{PROJECT_NAME}.report.variant_contrib.tsv",
        report_cons=f"data/final/{PROJECT_NAME}.report.consequence_matrix.tsv",
        report_gene_carr=f"data/final/{PROJECT_NAME}.report.gene_carriers.tsv",
        report_chek2_carr=f"data/final/{PROJECT_NAME}.report.chek2_carriers.tsv",
    output:
        f"data/final/{PROJECT_NAME}.regenie_report.html"
    params:
        project=PROJECT_NAME,
    shell:
        "R -e \"rmarkdown::render('report_regenie.Rmd',"
        "params=list(project_name='{params.project}'),"
        "output_file='{output}')\""
