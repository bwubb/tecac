#include: "select_variants3.smk"
#include: "build_files.smk"
#include: "annotation_files2.smk"
#include: "missingness_bias.smk"

import os
import datetime

os.makedirs('logs/lsf',exist_ok=True)
DATE=datetime.date.today().strftime('%Y%m%d')


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

rule qc_report:
    input:
        f"exwas_qc_report.{DATE}.html",

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


###### FINAL OUTPUT #########
# QC report does not list diff-miss files as Snakemake inputs; run missingness_bias.smk if you want those sections populated.

# Generate QC report
rule generate_qc_report:
    input:
        sex_report="data/qc/reports/sex_inference.txt",
        miss_report="data/qc/reports/missingness_report.txt",
        het_report="data/qc/reports/het_report.txt",
        het_exclusions="data/qc/exclusions/het_exclusions.txt",
        pca_eigenvec="data/preprocess/build_pca_clean.eigenvec",
        pca_eigenval="data/preprocess/build_pca_clean.eigenval",
        prefilter_stats=expand("data/qc/reports/chr{CHR}.prefilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postfilter_stats=expand("data/qc/reports/chr{CHR}.postfilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        plink_filter_metrics=expand("data/qc/reports/chr{CHR}.plink_filter_metrics.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        postplink_counts=expand("data/qc/reports/chr{CHR}.postplink.variant_count.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_afreq=expand("data/plink/chr{CHR}.afreq",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_vmiss=expand("data/plink/chr{CHR}.vmiss",CHR=CHROMOSOMES_AUTOSOMAL),
        variant_hardy=expand("data/plink/chr{CHR}.hardy",CHR=CHROMOSOMES_AUTOSOMAL),
        postmnp_stats=expand("data/qc/reports/chr{CHR}.postmnp.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        html="data/qc/reports/exwas_qc_report.html",
        dated=f"exwas_qc_report.{DATE}.html",
    params:
        project=PROJECT_NAME,
        data_dir="data",
        covariate_file=config.get('input',{}).get('covariates','')
    resources:
        lsf_err="logs/lsf/generate_qc_report.e",
        lsf_out="logs/lsf/generate_qc_report.o"
    shell:
        """
        Rscript -e "suppressWarnings(rmarkdown::render('exwas_qc_report2.Rmd', output_file='{output.html}', params=list(project='{params.project}', data_dir='{params.data_dir}', covariate_file='{params.covariate_file}')))"
        cp {output.html} {output.dated}
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
    resources:
        lsf_err="logs/lsf/final_annotation.chr{CHR}.e",
        lsf_out="logs/lsf/final_annotation.chr{CHR}.o"
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
        eigenvec="data/preprocess/build_pca_clean.eigenvec",
        eigenval="data/preprocess/build_pca_clean.eigenval"
    output:
        pgen="data/final/build.pgen",
        pvar="data/final/build.pvar",
        psam="data/final/build.psam",
        snplist="data/final/build.snplist",
        eigenvec="data/final/build.eigenvec",
        eigenval="data/final/build.eigenval"
    resources:
        lsf_err="logs/lsf/final_build.e",
        lsf_out="logs/lsf/final_build.o"
    shell:
        """
        rsync -av {input.pgen} {output.pgen}
        rsync -av {input.pvar} {output.pvar}
        rsync -av {input.psam} {output.psam}
        rsync -av {input.snplist} {output.snplist}
        rsync -av {input.eigenvec} {output.eigenvec}
        rsync -av {input.eigenval} {output.eigenval}
        """
    
# Splits passing samples into controls vs cases. Prefers covariate file (col4=STATUS, 1=control 2=case); else separate controls file.
rule final_lists:
    input:
        samples="data/qc/passing_samples.txt",
        exclusions1="data/qc/exclusions/sexcheck_fail.txt"
    output:
        samples="data/final/samples.txt",
        controls="data/final/controls.txt",
        cases="data/final/cases.txt"
    params:
        covariates=config.get('input',{}).get('covariates',''),
        controls_file=config.get('input',{}).get('controls','')
    resources:
        lsf_err="logs/lsf/final_lists.e",
        lsf_out="logs/lsf/final_lists.o"
    run:
        with open(input.exclusions1,'r') as e1:
            exclusions = e1.read().splitlines()
        # Optional (select_variants2-style); omit if not used
        het_miss_path = "data/qc/exclusions/het_missingness_failures.txt"
        if os.path.isfile(het_miss_path):
            with open(het_miss_path,'r') as e2:
                exclusions.extend(e2.read().splitlines())

        # Controls: from covariate file (col2=IID, col4=STATUS) if present, else separate controls file
        controls = set()
        if params.covariates and os.path.isfile(params.covariates):
            with open(params.covariates,'r') as f:
                lines = f.readlines()
            if len(lines) > 1:
                for line in lines[1:]:
                    parts = line.strip().split()
                    if len(parts) >= 4 and parts[3] in ('1', 'Control', 'control', 'CONTROL'):
                        controls.add(parts[1])  # IID (STATUS 1=control)
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



#NOTE I could exclude the mnp pairs from vep annotated files, and just annotete the mnp pairs + mnp and concat and sort...

#Negative Controls: synonymous snps, but we should make sure none of them appear in the mnp files.
#We are only keep/correcting mnps that have at least one case.
#So no snp that are in the mnp PASS files.
#But we should only pull out mnp pair that are validated.
#The non validated ones can stay as snps, but should be excluded from 
#the synonymous snps into regenie.