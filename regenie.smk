from datetime import datetime

#Raw WES VCF File: TECAC_WES.vcf.gz
#Covariate File: TECAC_WES.covar.txt
#Phenotype File: TECAC_WES.TECAC.txt
#Annotation File (VEP-generated): TECAC_WES.PV.annotation.txt

rule run_regenie:
    input:
        expand("preprocess/TECAC_WES.chr{chrom}.vcf.gz", chrom=list(range(1,23))+['X','Y']),
        expand("preprocess/TECAC_WES.chr{chrom}.pgen", chrom=list(range(1,23))+['X','Y']),
        "preprocess/TECAC_WES.all_chr.pgen",
        "TECAC_WES.regenie.annotation.txt",
        "TECAC_WES.regenie.set.txt",
        "TECAC_WES.regenie.mask.txt",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "step1_1.loco.gz",
        "step2_single_variant_STATUS.regenie",
        "step2_gene_based_STATUS.regenie"

rule run_report:
    input:
        expand("TECAC_{date}_regenie_report.html",date=datetime.now().strftime("%Y%m%d"))

#TECAC_WES to be replaced with config name
rule preprocess_vcf:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.bcf"
    output:
        "preprocess/TECAC_WES.chr{chrom}.vcf.gz"
    shell:
        "bcftools view -r chr{wildcards.chrom} {input} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Oz -o {output}"
#This might change if we include -f PASS

#--vcf-half-call <mode>

#The current VCF standard does not specify how '0/.' and similar GT values should be interpreted.
#By default (mode 'error'/'e'), PLINK 1.9 errors out and reports the line number of the anomaly.
#Should the half-call be intentional, though (this can be the case with Complete Genomics data),
#you can request the following other modes:
#
#    'haploid'/'h': Treat half-calls as haploid/homozygous (the PLINK 1 file format does not distinguish between the two). This maximizes similarity between the VCF and BCF2 parsers.
#    'missing'/'m': Treat half-calls as missing.
#    'reference'/'r': Treat the missing part as reference.


#This can be --bcf instead of --vcf if we want to use the bcf file.
rule convert_to_plink:
    input:
        "preprocess/TECAC_WES.chr{chrom}.vcf.gz"
    output:
        "preprocess/TECAC_WES.chr{chrom}.pgen"
    params:
        out_basename="preprocess/TECAC_WES.chr{chrom}",
        sex_info="sex_info.txt"
    shell:
        "plink2 --vcf {input} --double-id --make-pgen --update-sex {params.sex_info} --split-par b38 --vcf-half-call reference --out {params.out_basename}"

rule merge_plink:
    input:
        expand("preprocess/TECAC_WES.chr{chrom}.pgen", chrom=list(range(1,23))+['X','Y'])
    output:
        "preprocess/TECAC_WES.all_chr.pgen",
        "preprocess/TECAC_WES.all_chr.eigenvec"
    params:
        out_basename="preprocess/TECAC_WES.all_chr"
    shell:
        "rm -f file_list.txt && "
        "for chrom in {{2..22}} X Y; do "
        "echo \"preprocess/TECAC_WES.chr${{chrom}}.pgen preprocess/TECAC_WES.chr${{chrom}}.pvar preprocess/TECAC_WES.chr${{chrom}}.psam\" >> file_list.txt; "
        "done && "
        "plink2 --pfile preprocess/TECAC_WES.chr1 --pmerge-list file_list.txt --make-pgen --pca --out {params.out_basename}"

rule preprocess_regenie:
    input:
        vep_csv="data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vep.report.csv",
        eigenvec="preprocess/TECAC_WES.all_chr.eigenvec",
        sites="seq_sites.txt",
        controls="controls.txt"
    output:
        annotation="TECAC_WES.regenie.annotation.txt",
        set_file="TECAC_WES.regenie.set.txt",
        mask="TECAC_WES.regenie.mask.txt",
        covar="TECAC_WES.regenie.covar.txt",
        pheno="TECAC_WES.regenie.pheno.txt"
    params:
        output_prefix="TECAC_WES"
    shell:
        "python preprocess_regenie.py "
        "--vep-csv {input.vep_csv} "
        "--eigenvec {input.eigenvec} "
        "--sites {input.sites} "
        "--controls {input.controls} "
        "--remove-PCs " # Remove PCs from eigenvec file
        "-O {params.output_prefix}"

rule step1_variant_list:
    input:
        "preprocess/TECAC_WES.all_chr.pgen"
    output:
        "step1_variant_list.snplist"
    params:
        input_basename="preprocess/TECAC_WES.all_chr",
        output_basename="step1_variant_list"
    shell:
        "plink2 --pfile {params.input_basename} --chr 1-22 X Y --maf 0.01 --snps-only --geno 0.1 --hwe 1e-15 --write-snplist --out {params.output_basename}"

rule run_step1_regenie:
    input:
        "preprocess/TECAC_WES.all_chr.pgen",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "step1_variant_list.snplist"
    output:
        "step1_1.loco.gz",
        "step1_pred.list"
    params:
        input_basename="preprocess/TECAC_WES.all_chr"
    shell:
        "regenie --step 1 --pgen {params.input_basename} --covarFile {input[1]} --phenoFile {input[2]} --extract {input[3]} "
        "--bsize 1000 --gz --bt --lowmem --lowmem-prefix tmp_rg_ --out step1"

#Example broke this down by chrom, useing chrom pgen.
rule run_step2_single_variant:
    input:
        "preprocess/TECAC_WES.all_chr.pgen",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "step1_pred.list"
    output:
        "step2_single_variant_STATUS.regenie"
    params:
        input_basename="preprocess/TECAC_WES.all_chr"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc "
        "--out step2_single_variant"


rule run_step2_gene_based:
    input:
        "preprocess/TECAC_WES.all_chr.pgen",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "TECAC_WES.regenie.annotation.txt",
        "TECAC_WES.regenie.set.txt",
        "TECAC_WES.regenie.mask.txt",
        "step1_pred.list"
    output:
        "step2_gene_based_STATUS.regenie"
    params:
        input_basename="preprocess/TECAC_WES.all_chr"
    shell:
        "regenie --step 2 --pgen {params.input_basename} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato --out step2_gene_based"

rule report_regenie:
    input:
        "TECAC_WES.regenie.pheno.txt",
        "step2_single_variant_STATUS.regenie",
        "step2_gene_based_STATUS.regenie"
    output:
        "TECAC_{date}_regenie_report.html"
    shell:
        """
        Rscript -e "rmarkdown::render('report_regenie.Rmd',output_file='{output}')"
        """
