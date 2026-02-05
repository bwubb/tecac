from datetime import datetime

#Raw WES VCF File: TECAC_WES.vcf.gz
#Covariate File: TECAC_WES.covar.txt
#Phenotype File: TECAC_WES.TECAC.txt
#Annotation File (VEP-generated): TECAC_WES.PV.annotation.txt

wildcard_constraints:
    chrom='[1-22]'

rule run_regenie:
    input:
        expand("preprocess/TECAC_WES.chr{chrom}.vcf.gz",chrom=list(range(1,23))),
        expand("preprocess/TECAC_WES.chr{chrom}.pgen",chrom=list(range(1,23))),
        "preprocess/TECAC_WES.merge.step1.pgen",
        "TECAC_WES.regenie.annotation.txt",
        "TECAC_WES.regenie.set.txt",
        "TECAC_WES.regenie.mask.txt",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "preprocess/TECAC_WES.merge.step1.snplist"

rule run_report:
    input:
        expand("TECAC_{date}_regenie_report.html",date=datetime.now().strftime("%Y%m%d"))

#TECAC_WES to be replaced with config name
rule preprocess_vcf:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.bcf"
    output:
        "preprocess/TECAC_WES.chr{chrom}.bcf"
    shell:
        "bcftools view -r chr{wildcards.chrom} {input} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ob -o {output}"
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
rule preprocess_bcf:
    input:
        "preprocess/TECAC_WES.chr{chrom}.bcf"
    output:
        "preprocess/TECAC_WES.chr{chrom}.pgen"
    params:
        out_basename="preprocess/TECAC_WES.chr{chrom}",
        sex_info="sex_info.txt"
    shell:
        "plink2 --bcf {input} --double-id --make-pgen --update-sex {params.sex_info} --split-par b38 --vcf-half-call reference --out {params.out_basename}"

rule preprocess_variant_list:
    input:
        "preprocess/TECAC_WES.chr{chrom}.pgen"
    output:
        "preprocess/TECAC_WES.chr{chrom}.step1.pgen",
        "preprocess/TECAC_WES.chr{chrom}.step1.snplist"
    params:
        pfile="preprocess/TECAC_WES.chr{chrom}",
        out="preprocess/TECAC_WES.chr{chrom}.step1"
    shell:
        "plink2 --pfile {params.pfile} --double-id --indep-pairwise 500 50 0.4 --maf 0.05 --snps-only --geno 0.1 --hwe 1e-6 --make-pgen --write-snplist --out {params.out}"

rule merge_chr_pgen:
    input:
        expand("preprocess/TECAC_WES.chr{chrom}.step1.pgen",chrom=list(range(1,23)))
    output:
        "preprocess/TECAC_WES.merge.step1.pgen"
    params:
        out="preprocess/TECAC_WES.merge.step1"
    shell:
        "plink2 --pfile preprocess/TECAC_WES.chr1.step1 --pmerge-list file_list.txt --make-pgen --out {params.out}"

rule merge_chr_snplist:
    input:
        expand("preprocess/TECAC_WES.chr{chrom}.step1.snplist",chrom=list(range(1,23)))
    output:
        "preprocess/TECAC_WES.merge.step1.snplist"
    run:
        with open(output[0],'w') as f:
            for i in range(1,23):
                with open(f"preprocess/TECAC_WES.chr{i}.step1.snplist",'r') as f2:
                    f.write(f2.read())

rule preprocess_regenie:
    input:
        vep_csv="data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vep.report.csv",
        sites="seq_sites.txt",
        controls="controls.txt",
        covariates="covariates.txt"
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
        "--sites {input.sites} "
        "--controls {input.controls} "
        "-O {params.output_prefix}"


rule run_step1_regenie:
    input:
        "preprocess/TECAC_WES.merge.step1.pgen",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "preprocess/TECAC_WES.merge.step1.snplist"
    output:
        "TECAC_WES.step1_1.loco.gz",
        "TECAC_WES.step1_pred.list"
    params:
        pgen="preprocess/TECAC_WES.merge.step1",
        out="TECAC_WES.step1"
    shell:
        "regenie --step 1 --pgen {params.pgen} --covarFile {input[1]} --phenoFile {input[2]} --extract {input[3]} "
        "--bsize 1000 --gz --bt --lowmem --lowmem-prefix tmp_rg_ --out {params.out}"

#Example broke this down by chrom, useing chrom pgen.
rule run_step2_single_variant:
    input:
        "preprocess/TECAC_WES.chr{chrom}.pgen",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "TECAC_WES.step1_pred.list"
    output:
        "TECAC_WES.chr{chrom}.step2_single_variant_STATUS.regenie"
    params:
        pgen="preprocess/TECAC_WES.chr{chrom}",
        out="TECAC_WES.chr{chrom}.step2_single_variant"
    shell:
        "regenie --step 2 --pgen {params.pgen} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[3]} --bsize 400 --af-cc "
        "--out step2_single_variant"


rule run_step2_gene_based:
    input:
        "preprocess/TECAC_WES.chr{chrom}.pgen",
        "TECAC_WES.regenie.covar.txt",
        "TECAC_WES.regenie.pheno.txt",
        "TECAC_WES.regenie.annotation.txt",
        "TECAC_WES.regenie.set.txt",
        "TECAC_WES.regenie.mask.txt",
        "TECAC_WES.step1_pred.list"
    output:
        "TECAC_WES.chr{chrom}.step2_gene_based_STATUS.regenie"
    params:
        pgen="preprocess/TECAC_WES.chr{chrom}",
        out="TECAC_WES.chr{chrom}.step2_gene_based"
    shell:
        "regenie --step 2 --pgen {params.pgen} --phenoFile {input[2]} --covarFile {input[1]} --bt "
        "--firth --approx --pThresh 0.999 --firth-se --pred {input[6]} --anno-file {input[3]} "
        "--set-list {input[4]} --mask-def {input[5]} --build-mask 'max' --write-mask-snplist "
        "--check-burden-files --af-cc --bsize 200 --vc-tests skat,skato --out {params.out}"

rule report_regenie:
    input:
        "TECAC_WES.regenie.pheno.txt",
        "TECAC_WES.chr{chrom}.step2_single_variant_STATUS.regenie",
        "TECAC_WES.chr{chrom}.step2_gene_based_STATUS.regenie"
    output:
        "TECAC_{date}_regenie_report.html"
    shell:
        """
        Rscript -e "rmarkdown::render('report_regenie.Rmd',output_file='{output}',params=list(pheno_file='TECAC_WES.regenie.pheno.txt'))"
        """
