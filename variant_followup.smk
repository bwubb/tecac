
import csv

def _load_gene_chr_map(path):
    gene_to_chr={}
    with open(path,'r',newline='') as f:
        r=csv.reader(f)
        next(r,None)
        for row in r:
            if len(row)<2:
                continue
            vid=row[0].strip()
            gene=row[1].strip()
            if not vid or not gene:
                continue
            chrom=vid.split('_',1)[0].replace('chr','')
            if gene in gene_to_chr and gene_to_chr[gene]!=chrom:
                raise ValueError(f"Gene {gene} maps to multiple chromosomes: {gene_to_chr[gene]}, {chrom}")
            gene_to_chr[gene]=chrom
    return gene_to_chr

_GENE_CHR_MAP=None
def _bcf_for_gene(wildcards):
    global _GENE_CHR_MAP
    if _GENE_CHR_MAP is None:
        _GENE_CHR_MAP=_load_gene_chr_map("data/regenie/pathogenic_vus.csv")
    gene=wildcards.GENE
    if gene not in _GENE_CHR_MAP:
        raise ValueError(f"GENE {gene} not found in data/regenie/pathogenic_vus.csv")
    chrom=_GENE_CHR_MAP[gene]
    return f"data/bcftools/chr{chrom}.qc_filter1.bcf"

#We need to load some sort of file to run followup on.
#Which will come from regenie/meta analysis.

#Save for production use
#with open('','r') as f:
#    GENES=file.read().splitlines()

GENES=["CHEK2"]#testing

rule variant_followup:
    input:
        "data/variant/variant-followup-report.html",
        expand("data/variant/variants-pathogenic-summary.{GENE}.tsv",GENE=GENES),
        expand("data/variant/variants-vus-summary.{GENE}.tsv",GENE=GENES)
    

rule awk_gene_variants:
    input:
        "data/regenie/pathogenic_vus.csv"
    output:
        pathogenic="data/variant/variants-pathogenic.{GENE}.txt",
        vus="data/variant/variants-vus.{GENE}.txt"
    shell:
        """
        awk -F',' 'NR>1 {{for(i=1;i<=NF;i++) gsub(/^"|"$/,"",$i); if($2=="{wildcards.GENE}" && $3=="1") print $1}}' {input} > {output.pathogenic}
        awk -F',' 'NR>1 {{for(i=1;i<=NF;i++) gsub(/^"|"$/,"",$i); if($2=="{wildcards.GENE}" && $3=="2") print $1}}' {input} > {output.vus}
        """

#qc_filter1 is fine to use. It has the IDS and the passing samples are from post het_miss sample filtering.
rule query_variant_genotypes:
    input:
        bcf=_bcf_for_gene,
        samples="data/qc/passing_samples.txt",
        pathogenic="data/variant/variants-pathogenic.{GENE}.txt",
        vus="data/variant/variants-vus.{GENE}.txt"
    output:
        header="data/variant/samples.{GENE}.tsv",
        pathogenic_gt="data/variant/variants-pathogenic-genotypes.{GENE}.tsv",
        vus_gt="data/variant/variants-vus-genotypes.{GENE}.tsv"
    shell:
        """
        bcftools view -h -S {input.samples} {input.bcf} | awk -F'\\t' '$1=="#CHROM" {{printf "ID\\tTYPE"; for(i=10;i<=NF;i++) printf "\\t%s",$i; printf "\\n"}}' > {output.header}
        
        cat {output.header} > {output.pathogenic_gt}
        bcftools query -S {input.samples} -i 'ID=@{input.pathogenic}' -f '%ID\\t%TYPE[\\t%GT]\\n' {input.bcf} >> {output.pathogenic_gt}
        
        cat {output.header} > {output.vus_gt}
        bcftools query -S {input.samples} -i 'ID=@{input.vus}' -f '%ID\\t%TYPE[\\t%GT]\\n' {input.bcf} >> {output.vus_gt}
        """

rule parse_variant_genotypes:
    input:
        covariates=config.get('input',{}).get('covariates',''),
        pathogenic_gt="data/variant/variants-pathogenic-genotypes.{GENE}.tsv",
        vus_gt="data/variant/variants-vus-genotypes.{GENE}.tsv"
    output:
        pathogenic_report="data/variant/variants-pathogenic-summary.{GENE}.tsv",
        vus_report="data/variant/variants-vus-summary.{GENE}.tsv"
    shell:
        """
        python parse_variant_genotypes.py -c {input.covariates} -p {input.pathogenic_gt} -v {input.vus_gt} -P {output.pathogenic_report} -V {output.vus_report}
        """

rule variant_followup_report:
    input:
        #Gene list file
        expand("data/variant/variants-pathogenic-summary.{GENE}.tsv",GENE=GENES),
        expand("data/variant/variants-vus-summary.{GENE}.tsv",GENE=GENES)
    output:
        report="data/variant/variant-followup-report.html"
    shell:
        """
        Rscript -e "rmarkdown::render('variant_followup.Rmd', output_file='{output.report}')"
        """
    
