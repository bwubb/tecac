




rule gene_variant_stats:
    input:
    output:
        "data/work/bcftools/{gene}/variant_coverage.stats"
    params:
        region=lambda x: GENE[wildcards.gene]
    shell:
        """
        bcftools query -r {params.region} -f '%CHROM,%POS,%REF,%ALT,%INFO/MAF,%INFO/NS[,%SAMPLE:%DP]' {input} > {output}
        """

rule gene_coverage_qc:
    input:
    output:
    params:
    shell:

rule subset_pgen:
    input:
    output:
    params:
    shell:
        """
        plink --pfile {input} --chr {params.chr} --from-bp {params.from_bp} --to-bp {params.to_bp} --make-pgen --out {output}
        """

rule subset_regenie_files:
    input:
    output:

    params:
    shell:
        """
        grep {wildcards.gene} {input[0]} > {output[0]}
        grep {wildcards.gene} {input[1]} > {output[1]}
        """
