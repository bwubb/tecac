with open('qcsamples.list','r') as i:
    QCSAMPLES=i.read().splitlines()
    for sample in QCSAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)
        os.makedirs(f'tmp/{sample}',exist_ok=True)


rule merge_cohort:
    input:
        expand("data/cohort/chr{chr}.output.vcf.gz",chr=list(range(1,23))+['X','Y'])

rule split_vcf:
    input:
        "data/work/{sample}/deep_variant.output.vcf.gz"
    output:
        "data/work/{sample}/chr{chr}.output.vcf.gz",
        "data/work/{sample}/chr{chr}.output.vcf.gz.csi"
    shell:
        """
        bcftools view -f PASS --type snps,indels,mnps -O z -o {output[0]} {input} {wildcards.chr}
        bcftools index {output}
        """

rule merge_vcf:
    input:
        expand("data/work/{sample}/chr{chr}.output.vcf.gz",sample=QCSAMPLES,chr=f"{wildcards.chr}")
    output:
        "data/cohort/chr{chr}.output.vcf.gz"
    shell:
        """
        bcftools merge -m none -O z -o {output} {input}
        """
