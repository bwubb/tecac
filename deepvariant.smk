

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)
        os.makedirs(f'tmp/{sample}',exist_ok=True)

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def bam_input(wildcards):
    return BAMS[wildcards.sample]

rule all_deepvariant:
    input:
        expand("data/work/{sample}/deep_variant.output.vcf.gz",sample=SAMPLES)

rule stats_filter:
    input:
        "vcf_stats/stats_filter.txt"

rule run_deepvariant:
    input:
        bam=bam_input
    output:
        vcf="data/work/{sample}/deep_variant.output.vcf.gz",
        g_vcf="data/work/{sample}/deep_variant.output.g.vcf.gz"
    params:
        bed=config['resources']['targets_bed'],
        ref=config['reference']['fasta'],
        tmp="tmp/{sample}",
        outdir="data/work/{sample}"
    threads:
        4
    shell:
        """
        mkdir -p {params.outdir}

        singularity run -B /usr/lib/locale/:/usr/lib/locale/ \
          docker://google/deepvariant:"1.4.0" \
          /opt/deepvariant/bin/run_deepvariant \
          --model_type=WES \
          --ref={params.ref} \
          --reads={input.bam} \
          --regions={params.bed} \
          --output_vcf={output.vcf} \
          --output_gvcf={output.g_vcf} \
          --intermediate_results_dir={params.tmp} \
          --num_shards={threads}
        """

rule run_bcftools_stats:
    input:
        "data/work/{sample}/deep_variant.output.vcf.gz"
    output:
        "data/work/{sample}/deep_variant.pass.bcftools_stats.output"
    shell:
        """
        bcftools view -i \"%FILTER='PASS'\" -s {wildcards.sample} {input} | bcftools stats -s- > {output}
        """
        #There appears a bug in bcftools stats.
        #If I use the -i EXPRESSION option it will not count the PSI and PSC values.

rule run_qc_stats:
    input:
        expand("data/work/{sample}/deep_variant.pass.bcftools_stats.output",sample=SAMPLES)
    output:
        "vcf_stats/qc_stats.json",
        "vcf_stats/qc_stats.csv",
        "vcf_stats/stats_filter.txt"
    shell:
        """
        python qc_stats.py {input}
        """
