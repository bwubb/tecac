


#new realignment config
with open(config.get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config['bam_input'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def bam_input(wildcards):
    return BAMS[wildcards.sample]

rule revert_all:
    input:
        expand("FASTQ/{sample}_{lib}_{run}_{lane}_{index}_R1.fastq.gz",sample=SAMPLES,lib="S07604514",run="CNIO-WES",lane="1",index="NNNNNNNN"),
        expand("FASTQ/{sample}_{lib}_{run}_{lane}_{index}_R2.fastq.gz",sample=SAMPLES,lib="S07604514",run="CNIO-WES",lane="1",index="NNNNNNNN")

rule samtools_query_sort:
    input:
        bam_input
    output:
        temp("bam_input/work/{sample}/samtools/qsort.bam")
    shell:
        "samtools sort -n {input} -o {output}"

rule samtools_bam2fq:
    input:
        "bam_input/work/{sample}/samtools/qsort.bam"
    output:
        R1="FASTQ/{sample}_{lib}_{run}_{lane}_{index}_R1.fastq.gz",
        R2="FASTQ/{sample}_{lib}_{run}_{lane}_{index}_R2.fastq.gz"
    shell:
        "samtools fastq -1 {output.R1} -2 {output.R2} {input}"
