


#new realignment config
with open(config.get('sample_list','sample.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('bam_input','bam.input'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def bam_input(wildcards):
    return BAMS[wildcards.sample]

rule revert_all:
    input:
        expand("FASTQ/{sample}_{lib}_{run}_{lane}_{index}_R1.fastq.gz",sample=SAMPLES,lib="TruSeq",run="ICR",lane="1",index="NNNNNNNN")

rule samtools_query_sort:
    input:
        bam_input
    output:
        temp("bam_input/work/{sample}/samtools/qsort.bam")
    shell:
        "samtools sort -n -o {output} {input}"

rule samtools_bam2fq:
    input:
        "bam_input/work/{sample}/samtools/qsort.bam"
    output:
        R1="bam_input/work/{sample}/samtools/R1.fastq.gz",
        R2="bam_input/work/{sample}/samtools/R2.fastq.gz"
    shell:
        "samtools fastq -1 {output.R1} -2 {output.R2} {input}"

rule repair_fastq:
    input:
        R1="bam_input/work/{sample}/samtools/R1.fastq.gz",
        R2="bam_input/work/{sample}/samtools/R2.fastq.gz"
    output:
        R1="FASTQ/{sample}_{lib}_{run}_{lane}_{index}_R1.fastq.gz",
        R2="FASTQ/{sample}_{lib}_{run}_{lane}_{index}_R2.fastq.gz"
    params:
        memory="61440m"
    shell:
        "repair.sh -Xmx{params.memory} in1={input.R1} in2={input.R2} out1={output.R1} out2={output.R2}"
