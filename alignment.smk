import os
import yaml

### FUNCTIONS ###

def map_input(wildcards):
    inputs=[]
    for RUN,_run in list(FILES[wildcards.sample].items()):
        run,lane,index=_run['PU'].split('-',2)
        inputs.append(f'bam_input/work/{wildcards.sample}/{wildcards.reference}/{run}/{lane}/{index}/5.markdup.bam')
    assert len(inputs)>0
    return sorted(inputs)

def get_fastqs(wildcards):
    return {'R1':'FASTQ/'+FILES[wildcards.sample][f'{wildcards.run}-{wildcards.lane}']['files'][0],'R2':'FASTQ/'+FILES[wildcards.sample][f'{wildcards.run}-{wildcards.lane}']['files'][1]}
    #fastq.yaml does not currently include file paths.

### ### PYTHON ### ###

with open(config['project']['fastq_config']) as file:
    FILES=yaml.load(file,Loader=yaml.BaseLoader)
    #SAMPLES=sorted(list(FILES.keys()))

#This allows me to subset the sample list easier.
with open(config['project']['sample_list']) as file:
    SAMPLES=file.read().splitlines()
    for sample in SAMPLES:
            os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

### ### ### RULES ### ### ###

localrules:aln_all,bam_table,write_bam_table

rule aln_all:
    input:
        expand("bam_input/final/{sample}/{sample}.{reference}.bam",sample=SAMPLES,reference=config['reference']['key'])

rule bam_table:
    input:
        "bam.table"

rule bwa_mem:
    input:
        unpack(get_fastqs)
    output:
        temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam"),
    params:
        fasta=config['reference']['fasta']
    threads:
        4
    shell:
        #fastq config has no path.
        """
        bwa mem -M -t {threads} {params.fasta} {input.R1} {input.R2} | samtools view -bS -o {output}
        """

rule samtools_readgroup:
    input:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam"
    output:
        rg=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/2.readgroup.bam"),
        fm=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/3.fixmate.bam"),
        qs=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/4.qsort.bam"),
    params:
        LB=config['resources']['library_key']
    threads:
        4
    shell:
        """
        samtools addreplacerg -@ {threads} -r 'ID:{wildcards.run}.{wildcards.lane}' -r 'PU:{wildcards.run}.{wildcards.lane}.{wildcards.index}' -r 'PL:illumina' -r 'LB:{params.LB}' -r 'SM:{wildcards.sample}' -o {output.rg} {input}
        samtools fixmate -m -@ {threads} {output.rg} {output.fm}
        samtools sort -@ {threads} -o {output.qs} {output.fm}
        """

rule samtools_markdup:
    input:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/4.qsort.bam"
    output:
        temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/5.markdup.bam")
    params:
        stats="bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/5.stats.txt"
    threads:
        4
    shell:
        """
        samtools markdup -s -f {params.stats} -@ {threads} {input} {output}
        """

rule input_ready:
    input:
        map_input
    output:
        temp("bam_input/work/{sample}/{reference}/input.bam")
    threads:
        4
    shell:
        """
        samtools merge -f -@ {threads} {output} {input}
        samtools index {output}
        """

rule ValidateSamFile:
    input:
        "bam_input/work/{sample}/{reference}/input.bam"
    output:
        "metrics/{reference}/{sample}/validation_data.table"
    params:
        memory="10240m"
    shell:
        """
        set +e
        exitcode=$?
        java -Xmx{params.memory} -jar $HOME/software/picard/2.20.7/picard.jar ValidateSamFile I={input} O={output} MODE=SUMMARY
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

rule ready_bam:
    input:
        bam="bam_input/work/{sample}/{reference}/input.bam",
        table="metrics/{reference}/{sample}/validation_data.table"
    output:
        "bam_input/final/{sample}/{sample}.{reference}.bam"
    shell:
        """
        rsync -v {input.bam} {output}
        samtools index {output}
        """

rule write_bam_table:
    input:
        expand("bam_input/final/{sample}/{sample}.{reference}.bam",sample=SAMPLES,reference=config['reference']['key'])
    output:
        "bam.table"
    params:
        ref=config['reference']['key']
    run:
        with open(output[0],'w') as file:
            for sample in SAMPLES:
                file.write(f"{sample}\tbam_input/final/{sample}/{sample}.{params.ref}.bam\n")
