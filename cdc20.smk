import os
import csv


with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

localrules:merge_csv

rule all_cdc20:
    input:
        expand("{project}.deep_variant.CDC20.PASS.vep.report.csv",project=config['project']['name'])

rule view_cdc20:
    input:
        "data/work/{sample}/deep_variant.output.vcf.gz"
    output:
        u="analysis/work/{sample}/CDC20/deep_variant.CDC20.vcf.gz",
        p="analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vcf.gz",

    shell:
        """
        tabix -fp vcf {input}
        bcftools view -Oz -o {output.u} {input} 1:43358980-43363203
        bcftools view -i \"%FILTER='PASS'\" -Oz -o {output.p} {output.u}
        """

rule annotate_vep:
    input:
        "analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vcf.gz"
    output:
        "analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vep.vcf.gz"
    params:
        in_vcf=temp('analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vcf'),
        out_vcf='analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vep.vcf',
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='$HOME/.vep/Plugins/loftee',#check
        utr=config['resources']['utrannotator'],
    shell:
        """
        bcftools view -O v -o {params.in_vcf} {input}

        vep -i {params.in_vcf} -o {params.out_vcf} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf \
        --everything \
        --canonical \
        --assembly {params.assembly} \
        --species homo_sapiens \
        --fasta {params.fa} \
        --plugin NMD \
        --plugin Downstream \
        --plugin REVEL,{params.revel} \
        --plugin SpliceAI,snv={params.splice_snv},indel={params.splice_indel} \
        --plugin gnomADc,{params.gnomAD} \
        --plugin LoF,loftee_path:{params.loftee} \
        --plugin UTRannotator,{params.utr} \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

        bgzip {params.out_vcf} && tabix -fp vcf {output}
        """

rule vep_report_cdc20:
    input:
        "analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vep.vcf.gz"
    output:
        "analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vep.report.csv"
    shell:
        """
        python vep_vcf_parser.py -i {input} -o {output} --mode single,{wildcards.sample} everything
        """

rule merge_csv:
    input:
        expand("analysis/work/{sample}/CDC20/deep_variant.CDC20.PASS.vep.report.csv",sample=SAMPLES)
    output:
        "{project}.deep_variant.CDC20.PASS.vep.report.csv"
    run:
        for n,f in enumerate(input):
            print(f"File {n+1} of {len(input)}")
            with open(f,'r') as infile:
                reader=csv.DictReader(infile,delimiter=',')
                if n==0:
                    outfile=open(f"{output[0]}",'w')
                    writer=csv.DictWriter(outfile,delimiter=',',fieldnames=reader.fieldnames)
                    writer.writeheader()
                for row in reader:
                    writer.writerow(row)
        else:
            outfile.close()
