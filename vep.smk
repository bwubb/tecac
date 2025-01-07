
##INIT
with open('pass_samples.list','r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

rule all:
    input: expand("data/work/{sample}/deep_variant.output.vep.report.csv",sample=SAMPLES)

#rule annotate_only:
#    input: expand("data/work/{sample}/deep_variant.output.vep.vcf.gz",sample=SAMPLES)

rule bcftools_norm:
    input:
        "{sample}.deep_variant.output.vcf.gz"
    output:
        "{sample}.deep_variant.output.norm.vcf.gz"
    shell:
        "bcftools norm -m -both {input} | bcftools view -a -O z -o {output}"

rule deepvariant_vep:
    input:
        "data/work/{sample}/deep_variant.output.norm.vcf.gz"
    output:
        "data/work/{sample}/deep_variant.output.norm.vep.vcf.gz"
    params:
        out_vcf="data/work/{sample}/deep_variant.output.norm.vep.vcf",
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='/home/bwubb/.vep/Plugins/loftee',#check
        human_ancestor_fa='/home/bwubb/.vep/Plugins/loftee/human_ancestor.fa.gz',
        conservation_file='/home/bwubb/.vep/Plugins/loftee/loftee.sql',
        gerp_bigwig='/home/bwubb/.vep/Plugins/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        utr=config['resources']['utrannotator'],
        alphamissense='/home/bwubb/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz',
        mavedb='/home/bwubb/.vep/mavedb/MaveDB_variants.tsv.gz'
    shell:
        """
        singularity run -H $PWD:/home \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep111_021324.sif vep \
        --dir /opt/vep/.vep \
        -i {input} \
        -o {params.out_vcf} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta {params.fa} \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
        --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
        --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz

        bgzip {params.out_vcf}
        tabix -fp vcf {output}
        """
#--plugin Downstream \
#--plugin LoF,loftee_path:{params.loftee},human_ancestor_fa:{params.human_ancestor_fa},conservation_file:{params.conservation_file},gerp_bigwig:{params.gerp_bigwig} \

rule vep_report:
    input:
        vcf="data/work/{sample}/deep_variant.output.norm.vep.vcf.gz"
    output:
        csv="data/work/{sample}/deep_variant.output.norm.vep.report.csv"
    shell:
        """
        python vep_vcf_parser.py -i {input.vcf} -o {output.csv} --mode single,{wildcards.sample} everything
        """
