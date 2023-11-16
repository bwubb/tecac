
##INIT
with open('pass_samples.list','r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

rule all:
    input: expand("data/work/{sample}/deep_variant.output.vep.report.csv",sample=SAMPLES)

rule annotate_only:
    input: expand("data/work/{sample}/deep_variant.output.vep.vcf.gz",sample=SAMPLES)

rule deepvariant_vep:
    input:
        "data/work/{sample}/deep_variant.output.vcf.gz"
    output:
        "data/work/{sample}/deep_variant.output.vep.vcf.gz"
    params:
        out_vcf="data/work/{sample}/deep_variant.output.vep.vcf",
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='$HOME/.vep/Plugins/loftee',#check
        human_ancestor_fa='$HOME/.vep/Plugins/loftee/human_ancestor.fa.gz',
        conservation_file='$HOME/.vep/Plugins/loftee/loftee.sql',
        gerp_bigwig='$HOME/.vep/Plugins/loftee/gerp_conservation_scores.homo_sapiens.GRCh38.bw',
        utr=config['resources']['utrannotator'],
        alphamissense='$HOME/.vep/alphamissense/AlphaMissense_hg38.tsv.gz',
        mavedb='$HOME/.vep/mavedb/MaveDB_variants.tsv.gz'
    shell:
        """
        vep -i {input} -o {params.out_vcf} \
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
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,{params.revel} \
        --plugin SpliceAI,snv={params.splice_snv},indel={params.splice_indel} \
        --plugin gnomADc,{params.gnomAD} \
        --plugin UTRannotator,{params.utr} \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
        --plugin AlphaMissense,file={params.alphamissense} \
        --plugin MaveDB,file={params.mavedb}

        bgzip {params.out_vcf}
        tabix -fp vcf {output}
        """
#--plugin Downstream \
#--plugin LoF,loftee_path:{params.loftee},human_ancestor_fa:{params.human_ancestor_fa},conservation_file:{params.conservation_file},gerp_bigwig:{params.gerp_bigwig} \

rule vep_report:
    input:
        vcf="data/work/{sample}/deep_variant.output.vep.vcf.gz"
    output:
        csv="data/work/{sample}/deep_variant.output.vep.report.csv"
    shell:
        """
        python vep_vcf_parser.py -i {input.vcf} -o {output.csv} --mode single,{wildcards.sample} everything
        """
