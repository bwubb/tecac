

rule all:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.bcf",
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vep.report.csv",
        "vcf_stats/qc_stats.csv",

rule norm_output:
    input:
        "glnexus.bcf"
    output:
        "data/work/glnexus/glnexus.norm.bcf"
    shell:
        "bcftools norm -m-both {input} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -W=csi -Ob -o {output}"

rule split_sex_chr:
    input:
        "data/work/glnexus/glnexus.norm.bcf"
    output:
        "data/work/glnexus/glnexus.norm.autosomes.bcf",
        "data/work/glnexus/glnexus.norm.sex_chr.bcf",
        "data/work/glnexus/glnexus.norm.par.bcf"
    params:
        par="$HOME/resources/Bed_files/par.bed"
    shell:
        """
        # 1. Autosomes (exclude chrX and chrY)
        bcftools view -t ^chrX,chrY {input} -W=csi -Ob -o {output[0]}

        # 2. Non-PAR chrX and chrY (haploid in males)
        bcftools view -T ^{params.par} -r chrX,chrY {input} -W=csi -Ob -o {output[1]}

        # 3. PAR regions (diploid in males)
        bcftools view -R {params.par} {input} -W=csi -Ob -o {output[2]}
        """

rule split_autosomes_by_chrom:
    input:
        "data/work/glnexus/glnexus.norm.autosomes.bcf"
    output:
        "data/work/glnexus/glnexus.norm.chr{chr}.split.bcf"
    shell:
        """
        bcftools view -r chr{wildcards.chr} -W=csi -Ob -o {output} {input}
        """

rule set_chromosome_gt:
    input:
        "data/work/glnexus/glnexus.norm.chr{chr}.split.bcf"
    output:
        "data/work/glnexus/glnexus.norm.chr{chr}.ref.bcf",
        "data/work/glnexus/glnexus.norm.chr{chr}.ref.het.bcf",
        "data/work/glnexus/glnexus.norm.chr{chr}.ref.het.hom.bcf",
        "data/work/glnexus/glnexus.norm.chr{chr}.ref.het.hom.gt.bcf"
    shell:
        """
        # 0/0
        bcftools +setGT {input} -Ob -o {output[0]} -- -t q -n c:'0/0' -i '(AD[:0] > 0 & AD[:1] == 0)'

        # 0/1
        bcftools +setGT {output[0]} -Ob -o {output[1]} -- -t q -n c:'0/1' -i '(AD[:0] > 0 & AD[:1] > 0)'

        # 1/1
        bcftools +setGT {output[1]} -Ob -o {output[2]} -- -t q -n c:'1/1' -i '(AD[:0] == 0 & AD[:1] > 0)'

        # Filter weak hets
        bcftools +setGT {output[2]} -W=csi -Ob -o {output[3]} -- -t 'b:AD<1e-3' -n c:'0/0'
        """

# Concatenate all autosome chromosomes back together
rule concat_autosomes:
    input:
        expand("data/work/glnexus/glnexus.norm.chr{chr}.ref.het.hom.gt.bcf",chr=list(range(1,23)))
    output:
        "data/work/glnexus/glnexus.norm.autosomes.ref.het.hom.gt.bcf"
    params:
        bcf=temp("data/work/glnexus/temp1.bcf"),
        tmpdir="data/work/glnexus/"
    shell:
        """
        bcftools concat -a {input} -Ob -o {params.bcf}
        bcftools sort -T {params.tmpdir} -W=csi -Ob -o {output} {params.bcf}
        """

rule set_par_gt:
    input:
        "data/work/glnexus/glnexus.norm.par.bcf"
    output:
        "data/work/glnexus/glnexus.norm.par.ref.bcf",
        "data/work/glnexus/glnexus.norm.par.ref.het.bcf",
        "data/work/glnexus/glnexus.norm.par.ref.het.hom.bcf",
        "data/work/glnexus/glnexus.norm.par.ref.het.hom.gt.bcf"
    shell:
        """
        # 0/0
        bcftools +setGT {input} -Ob -o {output[0]} -- -t q -n c:'0/0' -i '(AD[:0] > 0 & AD[:1] == 0)'

        # 0/1
        bcftools +setGT {output[0]} -Ob -o {output[1]} -- -t q -n c:'0/1' -i '(AD[:0] > 0 & AD[:1] > 0)'

        # 1/1
        bcftools +setGT {output[1]} -Ob -o {output[2]} -- -t q -n c:'1/1' -i '(AD[:0] == 0 & AD[:1] > 0)'

        # Filter weak hets
        bcftools +setGT {output[2]} -W=csi -Ob -o {output[3]} -- -t 'b:AD<1e-3' -n c:'0/0'
        """

#This is for plink specification.
rule set_sex_chr_gt:
    input:
        "data/work/glnexus/glnexus.norm.sex_chr.bcf"
    output:
        "data/work/glnexus/glnexus.norm.sex_chr.ref.bcf",
        "data/work/glnexus/glnexus.norm.sex_chr.ref.alt.bcf"
    shell:
        """
        # REF allele only → haploid 0
        bcftools +setGT {input} -Ob -o {output[0]} -- -t q -n c:'0' -i '(AD[:0] > 0 & AD[:1] == 0)'

        # ALT allele only → haploid 1
        bcftools +setGT {output[0]} -W=csi -Ob -o {output[1]} -- -t q -n c:'1' -i '(AD[:0] == 0 & AD[:1] > 0)'
        """

rule merge_genotypes:
    input:
        "data/work/glnexus/glnexus.norm.autosomes.ref.het.hom.gt.bcf",
        "data/work/glnexus/glnexus.norm.sex_chr.ref.alt.bcf",
        "data/work/glnexus/glnexus.norm.par.ref.het.hom.gt.bcf"
    output:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.bcf"
    params:
        bcf=temp("data/work/glnexus/temp2.bcf"),
        tmpdir="data/work/glnexus/"
    shell:
        """
        bcftools concat -a {input} -W=csi -Ob -o {params.bcf}
        bcftools sort -T {params.tmpdir} -W=csi -Ob -o {output} {params.bcf}
        """

rule fill_tags:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.bcf"
    output:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.bcf"
    shell:
        "bcftools +fill-tags {input} -W=csi -Ob -o {output} -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE"

#FORMAT/VAF     Number:A  Type:Float    ..  The fraction of reads with the alternate allele, requires FORMAT/AD or ADF+ADR
#INFO/AC        Number:A  Type:Integer  ..  Allele count in genotypes
#INFO/AC_Hom    Number:A  Type:Integer  ..  Allele counts in homozygous genotypes
#INFO/AC_Het    Number:A  Type:Integer  ..  Allele counts in heterozygous genotypes
#INFO/AC_Hemi   Number:A  Type:Integer  ..  Allele counts in hemizygous genotypes
#INFO/AF        Number:A  Type:Float    ..  Allele frequency
#INFO/AN        Number:1  Type:Integer  ..  Total number of alleles in called genotypes
#INFO/ExcHet    Number:A  Type:Float    ..  Test excess heterozygosity; 1=good, 0=bad
#INFO/END       Number:1  Type:Integer  ..  End position of the variant
#INFO/F_MISSING Number:1  Type:Float    ..  Fraction of missing genotypes (all samples, experimental)
#INFO/HWE       Number:A  Type:Float    ..  HWE test (PMID:15789306); 1=good, 0=bad
#INFO/MAF       Number:A  Type:Float    ..  Minor Allele frequency
#INFO/NS        Number:1  Type:Integer  ..  Number of samples with data
#INFO/TYPE      Number:.  Type:String   ..  The record type (REF,SNP,MNP,INDEL,etc)

rule variant_qc:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.bcf"
    output:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.bcf"
    shell:
        """
        bcftools filter -s ExcessHet -e 'INFO/ExcHet=0' -m + {input} |
        bcftools filter -s LowCallRate -e 'INFO/F_MISSING > 0.05' -m + |
        bcftools filter -s LowSampleCount -e 'INFO/NS < 1000' -m + |
        bcftools filter -s NoHQHet -e 'COUNT(FORMAT/GT="0/1" && FORMAT/DP>=10 && FORMAT/GQ>=20 && FORMAT/VAF>0.2)=0' -m + -W=csi -Ob -o {output}
        """

rule vcf_Stats:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.bcf"
    output:
        "vcf_stats/glnexus.norm.all.ref.het.hom.gt.tags.qc.stats"
    shell:
        "bcftools stats -f PASS -s - {input} > {output}"
        # added -f PASS
        #Maybe I should add a bed file for targets here? Requires testing.

rule qc_stats:
    input:
        "vcf_stats/glnexus.norm.all.ref.het.hom.gt.tags.qc.stats"
    output:
        "vcf_stats/qc_stats.csv",
        "vcf_stats/qc_report.html"
    params:
        info="meta4qc.csv",
        bed="xgen_plus_spikein.Covered.slop10.hg38.bed"#Hardcoded for now.
    shell:
        "python qc_stats.py --metadata {params.info} --bed {params.bed} {input}"

rule no_sample:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.bcf"
    output:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vcf"
    shell:
        "bcftools view -G -Ov -o {output} {input}"

#Filter PASS?
rule annotate_variants:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vcf"
    output:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vep.vcf"
    shell:
        """
        singularity run -H $PWD:/home \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i {input} \
        -o {output} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta /opt/vep/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.fasta \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
        --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
        --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar_20250106.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
        """

rule parse_vep:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vep.vcf"
    output:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vep.report.csv"
    params:
        blacklist="blacklist.txt"
    shell:
        "python vep_vcf_parser.py -i {input} -o {output} -b {params.blacklist} -m no_sample"

rule qc_variants:
    input:
        "data/work/glnexus/glnexus.norm.all.ref.het.hom.gt.tags.qc.no_sample.vep.report.csv"
    output:
        "variant_qc_report.html"
    shell:
        "python qc_variants.py {input} --json"
