import os

CHROMOSOMES_AUTOSOMAL=list(range(1,23))
CHROMOSOMES_ALL=list(range(1,23))+['X','Y']

rule all:
    input:
        expand("data/work/glnexus/glnexus.norm.chr{chr}.ref.het.hom.gt.tags.bcf",chr=CHROMOSOMES_AUTOSOMAL),
        "data/work/glnexus/glnexus.norm.chrX.ref.het.hom.gt.tags.bcf",
        "data/work/glnexus/glnexus.norm.chrY.ref.hom.tags.bcf",
        expand("vcf_stats/glnexus.norm.chr{chr}.ref.het.hom.gt.tags.stats",chr=CHROMOSOMES_AUTOSOMAL),
        "vcf_stats/glnexus.norm.chrX.ref.het.hom.gt.tags.stats",
        "vcf_stats/glnexus.norm.chrY.ref.hom.tags.stats"

rule run_glnexus:
    input:
        "gvcf_list.txt"
    output:
        "data/work/glnexus/glnexus.bcf",
        "data/work/glnexus/glnexus.bcf.csi"
    shell:
        """
        glnexus -t 16 -b {input} -O {output[0]}
        bcftools index {output[0]}
        """

rule norm_output:
    input:
        "data/work/glnexus/glnexus.bcf"
    output:
        "data/work/glnexus/glnexus.norm.bcf",
    params:
        targets=config["targets"]["bed"]
    shell:
        "bcftools norm -R {params.targets} -m-both {input} | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -W=csi -Ob -o {output}"

rule split_chromosomes:
    input:
        "data/work/glnexus/glnexus.norm.bcf"
    output:
        "data/work/glnexus/glnexus.norm.autosomes.bcf",
        "data/work/glnexus/glnexus.norm.chrX.bcf",
        "data/work/glnexus/glnexus.norm.chrY.bcf"
    params:
        par="$HOME/resources/Bed_files/par.bed"
    shell:
        """
        # 1. Autosomes (chr1-22)
        bcftools view -t ^chrX,chrY {input} -W=csi -Ob -o {output[0]}

        # 2. Non-PAR chrX (exclude PAR regions)
        bcftools view -T ^{params.par} -r chrX {input} -W=csi -Ob -o {output[1]}

        # 3. Non-PAR chrY (exclude PAR regions)
        bcftools view -T ^{params.par} -r chrY {input} -W=csi -Ob -o {output[2]}
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

# Fill tags for autosomal chromosomes
rule fill_tags_autosomes:
    input:
        "data/work/glnexus/glnexus.norm.chr{chr}.ref.het.hom.gt.bcf"
    output:
        "data/work/glnexus/glnexus.norm.chr{chr}.ref.het.hom.gt.tags.bcf"
    shell:
        "bcftools +fill-tags {input} -W=csi -Ob -o {output} -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE"


# Diploid genotype setting for chrX (non-PAR, handles both males and females)
rule set_chrX_gt:
    input:
        "data/work/glnexus/glnexus.norm.chrX.bcf"
    output:
        "data/work/glnexus/glnexus.norm.chrX.ref.bcf",
        "data/work/glnexus/glnexus.norm.chrX.ref.het.bcf",
        "data/work/glnexus/glnexus.norm.chrX.ref.het.hom.bcf",
        "data/work/glnexus/glnexus.norm.chrX.ref.het.hom.gt.bcf"
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

# Haploid genotype setting for chrY (non-PAR, males only)
rule set_chrY_gt:
    input:
        "data/work/glnexus/glnexus.norm.chrY.bcf"
    output:
        "data/work/glnexus/glnexus.norm.chrY.ref.bcf",
        "data/work/glnexus/glnexus.norm.chrY.ref.hom.bcf"
    shell:
        """
        # REF allele only → haploid 0
        bcftools +setGT {input} -Ob -o {output[0]} -- -t q -n c:'0' -i '(AD[:0] > 0 & AD[:1] == 0)'

        # ALT allele only → haploid 1
        bcftools +setGT {output[0]} -W=csi -Ob -o {output[1]} -- -t q -n c:'1' -i '(AD[:0] == 0 & AD[:1] > 0)'
        """

# Fill tags for chrX
rule fill_tags_chrX:
    input:
        "data/work/glnexus/glnexus.norm.chrX.ref.het.hom.gt.bcf"
    output:
        "data/work/glnexus/glnexus.norm.chrX.ref.het.hom.gt.tags.bcf"
    shell:
        "bcftools +fill-tags {input} -W=csi -Ob -o {output} -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE"

# Fill tags for chrY
rule fill_tags_chrY:
    input:
        "data/work/glnexus/glnexus.norm.chrY.ref.hom.bcf"
    output:
        "data/work/glnexus/glnexus.norm.chrY.ref.hom.tags.bcf"
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


# Stats for autosomal chromosomes
rule vcf_stats_autosomes:
    input:
        "data/work/glnexus/glnexus.norm.chr{chr}.ref.het.hom.gt.tags.bcf"
    output:
        "vcf_stats/glnexus.norm.chr{chr}.ref.het.hom.gt.tags.stats"
    shell:
        "bcftools stats -s - {input} > {output}"

# Stats for chrX
rule vcf_stats_chrX:
    input:
        "data/work/glnexus/glnexus.norm.chrX.ref.het.hom.gt.tags.bcf"
    output:
        "vcf_stats/glnexus.norm.chrX.ref.het.hom.gt.tags.stats"
    shell:
        "bcftools stats -s - {input} > {output}"

# Stats for chrY
rule vcf_stats_chrY:
    input:
        "data/work/glnexus/glnexus.norm.chrY.ref.hom.tags.bcf"
    output:
        "vcf_stats/glnexus.norm.chrY.ref.hom.tags.stats"
    shell:
        "bcftools stats -s - {input} > {output}"

