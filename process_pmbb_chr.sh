#!/bin/bash

# Emergency script to process PMBB BCF files
# Usage: bash process_pmbb_chr.sh $i (where i is chromosome number)

set -e  # Exit on any error

CHR=$1
if [ -z "$CHR" ]; then
    echo "Usage: bash process_pmbb_chr.sh <chromosome_number>"
    echo "Example: bash process_pmbb_chr.sh 9"
    exit 1
fi

# Input file pattern
INPUT_BCF="PMBB-Release-2024-3.0_genetic_exome_chr${CHR}_NF.bcf"

# Check if input file exists
if [ ! -f "$INPUT_BCF" ]; then
    echo "Error: Input file $INPUT_BCF not found"
    exit 1
fi

echo "Processing chromosome $CHR: $INPUT_BCF"

# Create output directory
mkdir -p data/work/pmbb_chr${CHR}

# Step 1: Normalize
echo "Step 1: Normalizing $INPUT_BCF"
bcftools norm -m-both "$INPUT_BCF" | \
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -W=csi -Ob -o "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.bcf"

# Step 2: Fill tags
echo "Step 2: Filling tags"
bcftools +fill-tags "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.bcf" -W=csi -Ob -o "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.bcf" -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE

# Step 3: Variant QC
echo "Step 3: Variant QC filtering"
bcftools filter -s ExcessHet -e 'INFO/ExcHet=0' -m + "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.bcf" | \
bcftools filter -s LowCallRate -e 'INFO/F_MISSING > 0.05' -m + | \
bcftools filter -s LowSampleCount -e 'INFO/NS < 1000' -m + | \
bcftools filter -s NoHQHet -e 'COUNT(FORMAT/GT="0/1" && FORMAT/DP>=10 && FORMAT/GQ>=20 && FORMAT/VAF>0.2)=0' -m + -W=csi -Ob -o "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.qc.bcf"

# Step 4: Remove sample data (no_sample)
echo "Step 4: Removing sample data"
bcftools view -G -Ov -o "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.qc.no_sample.vcf" "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.qc.bcf"

# Step 5: Annotate variants with VEP
echo "Step 5: VEP annotation"
singularity run -H $PWD:/home \
--bind /home/bwubb/resources:/opt/vep/resources \
--bind /home/bwubb/.vep:/opt/vep/.vep \
/appl/containers/vep112.sif vep \
--dir /opt/vep/.vep \
-i "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.qc.no_sample.vcf" \
-o "data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.qc.no_sample.vep.vcf" \
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
--custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar_20250831.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
--plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
--plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz

echo "Completed processing chromosome $CHR"
echo "Final output: data/work/pmbb_chr${CHR}/pmbb_chr${CHR}.norm.tags.qc.no_sample.vep.vcf"
