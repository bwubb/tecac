#!/bin/bash
#This makes input for denisty plots

i=$1
o=$2
chr=$3

bcftools view -f PASS --type snps,indels,mnps $i $chr | bcftools query -i 'GT="het"' -f '%CHROM:%POS:%ALT [%VAF %SAMPLE]\n' > $o
