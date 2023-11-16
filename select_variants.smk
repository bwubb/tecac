







rule



#We will deal with the biallelic sites first.
#Need to check annotation code for complex.

#bcftools select biallelic sites only
#bcftools view --max-alleles 2 --exclude-types indels input.vcf.gz
#A typical command to filter out anything but biallelic SNPs, as stated in the bcftools manual, is the following:
#bcftools view -m2 -M2 -v snps input.vcf.gz
