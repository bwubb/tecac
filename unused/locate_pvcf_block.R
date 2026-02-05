#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)


args <- commandArgs(trailingOnly = TRUE)

# First argument is the path to the small_genes.bed file
genes_file <- args[1]

genes <- read.delim(genes_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(genes) <- c("chr", "start", "end", "gene_name")

blocks <- read.delim("pvcf_blocks.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(blocks) <- c("line_number", "chr", "block_number", "block_start", "block_end")

genes$chr <- as.character(genes$chr)
blocks$chr <- as.character(blocks$chr)

# Cross join genes with blocks and then filter to find overlaps
gene_blocks <- suppressWarnings(
  genes %>%
    inner_join(blocks, by = "chr") %>%
    filter(block_start <= end, block_end >= start) %>%
    mutate(filename = paste0("ukb23157_c", chr, "_b", block_number, "_v1.vcf.gz")) %>%
    select(chr, gene_name, formatted_block) %>%
    distinct()
)

print(gene_blocks)

#write.csv(gene_blocks, "genes_with_file_blocks.csv", row.names = FALSE)
