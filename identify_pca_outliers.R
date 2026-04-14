#!/usr/bin/env Rscript
# Remove samples with PC4 > 0.3 from PLINK --pca eigenvec; write FID IID for plink --remove.

PC_COL <- "PC4"
PC_THRESHOLD <- 0.3

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(NA_character_)
  args[idx + 1]
}

eigenvec_path <- get_arg("--eigenvec")
outliers_path <- get_arg("--outliers")
plot_path <- get_arg("--plot")
report_path <- get_arg("--report")

if (!nzchar(eigenvec_path) || !nzchar(outliers_path) || !nzchar(plot_path) || !nzchar(report_path)) {
  stop("Usage: Rscript identify_pca_outliers.R --eigenvec <eigenvec> --outliers <out.txt> --plot <png> --report <pca_outlier_report.txt>")
}

d <- read.table(
  eigenvec_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)

if (!all(c("PC1", "PC2", PC_COL) %in% names(d))) {
  stop("Input eigenvec must contain PC1, PC2, and ", PC_COL, ": ", eigenvec_path)
}

fid_col <- if ("#FID" %in% names(d)) "#FID" else names(d)[1]
iid_col <- if ("IID" %in% names(d)) "IID" else if (length(names(d)) >= 2) names(d)[2] else names(d)[1]

d$PC1 <- suppressWarnings(as.numeric(d$PC1))
d$PC2 <- suppressWarnings(as.numeric(d$PC2))
pcv <- suppressWarnings(as.numeric(d[[PC_COL]]))

out_idx <- which(is.finite(pcv) & pcv > PC_THRESHOLD)
out <- d[out_idx, , drop = FALSE]

n_in <- nrow(d)
n_out <- length(out_idx)
n_after <- n_in - n_out
rule_str <- paste0(PC_COL, ">", PC_THRESHOLD)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(
  data.frame(
    n_in_pca = n_in,
    n_pca_outliers = n_out,
    rule = rule_str,
    n_after_pca = n_after,
    stringsAsFactors = FALSE
  ),
  file = report_path,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

dir.create(dirname(outliers_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)

if (nrow(out) > 0) {
  write.table(
    out[, c(fid_col, iid_col), drop = FALSE],
    file = outliers_path,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
} else {
  file.create(outliers_path)
}

rule_lab <- paste0(PC_COL, " > ", PC_THRESHOLD)

png(plot_path, width = 1200, height = 900, res = 150)
plot(
  d$PC1, d$PC2,
  pch = 16, cex = 0.55, col = rgb(0, 0, 0, 0.35),
  xlab = "PC1", ylab = "PC2",
  main = paste0("PCA outliers (", rule_lab, ")")
)
if (nrow(out) > 0) {
  points(out$PC1, out$PC2, pch = 16, cex = 0.9, col = "red")
}
dev.off()

plot_pc <- dirname(plot_path)
plot_base <- tools::file_path_sans_ext(basename(plot_path))
plot_pcgt <- file.path(plot_pc, paste0(plot_base, "_", PC_COL, ".png"))
png(plot_pcgt, width = 1200, height = 900, res = 150)
plot(
  d$PC1, pcv,
  pch = 16, cex = 0.55, col = rgb(0, 0, 0, 0.35),
  xlab = "PC1", ylab = PC_COL,
  main = paste0("PCA: ", rule_lab)
)
abline(h = PC_THRESHOLD, lty = 2, col = "blue")
if (length(out_idx) > 0) {
  points(d$PC1[out_idx], pcv[out_idx], pch = 16, cex = 0.9, col = "red")
}
dev.off()
