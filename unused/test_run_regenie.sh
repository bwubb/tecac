#!/usr/bin/env bash
set -euo pipefail

# Simple REGENIE test run script:
# - Step 1 once
# - Step 2 single-variant for chr1..22
# - Step 2 gene-based for chr1..22
# Uses freeze_qc_candidate annotation/set files and DM count 10.

DM_COUNT=10
DM_COVAR_FLAGS="DM{1:${DM_COUNT}}"

# Step 1
#regenie \
#  --step 1 \
#  --pgen data/preprocess/build \
#  --covarFile data/regenie/regenie.covar.txt \
#  --covarCol FREEZE \
#  --covarCol "${DM_COVAR_FLAGS}" \
#  --phenoFile data/regenie/regenie.pheno.txt \
#  --phenoCol STATUS \
#  --extract data/preprocess/build.snplist \
#  --bsize 1000 \
#  --gz \
#  --bt \
#  --lowmem \
#  --lowmem-prefix data/regenie/tmp_rg_ \
#  --out data/regenie/step1

# Step 2 (single variant) and Step 2 (gene-based) for chr1..22
for CHR in $(seq 1 22); do
  regenie \
    --step 2 \
    --pgen data/preprocess/chr${CHR}.annotation \
    --covarFile data/regenie/regenie.covar.txt \
    --covarCol FREEZE \
    --covarCol "${DM_COVAR_FLAGS}" \
    --phenoFile data/regenie/regenie.pheno.txt \
    --phenoCol STATUS \
    --bt \
    --firth \
    --approx \
    --pThresh 0.999 \
    --firth-se \
    --pred data/regenie/step1_pred.list \
    --bsize 400 \
    --af-cc \
    --minMAC 10 \
    --out data/regenie/chr${CHR}.step2_single_variant_STATUS

  regenie \
    --step 2 \
    --pgen data/preprocess/chr${CHR}.annotation \
    --phenoFile data/regenie/regenie.pheno.txt \
    --phenoCol STATUS \
    --covarFile data/regenie/regenie.covar.txt \
    --covarCol FREEZE \
    --covarCol "${DM_COVAR_FLAGS}" \
    --bt \
    --firth \
    --approx \
    --pThresh 0.999 \
    --firth-se \
    --pred data/regenie/step1_pred.list \
    --anno-file data/regenie/regenie.annotation.freeze_qc_candidate.txt \
    --set-list data/regenie/regenie.set.freeze_qc_candidate.txt \
    --mask-def data/regenie/regenie.mask.txt \
    --build-mask max \
    --write-mask-snplist \
    --aaf-bins 0.01,0.001,0.0001 \
    --strict-check-burden \
    --check-burden-files \
    --af-cc \
    --bsize 200 \
    --vc-tests skat,skato \
    --out data/regenie/chr${CHR}.step2_gene_based_STATUS.freeze_qc_candidate
  # REGENIE appends _${phenoCol} to --out (--phenoCol STATUS → ...freeze_qc_candidate_STATUS.regenie).
  # After copying to data/final/ with a project prefix: *.chrN.step2_gene_based_STATUS.freeze_qc_candidate_STATUS.regenie
done

