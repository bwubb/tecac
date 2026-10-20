# TECAC exome association

QC for a TECAC/PMBB exome callset, an optional MNP genotype repair, then REGENIE single-variant and gene-based tests.

Run Snakemake from this directory. There are two entry points:

```bash
snakemake -s qc_pipeline.smk <target> --configfile config.yml
snakemake -s regenie.smk <target> --configfile config.yml
```

Copy `config.example.yml` to `config.yml` and fill in paths. The example is the key list the live Snakefiles read. `config.yml` is local and is not part of the repo.

Files stay in this directory. The prefix is the block:

- `qc_` — sample QC and PCA helpers called from `build_files.smk`
- `mnp_` — the MNP branch called from `mnp.smk`
- `regenie_` — mask prep and report tables called from `regenie.smk`
- no prefix — not one block. `build_covariates.py` is prep. `vep_vcf_parser.py` is called by both annotation and MNP.

`rare_variant_qc.py` matches `rare_variant_qc.smk` and `rare_variant_qc.Rmd`. Snakemake rule names are unchanged, so `check_mnp_reads` is still the target you request. The shell line inside that rule calls `mnp_check_reads.py`.

## QC

`select_variants` is per-chromosome site QC and sex inference from chrX heterozygosity. It subsets each input BCF to `input.samples` and `input.targets`, splits multiallelics, and writes a site-QC BCF plus a PLINK pgen.

`build_files` is sample QC on the common SNPs: missingness (`qc_aggregate_missingness.py`), heterozygosity (`qc_het_outliers.py`), an LD-pruned SNP list, and PCA. The genotype set it writes is `data/preprocess/build` (het outliers removed; full post-QC SNP set kept for association). Neither target needs MNP calls or CRAMs.

Relatedness SNP filters (LD prune, PCA, and later REGENIE step 1) are separate from `build`. They drop chromosome 8, intervals in `long_ld_regions.bed`, and — when `qc.af_delta_thr` is set — IDs with `|AF_TECAC − AF_PMBB|` above that threshold (`data/qc/exclusions/af_delta_exclude.snplist`). Those SNPs stay in `build` and in the ExWAS annotation set. Extra long-LD windows were chosen from PCA allele-loading spikes (tight multi-SNP windows on `|PC|`, especially PC2/PC3), then added to the BED by hand. The AF table and exclude list come from `build_af_tecac_vs_pmbb` / `make_af_delta_exclude_snplist`. PCA output is `data/preprocess/build_pca.*` (no second pass that drops samples).

```bash
snakemake -s qc_pipeline.smk select_variants --configfile config.yml
snakemake -s qc_pipeline.smk build_files --configfile config.yml
```

`build_covariates.py` is not part of that run. Use it beforehand if you need the `FID IID FREEZE STATUS` file (`STATUS` 1 = control, 2 = case).

## MNP

MNP repair sits inside annotation (`annotation_files.smk` includes `mnp.smk`). `mnp.correct_genotypes` in `config.yml` turns the genotype rewrite on or off. Default is `true`.

The CRAM read-check built by Snakemake is **case carriers**. A pair is kept only when at least one case has both alleles and close VAFs. The PASS table still lists every cis-like carrier of that pair, cases and controls. `choose_mnp_read_check_samples` then writes a read-check target only for samples in `input.cases`. Controls in the PASS table are not sent to samtools. `mnp_carrier_stats.py` can write a control assignment list by hand; the Snakefile does not feed that list into `check_mnp_reads`.

With `mnp.per_sample: true` (the default), a genotype is rewritten only for a sample that was read-checked and called Strong_Cis. On the default path that is the case carriers who passed the CRAM check.

### Run the read-check and the rewrite

`mnp.correct_genotypes: true`

1. `mnp_before_read_check` — no CRAMs. Finds candidate pairs from the first VEP CSV, pulls het/hom-alt genotypes in those regions, and keeps pairs with at least one cis-like case.
   - `data/mnp/chr*.annotation.mnp.pairs.txt` — one candidate pair per line (`id1`, `id2`).
   - `data/mnp/chr*.mnp_sample_info.PASS.txt` — cis-like carriers, cases and controls, for pairs that cleared the case threshold.
   - `data/mnp/chr*.mnp_sample_info.FAIL.txt` — carriers who had only one allele, or both alleles with VAF too far apart, plus pairs dropped for having no cis-like case.
   - `data/mnp/mnp_read_check.assignments.tsv` — one row per **case** carrier to send to samtools.
   - `data/mnp/mnp_read_check.samples.txt` — the case IDs in that assignment file.
2. `mnp_carrier_stats.py` — optional, not a Snakemake target. Reads the PASS files and writes per-pair case and control counts. `--targets-out` is an `id1,id2` list; point `mnp.target_mnps` at it to limit which pairs become MNP records. `--assignments-out` is a control CRAM list. Using it means replacing or appending `mnp_read_check.assignments.tsv` yourself before the next step.
3. `check_mnp_reads` — samtools on every row in the assignment file. Needs `mnp.bam_table` and `resources/GRCh38_full_analysis_set_plus_decoy_hla.fa`.
   - `data/mnp/mnp_read_check.results.mapq20.csv` — one row per checked sample, including the Strong_Cis call.
   - `data/mnp/mnp_read_check.read_details.mapq20.tsv` — the reads behind that call.
   If that CSV is already on disk, later rules use it and do not call samtools again.
4. `mnp_check_and_blacklist` — keeps Strong_Cis rows, rewrites those genotypes into an MNP allele, VEPs the new records, and writes the M4 blacklist.
   - `data/mnp/chr*.mnp_sample_info.read_check_filtered.PASS.txt` — carriers whose genotypes will be rewritten.
   - `data/mnp/chr*.mnp_gt.plan.txt` — the rewrite plan (`mnp_id`, position, ref, alt, the two SNP ids, carriers).
   - `data/mnp/chr*.mnp.pass.gt.sorted.bcf` — rewritten genotypes.
   - `data/mnp/chr*.mnp.pass.no_sample.sorted.vep.report.csv` — VEP table for the new MNP records.
   - `data/mnp/mnp.validation.metrics.tsv` — counts from candidate pairs through the CRAM check.
   - `data/mnp/mnp_blacklist.txt` — SNP ids from the candidate pair list. REGENIE drops these from the M4 synonymous mask so a repaired pair is not also counted as two synonymous SNPs.

`annotation_files` then splices `mnp.pass.gt.sorted.bcf` into the passing-sample BCF (`adjust_mnp_gt`) and splices the MNP VEP records into the annotation VEP (`adjust_mnp_vep`).

### Run annotation, finalize, and REGENIE with no MNP rewrite

Set `mnp.correct_genotypes: false`. Do not request `mnp_before_read_check`, `check_mnp_reads`, or `mnp_check_and_blacklist`.

`annotation_files` copies the passing-sample BCF and the first VEP forward under the usual post-MNP filenames and does not call samtools. `finalize` still runs. It writes an empty `data/mnp/mnp_blacklist.txt`, so M4 is not stripped of MNP-pair SNP ids. `run_regenie` and `run_regenie_report` use that empty file.

```bash
snakemake -s qc_pipeline.smk annotation_files --configfile config.yml
snakemake -s qc_pipeline.smk finalize --configfile config.yml
snakemake -s regenie.smk run_regenie --configfile config.yml
snakemake -s regenie.smk run_regenie_report --configfile config.yml
```

The QC HTML report does not wait on `mnp.validation.metrics.tsv` when correction is off. The post-MNP variant-type counts are still written; they are counts of the uncorrected annotation BCF.

## Annotation and rare-variant QC

`annotation_files` subsets to samples that passed QC, runs VEP, and writes `data/preprocess/chr*.annotation`. With `mnp.correct_genotypes: true` it splices in the rewritten MNP genotypes. With `false` it copies the uncorrected BCF and VEP forward. `vep_vcf_parser.py` turns those VCFs into the CSVs REGENIE uses.

`rare_variant_qc` is controls-only Freeze-vs-Freeze missingness and carrier tests. It always writes an exclusion list. Those IDs are dropped from the annotation BCF only when `rare_variant_qc.apply_exclusions` is true.

`qc_report` renders `exwas_qc_report.Rmd`. `finalize` copies the annotation pgens, the `build` pgen, the PCA, and the case/control lists into `data/final/`.

```bash
snakemake -s qc_pipeline.smk annotation_files --configfile config.yml
snakemake -s qc_pipeline.smk qc_report --configfile config.yml
snakemake -s qc_pipeline.smk finalize --configfile config.yml
```

## REGENIE

Second entry point, after `finalize`. Step 1 starts from `build`, drops `long_ld_regions.bed`, then drops any AF-delta IDs. Step 2 is single-variant and gene-based (masks). `mask_gene_consequence_exclude.txt` removes CHEK2 missense from every mask.

```bash
snakemake -s regenie.smk run_regenie --configfile config.yml
snakemake -s regenie.smk run_regenie_report --configfile config.yml
```

`run_regenie_report` renders `report_regenie.Rmd` from tables written by `regenie_build_report_tables.py`.

`dm_count: 0` keeps FREEZE as a covariate and does not add ancestry PCs.

## Reports that are not targets

- `exwas_qc_analysis_plan.Rmd` — the written QC plan. Knit it by hand.
- `rare_variant_qc.Rmd` — rendered by the `rare_variant_qc` target.
- `mnp_validation_report.Rmd` — reads `data/mnp/mnp_carrier_stats.tsv` from `mnp_carrier_stats.py`.

## Long-range LD and AF-delta (PCA / GRM only)

`long_ld_regions.bed` is the hg38 list used for LD prune, PCA, and REGENIE step 1 (plink2 `--exclude range`, 0-based half-open). It started as HLA / INV17q21 and grew after PC allele-loading diagnostics in `exwas_qc_report.Rmd` (spike-body windows at `|PC| ≥ 2`, plus the 99.9th-percentile tables). Those are region filters for relatedness, not a variant blacklist on `build` or annotation.

`qc.af_delta_thr` (e.g. `0.05`) flags common SNPs whose allele frequency differs between TECAC (FREEZE 5) and PMBB (FREEZE 1–4). IDs already inside the BED are not an extra drop; the leftover IDs are the second plink2 `--exclude` after the BED. Omit the key or set it `≤ 0` for regions-only. `input.long_ld_bed` overrides the default BED path.

## What a clone does not include

Genotypes and working files under `data/`, your `config.yml`, the GRCh38 FASTA at `resources/GRCh38_full_analysis_set_plus_decoy_hla.fa`, CRAM/BAM paths in `mnp_bam.table`, and the cluster profile. The shell rules log under `logs/lsf`.

Tools those rules call: bcftools, plink2, VEP, REGENIE, samtools (MNP read check), Python (pandas; vcfpy for the VEP parser and MNP genotype rewrite), R (rmarkdown, tidyverse, yaml).
