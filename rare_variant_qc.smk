CHROMOSOMES_AUTOSOMAL=list(range(1,23))

MISSING_ABS_DIFF_THR=0.03
HIGH_MISSING_THR=0.05
MISSING_FISHER_P_THR=1e-6
CARRIER_FISHER_P_THR=1e-4



rule rare_variant_qc:
    input:
        "data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.html",
        "data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.all.tsv",
        "data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.flagged.tsv"

rule rare_variant_qc_extract_per_chr:
    input:
        bcf="data/bcftools/chr{CHR}.qc_filter1.het_miss.mnp.gt.bcf",
        samples="data/qc/passing_samples.txt",
        covariates=config.get("input",{}).get("covariates","")
    output:
        matrix="data/qc/rare_variant_qc/chr{CHR}.controls.freeze_gt.matrix.tsv",
        sample_order="data/qc/rare_variant_qc/chr{CHR}.controls.sample_order.txt",
        sample_freeze="data/qc/rare_variant_qc/chr{CHR}.controls.sample_freeze.tsv"
    shell:
        """
        mkdir -p data/qc/rare_variant_qc
        # covariates columns are fixed: FID IID FREEZE STATUS; controls are STATUS==1.
        awk 'BEGIN{{FS="[ \t]+";OFS="\\t"}} NR==FNR{{if(FNR>1 || $1!="IID")pass[$1]=1;next}} FNR==1{{next}} {{gsub(/\\r/,"",$2);gsub(/\\r/,"",$3);gsub(/\\r/,"",$4);if(($2 in pass) && ($4=="1") && (($3=="2") || ($3=="3")))print $2,$3}}' {input.samples} {input.covariates} > {output.sample_freeze}
        test -s {output.sample_freeze}
        echo "Controls written: $(wc -l < {output.sample_freeze})"
        echo "Freeze 2 controls: $(awk 'BEGIN{{FS=\"\\t\"}} $2==\"2\"{{c++}} END{{print c+0}}' {output.sample_freeze})"
        echo "Freeze 3 controls: $(awk 'BEGIN{{FS=\"\\t\"}} $2==\"3\"{{c++}} END{{print c+0}}' {output.sample_freeze})"
        awk 'BEGIN{{FS="\\t"}} $2=="2"{{c2++}} $2=="3"{{c3++}} END{{if((c2+0)==0 || (c3+0)==0){{print "ERROR: freeze split invalid in sample_freeze.tsv (freeze2=" (c2+0) ", freeze3=" (c3+0) ")" > "/dev/stderr"; exit 1}}}}' {output.sample_freeze}
        awk 'BEGIN{{FS="\\t"}} {{print $1}}' {output.sample_freeze} > {output.sample_order}.tmp
        test -s {output.sample_order}.tmp
        bcftools query -l -S {output.sample_order}.tmp {input.bcf} > {output.sample_order}
        test -s {output.sample_order}
        bcftools query -S {output.sample_order} -f '%ID\\t%TYPE[\\t%GT]\\n' {input.bcf} > {output.matrix}
        rm -f {output.sample_order}.tmp
        """

rule rare_variant_qc_per_chr:
    input:
        matrix="data/qc/rare_variant_qc/chr{CHR}.controls.freeze_gt.matrix.tsv",
        sample_order="data/qc/rare_variant_qc/chr{CHR}.controls.sample_order.txt",
        sample_freeze="data/qc/rare_variant_qc/chr{CHR}.controls.sample_freeze.tsv"
    output:
        "data/qc/rare_variant_qc/chr{CHR}.rare_variant_qc.tsv"
    shell:
        """
        python rare_variant_qc.py \
          --matrix {input.matrix} \
          --sample-order {input.sample_order} \
          --sample-freeze {input.sample_freeze} \
          --min-abs-missing-diff {MISSING_ABS_DIFF_THR} \
          --high-missing-threshold {HIGH_MISSING_THR} \
          --missing-fisher-p-threshold {MISSING_FISHER_P_THR} \
          --carrier-fisher-p-threshold {CARRIER_FISHER_P_THR} \
          --output {output}
        """

rule rare_variant_qc_report:
    input:
        expand("data/qc/rare_variant_qc/chr{CHR}.rare_variant_qc.tsv",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        html="data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.html",
        all="data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.all.tsv",
        flagged="data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.flagged.tsv"
    shell:
        """
        Rscript -e "rmarkdown::render('rare_variant_qc.Rmd',output_file='{output.html}')"
        """


