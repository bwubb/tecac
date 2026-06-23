#needs rare_variant_qc.py
#needs rare_variant_qc.Rmd

import datetime

CHROMOSOMES_AUTOSOMAL=list(range(1,23))
DATE=datetime.date.today().strftime('%Y%m%d')

MISSING_ABS_DIFF_THR=0.03
HIGH_MISSING_THR=0.05
MISSING_FISHER_P_THR=1e-6
CARRIER_FISHER_P_THR=1e-4


rule rare_variant_qc:
    input:
        f"data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.{DATE}.html",
        "data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.all.tsv",
        "data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.flagged.tsv",
        "data/qc/exclusions/rare_variant_exclusions.txt",

rule rare_variant_qc_extract_per_chr:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.het_miss.mnp.gt.bcf",
        samples="data/qc/passing_samples.txt",
    params:
        covariates=config.get("input",{}).get("covariates","covariates.txt")
    output:
        matrix="data/qc/rare_variant/chr{CHR}.controls.freeze_gt.matrix.tsv",
        sample_order="data/qc/rare_variant/chr{CHR}.controls.sample_order.txt",
        sample_freeze="data/qc/rare_variant/chr{CHR}.controls.sample_freeze.tsv"
    shell:
        """
        mkdir -p data/qc/rare_variant
        # covariates: FID IID FREEZE STATUS; controls STATUS==1.passing_samples: IID per line OR FID IID;register both FID and IID for membership (match covariates IID col2).
        awk 'BEGIN{{FS="[ \t]+";OFS="\\t"}} NR==FNR{{if(NR==1 && ((NF>=2 && $1 ~ /^#?FID$/ && $2 ~ /^#?IID$/) || (NF==1 && $1 ~ /^#?IID$/)))next;if(NF>=2){{pass[$1]=1;pass[$2]=1}}else if($1!="")pass[$1]=1;next}} FNR==1{{next}} {{gsub(/\\r/,"",$2);gsub(/\\r/,"",$3);gsub(/\\r/,"",$4);if(($2 in pass) && ($4=="1") && (($3=="2") || ($3=="3")))print $2,$3}}' {input.samples} {params.covariates} > {output.sample_freeze}
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
        matrix="data/qc/rare_variant/chr{CHR}.controls.freeze_gt.matrix.tsv",
        sample_order="data/qc/rare_variant/chr{CHR}.controls.sample_order.txt",
        sample_freeze="data/qc/rare_variant/chr{CHR}.controls.sample_freeze.tsv"
    output:
        tsv="data/qc/rare_variant/chr{CHR}.rv-qc.tsv",
        flagged_ids="data/qc/rare_variant/chr{CHR}.rv-qc.flagged.ids",
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
          --output {output.tsv} \
          --exclude-ids {output.flagged_ids}
        """

rule rare_variant_qc_report:
    input:
        expand("data/qc/rare_variant/chr{CHR}.rv-qc.tsv",CHR=CHROMOSOMES_AUTOSOMAL),
        pathogenic="data/preprocess/pathogenic_vus.csv",
    output:
        html=f"data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.{DATE}.html",
        all="data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.all.tsv",
        flagged="data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.flagged.tsv",
    shell:
        """
        Rscript -e "rmarkdown::render('rare_variant_qc.Rmd',output_file='{output.html}')"
        """

rule rare_variant_qc_exclusions:
    input:
        expand("data/qc/rare_variant/chr{CHR}.rv-qc.flagged.ids",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/qc/exclusions/rare_variant_exclusions.txt"
    shell:
        """
        mkdir -p data/qc/exclusions
        cat {input} | awk 'NF>0' | sort -u > {output}
        """
