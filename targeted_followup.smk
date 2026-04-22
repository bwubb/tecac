rule targeted_followup_all:
    input:
        "data/qc/targeted_followup/freeze_case_control_summary.tsv",
        "data/qc/targeted_followup/freeze_qc_m4_overlap.tsv",
        "data/regenie/regenie.annotation.freeze_qc_candidate.txt"

rule targeted_followup_variant_ids:
    input:
        pathogenic="data/regenie/pathogenic_vus.csv",
        flagged="data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.flagged.tsv"
    output:
        "data/qc/targeted_followup/target_variant_ids.txt"
    params:
        covariates=config.get("input",{}).get("covariates","")
    shell:
        """
        mkdir -p data/qc/targeted_followup
        python targeted_freeze_followup.py \
          --covariates {params.covariates} \
          --passing-samples data/qc/passing_samples.txt \
          --pathogenic-vus {input.pathogenic} \
          --genes CHEK2,HEATR3 \
          --flagged-tsv {input.flagged} \
          --top-flagged 50 \
          --output-ids {output} \
          --ids-only
        """

rule targeted_followup_summary:
    input:
        covariates=config.get("input",{}).get("covariates",""),
        passing="data/qc/passing_samples.txt",
        ids="data/qc/targeted_followup/target_variant_ids.txt"
    output:
        "data/qc/targeted_followup/freeze_case_control_summary.tsv"
    shell:
        """
        python targeted_freeze_followup.py \
          --covariates {input.covariates} \
          --passing-samples {input.passing} \
          --variant-list {input.ids} \
          --output {output}
        """

rule freeze_qc_m4_overlap:
    input:
        flagged="data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.flagged.tsv",
        pathogenic="data/regenie/pathogenic_vus.csv"
    output:
        "data/qc/targeted_followup/freeze_qc_m4_overlap.tsv"
    shell:
        """
        mkdir -p data/qc/targeted_followup
        python freeze_qc_m4_overlap.py \
          --flagged {input.flagged} \
          --pathogenic-vus {input.pathogenic} \
          --output {output}
        """

rule regenie_annotation_freeze_qc_candidate:
    input:
        annotation="data/regenie/regenie.annotation.txt",
        flagged="data/qc/rare_variant_qc/rare_variant_qc.freeze2_vs_freeze3.flagged.tsv"
    output:
        "data/regenie/regenie.annotation.freeze_qc_candidate.txt"
    shell:
        """
        awk -F'\\t' 'NR>1{{print $1}}' {input.flagged} > data/qc/targeted_followup/_exclude_ids.txt
        python filter_regenie_annotation.py \
          --annotation {input.annotation} \
          --exclude-ids data/qc/targeted_followup/_exclude_ids.txt \
          --output {output}
        """
