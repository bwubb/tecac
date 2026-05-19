#When building target_variant_ids: all IDs for these Gene values in pathogenic_vus.csv,plus top freeze-QC flagged.Use FOCUS_GENES=NONE for no gene pull (flagged-only).
FOCUS_GENES="CHEK2,HEATR3"
TOP_FLAGGED=50

rule targeted_followup_all:
    input:
        "data/qc/targeted_followup/freeze_case_control_summary.tsv",
        "data/qc/targeted_followup/freeze_case_control_summary_carrier_fisher.tsv",
        "data/qc/targeted_followup/freeze_qc_m4_overlap.tsv",
        "data/regenie/regenie.annotation.freeze_qc_candidate.txt",
        "data/regenie/regenie.set.freeze_qc_candidate.txt"

rule targeted_followup_variant_ids:
    input:
        pathogenic="data/regenie/pathogenic_vus.csv",
        flagged="data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.flagged.tsv"
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
          --genes {FOCUS_GENES} \
          --flagged-tsv {input.flagged} \
          --top-flagged {TOP_FLAGGED} \
          --output-ids {output} \
          --ids-only
        """

rule targeted_followup_summary:
    input:
        covariates=config.get("input",{}).get("covariates",""),
        passing="data/qc/passing_samples.txt",
        ids="data/qc/targeted_followup/target_variant_ids.txt"
    output:
        long="data/qc/targeted_followup/freeze_case_control_summary.tsv",
        fisher="data/qc/targeted_followup/freeze_case_control_summary_carrier_fisher.tsv"
    shell:
        """
        python targeted_freeze_followup.py \
          --covariates {input.covariates} \
          --passing-samples {input.passing} \
          --variant-list {input.ids} \
          --output {output.long} \
          --carrier-fisher-out {output.fisher}
        """

rule freeze_qc_m4_overlap:
    input:
        flagged="data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.flagged.tsv",
        pathogenic="data/regenie/pathogenic_vus.csv",
        synonymous="data/regenie/synonymous.csv"
    output:
        "data/qc/targeted_followup/freeze_qc_m4_overlap.tsv"
    shell:
        """
        mkdir -p data/qc/targeted_followup
        python targeted_freeze_qc_m4_overlap.py \
          --flagged {input.flagged} \
          --pathogenic-vus {input.pathogenic} \
          --synonymous {input.synonymous} \
          --output {output}
        """

rule regenie_freeze_qc_candidate_files:
    input:
        annotation="data/regenie/regenie.annotation.txt",
        set_file="data/regenie/regenie.set.txt",
        flagged="data/qc/rare_variant/rv-qc.freeze2_vs_freeze3.flagged.tsv"
    output:
        annotation="data/regenie/regenie.annotation.freeze_qc_candidate.txt",
        set_file="data/regenie/regenie.set.freeze_qc_candidate.txt"
    shell:
        """
        python filter_regenie_annotation.py \
          --annotation {input.annotation} \
          --annotation-out {output.annotation} \
          --set-file {input.set_file} \
          --set-out {output.set_file} \
          --flagged-tsv {input.flagged}
        """
