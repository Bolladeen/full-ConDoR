mamba activate mosaic-custom
falcon_repo_dir=/home/zhangh5/work/Tapestri_project/cn-calling
nclone_selection_script=${falcon_repo_dir}/scripts/select_optimal_nclones.py
falcon_plot_script=${falcon_repo_dir}/scripts/PLOT-unique_cn_profiles.py

SAMPLE_NAMES=(RA17_13 RA20_05 RA17_22 TP12 M07 RA15_16 M12 RA16_08 RA21_17 M13 RA16_29 M04 RA15_06 RA19_02 M11 RA19_10 TP6)
SAMPLE_NAMES=(RA15_06 RA16_08 RA17_13 RA17_22)
SAMPLE_NAMES=(M13)
# RA21_17 problematic
SAMPLE_NAMES=(RA16_29)
SAMPLE_NAMES=(RA19_21)

OUTPUT_DIR=/home/zhangh5/work/Tapestri_batch2/analysis/condor_wd/opt_nclones

for sample_i in ${SAMPLE_NAMES[@]}; do

    python ${nclone_selection_script} \
        -d ${sample_i} \
        -i /home/zhangh5/work/Tapestri_batch2/analysis/cn_clones/RA19_21/RA19_21-cohort-cn_calling-falcon_homdel_from_scratch_AJ/cn_call_with_homdel/solutions \
        --output_dir ${OUTPUT_DIR} \
        && \
    python ${falcon_plot_script} \
        --cohort_name ${sample_i} \
        --amp_gene_map_f ${falcon_repo_dir}/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt \
        --cn_clone_profiles_csv $(ls ${OUTPUT_DIR}/${sample_i}*clone_profiles.csv) \
        --sample_sc_clone_assignment_csv $(ls ${OUTPUT_DIR}/${sample_i}*assignments.csv) \
        --omit_clone_cleanup \
        --output_dir ${OUTPUT_DIR}

done
