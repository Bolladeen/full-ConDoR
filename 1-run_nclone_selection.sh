mamba activate mosaic
FALCON_REPO_DIR=/data/iacobuzc/haochen/Tapestri_project/tap_cn_calling
nclone_selection_script=${FALCON_REPO_DIR}/scripts/select_optimal_nclones.py
falcon_plot_script=${FALCON_REPO_DIR}/scripts/PLOT-unique_cn_profiles.py

CN_CLONES_DIR=/data/iacobuzc/haochen/Tapestri_batch2/analysis/cn_clones
OUTPUT_DIR=/data/iacobuzc/haochen/Tapestri_batch2/analysis/cn_clones/_falcon_debug
SAMPLE_NAMES=(BPA-5-RSX)

for sample_i in ${SAMPLE_NAMES[@]}; do
    # selected_cn_profiles=$(
    #     find ${CN_CLONES_DIR}/${sample_i}/*/*/selected_solution/ -name "*clone_profiles.csv" \
    #     -print -quit
    # )
    # selected_assignments=$(
    #     find ${CN_CLONES_DIR}/${sample_i}/*/*/selected_solution/ -name "*assignments.csv" \
    #     -print -quit
    # )
    python ${nclone_selection_script} \
        -d ${sample_i} \
        -i /data/iacobuzc/haochen/Tapestri_batch2/analysis/cn_clones/${sample_i}/${sample_i}-cohort-cn_calling-falcon_homdel_from_scratch_AJ/cn_call_with_homdel/cleaned_solutions \
        --output_dir ${OUTPUT_DIR}
        # && \
    # python ${falcon_plot_script} \
    #     --cohort_name ${sample_i} \
    #     --amp_gene_map_f ${FALCON_REPO_DIR}/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt \
    #     --cn_clone_profiles_csv ${selected_cn_profiles} \
    #     --sample_sc_clone_assignment_csv ${selected_assignments} \
    #     --output_dir ${OUTPUT_DIR}

done
