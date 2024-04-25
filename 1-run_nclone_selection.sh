mamba activate mosaic
FALCON_REPO_DIR=/data/iacobuzc/haochen/Tapestri_project/tap_cn_calling
nclone_selection_script=${FALCON_REPO_DIR}/scripts/select_optimal_nclones.py
falcon_plot_script=${FALCON_REPO_DIR}/scripts/PLOT-unique_cn_profiles.py

CN_CLONES_DIR=/data/iacobuzc/haochen/Tapestri_batch2/analysis/cn_clones
OUTPUT_DIR=/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/falcon_solutions
SAMPLE_NAMES=(BPA-1 BPA-2-IR BPA-3)

for sample_i in ${SAMPLE_NAMES[@]}; do
    selected_cn_profiles=$(
        find ${CN_CLONES_DIR}/${sample_i}/*/*/selected_solution/ -name "*clone_profiles.csv" \
        -print -quit
    )
    selected_assignments=$(
        find ${CN_CLONES_DIR}/${sample_i}/*/*/selected_solution/ -name "*assignments.csv" \
        -print -quit
    )
    # python ${nclone_selection_script} \
    #     -d ${sample_i} \
    #     -i /home/zhangh5/work/Tapestri_batch2/analysis/cn_clones/RA19_21/RA19_21-cohort-cn_calling-falcon_homdel_from_scratch_AJ/cn_call_with_homdel/solutions \
    #     --output_dir ${OUTPUT_DIR} \
    #     && \
    python ${falcon_plot_script} \
        --cohort_name ${sample_i} \
        --amp_gene_map_f ${FALCON_REPO_DIR}/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt \
        --cn_clone_profiles_csv ${selected_cn_profiles} \
        --sample_sc_clone_assignment_csv ${selected_assignments} \
        --output_dir ${OUTPUT_DIR}

done
