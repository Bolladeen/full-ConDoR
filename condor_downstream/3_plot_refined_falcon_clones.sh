mamba activate condor

SAMPLE_NAMES=$(ls /juno/work/iacobuzc/haochen/Tapestri_batch2/analysis/condor_downstream/HZ_ete_trees_refined_subclonal_snvs)

TAP_CN_CALLING_DIR=/home/zhangh5/work/Tapestri_project/cn-calling
PLOT_SCRIPT=${TAP_CN_CALLING_DIR}/scripts/PLOT-unique_cn_profiles.py \
OUTPUT_DIR=/home/zhangh5/work/Tapestri_batch2/analysis/condor_downstream/HZ_ete_trees_refined_subclonal_snvs
for patient_i in ${SAMPLE_NAMES[@]}; do
    echo "Processing ----- ${patient_i}"

    # mkdir -p ${OUTPUT_DIR}/${patient_i}/before_refine
    # python /Users/haochen/Desktop/Tapestri_analysis/copy_number/tap_cn_calling/scripts/PLOT-unique_cn_profiles.py \
    #     --cohort_name ${patient_i} \
    #     --amp_gene_map_f /Users/haochen/Desktop/Tapestri_analysis/copy_number/tap_cn_calling/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt \
    #     --cn_clone_profiles_csv /Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/all_cases_plot_info_2/${patient_i}/original_clone_profiles.csv \
    #     --sample_sc_clone_assignment_csv /Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/all_cases_plot_info_2/${patient_i}/original_cell_assignments.csv \
    #     --output_dir ${OUTPUT_DIR}/${patient_i}/before_refine

    mkdir -p ${OUTPUT_DIR}/${patient_i}/refined_cn_info
    python ${PLOT_SCRIPT} \
        --cohort_name ${patient_i} \
        --amp_gene_map_f ${TAP_CN_CALLING_DIR}/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt \
        --cn_clone_profiles_csv ${OUTPUT_DIR}/${patient_i}/${patient_i}_final_clone_cn_profiles.csv \
        --sample_sc_clone_assignment_csv ${OUTPUT_DIR}/${patient_i}/${patient_i}_final_sc_clone_assignment.csv \
        --output_dir ${OUTPUT_DIR}/${patient_i}/refined_cn_info \
        --omit_clone_cleanup
        # --plot_homdel True \
        # --homdel_genes_oi CDKN2A \
        # --cn_call_yaml /Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/pipeline_results_custom/falcon_results/RA20_05/RA20_05-3samples-falcon-AJ.yaml
done

