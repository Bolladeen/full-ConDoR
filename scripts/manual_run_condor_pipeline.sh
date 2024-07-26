# # ===== 1a.generate condor input =====
# mamba activate mosaic-custom

# script=/home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/scripts/generate_condor_input.py

# python ${script} \
#     -d M04 \
#     -l /home/zhangh5/work/Tapestri_batch2/batch2_data_compiled/SNV_FILLOUT_results/ \
#     -i /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/opt_nclones/M04-homdel-nclones=9_solution.cell_assignments.csv \
#     -snvs /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/M04-patient/sc_heatmaps/all_vars/mut_prev=0.01/all_vars-voi.hz.txt \
#     -v /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/character_vaf_matrix.csv \
#     -m /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/character_bin_matrix.csv \
#     -a /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/alt_readcounts.csv \
#     -t /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/total_readcounts.csv \
#     -g /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/germline_mutations.txt \
#     -s /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/somatic_mutations.txt


# patients=("M04" "M07" "M11" "M12" "M13" "RA15_06" "RA15_16" "RA16_08" "RA16_29" "RA17_13" "RA17_22" "RA19_02" "RA19_10" "RA19_21" "RA20_05" "RA21_17" "TP6" "TP12")


# ===== 1b. generate pre-condor heatmaps =====
conda activate mosaic
SC_HEATMAP_SCRIPT=/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/full-ConDoR/scripts/1b_generate_sc_heatmaps.py

python ${SC_HEATMAP_SCRIPT} \
    --patient_name BPA-4 \
    --patient_h5 /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled/fillout_h5/BPA-4.patient_wide.genotyped.h5 \
    --snv_f /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled/manual_annotated_snv_lists/BPA-4-all_vars-voi.hz_curated.txt\
    --clone_assignment /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled/falcon_solutions/BPA-4-RSX_homdel_nclones=9.sample_sc_clone_assignment.updated.csv \
    --output_dir /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_pipeline_new


# # ===== generate dot and png plots =====
# mkdir "/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/dot_files"

# # Loop through each dataset
# for dataset in "${datasets[@]}"; do
#   source="/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/condor_outputs/${dataset}/out_tree.dot"
#   # Execute scp command
#   dot -Tpng "$source" -o "/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/dot_files/${dataset}.out_tree.png"

#   # Check scp exit status
#   if [ $? -eq 0 ]; then
#   echo "Copied dot file for $dataset successfully."
#   else
#   echo "Failed to copy dot file for $dataset."
#   fi
# done

# # ===== generate CNV plots =====
# dest_dir="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/het_plots/"
# mkdir "$dest_dir"

# for dataset in "${datasets[@]}"; do
#   python "../../../tap_cn_calling/scripts/PLOT-unique_cn_profiles.py" \
#   --cohort_name ${dataset} \
#   --amp_gene_map_f "../../../tap_cn_calling/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt" \
#   --cn_clone_profiles_csv "../results/${dataset}/opt_nclones/optimal_clone_profiles.csv" \
#   --sample_sc_clone_assignment_csv "../results/${dataset}/opt_nclones/optimal_cell_assignments.csv" \
#   --output_dir "../results/het_plots/" \
#   --clone_cleanup False


# done

# echo "Script execution complete."

# ===== 4a. condor trees --> ETE trees =====
conda activate condor
ETE_TREE_SCRIPT=/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/full-ConDoR/condor_downstream/4a_make_ete_tree_with_subclonal_snvs.py
python ${ETE_TREE_SCRIPT} \
    --amp_gene_map /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/full-ConDoR/references/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.gene_pos_ado_added.csv \
    --sample_name_map /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/Tapestri_batch2_samples_MASTER_INTERNAL.xlsx \
    --patient_name M04 \
    --condor_tree_pickle /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_pipeline_OLD/condor_outputs/pickle_files/M04_self.solT_cell \
    --snv_ann_f /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_main_manuscript_analysis/data_compiled/manual_annotated_snv_lists/M04-patient-all_vars-voi.hz_curated.txt \
    --output_dir /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_pipeline_new/condor_downstream