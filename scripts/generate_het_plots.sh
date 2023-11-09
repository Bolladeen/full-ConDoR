#!/bin/bash
datasets=("M04" "M07" "M11" "M12" "M13" "RA15_06" "RA15_16" "RA16_08" "RA16_29" "RA17_13" "RA17_22" "RA19_02" "RA19_10" "RA19_21" "RA20_05" "RA21_17" "TP6" "TP12")

dest_dir="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/het_plots/"
mkdir "$dest_dir"


for dataset in "${datasets[@]}"; do
  python "../../../tap_cn_calling/scripts/PLOT-unique_cn_profiles.py" \
  --cohort_name ${dataset} \
  --amp_gene_map_f "../../../tap_cn_calling/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt" \
  --cn_clone_profiles_csv "../results/${dataset}/opt_nclones/optimal_clone_profiles.csv" \
  --sample_sc_clone_assignment_csv "../results/${dataset}/opt_nclones/optimal_cell_assignments.csv" \
  --output_dir "../results/het_plots/" \
  --clone_cleanup False


done

echo "Script execution complete."

