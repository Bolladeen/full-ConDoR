#!/bin/bash
datasets=("M04" "M07" "M11" "M12" "M13" "RA15_06" "RA15_16" "RA16_08" "RA16_29" "RA17_13" "RA17_22" "RA19_02" "RA19_10" "RA19_21" "RA20_05" "RA21_17" "TP6" "TP12")

source_dir="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results"
dest_dir="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/plot_info"
mkdir "$dest_dir"


for dataset in "${datasets[@]}"; do

  mkdir "${dest_dir}/${dataset}"
  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/pickle_files/${dataset}_self.solT_cell"
  dest="${dest_dir}/${dataset}/${dataset}_self.solT_cell"


  cp "$source" "$dest"

  echo "Copied $source to $dest"

  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/opt_nclones/optimal_clone_profiles.csv"
  dest="${dest_dir}/${dataset}/refined_clone_profiles.csv"

  cp "$source" "$dest"

  echo "Copied $source to $dest"

  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/opt_nclones/optimal_cell_assignments.csv"
  dest="${dest_dir}/${dataset}/refined_cell_assignments.csv"

  cp "$source" "$dest"

  echo "Copied $source to $dest"

  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/condor_inputs/character_vaf_matrix.csv"
  dest="${dest_dir}/${dataset}/character_vaf_matrix.csv"

  cp "$source" "$dest"

  echo "Copied $source to $dest"

  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/condor_inputs/somatic_mutations.txt"
  dest="${dest_dir}/${dataset}/somatic_mutations.txt"

  cp "$source" "$dest"

  echo "Copied $source to $dest"

  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/condor_inputs/germline_mutations.txt"
  dest="${dest_dir}/${dataset}/germline_mutations.txt"

  cp "$source" "$dest"

  echo "Copied $source to $dest"

  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/opt_nclones/optimal_clone_profiles.csv"
  dest="${dest_dir}/${dataset}/original_clone_profiles.csv"

  cp "$source" "$dest"

  echo "Copied $source to $dest"

  source="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/opt_nclones/optimal_cell_assignments.csv"
  dest="${dest_dir}/${dataset}/original_cell_assignments.csv"

  cp "$source" "$dest"

  echo "Copied $source to $dest"



done

echo "Script execution complete."
