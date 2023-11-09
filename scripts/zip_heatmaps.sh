#!/bin/bash

# List of dataset names
datasets=("M04" "M07" "M11" "M12" "RA15_06" "RA15_16" "RA17_13" "RA17_22" "RA19_02" "RA19_10" "RA20_05")

myList=()
destination_path="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/heatmaps.zip"

# Loop through each dataset
for dataset in "${datasets[@]}"; do
  source_vaf="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/heatmaps/vaf_heatmap.png"
  myList+=("$source_vaf")
  source_vaf="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/heatmaps/condor_solution_heatmap.png"
  myList+=("$source_vaf")
  source_vaf="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/${dataset}/condor_outputs/out_tree.dot"
  myList+=("$source_vaf")

echo "$myList"
zip -r "$destination_path" "${myList[@]}"

# Check scp exit status
if [ $? -eq 0 ]; then
echo "Copied heatmap for $dataset successfully."
else
echo "Failed to copy heatmap for $dataset."
fi
done

