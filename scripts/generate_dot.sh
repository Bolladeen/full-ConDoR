datasets=("M04" "M07" "M11" "M12" "M13" "RA15_06" "RA15_16" "RA16_08" "RA16_29" "RA17_13" "RA17_22" "RA19_02" "RA19_10" "RA19_21" "RA20_05" "RA21_17" "TP6" "TP12")


mkdir "/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/dot_files"

# Loop through each dataset
for dataset in "${datasets[@]}"; do
  source="/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/condor_outputs/${dataset}/out_tree.dot"
  # Execute scp command
  dot -Tpng "$source" -o "/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/dot_files/${dataset}.out_tree.png"

  # Check scp exit status
  if [ $? -eq 0 ]; then
  echo "Copied dot file for $dataset successfully."
  else
  echo "Failed to copy dot file for $dataset."
  fi
  done

