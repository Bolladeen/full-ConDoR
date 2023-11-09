#!/bin/bash

source_dir="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results"
dest_dir="/n/fs/ragr-research/users/aj7381/falcon/pipeline/results/clonal_refinement_files"

mkdir "$dest_dir"

# Loop through each directory in the source directory
for dir in "$source_dir"/*/; do
  dir_name=$(basename "$dir")

  # Exclude directories named 'dot_files' or 'heatmaps'
  if [ "$dir_name" != "dot_files" ] && [ "$dir_name" != "heatmaps" ]; then
  # Copy 'clonal_refinement' directory into destination
    if [ -d "$dir/clonal_refinement" ]; then
      new_dir_name="${dir_name}_clonal_refinement"
      cp -r "$dir/clonal_refinement" "$dest_dir/$new_dir_name"
      echo "Copied $dir_name/clonal_refinement to $new_dir_name"
    fi
  fi
done

echo "Script execution complete."
