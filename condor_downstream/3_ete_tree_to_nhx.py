# %% Import modules
import pickle as pkl
# import ete3
import pandas as pd
from pathlib import Path

# %% Import tree

refined_tree_dir = Path("/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined")

sample_names = [str(l).rsplit('/',1)[1].split("_HZ_ETE")[0] for l in refined_tree_dir.glob('*.pkl')]

for sample_i in sample_names:
    if not sample_i == "M12":
        continue

    tree_f = f"/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined/{sample_i}_HZ_ETE_tree.pkl"
    ete_tree = pkl.load(open(tree_f, "rb"))


    # # %% Parse for NHX and then ggtree
    # unique_features = set()
    for n in ete_tree.traverse():
        # for f in n.features:
        #     unique_features.add(f)
        if n.is_root():
            continue
        somatic_snv_gains = [n.somatic_snv_events[v][0] for v in n.somatic_snv_events if n.somatic_snv_events[v][1] == "GAIN"]
        # somatic_snv_gains = ["\n".join(somatic_snv_gains)]

        somatic_snv_lohs = [n.somatic_snv_events[v][0] for v in n.somatic_snv_events if n.somatic_snv_events[v][1] == "LOH"]
        # somatic_snv_lohs = ["\n".join(somatic_snv_lohs)]

        n.add_features(
            n_germline_snp_lohs = len(n.germline_snp_events),
            somatic_snv_gains = somatic_snv_gains,
            somatic_snv_lohs = somatic_snv_lohs,
        )
    # # %% Write NHX for ggtree in R
    nhx_f = f"/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined/{sample_i}_HZ_ETE_tree.nhx"

    with open(nhx_f, "w") as f:
        f.write(ete_tree.write(format=9, features=["dist","leaf_color","leaf_size","name", "n_germline_snp_lohs", "somatic_snv_gains", "somatic_snv_lohs"]))

# %%
