# %% Import modules
import pickle as pkl
# import ete3
import pandas as pd
from pathlib import Path

# %% Import tree
refined_tree_dir = Path("/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_downstream/HZ_ete_trees_refined_subclonal_snvs")
output_dir = refined_tree_dir.parent
output_dir = output_dir / "NHX_trees"
output_dir.mkdir(exist_ok=True, parents=True)
patient_names = set([str(l).rsplit('/',1)[1].split("_HZ_ETE")[0] for l in refined_tree_dir.glob('**/*.pkl')])

event_table = pd.DataFrame(
    index=list(patient_names), 
    columns=["somatic_snv_gains", "somatic_snv_lohs", "germline_snp_events"]
    )

for patient_i in patient_names:
    # if not patient_i == "M12":
    #     continue

    tree_f = refined_tree_dir / patient_i / f"{patient_i}_HZ_ETE_tree.refined.subclonal_snvs.pkl"
    ete_tree = pkl.load(open(tree_f, "rb"))

    snv_gains = []
    snv_lohs = []
    snp_lohs = []
    for n in ete_tree.traverse():
        # for f in n.features:
        #     unique_features.add(f)
        if n.is_root():
            continue
        print(f"processing {n.name}")
        somatic_snv_gains_formatted = [v.split("-")[0] for v in n.somatic_snv_events if n.somatic_snv_events[v][1] == "GAIN"]
        somatic_snv_gains = [n.somatic_snv_events[v][0] for v in n.somatic_snv_events if n.somatic_snv_events[v][1] == "GAIN"]
        somatic_snv_lohs_formatted = [v.split("-")[0] for v in n.somatic_snv_events if n.somatic_snv_events[v][1] == "LOH"]
        somatic_snv_lohs = [n.somatic_snv_events[v][0] for v in n.somatic_snv_events if n.somatic_snv_events[v][1] == "LOH"]
        germline_snp_events_formatted = [v.split("-")[0] for v in n.germline_snp_events]
        n.add_features(
            n_germline_snp_lohs = len(n.germline_snp_events),
            somatic_snv_gains = somatic_snv_gains,
            somatic_snv_lohs = somatic_snv_lohs,
        )
        snv_gains += somatic_snv_gains_formatted
        snv_lohs += somatic_snv_lohs_formatted
        snp_lohs += germline_snp_events_formatted
    # # %% Write NHX for ggtree in R
    nhx_f = output_dir / f"{patient_i}_HZ_ETE_tree.nhx"

    with open(nhx_f, "w") as f:
        f.write(ete_tree.write(format=8, features=["dist","leaf_color","clone_size","name", "n_germline_snp_lohs", "somatic_snv_gains", "somatic_snv_lohs"]))
    
    event_table.loc[patient_i, "somatic_snv_gains"] = ",".join(snv_gains)
    event_table.loc[patient_i, "somatic_snv_lohs"] = ",".join(snv_lohs)
    event_table.loc[patient_i, "germline_snp_events"] = ",".join(snp_lohs)

# %%
event_table.to_csv(
    refined_tree_dir / "all_pts-event_table.csv"
)
