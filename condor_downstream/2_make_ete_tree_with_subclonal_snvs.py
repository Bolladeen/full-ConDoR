# %% import and setup
import sys, os
from copy import deepcopy
import pandas as pd
import numpy as np
import pickle
import ete3
# set cwd to script dir
os.chdir(os.path.dirname(os.path.realpath(__file__)))
print(f"cwd: {os.getcwd()}")

# per:
# - https://github.com/etetoolkit/ete/issues/101
# - https://github.com/ContinuumIO/anaconda-issues/issues/1806
os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["XDG_RUNTIME_DIR"] = "/juno/work/iacobuzc/tmp"
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from pathlib import Path

# get logger
# make sure to record all logging from modules
logging.basicConfig(stream=sys.stdout, level=logging.INFO)

# del sys.modules['src.tree_io']

from src.tree_io import nx_to_ete_tree, add_info_to_ete_tree
from src.plotting_utils import rgb_string_to_hex, __hgvsp_3to1
# clone_color_sequence = sns.cubehelix_palette(n_colors = len(ete_tree.get_leaves()),start=.0, rot=2, light=0.8, dark=0.5, hue=3.5).as_hex()
color_sequence = [rgb_string_to_hex(x) for x in px.colors.qualitative.Pastel]

from src.refine_condor_tree import (
    rm_blacklisted_events,
    find_diploid_clone,
    merge_clones_into_diploid,
    adjust_clone_number_and_color,
)
from src.sankoff_fixed_tree import get_sankoff_min_parsimony_tree

amplicon_gene_mapping_f = "../references/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.gene_pos_ado_added.csv"
# ----- read in sample anatomical name mapping file -----
sample_name_map = pd.read_excel('../references/Tapestri_batch2_samples_MASTER.xlsx').set_index('sample')['HZ_official_site_name'].to_dict()

# %% Read in data
pickle_dir = Path("../../condor_pipeline/condor_outputs/pickle_files")
manual_snv_dir = Path(
    "/home/zhangh5/work/Tapestri_batch2/batch2_data_compiled/manual_annotated_snv_lists"
)
output_dir = Path("../../condor_downstream/HZ_ete_trees_refined_subclonal_snvs")
output_dir.mkdir(exist_ok=True, parents=True)

patient_names = [str(f.stem).split("_self")[0] for f in pickle_dir.glob("*")]
# patient_name = "TP6"
for patient_name in patient_names:

    f = str(pickle_dir / f"{patient_name}_self.solT_cell")

    with open(f, "rb") as f:
        nx_tree = pickle.load(f)

    som_event_dict = {"0": "MISSING", "1": "GAIN", "2": "LOSS", "3": "LOH"}
    germ_event_dict = {"0": "MISSING", "2": "LOSS", "3": "LOH"}

    # %% Read in SNV annotation info
    manual_snv_f = manual_snv_dir / f"{patient_name}-patient-all_vars-voi.hz_curated.txt"
    manual_snv = pd.read_csv(manual_snv_f, sep="\t", index_col=0, comment="#")
    manual_snv["annotation"].fillna("", inplace=True)
    # manual_snv["var_formatted"] = manual_snv.index.str.replace(":", "_").str.replace(
    #     "/", "_"
    # )
    snv_annotation_map = manual_snv["HGVSp"].to_dict()

    for k, v in snv_annotation_map.items():
        snv_annotation_map[k] = __hgvsp_3to1(v)

    germline_events_set = set(
        manual_snv[manual_snv["annotation"].str.contains("germline")].index
    )
    # somatic_events = set(manual_snv[manual_snv['annotation'].str.contains('somatic')]['var_formatted'])
    other_events_set = set(
        manual_snv[~manual_snv["annotation"].str.contains("germline")].index
    )


    # %% Make ETE tree
    subtrees = {node: ete3.Tree(name=node) for node in nx_tree.nodes()}
    [*map(lambda edge: subtrees[edge[0]].add_child(subtrees[edge[1]]), nx_tree.edges)]
    ete_tree = subtrees["root"]

    """
    # sanity check
    for node in nx_tree.nodes():
        node_attributes = nx_tree.nodes()[node]
        if "cell_attachment" in node_attributes:
            print(node)
    """

    import logging

    logger = logging.getLogger(__name__)

    # ===== Params =====
    minimum_clone_size = 10
    internal_subclonal_snv = True

    # ==================
    subclonal_trees = {}
    internal_subclonal_snvs = {}

    internal_nodes_count = 0
    events_cache = []

    # convert edge events
    for node in ete_tree.traverse("preorder"):
        nx_node_attributes = nx_tree.nodes()[node.name]
        # 0. root
        if node.name == "root":
            continue
        # 1. Leafs
        elif "cn_profile" in nx_node_attributes:
            node.name = str(node.name)  # convert integer to string
            if "cell_attachment" in nx_node_attributes:
                cn_clone_cells = nx_node_attributes["cell_attachment"]
                cn_clone_size = len(cn_clone_cells)
            else:
                cn_clone_cells = []
                cn_clone_size = 0
            # (1A) with subclonal SNVs
            if len(node.children) > 0:
                logger.info(
                    f"[Leaf node] {node.name} with subclonal SNV tree, events: {node.children}"
                )
                subclonal_trees[node.name] = node.copy()
                # get all children's cells and remove all children from the node
                for subclone in node.get_descendants():
                    subclone_cells = nx_tree.nodes()[subclone.name]["cell_attachment"]
                    subclone_size = len(subclone_cells)
                    cn_clone_cells = list(set(cn_clone_cells) | set(subclone_cells))
                    cn_clone_size += subclone_size
                node.children = []
            else:
                logger.info(f"[Leaf node] {node.name} with no subclonal SNV tree")
            node.add_features(
                events=events_cache,
                cn_profile=pd.Series(nx_node_attributes["cn_profile"]),
                cell_attachment=cn_clone_cells,
                clone_size=cn_clone_size,
            )
            events_cache = []

        # 2. edge events
        elif len(node.children) == 1:
            if internal_subclonal_snv and ("cell_attachment" in nx_node_attributes):
                logger.info(f"[Edge event node] {node.name} with {len(nx_node_attributes['cell_attachment'])} cells attached.")
                # (i) have > minimum_clone_size cells attached
                if len(nx_node_attributes["cell_attachment"]) < minimum_clone_size:
                    logger.warning(
                        f">> node {node.name} has < {minimum_clone_size} cells attached"
                    )
                else:
                    # (ii) have only one child, which is a gain event
                    # if not (node.children[0].name.endswith("_1")):
                    #     logger.warning(f">> node {node.name} is the last internal subclonal event in the chain")
                    # else:
                    logger.info(f">> node {node.name} is a valid internal subclonal event")
                    event_name = node.name.split("_")[0]
                    chrom = event_name.split(":")[0]
                    genomic_coord = event_name.split(":")[1]
                    ref = event_name.split(":")[2].split("/")[0]
                    alt = event_name.split(":")[2].split("/")[1]
                    tapestri_condensed_format_event = (
                        f"{chrom}:{genomic_coord}:{ref}/{alt}-GAIN"
                    )
                    internal_subclonal_snvs[tapestri_condensed_format_event] = dict(
                        events=events_cache,
                        cell_attachment=nx_node_attributes["cell_attachment"],
                        clone_size=len(nx_node_attributes["cell_attachment"]),
                    )
            # otherwise, it's a normal edge event
            events_cache += [node.name]
            logger.debug(f"dropping edge event node {node.name}")
            node.delete()  # here we use delete so that children can be reconnected
            # nodes_to_delete += [node]

        # 3. internal node
        elif len(node.children) > 1:
            internal_nodes_count += 1
            events_cache += [node.name]
            node.name = f"IN_{internal_nodes_count}"
            logger.info(f"[Internal node] {node.name} -- events: {events_cache} ")
            node.add_features(events=events_cache)
            events_cache = []

        # 4. others???
        else:
            raise NotImplementedError(f"node {node.name} seems weird!")

    # # delete nodes
    # for node in nodes_to_delete:
    #     node.delete()


    # %% Add info to tree
    for node in ete_tree.traverse("preorder"):
        # add edge events
        if "events" in node.features:
            germline_snp_events = {}
            somatic_snv_events = {}
            events = []
            for e in node.events:
                if e == "root":
                    continue
                if not e.startswith("chr"):
                    # @HZ 2023-10-31 consider homdels which are not in mutation format
                    logging.warning(f"Encountered non-mutation event: {e}")
                    events += [e]
                    continue
                else:
                    event_type = str(e.split("_")[-1])
                    event_name = e.split("_")[0]
                    chrom = event_name.split(":")[0]
                    genomic_coord = event_name.split(":")[1]
                    ref = event_name.split(":")[2].split("/")[0]
                    alt = event_name.split(":")[2].split("/")[1]
                    tapestri_condensed_format_event = (
                        f"{chrom}:{genomic_coord}:{ref}/{alt}-{event_type}"
                    )
                    events.append(tapestri_condensed_format_event)
                    try:
                        event_name_annotated = snv_annotation_map[event_name]
                    except KeyError:
                        raise KeyError(
                            f"for node {node.name}, event {event_name} not found in annotation map"
                        )
                    # determine color and type of event
                    if event_name in germline_events_set:
                        try:
                            event_type = germ_event_dict[event_type]
                            event_color = "Green"
                        except:
                            event_type = "UNKNOWN"
                            event_color = "Black"
                        germline_snp_events[tapestri_condensed_format_event] = (
                            event_name_annotated,
                            event_type,
                            event_color,
                        )
                    else:
                        # print(f"{event_name} -- {event_type} -- ")
                        event_type = som_event_dict[event_type]
                        # print(f"{event_name} -- {event_type} -- ")
                        if event_type == "LOSS" or event_type == "LOH":
                            event_color = "Purple"
                        else:
                            event_color = "Red"
                        somatic_snv_events[tapestri_condensed_format_event] = (
                            event_name_annotated,
                            event_type,
                            event_color,
                        )

            node.add_features(
                events_orig=node.events,
                germline_snp_events=germline_snp_events,
                somatic_snv_events=somatic_snv_events,
                events=events,
            )
            logger.info(f"added events info to node {node.name}")
            logger.debug(
                f"node {node.name} -- {len(germline_snp_events)} germline SNP events; {len(somatic_snv_events)} somatic SNV events"
            )

    # %% Find diploid clone, merge clones into diploid, adjust clone number and color
    try:
        diploid_clone = [l for l in ete_tree.get_leaves() if (l.cn_profile == 2).all()][0]
        diploid_clone_name = diploid_clone.name
        n_colors = len(ete_tree.get_leaves())
    except IndexError:
        logging.warning("No diploid clone found. Creating a dummy one!")
        # create a dummy diploid node that is not attached to anything
        diploid_clone = ete3.TreeNode()
        diploid_clone.name = "dummy_diploid"
        diploid_clone_name = diploid_clone.name
        diploid_clone.clone_size = 0
        ete_tree.add_child(diploid_clone)
        n_colors = len(ete_tree.get_leaves()) + 1

    if not diploid_clone_name == "0":
        logging.warning("diploid clone is not 0")
    leaves_to_merge_to_diploid = merge_clones_into_diploid(
        ete_tree, diploid_clone, rm_clones_with_same_snv_events_as_diploid=True
    )
    n_colors -= len(leaves_to_merge_to_diploid)
    clone_rename_map, clone_palette = adjust_clone_number_and_color(
        ete_tree,
        diploid_clone_name,
        # color_sequence=sns.cubehelix_palette(
        #     n_colors=n_colors, start=0.0, rot=2, light=0.8, dark=0.5, hue=3.5
        # ).as_hex(),
        color_sequence = [rgb_string_to_hex(x) for x in px.colors.qualitative.Pastel],
    )

    for c in leaves_to_merge_to_diploid:
        clone_rename_map[c] = clone_rename_map[diploid_clone_name]
        print(f"[INFO] renaming clone {c} to diploid")

    for f in diploid_clone.features:
        if not f in ["dist", "support"]:
            ete_tree.__dict__[f] = diploid_clone.__dict__[f]

    diploid_clone.detach()

    # clean up the internal nodes that only have one child, with no events attached
    # this is necessary as the sankoff algorithm below requires a bifurcating tree (every node except for the root must have 2 children)
    for n in ete_tree.get_descendants():
        if len(n.get_children()) == 1:
            if n.somatic_snv_events == {} and n.germline_snp_events == {}:
                logging.warning(f"removing internal node {n.name} with no events attached")
                n.delete()
            else:
                logging.warning(
                    f"internal node {n.name} has events attached, movin it to its child"
                )
                n.get_children()[0].somatic_snv_events.update(n.somatic_snv_events)
                n.get_children()[0].germline_snp_events.update(n.germline_snp_events)
                n.delete()


    # %% Get internal node CN profiles and branch lengths
    amplicon_gene_mapping = pd.read_csv(
        amplicon_gene_mapping_f,
    )
    # we are only interested in genes that have >= 3 amplicons covering them
    genes_of_interest = (
        amplicon_gene_mapping.groupby("gene_name")
        .filter(lambda x: len(x) >= 3)["gene_name"]
        .unique()
    )
    amplicons_of_interest = amplicon_gene_mapping.loc[
        amplicon_gene_mapping["gene_name"].isin(genes_of_interest), "amplicon_number"
    ].unique()

    get_sankoff_min_parsimony_tree(
        ete_tree=ete_tree,
        amplicons_of_interest=amplicons_of_interest,
    )

    """
    del sys.modules['src.sankoff_fixed_tree']
    from src.sankoff_fixed_tree import get_sankoff_min_parsimony_tree
    """

    # %% Add subclonal SNVs
    # remember to rename the CN clones in the subclonal trees
    ete_tree_with_subclonal_snvs = deepcopy(ete_tree)

    subclonal_trees_updated = {}
    for cn_clone_name, cn_clone_node in subclonal_trees.items():
        cn_clone_node.name = clone_rename_map[cn_clone_name]
        print(f"""
        Orig clone name: {cn_clone_name}
        New clone name: {cn_clone_node.name}
            """)
        subclonal_trees_updated[cn_clone_node.name] = cn_clone_node
        
    minimum_clone_size = 10
    sub_node_count = 0
    for _, subtree in subclonal_trees_updated.items():
        # move all the subtree strcutures to the corresponding node in the main tree
        node = ete_tree_with_subclonal_snvs & subtree.name
        # print(node.name)
        for sub_node in subclonal_trees_updated[subtree.name].get_descendants():
            nx_node_attributes = nx_tree.nodes()[sub_node.name]
            cells_attached = nx_node_attributes["cell_attachment"]
            clone_size = len(cells_attached)
            if clone_size < minimum_clone_size:
                logger.warning(
                    f"sub_node {sub_node.name} has {clone_size} cells attached, below threshold {minimum_clone_size}; skipping..."
                )
                sub_node.detach()
                continue
            else:
                sub_node_count += 1

                event_name = sub_node.name.rsplit("_", 1)[0]
                loc = event_name.rsplit(":", 1)[0]
                chrom = loc.split(":")[0]
                genomic_coord = loc.split(":")[1]
                ref = event_name.split(":")[2].split("/")[0]
                alt = event_name.split(":")[2].split("/")[1]
                tapestri_condensed_format_event = f"{chrom}:{genomic_coord}:{ref}/{alt}-GAIN"
                event_name_annotated = snv_annotation_map[event_name]

                sub_node.orig_name = sub_node.name
                sub_node.name = f"SUB_NODE_{sub_node_count}"
                sub_node.add_features(
                    germline_snp_events={},
                    somatic_snv_events={tapestri_condensed_format_event: (event_name_annotated, "GAIN", "Red")},
                    cell_attachment=nx_node_attributes["cell_attachment"],
                    clone_size=len(nx_node_attributes["cell_attachment"]),
                    clone_color = node.clone_color
                )
                node.__dict__["clone_size"] -= len(nx_node_attributes["cell_attachment"])

        for child in subclonal_trees_updated[subtree.name].get_children():
            if not child in node.get_children():
                node.add_child(child)
                logger.info(f"added subclonal SNV tree to node {node.name}")
            else:
                logger.warning(f"child {child.name} already exists in node {node.name}")

    # %% Plot and save the tree

    from src.ete_tree_style import ete_layout, add_edge_mutation_events, beautify_tree

    beautify_tree(ete_tree)
    beautify_tree(ete_tree_with_subclonal_snvs)

    def _layout_fn(node):
        ete_layout(node)
        add_edge_mutation_events(node)

    ts = ete3.TreeStyle()
    ts.layout_fn = _layout_fn
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    # ts.rotation = 90
    # tree.show(tree_style=ts)
    # tree.render("%%inline", tree_style=ts)
    (output_dir / patient_name).mkdir(exist_ok=True, parents=True)
    fig_file_name = str(output_dir / patient_name / f"{patient_name}_ETE_tree.refined.png")
    _ = ete_tree.render(fig_file_name, tree_style=ts, w=int(1800), units="mm")

    fig_file_name = str(output_dir / patient_name / f"{patient_name}_ETE_tree.refined.subclonal_snvs.png")
    _ = ete_tree_with_subclonal_snvs.render(fig_file_name, tree_style=ts, w=int(1800), units="mm")


    tree_pickle_file_name = str(
        output_dir / patient_name / f"{patient_name}_HZ_ETE_tree.pkl"
    )
    with open(tree_pickle_file_name, "wb") as f:
        pickle.dump(ete_tree, f)

    tree_pickle_file_name = str(
        output_dir / patient_name / f"{patient_name}_HZ_ETE_tree.refined.subclonal_snvs.pkl"
    )
    with open(tree_pickle_file_name, "wb") as f:
        pickle.dump(ete_tree_with_subclonal_snvs, f)

    # %% Plot clone composition

    # for each clone, get its constitutient single cells

    sc_clone_assignment_df = pd.DataFrame(columns=["sample_name", "cell_barcode", "final_clone_id"])
    cn_clone_profiles = {}
    for l in ete_tree.get_leaves():
        cells = pd.Series(l.cell_attachment)
        cn_clone_profiles[str(l.name)] = pd.Series(l.cn_profile)
        sc_clone_assignment_df = pd.concat([
            sc_clone_assignment_df, 
            pd.DataFrame({
                "sample_name": cells.str.split("-",n=1).str[1],
                "cell_barcode": cells.str.split("-",n=1).str[0],
                "final_clone_id": l.name,
                })
        ])
    cn_clone_profiles_df = pd.DataFrame(cn_clone_profiles).T

    sample_compo_stat = sc_clone_assignment_df.groupby('sample_name')['final_clone_id'].value_counts()

    sample_compo_stat.name = 'cell_count'
    sample_compo_stat = sample_compo_stat.reset_index()
    sample_compo_stat = sample_compo_stat.pivot(index='sample_name', columns='final_clone_id',values='cell_count').fillna(0)

    sample_compo_stat.to_csv(output_dir / patient_name / f"{patient_name}_clone_compo.csv", index=True)

    sample_compo_stat.columns = sample_compo_stat.columns.astype(int)
    # filter out diploid clone:
    sample_compo_stat = sample_compo_stat.loc[:, sample_compo_stat.columns != clone_rename_map[diploid_clone_name]]

    # rename samples
    sample_compo_stat.index = sample_compo_stat.index.map(sample_name_map)

    # drop the samples that have fewer than 10 tumor cells:
    sample_compo_stat = sample_compo_stat.loc[sample_compo_stat.sum(axis=1) >= 10]

    # convert keys of clone_palette to integers
    clone_palette = {int(k): v for k, v in clone_palette.items()}

    # set the figure size
    sns.set(style="whitegrid")
    f1, ax = plt.subplots(figsize=(len(sample_compo_stat.columns)*1.5, 5), dpi=200)
    sample_compo_stat.plot(
        ax = ax,
        kind='barh', 
        stacked=True, 
        color = clone_palette,
        )
    # Add a legend and informative axis label
    ax.legend(ncol=1, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.set(
        ylabel="",
        xlabel="Number of Single Cells"
        )
    sns.despine(left=True, bottom=True)
    # # save figure
    # output_dir = Path('/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ETE_trees')
    (output_dir / patient_name / f"{patient_name}_clone_compo").mkdir(parents=True, exist_ok=True)
    fig_file_name = str(output_dir / patient_name / f"{patient_name}_clone_compo" / f"{patient_name}_clone_compo-refined.png")
    f1.savefig(fig_file_name, bbox_inches='tight', dpi=300)

    # normalize the clone sizes
    sample_compo_stat_norm = sample_compo_stat.div(sample_compo_stat.sum(axis=1), axis=0)
    f2, ax = plt.subplots(figsize=(len(sample_compo_stat.columns)*1.5, 5), dpi=200)
    sample_compo_stat_norm.plot(
        ax = ax,
        kind='barh',
        stacked=True,
        color = clone_palette,
        )
    # Add a legend and informative axis label
    ax.legend(ncol=1, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.set(
        ylabel="",
        xlabel="Fraction of Single Cells"
        )
    sns.despine(left=True, bottom=True)
    # # save figure
    # output_dir = Path('/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ETE_trees')
    (output_dir / patient_name / f"{patient_name}_clone_compo").mkdir(parents=True, exist_ok=True)
    fig_file_name = str(output_dir / patient_name / f"{patient_name}_clone_compo" / f"{patient_name}_clone_compo-refined-normed.png")
    f2.savefig(fig_file_name, bbox_inches='tight', dpi=300)

    # Clone profiles to file
    cn_clone_profiles_df.sort_index().to_csv(output_dir / patient_name / f"{patient_name}_final_clone_cn_profiles.csv", index=True)
    sc_clone_assignment_df.to_csv(output_dir / patient_name / f"{patient_name}_final_sc_clone_assignment.csv", index=False, )



# %% Across all trees, how many internal subclonal cells are here?
    
pickle_dir = Path("../../condor_pipeline/condor_outputs/pickle_files")
for f in pickle_dir.glob("*solT_cell"):
    print(f"Processing {f.stem}")
    with open(f, "rb") as f:
        nx_tree = pickle.load(f)
        for node in nx_tree.nodes():
            node_attributes = nx_tree.nodes()[node]
            if "cell_attachment" in node_attributes and str(node).endswith("_1"):
                print(f" {node} -- size {len(node_attributes['cell_attachment'])}")

# # %% Add internal subclonal SNVs
# for internal_event_name in internal_subclonal_snvs


# %%
