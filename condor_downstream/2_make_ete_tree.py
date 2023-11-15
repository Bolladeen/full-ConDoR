# %%
import sys, os
import pandas as pd
import numpy as np
import pickle
import ete3
# per: 
# - https://github.com/etetoolkit/ete/issues/101
# - https://github.com/ContinuumIO/anaconda-issues/issues/1806
os.environ["QT_QPA_PLATFORM"]="offscreen"
os.environ["XDG_RUNTIME_DIR"]="/juno/work/iacobuzc/tmp"
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from pathlib import Path
# get logger
# make sure to record all logging from modules
logging.basicConfig(stream=sys.stdout, level=logging.INFO)

# clone_color_sequence = sns.cubehelix_palette(n_colors = len(ete_tree.get_leaves()),start=.0, rot=2, light=0.8, dark=0.5, hue=3.5).as_hex()

# del sys.modules['src.tree_io']

from src.tree_io import nx_to_ete_tree, add_info_to_ete_tree
from src.plotting_utils import rgb_string_to_hex
from src.refine_condor_tree import rm_blacklisted_events, find_diploid_clone, merge_clones_into_diploid, adjust_clone_number_and_color
from src.sankoff_fixed_tree import get_sankoff_min_parsimony_tree

import re
from Bio.SeqUtils import seq1

# %% Set up 
# get ete3 package path

# plot_info_dir = Path('/home/zhangh5/work/Tapestri_batch2/analysis/condor_downstream/AJ-all_cases_plot_info-before_refinement-oct03/plot_info-oct-2')

opt_nclones_dir = Path("/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/opt_nclones")
manual_snv_dir = Path('/n/fs/ragr-research/users/aj7381/falcon/transfer_to_RRG/manual_annotated_snv_lists')
condor_input_dir = Path("/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/condor_inputs")
pickle_dir = Path("/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/results/pickle_files")


output_dir = Path('/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/condor_downstream/AJ_ete_trees_refined')

# # For RA15_06
# plot_info_dir = Path("/home/zhangh5/work/Tapestri_batch2/analysis/per_case/M04")
# output_dir = plot_info_dir

output_dir.mkdir(exist_ok=True, parents=True)

patient_names = [str(l.stem).split('_self',1)[0] for l in pickle_dir.glob('*')]

# patient_names = [""]
# %%
# patient_names = ["M04"]
for patient_name in patient_names:
# patient_name = 'RA16_08'
    (output_dir / patient_name).mkdir(exist_ok=True, parents=True)
    if (output_dir / patient_name / f"{patient_name}_ETE_tree.refined.png").is_file():
        continue
    # if patient_name == 'TP6' or patient_name == 'TP12':
    #     continue
    print(f'[INFO] processing {patient_name}')
    f = pickle_dir / f'{patient_name}_self.solT_cell'

    with open(f, 'rb') as f:
        G = pickle.load(f)

    som_event_dict = {"0": "MISSING", "1": "GAIN", "2": "LOSS", "3": "LOH"}
    germ_event_dict = {"0": "MISSING", "2": "LOSS", "3": "LOH"} 


    # %% Load accompanying CN clone data for the tree
    # 
    # 1. CN clone profiles
    # 2. CN clone assignment for each single cell
    # 3. VAF matrix

    falcon_cn_profiles_f = opt_nclones_dir / f'{patient_name}' / f'optimal_clone_profiles.csv'
    falcon_cn_profiles = pd.read_csv(falcon_cn_profiles_f, index_col=0)

    cn_clone_assignment_f = opt_nclones_dir /f'{patient_name}' / 'optimal_cell_assignments.csv'
    cn_clone_assignment_df = pd.read_csv(cn_clone_assignment_f, index_col=0)
    cn_clone_assignment_df.index.name = 'sample_name'
    # get rid of sample RA18_18-11_1
    cn_clone_assignment_df = cn_clone_assignment_df.loc[cn_clone_assignment_df.index != 'RA18_18-11_1'].reset_index().rename(columns={'clone_id': 'falcon_clone_id'})

    condor_clone_assignment_f = condor_input_dir / patient_name / 'character_vaf_matrix.csv'
    condor_clone_assignment_df = pd.read_csv(condor_clone_assignment_f, index_col=0)['cluster_id'].reset_index().rename(columns={'cluster_id': 'condor_clone_id'})
    condor_clone_assignment_df['sample_name'] = condor_clone_assignment_df['index'].str.split('-',n=1).str[1]
    condor_clone_assignment_df['cell_barcode'] = condor_clone_assignment_df['index'].str.split('-',n=1).str[0]

    # merge above two clone assignment tables on sample_name and cell_barcode
    consensus_clone_assignment_df = pd.merge(cn_clone_assignment_df, condor_clone_assignment_df, on=['sample_name', 'cell_barcode'], how='inner')
    if not consensus_clone_assignment_df.shape[0] == condor_clone_assignment_df.shape[0]:
        raise ValueError('consensus clone assignment table has different number of rows than condor clone assignment table')

    falcon_to_condor_clone_map = consensus_clone_assignment_df[['falcon_clone_id', 'condor_clone_id']].drop_duplicates().set_index('falcon_clone_id')['condor_clone_id'].to_dict()
    # condor_to_falcon_clone_map = consensus_clone_assignment_df[['falcon_clone_id', 'condor_clone_id']].drop_duplicates().set_index('condor_clone_id')['falcon_clone_id'].to_dict()
    falcon_cn_profiles = falcon_cn_profiles.loc[consensus_clone_assignment_df['falcon_clone_id'].unique()]
    condor_cn_clone_profiles = falcon_cn_profiles.rename(index=falcon_to_condor_clone_map).sort_index()

    cn_clone_palette = {str(k): rgb_string_to_hex(px.colors.qualitative.Set3[i]) for i, k in enumerate(condor_cn_clone_profiles.index)}
    cn_clone_sizes = consensus_clone_assignment_df['condor_clone_id'].value_counts().to_dict()
    # make the biggest clone 800, and the rest scaled accordingly
    cn_clone_sizes_for_plotting = {k: v / max(cn_clone_sizes.values()) * 800 for k, v in cn_clone_sizes.items()}

    # read amplicon gene-name mapping file
    amplicon_gene_mapping_f =   '/n/fs/ragr-research/users/aj7381/CONDOR/full-ConDoR/references/amplicon_panel.csv'
    amplicon_gene_mapping = pd.read_csv(amplicon_gene_mapping_f, sep='\t',)
    amplicon_gene_mapping = amplicon_gene_mapping.loc[amplicon_gene_mapping['amplicon_number'].isin(condor_cn_clone_profiles.columns)]
    # we are only interested in genes that have >= 3 amplicons covering them
    genes_of_interest = amplicon_gene_mapping.groupby('gene_name').filter(lambda x: len(x) >= 3)['gene_name'].unique()
    amplicons_of_interest = amplicon_gene_mapping.loc[amplicon_gene_mapping['gene_name'].isin(genes_of_interest), 'amplicon_number'].unique()

    # %% Read in SNV annotation info
    manual_snv_f = manual_snv_dir / f'{patient_name}-patient-all_vars-voi.hz.txt'
    manual_snv = pd.read_csv(manual_snv_f, sep='\t', index_col=0, comment = '#')
    manual_snv['annotation'].fillna('', inplace=True)
    manual_snv['var_formatted'] = manual_snv.index.str.replace(':', '_').str.replace('/', '_')
    snv_annotation_map = manual_snv.set_index('var_formatted')['HGVSp'].to_dict()

    # for each value, parse out the amino acids (3 letters after every capital letter)

    def __hgvsp_3to1(hgvsp_string):
        """
        1. identify all the 3-letter amino acid codes (3 letters after every capital letter)
        2. replace those 3-letter amino acid codes with 1-letter amino acid codes using Bio.SeqUtils.se1(). Keep other characters intact
        """
        pattern = re.compile(r'[A-Z][a-z][a-z]')
        return re.sub(pattern, lambda x: seq1(x.group()), hgvsp_string)

    for k, v in snv_annotation_map.items():
        snv_annotation_map[k] = __hgvsp_3to1(v)



    germline_events_set = set(manual_snv[manual_snv['annotation'].str.contains('germline')]['var_formatted'])
    # somatic_events = set(manual_snv[manual_snv['annotation'].str.contains('somatic')]['var_formatted'])
    other_events_set = set(manual_snv[~manual_snv['annotation'].str.contains('germline')]['var_formatted'])
    # %% Make ETE tree

    ete_tree = nx_to_ete_tree(G)
    for n in ete_tree.traverse():
        if 'events' not in n.features:
            n.delete()

    add_info_to_ete_tree(
        ete_tree = ete_tree, 
        cn_profiles = condor_cn_clone_profiles, 
        cn_clone_sizes = cn_clone_sizes,
        cn_clone_palette = cn_clone_palette,
        snv_annotation_map = snv_annotation_map,
        germline_events_set = germline_events_set, 
        som_event_dict = {"0": "MISSING", "1": "GAIN", "2": "LOSS", "3": "LOH"},
        germ_event_dict = {"0": "MISSING", "2": "LOSS", "3": "LOH"},
        )

    # %% Refine the tree
    # rm_blacklisted_events(
    #     ete_tree,
    #     germline_snp_blacklist = snp_blacklist,
    #     somatic_snv_blacklist = snv_blacklist,
    # )

    # diploid_clone = find_diploid_clone(ete_tree)
    try:
        diploid_clone = [l for l in ete_tree.get_leaves() if (l.cn_profile == 2).all()][0]
        diploid_clone_name = diploid_clone.name
        n_colors = len(ete_tree.get_leaves())
    except IndexError:
        logging.warning('No diploid clone found. Creating a dummy one!')
        # create a dummy diploid node that is not attached to anything
        diploid_clone = ete3.TreeNode()
        diploid_clone.name = "dummy_diploid"
        diploid_clone_name = diploid_clone.name
        diploid_clone.leaf_size = 0 
        ete_tree.add_child(diploid_clone)
        n_colors = len(ete_tree.get_leaves()) + 1
        condor_cn_clone_profiles.loc['dummy_diploid'] = 2

    if not diploid_clone_name == '0':
        logging.warning('diploid clone is not 0')
    leaves_to_merge_to_diploid = merge_clones_into_diploid(
        ete_tree, diploid_clone, 
        rm_clones_with_same_snv_events_as_diploid = True
        )
    n_colors -= len(leaves_to_merge_to_diploid)
    clone_rename_map, clone_palette = adjust_clone_number_and_color(
        ete_tree, 
        diploid_clone_name,
        color_sequence = sns.cubehelix_palette(n_colors = n_colors, start=.0, rot=2, light=0.8, dark=0.5, hue=3.5).as_hex(),
        )

    for c in leaves_to_merge_to_diploid:
        clone_rename_map[c] = clone_rename_map[diploid_clone_name]
        print(f"[INFO] renaming clone {c} to diploid")

    # remember to rename the condor_cn_clone_profiles too
    condor_cn_clone_profiles.index = condor_cn_clone_profiles.index.astype(str)
    # drop the clones that are merged into diploid
    condor_cn_clone_profiles = condor_cn_clone_profiles.loc[~condor_cn_clone_profiles.index.isin(leaves_to_merge_to_diploid)]
    condor_cn_clone_profiles.rename(index=clone_rename_map, inplace=True)

    # # @HZ for organoids (no normal diploid clone)
    # # @HZ 2023-10-08 since we are removing the diploid clone below
    # # if the diploid clone is < min_leaf_size, remove it from the tree
    # min_leaf_size = 0.001
    # total_cell_num = consensus_clone_assignment_df.shape[0]
    # if diploid_clone.leaf_size < min_leaf_size * total_cell_num:
    #     from src.refine_condor_tree import __remove_leaf
    #     logging.warning(f"diploid clone {diploid_clone.name} has a clone size of {diploid_clone.leaf_size}, smaller than {min_leaf_size * 100}% of total cell number, removing")
    #     __remove_leaf(diploid_clone)
        
    # %% move the diploid clone to root; clean up the internal nodes that only have one child
    # move diploid clone's features to the root
    for f in diploid_clone.features:
        if not f in ['dist','support']:
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
                logging.warning(f"internal node {n.name} has events attached, movin it to its child")
                n.get_children()[0].somatic_snv_events.update(n.somatic_snv_events)
                n.get_children()[0].germline_snp_events.update(n.germline_snp_events)
                n.delete()


    # %% Get internal node CN profiles and branch lengths
    get_sankoff_min_parsimony_tree(
        tree = ete_tree,
        cn_profiles = condor_cn_clone_profiles,
        amplicons_of_interest = amplicons_of_interest,
        )

    # %% Plot tree

    from src.ete_tree_style import ete_layout, add_edge_mutation_events, beautify_tree
    beautify_tree(ete_tree)

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
    fig_file_name = str(output_dir / patient_name / f"{patient_name}_ETE_tree.refined.png")
    _ = ete_tree.render(fig_file_name, tree_style=ts, w=int(1800), units="mm")  

    tree_pickle_file_name = str(output_dir / patient_name / f"{patient_name}_HZ_ETE_tree.pkl")
    with open(tree_pickle_file_name, 'wb') as f:
        pickle.dump(ete_tree, f) 

    # %% Plot clone composition

    # ----- read in sample anatomical name mapping file -----
    sample_name_map = pd.read_excel('/home/zhangh5/work/publications/Tapestri_main_2023/Tapestri_batch2_samples_MASTER.xlsx').set_index('sample')['HZ_official_site_name'].to_dict()

    consensus_clone_assignment_df['tree_renamed_clone_id'] = consensus_clone_assignment_df['condor_clone_id'].astype(str).map(clone_rename_map).astype(str)
    sample_compo_stat = consensus_clone_assignment_df.groupby('sample_name')['tree_renamed_clone_id'].value_counts()

    # print(sample_compo_stat)

    sample_compo_stat.name = 'cell_count'
    sample_compo_stat = sample_compo_stat.reset_index()
    sample_compo_stat = sample_compo_stat.pivot(index='sample_name', columns='tree_renamed_clone_id',values='cell_count').fillna(0)

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

    # save renamed clone profiles file
    # remove clones merged into diploid
    # condor_cn_clone_profiles.index = condor_cn_clone_profiles.index.astype(str)
    condor_cn_clone_profiles = condor_cn_clone_profiles.loc[~condor_cn_clone_profiles.index.isin(leaves_to_merge_to_diploid)]
    condor_cn_clone_profiles.index = condor_cn_clone_profiles.index.astype(int)
    condor_cn_clone_profiles.sort_index().to_csv(output_dir / patient_name / f"{patient_name}_condor_clone_cn_profiles.csv", index=True)

    # need to rename the tree_renamed_clone_id for PLOT script of FALCON 
    consensus_clone_assignment_df[['sample_name', 'cell_barcode', 'falcon_clone_id', 'condor_clone_id', 'tree_renamed_clone_id']].to_csv(
        output_dir / patient_name / f"{patient_name}_consensus_clone_assignment.csv", 
        index=False, 
        columns=['sample_name', 'cell_barcode', 'falcon_clone_id', 'condor_clone_id', 'tree_renamed_clone_id']
        )
# %%

# for sample_i in ['M13', 'RA15_06', 'RA16_08', 'RA17_13', 'RA17_22']:
#     f = f'/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined/{sample_i}_consensus_clone_assignment.csv'
#     df = pd.read_csv(f)
#     df[['sample_name', 'cell_barcode','tree_renamed_clone_id']].rename(columns={'tree_renamed_clone_id': 'clone_id'}).to_csv(
#         f'/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined/{sample_i}_final_clone_assignment.csv',
#         index=False
#     )

# # %% deal with suspicious homdel amplicons

# amp_ado_f = '/home/zhangh5/work/Tapestri_project/cn-calling/reference/amp_blacklist/homdel_amps_count.csv'
# amp_ado_df = pd.read_csv(amp_ado_f, index_col=0)
# amp_ado_df['gene_name'] = amp_ado_df.index.str.split('-').str[0]

# amp_homdel_blacklist = amp_ado_df.loc[
#     (amp_ado_df['count'] > 1) & (amp_ado_df['ado_in_8_normals'] > 0.1),
#     'amplicon_number'
# ]

# amplicon_gene_mapping_f = '/home/zhangh5/work/Tapestri_project/cn-calling/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt'
# amplicon_gene_mapping = pd.read_csv(amplicon_gene_mapping_f, sep='\t',)

# patient_names = ["RA15_06"]

# # for amplicons in amp_homdel_blacklist, if they have homdel, assign its state to the gene's state (mode of other amplicons belonging to the same gene); if still 0, assign it np.nan
# for sample_i in patient_names:
#     f = output_dir / patient_name / f"{sample_i}_condor_clone_cn_profiles.csv"
#     cn_clone_profiles_df = pd.read_csv(f, index_col=0)
#     for amp in amp_homdel_blacklist:
#         for clone_i in cn_clone_profiles_df.index:
#             if cn_clone_profiles_df.loc[clone_i, amp] == 0:
#                 gene = amplicon_gene_mapping.loc[amplicon_gene_mapping['amplicon_number'] == amp, 'gene_name'].iloc[0]
#                 other_amplicons_in_gene = amplicon_gene_mapping.loc[amplicon_gene_mapping['gene_name'] == gene, 'amplicon_number'].tolist()
#                 other_amplicons_in_gene = [x for x in other_amplicons_in_gene if x in cn_clone_profiles_df.columns]

#                 inferred_gene_cn_val = cn_clone_profiles_df.loc[clone_i, other_amplicons_in_gene].mode()[0]
#                 if inferred_gene_cn_val > 0:
#                     cn_clone_profiles_df.loc[clone_i, amp] = inferred_gene_cn_val
#                     print(f'[INFO] {sample_i} clone {clone_i} amp {amp} inferred to be {inferred_gene_cn_val}')
#                 else:
#                     cn_clone_profiles_df.loc[clone_i, amp] = np.nan
#                     print(f'[INFO] {sample_i} clone {clone_i} amp {amp} inferred to be nan')
#     cn_clone_profiles_df.to_csv(f, index=True)

# # # %%

# # %%
