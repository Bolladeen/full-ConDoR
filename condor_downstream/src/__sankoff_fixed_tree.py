# %%
"""
@Akhil Jakatdar
2023-09

Solves a small parsimony problem using the Sankoff algorithm on a fixed tree topology
Givens: 
    n - # of clones
    k - # of cn_states
    cn_profiles - dataframe of leaf profiles
    T - fixed tree phylogeny
Return:
    Labels(T) - copy number profile labels of internal nodes
"""

import numpy as np
import pandas as pd
import pickle as pkl
import ete3
from pathlib import Path
from ete_tree_style import ete_layout, add_edge_events

output_dir = Path('/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined')

amplicon_gene_mapping_f = '/Users/haochen/Desktop/Tapestri_analysis/copy_number/tap_cn_calling/reference/panel3559_analysis.gene_cytoband.formatted.manual_adjusted.txt'
amplicon_gene_mapping = pd.read_csv(amplicon_gene_mapping_f, sep='\t',)

# %%
sampel_name = 'RA17_13'

tree = pkl.load(open('/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined/RA17_13_HZ_ETE_tree.pkl', 'rb'))

# def get_sankoff_min_parsimony_tree(tree):

# first, fetch all cn_profiles
cn_profiles = {}

for l in tree.get_leaves():
    l_cn_profile = l.cn_profile
    cn_profiles[int(l.name)] = l_cn_profile
cn_profiles = pd.DataFrame(cn_profiles).T

# ----- only select amplicons that have >= 3 amplicons covering them -----
amplicon_gene_mapping = amplicon_gene_mapping.loc[amplicon_gene_mapping['amplicon_number'].isin(cn_profiles.columns)]
# we are only interested in genes that have >= 3 amplicons covering them
genes_of_interest = amplicon_gene_mapping.groupby('gene_name').filter(lambda x: len(x) >= 3)['gene_name'].unique()
amplicons_of_interest = amplicon_gene_mapping.loc[amplicon_gene_mapping['gene_name'].isin(genes_of_interest), 'amplicon_number'].unique()
print(f"[INFO] Number of amplicons on genes with >= 3 amplicons: {len(amplicons_of_interest)}")

cn_profiles = cn_profiles[amplicons_of_interest]
max_cn_state = int(max(cn_profiles.max()))


# %% Traverse the tree in post-order (bottom-up)

# get number of nodes in tree
num_nodes = len(tree.get_descendants()) + 1 # consider root node
DP = np.zeros((num_nodes, int(max_cn_state) + 1, cn_profiles.shape[1])) # T.size represents the number of nodes, cn_profiles.shape represents the number of amplicons
# fill DP with infinity
DP[:] = np.inf

node_name_num_map = {}
count = 0

for node in tree.traverse('postorder'):

    # leaf
    if node.is_leaf():
        cn_profiles_2d = np.array([cn_profiles.loc[int(node.name)] == cn_state for cn_state in range(int(max_cn_state) + 1)])
        DP[count, cn_profiles_2d] = 0
        node_name_num_map[node.name] = count
        count += 1

    # internal node or root
    else:
        # DP[count, :, :] = np.inf
        print(f"[INFO] Processing node {node.name}")
        for cn_state in range(int(max_cn_state) + 1):
            for amplicon in range(cn_profiles.shape[1]):

                cost = DP[node_name_num_map[node.children[0].name], np.arange(int(max_cn_state) + 1), amplicon] + DP[node_name_num_map[node.children[1].name], cn_state, amplicon] + abs(np.arange(int(max_cn_state) + 1) - cn_state)

                DP[count, cn_state, amplicon] = min(cost)
        node_name_num_map[node.name] = count
        count += 1

#find the minimum value for each amplicon
min_parsimony_score = DP[node_name_num_map[tree.get_tree_root().name]].min(axis=0)

# %% Backtrack to find the copy number states at internal nodes
for node in tree.traverse('preorder'):
    if node.is_leaf():
        parent = node.up
        distance = (abs(cn_profiles.loc[int(node.name), amplicons_of_interest] - parent.cn_profile)).sum()
        node.dist = distance
    else:
        cn_profile = DP[node_name_num_map[node.name]].min(axis=0)
        node.add_features(
            cn_profile = cn_profile
        )
        if node.is_root():
            node.dist = 0
        else:
            parent = node.up
            distance = (abs(cn_profile - parent.cn_profile)).sum()
            node.dist = distance
            
# %%
def _layout_func(node):
    ete_layout(node)
    add_edge_events(node)
ts = ete3.TreeStyle()
ts.layout_fn = _layout_func
# ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = False
# tree.show(tree_style=ts)
fig_file_name = str(output_dir / f"{sample_name}_HZ_ETE_tree.sankoffed.png")
_ = tree.render(fig_file_name, tree_style=ts, w=1800, units="mm") 

# %%
