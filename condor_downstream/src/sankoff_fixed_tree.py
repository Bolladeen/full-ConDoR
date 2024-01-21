# import everything
import numpy as np
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def get_sankoff_min_parsimony_tree(ete_tree, amplicons_of_interest=None):
    """
    Label internal vertices' copy number states, calculate branch lengths using the Sankoff algorithm

    Parameters
    ----------
    ete_tree: ete3.Tree
        The leaves should all have a cn_profile attribute
    
    amplicons_of_interest: list of str
        A list of amplicon names that we are interested in
    """

    # sanity check:
    for l in ete_tree.get_leaves():
        if not hasattr(l, "cn_profile"):
            raise ValueError(f"Leaf {l.name} does not have a cn_profile attribute!")

    # construct a cn_profiles dataframe
    cn_profiles = pd.DataFrame(
        [l.cn_profile for l in ete_tree.get_leaves()],
        index = [l.name for l in ete_tree.get_leaves()]
    )

    if amplicons_of_interest is None:
        amplicons_of_interest = cn_profiles.columns
    else:
        logger.warning(f"""
            amplicons of interest: {
            [a for a in amplicons_of_interest if not a in cn_profiles.columns]
            } not in the tree, ignoring
        """)
        amplicons_of_interest = [a for a in amplicons_of_interest if a in cn_profiles.columns]
    # aoi_index = [i for i, amplicon in enumerate(cn_profiles.columns) if amplicon in amplicons_of_interest]
    max_cn_state = int(max(cn_profiles[amplicons_of_interest].max()))

    num_nodes = len(ete_tree.get_descendants()) # don't consider root node
    DP = np.zeros((num_nodes, int(max_cn_state) + 1, cn_profiles.shape[1]))
    DP[:] = np.inf

    node_name_num_map = {}
    count = 0

    for node in ete_tree.traverse('postorder'):
        
        if node.is_root():
            continue
        # leaf
        elif node.is_leaf():
            cn_profiles_2d = np.array([node.cn_profile == cn_state for cn_state in range(int(max_cn_state) + 1)])
            DP[count, cn_profiles_2d] = 0
            node_name_num_map[node.name] = count
            count += 1

        # internal node
        else:
            # DP[count, :, :] = np.inf
            print(f"[INFO] Processing node {node.name}")
            for cn_state in range(int(max_cn_state) + 1):
                for amplicon in range(cn_profiles.shape[1]):

                    cost = DP[node_name_num_map[node.children[0].name], np.arange(int(max_cn_state) + 1), amplicon] + DP[node_name_num_map[node.children[1].name], cn_state, amplicon] + abs(np.arange(int(max_cn_state) + 1) - cn_state)

                    DP[count, cn_state, amplicon] = min(cost)
            node_name_num_map[node.name] = count
            count += 1

    # find the minimum value for each amplicon
    # min_parsimony_score = DP[node_name_num_map[ete_tree.get_tree_root().name]].min(axis=0)

    # backtrack to find the copy number states at internal nodes
    for node in ete_tree.traverse('preorder'):
        if node.is_root():
            if not "cn_profile" in node.features:
                cn_profile = pd.Series(
                    2, index = cn_profiles.columns
                )
                node.add_features(
                    cn_profile = cn_profile
                )
            elif not (node.cn_profile == 2).all():
                raise ValueError("Root node must be diploid!")
        elif node.is_leaf():
            cn_profile = node.cn_profile
            parent = node.up
            distance = (abs(node.cn_profile[amplicons_of_interest] - parent.cn_profile[amplicons_of_interest])).sum()
            node.dist = distance
        else:
            cn_profile = pd.Series(
                DP[node_name_num_map[node.name]].argmin(axis=0),
                index=cn_profiles.columns
            )
            node.add_features(
                cn_profile = cn_profile
            )
            parent = node.up
            distance = (abs(cn_profile[amplicons_of_interest] - parent.cn_profile[amplicons_of_interest])).sum()
            node.dist = distance

def update_edge_dist(ete_tree, snv_norm=3, cnv_norm=200):
    """
    Update edge distances considering both SNV events and CNV events
    """
    
            
