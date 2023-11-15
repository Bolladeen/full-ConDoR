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
# %%
pickle_dir = Path("/home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/results/pickle_files")
manual_snv_dir = Path("/home/zhangh5/work/Tapestri_batch2/batch2_data_compiled/manual_annotated_snv_lists")
output_dir = Path("/home/zhangh5/work/Tapestri_batch2/analysis/condor_downstream/HZ_ete_trees_refined_subclonal_snvs")
output_dir.mkdir(exist_ok=True, parents=True)

patient_name = "RA17_13"

f = str(pickle_dir / f'{patient_name}_self.solT_cell')

with open(f, 'rb') as f:
    nx_tree = pickle.load(f)

som_event_dict = {"0": "MISSING", "1": "GAIN", "2": "LOSS", "3": "LOH"}
germ_event_dict = {"0": "MISSING", "2": "LOSS", "3": "LOH"} 

# %% Read in SNV annotation info
manual_snv_f = manual_snv_dir / f'{patient_name}-patient-all_vars-voi.hz_curated.txt'
manual_snv = pd.read_csv(manual_snv_f, sep='\t', index_col=0, comment = '#')
manual_snv['annotation'].fillna('', inplace=True)
manual_snv['var_formatted'] = manual_snv.index.str.replace(':', '_').str.replace('/', '_')
snv_annotation_map = manual_snv.set_index('var_formatted')['HGVSp'].to_dict()

def __hgvsp_3to1(hgvsp_string):
    """
    1. identify all the 3-letter amino acid codes (3 letters after every capital letter)
    2. replace those 3-letter amino acid codes with 1-letter amino acid codes using Bio.SeqUtils.se1(). Keep other characters intact
    """
    pattern = re.compile(r'[A-Z][a-z][a-z]')
    return re.sub(pattern, lambda x: seq1(x.group()), hgvsp_string)

for k, v in snv_annotation_map.items():
    snv_annotation_map[k] = __hgvsp_3to1(v)

# %% Make ETE tree

subtrees = {node:ete3.Tree(name=node) for node in nx_tree.nodes()}
[*map(lambda edge:subtrees[edge[0]].add_child(subtrees[edge[1]]), 
nx_tree.edges)]
ete_tree = subtrees['root']

"""
# sanity check
for node in nx_tree.nodes():
    node_attributes = nx_tree.nodes()[node]
    if "cell_attachment" in node_attributes:
        print(node)
"""

# add two features from nx_tree:
# - cn_profiles
# - cell_attachments

minimum_clone_size = 5
internal_nodes_count = 0 
events_cache = []

for node in ete_tree.traverse():
    nx_node_attributes = nx_tree.nodes()[node.name]
    if str(node.name).startswith("root_") or str(node.name).startswith("subtree_"):
        node.delete()
    if nx_node_attributes != {}:
        for attribute in nx_tree.nodes()[node.name]:
            node.add_feature(attribute, nx_tree.nodes()[node.name][attribute])
    # if none of the attributes are added, then merge it
    else:
        


# %% convert edge events

minimum_clone_size = 5
internal_nodes_count = 0 
events_cache = []
for node in ete_tree.traverse('preorder'):
    if "cell_attachment" not in node.features





    # 1. internal node
    elif len(node.children) > 1: 
        internal_nodes_count += 1
        events_cache += [node.name]
        logger.info(f"internal node {node.name} -- events: {events_cache} ")
        node.add_features(
            events = events_cache
        )
        events_cache = []
        node.name = f"IN_{internal_nodes_count}"
    # 2. edge event
    elif len(node.children) == 1:
        events_cache += [node.name]
        # node.add_features(
        #     events = None
        # )
        logger.debug(f"dropping node {node.name}")
        node.delete() # here we use delete so that children can be reconnected
    # 3. leaf node
    else:
        node.name = str(node.name)
        node.add_features(
            events = events_cache,
        )
        events_cache = []

