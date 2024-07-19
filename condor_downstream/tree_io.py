import pandas as pd
import numpy as np
import ete3
import plotly.express as px
import logging

logger = logging.getLogger(__name__)

def nx_to_ete_tree(nx_tree):
    """
    Given a networkx tree, convert it to ete3 tree. Need to convert all edge events from nodes to connecting nodes' attributes
    """
    subtrees = {node:ete3.Tree(name=node) for node in nx_tree.nodes()}
    [*map(lambda edge:subtrees[edge[0]].add_child(subtrees[edge[1]]), 
    nx_tree.edges)]
    ete_tree = subtrees['root']

    # convert edge events
    internal_nodes_count = 0 
    events_cache = []
    for node in ete_tree.traverse('preorder'):
        # 1. internal node
        if len(node.children) > 1: 
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
    return ete_tree

# @HZ deprecated
def add_info_to_ete_tree(
        ete_tree, 
        cn_profiles, 
        cn_clone_sizes,
        cn_clone_palette,
        snv_annotation_map,
        germline_events_set, 
        # somatic_events_set,
        som_event_dict = {"0": "MISSING", "1": "GAIN", "2": "LOSS", "3": "LOH"},germ_event_dict = {"0": "MISSING", "2": "LOSS", "3": "LOH"},
        ):
    """
    We are adding:
    1. falcon clone profile 
    2. falcon clone size
    3. falcon clone color
    4. germline and somatic SNV events
    """
    # cn_clone_sizes_for_plotting = {k: v / max(cn_clone_sizes.values()) * 800 for k, v in cn_clone_sizes.items()}
    
    for leaf in ete_tree.iter_leaves():

        leaf_cn_profile = cn_profiles.loc[int(leaf.name), :]
        clone_size = cn_clone_sizes[int(leaf.name)]
        # clone_size_for_plotting = cn_clone_sizes_for_plotting[int(leaf.name)]
        leaf_color = cn_clone_palette[leaf.name]
        leaf.add_features(
            cn_profile = leaf_cn_profile,
            clone_size = clone_size,
            # clone_size_for_plotting = clone_size_for_plotting,
            leaf_color = leaf_color,
        )
        logger.info(f"added cn_clone info to leaf {leaf.name}")
        logger.debug(f"leaf {leaf.name} -- size: {clone_size}, color: {leaf_color}")

    for node in ete_tree.traverse('preorder'):
        # add edge events
        if 'events' in node.features:
            germline_snp_events = {}
            somatic_snv_events = {}
            events = []
            for e in node.events:
                if e == 'root':
                    continue
                if not e.startswith("chr"):
                    # @HZ 2023-10-31 consider homdels which are not in mutation format
                    logging.warning(f"Encountered non-mutation event: {e}")
                    events += [e]
                    continue
                else:
                    event_type = e.split('_')[-1]
                    event_name = e.rsplit('_', 1)[0]
                    loc = event_name.rsplit('_', 2)[0]
                    chr = loc.split('_')[0]
                    genomic_coord = loc.split('_')[1]
                    ref = event_name.rsplit('_', 2)[2]
                    alt = event_name.rsplit('_', 2)[1]
                    event_name = f"{loc}_{ref}_{alt}"
                    tapestri_condensed_format_event = f"{chr}:{genomic_coord}:{ref}/{alt}-{event_type}"
                    events.append(tapestri_condensed_format_event)
                    try:
                        event_name_annotated = snv_annotation_map[event_name]
                    except KeyError:
                        raise KeyError(f"for node {node.name}, event {event_name} not found in annotation map")
                    # determine color and type of event
                    if event_name in germline_events_set:
                        try:
                            event_type = germ_event_dict[event_type]
                            event_color = 'Green'
                        except:
                            event_type = 'UNKNOWN'
                            event_color = 'Black'
                        germline_snp_events[tapestri_condensed_format_event] = (event_name_annotated, event_type, event_color)

                    else:
                        event_type = som_event_dict[event_type]
                        if event_type == 'LOSS' or event_type == 'LOH':
                            event_color = 'Purple'
                        else:
                            event_color = 'Red'
                        somatic_snv_events[tapestri_condensed_format_event] = (event_name_annotated, event_type, event_color)
                
            node.add_features(
                events_orig = node.events,
                germline_snp_events = germline_snp_events,
                somatic_snv_events = somatic_snv_events,
                )
            node.add_features(
                events = events,
            )
            logger.info(f"added events info to node {node.name}")
            logger.debug(f"node {node.name} -- {len(germline_snp_events)} germline SNP events; {len(somatic_snv_events)} somatic SNV events")

