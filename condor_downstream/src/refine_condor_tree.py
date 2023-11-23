import logging
import plotly.express as px

logger = logging.getLogger(__name__)

# 1. Remove blacklisted events
def rm_blacklisted_events(
        ete_tree,
        germline_snp_blacklist,
        somatic_snv_blacklist,
    ):
    for node in ete_tree.traverse("preorder"):
        # Do some analysis on node
        if 'somatic_snv_events' not in node.features:
            continue
        else:
            # need to parse out the event name
            somatic_snv_events = set([x.split('-')[0] for x in node.somatic_snv_events])
            germline_snp_events = set([x.split('-')[0] for x in node.germline_snp_events])
            somatic_snvs_to_be_removed = []
            germline_snps_to_be_removed = []
            
            for x in [x for x in somatic_snv_events if x in somatic_snv_blacklist]:
                for y in node.somatic_snv_events:
                    if y.startswith(x):
                        somatic_snvs_to_be_removed.append(y)
            somatic_snvs_to_be_removed = set(somatic_snvs_to_be_removed)
            if len(somatic_snvs_to_be_removed) > 0:
                logger.warning(f"removing somatic SNV: {somatic_snvs_to_be_removed}")
                for y in somatic_snvs_to_be_removed:
                    del node.somatic_snv_events[y]
            
            for x in [x for x in germline_snp_events if x in germline_snp_blacklist]:
                for y in node.germline_snp_events:
                    if y.startswith(x):
                        germline_snps_to_be_removed.append(y)
            germline_snps_to_be_removed = set(germline_snps_to_be_removed)
            if len(germline_snps_to_be_removed) > 0:
                logger.warning(f"removing germline SNP: {germline_snps_to_be_removed}")
                for y in germline_snps_to_be_removed:
                    del node.germline_snp_events[y]
        
# 2. Find diploid clone     
def find_diploid_clone(ete_tree):
    leaf_levels = {}
    node_num_events = {}
    for node in ete_tree.traverse("preorder"):
        if node.is_root():
            logger.info(f"original root node: {node.name}")
            # original_root_node = node
            continue
        node.add_feature('level', len(node.get_ancestors()))
        # print(f"node: {node.name}, level: {node.level}")
        if node.is_leaf():
            leaf_levels[node.name] = node.level
        if 'somatic_snv_events' in node.features:
            node_num_events[node.name] = len(node.somatic_snv_events) + len(node.germline_snp_events)
        else:
            node_num_events[node.name] = 0

    # diploid clone should be the one with the least level and number of events
    root_candidates = [n for n in leaf_levels if leaf_levels[n] == min(leaf_levels.values())]
    if len(root_candidates) > 1:
        root_name = min(root_candidates, key=lambda x: node_num_events[x])
    else:
        root_name = root_candidates[0]
    logger.info(f"diploid clone should be {root_name}")
    return (ete_tree & root_name)

# 3. Merge clones into diploid (root)
def merge_clones_into_diploid(
    ete_tree, 
    diploid_clone, 
    min_leaf_size = 0.001,
    rm_clones_with_same_snv_events_as_diploid = False,
    rm_clones_with_no_somatic_snv_events = False,
    ):
    """
    (A) if the clone is smaller than (by default) 0.1% of total cell number #@HZ to-do: make this a parameter
    (B) if no event exists between a leaf and the root, merge the leaf into the root
    """
    leaves_to_merge_to_diploid = []
    total_cell_num = 0
    for leaf in ete_tree.get_leaves():
        total_cell_num += leaf.clone_size

    pre_diploid_events = []

    for n in diploid_clone.get_ancestors():
        if 'somatic_snv_events' in n.features:
            pre_diploid_events += list(n.somatic_snv_events.keys()) + list(n.germline_snp_events.keys())

    for leaf in ete_tree.get_leaves():
        logger.info(f"leaf: {leaf.name}, clone size: {leaf.clone_size}")

        if leaf.name == '':
            logger.info(f"leaf {leaf.name} has no name, removing")
            leaves_to_merge_to_diploid.append(leaf.name)
        elif leaf is diploid_clone:
            logger.info(f"leaf {leaf.name} is the diploid, skipping")
            continue
        elif leaf.clone_size < 0.001 * total_cell_num:
            logger.info(f"leaf {leaf.name} has a clone size of {leaf.clone_size}, smaller than {min_leaf_size*100}% of total cell number, merging into diploid")
            logger.warning(f"merging leaf {leaf.name} into diploid_clone {diploid_clone.name}")
            diploid_clone.clone_size += leaf.clone_size
            leaves_to_merge_to_diploid.append(leaf.name)
            # @HZ 2023-09-21: it's problematic to delete leaves here, as the internal node events will always get removed if there's <2 leaves left!
            # leaf.delete()
            continue
        elif 'somatic_snv_events' in leaf.features and (len(leaf.somatic_snv_events) > 0 or len(leaf.germline_snp_events) > 0):
            # logger.info(f"leaf {leaf.name} has both somatic SNV and germline SNP events, skipping")
            continue
        else:
            # get all events from the leaf to the root
            events = []
            somatic_snv_events = []
            for n in leaf.get_ancestors():
                if 'somatic_snv_events' in n.features:
                    somatic_snv_events += list(n.somatic_snv_events.keys())
                    events += list(n.somatic_snv_events.keys()) + list(n.germline_snp_events.keys())
            # print(event_count)
            if len(events) == 0:
                logger.info(f"leaf {leaf.name} has no events to the root, merging into diploid")
                logger.warning(f"merging leaf {leaf.name} into diploid {diploid_clone.name}")
                diploid_clone.clone_size += leaf.clone_size
                leaves_to_merge_to_diploid.append(leaf.name)
                # @HZ 2023-09-21: this is ok, because there's no internal node event to consider
                leaf.delete()
            elif events == pre_diploid_events and rm_clones_with_same_snv_events_as_diploid:
                logger.info(f"leaf {leaf.name} has the same events as the diploid, merging into diploid")
                logger.warning(f"merging leaf {leaf.name} into diploid {diploid_clone.name}")
                diploid_clone.clone_size += leaf.clone_size
                leaves_to_merge_to_diploid.append(leaf.name)
            elif len(somatic_snv_events) == 0 and rm_clones_with_no_somatic_snv_events:
                logger.info(f"leaf {leaf.name} has no somatic SNV events, merging into diploid")
                logger.warning(f"merging leaf {leaf.name} into diploid {diploid_clone.name}")
                diploid_clone.clone_size += leaf.clone_size
                leaves_to_merge_to_diploid.append(leaf.name)
            else:
                logger.info(f"leaf {leaf.name} has {len(events)} events to the root")

    # do one more iteration to remove the smaller clones and preserve internal node events
    for leaf in ete_tree.get_leaves():
        if not leaf.name in leaves_to_merge_to_diploid:
            continue
        else:
            __remove_leaf(leaf)
    return leaves_to_merge_to_diploid

def __remove_leaf(leaf):
    """
    remove a leaf from the tree, transferring its events to the other child of its parent if needed
    """
    parent = leaf.up
    if len(parent.get_children()) <= 2:
        # if the parent has only 2 or fewer children, merge the parent's events to the other child
        other_child = [x for x in parent.get_children() if x != leaf][0]
        logger.info(f"merging parent {parent.name}'s event to child {other_child.name}")
        if 'somatic_snv_events' in parent.features:
            if 'somatic_snv_events' in other_child.features:
                other_child.somatic_snv_events.update(parent.somatic_snv_events)
            else:
                other_child.add_features(somatic_snv_events=parent.somatic_snv_events)
        if 'germline_snp_events' in parent.features:
            if 'germline_snp_events' in other_child.features:
                other_child.germline_snp_events.update(parent.germline_snp_events)
            else:
                other_child.add_features(germline_snp_events=parent.germline_snp_events)
    leaf.delete()


def rgb_string_to_hex(rgb):
    """
    e.g. input: 'rgb(141,211,199)'
    """
    rgb = tuple(map(int, rgb[4:-1].split(',')))
    return '#%02x%02x%02x' % rgb


# 4. Adjust clone number
def adjust_clone_number_and_color(ete_tree, diploid_clone_name, color_sequence = px.colors.qualitative.Set3):

    if color_sequence[0].startswith('rgb'):
        color_sequence = [rgb_string_to_hex(x) for x in color_sequence]
    # visit the leaves by level, remove the diploid clone if it is present
    leaf_levels = {leaf.name: len(leaf.get_ancestors()) for leaf in ete_tree.get_leaves() if leaf.name != diploid_clone_name}

    leaf_levels = {k: v for k, v in sorted(leaf_levels.items(), key=lambda item: item[1])}
    clone_rename_map = {orig_clone_name: idx + 1 for idx, orig_clone_name in enumerate(leaf_levels.keys())}
    clone_rename_map = {diploid_clone_name: 0, **clone_rename_map}
    clone_palette = {idx: color_sequence[idx] for idx in clone_rename_map.values()}

    for leaf in ete_tree.get_leaves():
        clone_color = color_sequence[clone_rename_map[leaf.name]]
        if clone_color.startswith('rgb'):
            clone_color = rgb_string_to_hex(clone_color)
        leaf.add_features(
            clone_color = clone_color
        )
        print(f"renaming clone {leaf.name} to {clone_rename_map[leaf.name]}")
        leaf.name = str(clone_rename_map[leaf.name])
    
    return clone_rename_map, clone_palette
