import ete3

#  Create the Layout function
def ete_layout(node):
    """
    Formatting of tree nodes while tree is rendered
    :param node: ete node
    """

    nstyle = ete3.NodeStyle()
    nstyle["shape"] = "sphere"

    nstyle["vt_line_width"] = 1     # line width of vertical branches
    nstyle["hz_line_width"] = 1     # line width of horizontal branches

    node.set_style(nstyle)

    if node.is_root():
        if not (node.cn_profile == 2).all():
            raise ValueError("Root node must be diploid!")
        
        root_text = "Diploid"
        tf_root = ete3.TextFace(root_text, ftype='Arial', fsize=10, fgcolor="White")
        tf_root.background.color = "Green"
        tf_root.margin_top = 10
        tf_root.margin_right = 10
        tf_root.margin_left = 10
        tf_root.margin_bottom = 10
        ete3.faces.add_face_to_node(tf_root, node, column=0, position='branch-right')
        # return
    # if node.is_leaf():
    if 'clone_size_for_plotting' in node.features:

        # nstyle['size'] = node.leaf_size / 10
        # nstyle['fgcolor'] = node.clone_color
        nstyle['size'] = 0
        if not node.name.startswith("SUB"):
            # a. rectangle face
            clone_face = ete3.AttrFace("name", ftype='Arial', fsize=10, fgcolor=node.clone_color, text_prefix='clone-', text_suffix='')
            clone_face.background.color = "White"
            # clone_face.border.type = 0
            # clone_face.border.width = 0
            # clone_face.margin_bottom = 0
            # clone_face.margin_top = 0
            ete3.faces.add_face_to_node(clone_face, node, column=1, position="branch-right") 

        # b. circle face
        clone_face = ete3.CircleFace(
            radius = node.clone_size_for_plotting / 20, 
            color = node.clone_color,
            style = 'circle',
            # label = {
            #     'text': f"clone-{node.name}",
            #     'font': 'Arial', 
            #     'fontsize': 5,
            #     'color': 'Black',
            # }
        )
        ete3.faces.add_face_to_node(clone_face, node, column=0, position="branch-right") 

def add_edge_mutation_events(node):
    # add mutation edge events
    if 'somatic_snv_events' in node.features or 'germline_snp_events' in node.features:
        for germline_snp in node.germline_snp_events:
            event_name_annotated, event_type, event_color = node.germline_snp_events[germline_snp]
            
            if event_type == 'LOSS':
                event_type = 'LOH'
                
            # add face to branch top
            event_face = ete3.TextFace(f"{event_name_annotated} {event_type}", fsize=8, ftype='Arial', fgcolor=event_color)
            # event_face.rotation = -30
            ete3.faces.add_face_to_node(event_face, node, column=2, position="branch-top")
        for somatic_snv in node.somatic_snv_events:
            event_name_annotated, event_type, event_color = node.somatic_snv_events[somatic_snv]

            if event_type == 'GAIN':
                event_type = 'mut'

            # add face to branch bottom
            event_face = ete3.TextFace(f"{event_name_annotated} {event_type}", fsize=8, ftype='Arial', fgcolor=event_color)
            # event_face.rotation = -30
            ete3.faces.add_face_to_node(event_face, node, column=2, position="branch-bottom")



# function to add more features for plotting
def beautify_tree(
    ete_tree,
    cn_clone_size_scale_factor = 800,
    ):
    # get all leaves' clone sizes
    cn_clone_sizes = {}
    for clone in ete_tree.get_descendants():
        if "clone_size" in clone.features:
            cn_clone_sizes[clone.name] = clone.clone_size 
    cn_clone_sizes_for_plotting = {k: v / max(cn_clone_sizes.values()) * cn_clone_size_scale_factor for k, v in cn_clone_sizes.items()}

    for clone in ete_tree.get_descendants():
        if "clone_size" in clone.features:
            clone.add_features(
                clone_size_for_plotting = cn_clone_sizes_for_plotting[clone.name],
            )
