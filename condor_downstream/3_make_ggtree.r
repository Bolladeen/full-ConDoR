# BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)

# Load the necessary libraries
library(ggtree)
library(dplyr)
library(ggrepel)
library(stringr)
library(patchwork)

# Read in NHX tree
nhx_tree_f <- "/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined/RA17_13_HZ_ETE_tree.nhx"

# Create a function to read in the NHX tree file and process it for plotting
process_nhx_tree <- function(nhx_tree_f) {
  nhx_tree <- read.nhx(nhx_tree_f)
  tree_tbl <- as_tibble(nhx_tree)

  # fill in the internal nodes' (the nodes with label as `NA`) names with IN-x, where x is the number of the internal node
  tree_tbl$label[is.na(tree_tbl$label)] <- paste0("IN-", 1:sum(is.na(tree_tbl$label)))
  
  # convert empty strings to NA
  tree_tbl[tree_tbl == ""] <- NA

  # Create n stars to represent the number of somatic_snv_gains, and store in the n_somatic_snv_gain_or_loh_geom column, ignore when n_somatic_snv_gain_or_loh is NA

  # count the number of somatic SNV gains, ignore NAs
  tree_tbl$n_somatic_snv_gains <- sapply(str_split(tree_tbl$somatic_snv_gains, "\\|"), function(x) sum(x != "NA"))
  tree_tbl$n_somatic_snv_gains[is.na(tree_tbl$n_somatic_snv_gains)] <- 0

  # count the number of somatic SNV LOHs, assign 0 to NA values
  tree_tbl$n_somatic_snv_lohs <- sapply(str_split(tree_tbl$somatic_snv_lohs, "\\|"), function(x) sum(x != "NA"))
  tree_tbl$n_somatic_snv_lohs[is.na(tree_tbl$n_somatic_snv_lohs)] <- 0

  # create stars to represent the number of somatic SNV gains and LOHs, ignore NAs # type out a triangle here --> ▲
  tree_tbl$n_somatic_snv_gains_geom <- sapply(tree_tbl$n_somatic_snv_gains, function(x) paste0(rep("▲", x), collapse = ""))
  tree_tbl$n_somatic_snv_lohs_geom <- sapply(tree_tbl$n_somatic_snv_lohs, function(x) paste0(rep("▲", x), collapse = ""))

  # Split the labels in somatic_snv_gains by the "|" character and concatenate with a newline character
  tree_tbl$somatic_snv_gains <- sapply(str_split(tree_tbl$somatic_snv_gains, "\\|"), function(x) paste(x, collapse = "\n"))
  # substitute the NA values in somatic_snv_gains with an empty string
  tree_tbl$somatic_snv_gains[tree_tbl$somatic_snv_gains == "NA"]  <- ""
  tree_tbl$somatic_snv_lohs <- sapply(str_split(tree_tbl$somatic_snv_lohs, "\\|"), function(x) paste(x, collapse = "\n"))
  tree_tbl$somatic_snv_lohs[tree_tbl$somatic_snv_lohs == "NA"]  <- ""

  tree_tbl$n_germline_snp_lohs[is.na(tree_tbl$n_germline_snp_lohs)] <- 0

  # Create n squares to represent the number of germline_snp_lohs, and store in the n_germline_snp_lohs_geom column. Create a new line every 3 squares
  tree_tbl$n_germline_snp_lohs_geom <- sapply(tree_tbl$n_germline_snp_lohs, function(x) paste0(rep("■", x), collapse = ""))

  nhx_tree <- as.treedata(tree_tbl)
  return(nhx_tree)
}


tree <- process_nhx_tree(nhx_tree_f)

# plot the tree
ggtree(tree, layout = "rectangular") +
  geom_tiplab() +
  geom_tree(aes(color = node), size = 0.5)

# Plot the tree with the added tags and annotations
ggtree(tree, layout = "rectangular", branch.length = "dist") +
  geom_tippoint() + 
  geom_nodepoint() + 
  # geom_tiplab(size=5, color="black", hjust = -2) + 
  # add in leaf size: if not NA, size = leaf_size; if NA, size = 5 (default)
  # only color the tips 
  geom_point(aes(size = leaf_size, color=isTip)) +
  theme_tree2() + 
  geom_label(aes(x=branch, label=n_somatic_snv_gains_geom), vjust=-1, color="#ff0000", size=5, label.size=NA, fill=NA) + 
  geom_label(aes(x=branch, label=n_somatic_snv_lohs_geom), vjust=2, color="#9500ffff", size=5, label.size=NA, fill=NA) + 
  geom_label(aes(x=branch, label=n_germline_snp_lohs_geom), vjust=1, color="#9500ffff", size=5, label.size=NA, fill=NA)

# ==========
# read in all the tree files
tree_files <- list.files("/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/HZ_ete_trees_refined/", pattern = "*.nhx", full.names = TRUE)

# create a list to store the plots
tree_objs <- list()
tree_plots <- list()

# plot each patient's tree, store in the list
for (tree_f in tree_files) {
  
  sample_name <- strsplit(basename(tree_f), "_HZ_ETE")[[1]][1]
  if (sample_name == "M12" | sample_name == "RA19_21") {
    next
  }

  tree <- process_nhx_tree(tree_f)
  tree_objs[[sample_name]] <- tree
  # tree_plots[[tree_f]] <- ggtree(tree, layout = "rectangular", branch.length = "dist") +
  #   geom_tippoint() + 
  #   geom_nodepoint() + 
  #   # geom_tiplab(size=5, color="black", hjust = -2) + 
  #   # add in leaf size: if not NA, size = leaf_size; if NA, size = 5 (default)
  #   geom_point(aes(size = leaf_size, color=isTip)) +
  #   theme_tree2() + 
  #   geom_label_repel(aes(x=branch, label=somatic_snv_gains), vjust=-0.5, color="steelblue", size=3, label.size=NA, fill=NA) + 
  #   geom_label_repel(aes(x=branch, label=somatic_snv_lohs), vjust=0.5, color="purple", size=3, label.size=NA, fill=NA) + 
  #   geom_label_repel(aes(x=branch, label=n_germline_snp_lohs_geom), vjust=0, color="#b51111", size=5, label.size=NA, fill=NA)
}

# 
class(tree_objs) = "multiPhylo"
all_trees_plot <- ggtree(tree_objs, layout = "roundrect", branch.length = "dist") + 
  facet_wrap(~.id, scales="fixed", ncol=1, strip.position = "left") + 
  geom_tippoint() + 
  geom_nodepoint() + 
  geom_tippoint(aes(size=leaf_size), color="#7e7e7e") + guides(colour=FALSE, size=FALSE) + 
  geom_label(aes(x=branch, label=n_somatic_snv_gains_geom), vjust=-1, color="#ff0000", size=3, label.size=NA, fill=NA) + 
  geom_label(aes(x=branch, label=n_somatic_snv_lohs_geom), vjust=2, color="#9500ffff", size=3, label.size=NA, fill=NA) + 
  geom_label(aes(x=branch, label=n_germline_snp_lohs_geom), vjust=1, color="#9500ffff", size=3, label.size=NA, fill=NA) + 
  coord_fixed(ratio = 16) + 
  theme_tree2()

# save to PDF
ggsave(
  "/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/ggtree_plots/all_trees_plot.pdf", 
  all_trees_plot, 
  device=cairo_pdf, family="Arial Unicode MS",
  width = 10, height = 20)

# save a PNG with dpi=300
ggsave(
  "/Users/haochen/Desktop/Tapestri_analysis/Tapestri_data_batch2/falcon+condor/ggtree_plots/all_trees_plot.png", 
  all_trees_plot, 
  width = 10, height = 20, dpi=300)


# use Patchwork to combine all the plots
combined_plots <- wrap_plots(tree_plots, ncol = 2)

# Display the combined plot
combined_plots




