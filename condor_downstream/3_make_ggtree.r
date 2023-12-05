# BiocManager::install("YuLab-SMU/treedataverse")
library(treedataverse)
set.seed(42)
# Load the necessary libraries
library(ggtree)
library(dplyr)
library(ggrepel)
library(stringr)
library(patchwork)

# Read in NHX tree
NHX_tree_dir <- "/juno/work/iacobuzc/haochen/Tapestri_batch2/analysis/condor_downstream/NHX_trees"
nhx_tree_f <- file.path(NHX_tree_dir, "RA21_17_HZ_ETE_tree.nhx")

NONFUNC_SO <- c('2kb_upstream_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'intron_variant', 'synonymous_variant')

# Create a function to read in the NHX tree file and process it for plotting
process_nhx_tree <- function(nhx_tree_f) {
  nhx_tree <- read.nhx(nhx_tree_f)
  tree_tbl <- as_tibble(nhx_tree)

  # fill in the internal nodes' (the nodes with label as `NA`) names with IN-x, where x is the number of the internal node
  tree_tbl$label[is.na(tree_tbl$label)] <- paste0("IN-", 1:sum(is.na(tree_tbl$label)))

  # filter out all the somatic SNV gain/LOH elements that have any of the NONFUNC_SO as a substring, and filter out any NA or empty string
  tree_tbl$somatic_snv_gains <- sapply(
    str_split(tree_tbl$somatic_snv_gains, "\\|"), 
    function(x) paste(x[!str_detect(x, paste(NONFUNC_SO, collapse = "|"))], collapse = "|")
  )
  # convert empty strings to NA
  tree_tbl[tree_tbl == ""] <- NA
  
  # count the number of somatic SNV gains, ignore NAs
  tree_tbl$n_somatic_snv_gains <- sapply(str_split(tree_tbl$somatic_snv_gains, "\\|"), function(x) sum(x != "NA"))
  tree_tbl$n_somatic_snv_gains[is.na(tree_tbl$n_somatic_snv_gains)] <- 0

  # count the number of somatic SNV LOHs, assign 0 to NA values
  tree_tbl$n_somatic_snv_lohs <- sapply(str_split(tree_tbl$somatic_snv_lohs, "\\|"), function(x) sum(x != "NA"))
  tree_tbl$n_somatic_snv_lohs[is.na(tree_tbl$n_somatic_snv_lohs)] <- 0

  # create triangles to represent the number of somatic SNV gains and LOHs, ignore NAs # type out a triangle here --> ▲
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
  # add in leaf size: if not NA, size = clone_size; if NA, size = 5 (default)
  # only color the tips 
  geom_point(aes(size = clone_size, color=isTip)) +
  theme_tree2() + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_gains_geom, paste0("(.{8})"), "\\1\n")), 
    color="#ff0000", vjust=0, size=5, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#9500ffff", vjust=0.6, size=5, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_germline_snp_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#9500ffff", vjust=1.2, size=3, label.size=NA, fill=NA) + 
  theme_tree2()


# ==========
# read in all the tree files
tree_files <- list.files(NHX_tree_dir, pattern = "*.nhx", full.names = TRUE)

# create a list to store the plots
tree_objs <- list()
tree_plots <- list()

# plot each patient's tree, store in the list
for (tree_f in tree_files) {
  
  sample_name <- strsplit(basename(tree_f), "_HZ_ETE")[[1]][1]
  # if (sample_name == "M12" | sample_name == "RA19_21") {
  #   next
  # }

  tree <- process_nhx_tree(tree_f)
  tree_objs[[sample_name]] <- tree
  # tree_plots[[tree_f]] <- ggtree(tree, layout = "rectangular", branch.length = "dist") +
  #   geom_tippoint() + 
  #   geom_nodepoint() + 
  #   # geom_tiplab(size=5, color="black", hjust = -2) + 
  #   # add in leaf size: if not NA, size = clone_size; if NA, size = 5 (default)
  #   geom_point(aes(size = clone_size, color=isTip)) +
  #   theme_tree2() + 
  #   geom_label_repel(aes(x=branch, label=somatic_snv_gains), vjust=-0.5, color="steelblue", size=3, label.size=NA, fill=NA) + 
  #   geom_label_repel(aes(x=branch, label=somatic_snv_lohs), vjust=0.5, color="purple", size=3, label.size=NA, fill=NA) + 
  #   geom_label_repel(aes(x=branch, label=n_germline_snp_lohs_geom), vjust=0, color="#b51111", size=5, label.size=NA, fill=NA)
}

# 
class(tree_objs) = "multiPhylo"
# add yscale which is the number of clones
all_trees_plot <- ggtree(tree_objs, layout = "roundrect", branch.length = "dist") +
  facet_wrap(~.id, scales="free_y", ncol=1, strip.position = "left") + 
  geom_nodepoint() + 
  geom_tippoint(aes(size=clone_size), color="#7e7e7e") + guides(colour=FALSE, size=FALSE) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_gains_geom, paste0("(.{8})"), "\\1\n")), 
    color="#ff0000", vjust=0, size=5, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#9500ffff", vjust=0.5, size=5, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_germline_snp_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#008000", vjust=1.2, size=3, label.size=NA, fill=NA) + 
  theme_tree2() + 
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3)))
  # layout_dendrogram()

# save to PDF
ggsave(
  "/juno/work/iacobuzc/haochen/Tapestri_batch2/analysis/condor_downstream/all_trees_plot.pdf", 
  all_trees_plot, 
  device=cairo_pdf, family="Arial Unicode MS",
  width = 15, height = 20)

# save a PNG with dpi=400
ggsave(
  "/juno/work/iacobuzc/haochen/Tapestri_batch2/analysis/condor_downstream/all_trees_plot.png", 
  all_trees_plot, 
  width = 15, height = 20, dpi=400)


# use Patchwork to combine all the plots
combined_plots <- wrap_plots(tree_plots, ncol = 2)

# Display the combined plot
combined_plots




