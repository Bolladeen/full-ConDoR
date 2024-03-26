# BiocManager::install("YuLab-SMU/treedataverse")
# BiocManager::install("phytools")
library(treedataverse)
set.seed(42)
# Load the necessary libraries
library(ggtree)
library(phytools)
library(dplyr)
library(ggrepel)
library(stringr)
library(patchwork)

# Read in NHX tree
NHX_tree_dir <- "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_downstream/NHX_trees"
OUTPUT_DIR <- "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_downstream"
nhx_tree_f <- file.path(NHX_tree_dir, "M12_HZ_ETE_tree.nhx")

NONFUNC_SO <- c('2kb_upstream_variant', '3_prime_UTR_variant', '5_prime_UTR_variant', 'intron_variant', 'synonymous_variant')

# Create a function to read in the NHX tree file and process it for plotting
process_nhx_tree <- function(nhx_tree_f) {
  nhx_tree <- read.nhx(nhx_tree_f)
  tree_tbl <- as_tibble(nhx_tree)

  # see if the tree is rooted: that is the row whose `parent` and `node` are equal, but with name as NA 
  # <---- this might be a better way than the below temporary fix
  tree_tbl$root <- tree_tbl$parent == tree_tbl$node & is.na(tree_tbl$name)
  # if no root
  # ===== @HZ TODO TEMPORARY FIX FOR M12 =====
  if (grepl("M12", nhx_tree_f)) {
    p <- as.phylo(nhx_tree)
    t2 <- bind.tip(tree=p, tip.label = "Diploid", where=2)
    d2 <- as_tibble(t2)
    d2$parent <- c(3,2,2) # alter the parental relationship
    t2 <- as.treedata(d2)
    # make sure to add the Diploid clone to tree_tbl
    tree_tbl <- bind_rows(tree_tbl, tibble(parent=2, node=2, name="Diploid", label="Diploid"))
    tree_tbl$parent <- c(2,3,3)
    tree_tbl$node <- c(1,2,3)
  }

  # normalize clone_size by total number of cells
  tree_tbl$clone_size <- tree_tbl$clone_size / sum(tree_tbl$clone_size, na.rm=TRUE)

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
d <- as_tibble(tree)
# attach the other info 
ggtree(t2, layout = "rectangular", ladderize = FALSE) +
  geom_label(aes(label=label))

# plot the tree
ggtree(tree, layout = "rectangular", ladderize = FALSE) +
  geom_label(aes(label=label)) + 
  geom_point(aes(size=clone_size), color="#7e7e7e") + guides(colour=FALSE, size=FALSE)


# Plot the tree with the added tags and annotations
ggtree(tree, ladderize = FALSE, layout = "roundrect", branch.length = "dist") +
  geom_tippoint(aes(size=clone_size), color="#7e7e7e") + 
  geom_nodepoint(aes(size=clone_size), color="#7e7e7e") + guides(colour=FALSE, size=FALSE) + 
  # geom_tiplab(size=5, color="black", hjust = -2) + 
  # add in leaf size: if not NA, size = clone_size; if NA, size = 5 (default)
  # only color the tips 
  # geom_point(aes(size = clone_size, color=isTip)) +
  theme_tree2() + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_gains_geom, paste0("(.{8})"), "\\1\n")), 
    color="#ff0000", vjust=0, size=10, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_somatic_snv_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#9500ffff", vjust=0.6, size=10, label.size=NA, fill=NA) + 
  geom_label(
    aes(x=branch, label=str_replace_all(n_germline_snp_lohs_geom, paste0("(.{8})"), "\\1\n")), 
    color="#008000", vjust=1.2, size=8, label.size=NA, fill=NA) + 
  theme_tree2() +
  scale_size_area(max_size=20) + 
  coord_fixed(ratio=50)


# ==========
# read in all the tree files
tree_files <- list.files(NHX_tree_dir, pattern = "*.nhx", full.names = TRUE)
# # @HZ TODO exclude M12 for now
# tree_files <- tree_files[!str_detect(tree_files, "M12")]

# create a list to store the plots
tree_objs <- list()

# plot each patient's tree, store in the list
for (tree_f in tree_files) {
  
  sample_name <- strsplit(basename(tree_f), "_HZ_ETE")[[1]][1]
  # if (sample_name == "M12" | sample_name == "RA19_21") {
  #   next
  # }

  tree <- process_nhx_tree(tree_f)
  tree_objs[[sample_name]] <- tree
# 
class(tree_objs) = "multiPhylo"
# add yscale which is the number of clones
# all_trees_plot <- 
ggtree(tree_objs, layout = "roundrect", branch.length = "dist", ladderize = FALSE, ) +
  facet_wrap(~.id, scales="free_x", nrow=2, strip.position = "top") + 
  geom_tippoint(aes(size=clone_size*100), color="#7e7e7e") + 
  geom_nodepoint(aes(size=clone_size*100), color="#7e7e7e") + guides(colour=FALSE, size=FALSE) + 
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
  scale_size_area(max_size=10) + 
  # scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) + 
  layout_dendrogram() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# save to PDF
ggsave(
  paste0(OUTPUT_DIR, "/all_trees_plot.pdf"),
  all_trees_plot, 
  device=pdf, 
  # family="Arial Unicode MS",
  width = 15, height = 20)

# save a PNG with dpi=400
ggsave(
  paste0(OUTPUT_DIR, "/all_trees_plot.png"),
  all_trees_plot, 
  width = 20, height = 15, dpi=400)

# ================= save a table 


