# ============================================================================
# grf_dis_plot_clean.R
# ----------------------------------------------------------------------------
# Purpose:    Compute Generalised Robinson-Foulds (GRF) distances between
#             trees using TreeDist::TreeDistance and write 2D MDS coords TSV.
#             Plotting itself is done in plot_grf_distance.py.
# Usage:      Rscript grf_dis_plot_clean.R --tree all_trees.tree \
#                  --coords grf_mds_coords.tsv
# Paper ref:  Suppl Fig 2 (GRF panels)
# Software:   R (>= 4.2), TreeDist, ape, optparse
# ============================================================================
#!/usr/bin/env Rscript

library(optparse)
library(TreeDist)
library(ape)

option_list <- list(
  make_option(c("-t", "--tree"), type="character", default=NULL),
  make_option(c("-c", "--coords"), type="character", default="grf_mds_coords.tsv")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tree)) {
  stop("You must provide the raw .tree file using --tree (-t).\n", call.=FALSE)
}

cat("[1/4] Reading trees...\n")
trees <- read.tree(opt$tree)

cat("[2/4] Calculating Generalized RF distances...\n")
grf_dist <- TreeDistance(trees)
dist_mat <- as.matrix(grf_dist)

cat("[3/4] Performing MDS...\n")
mds_res <- cmdscale(dist_mat, k = 2)
mds_df <- as.data.frame(mds_res)
colnames(mds_df) <- c("MDS1", "MDS2")

# Create TreeIDs (Tree0, Tree1, Tree2...) to match your Python map logic
mds_df$TreeID <- paste0("Tree", 0:(length(trees)-1))

cat("[4/4] Saving coordinates...\n")
write.table(mds_df, opt$coords, sep="\t", row.names=FALSE, quote=FALSE)
cat("[OK] Coords successfully saved to", opt$coords, "\n")