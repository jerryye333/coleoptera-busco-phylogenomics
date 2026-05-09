# ============================================================================
# tree_main_with_support.R
# ----------------------------------------------------------------------------
# Purpose:    Plot the main Coleoptera phylogenomic tree (PMSF 60% retained)
#             with UFBoot/gCF node support and superfamily colour bands.
# Inputs:     rename_root_pmsf_60_tips.tree (PMSF 60% retained tree with tips
#             relabelled to Family-Subfamily-Species via phylostuff/phylabel.R;
#             provided on Zenodo with the alignments and gene trees);
#             Coleoptera.csv (master taxonomy)
# Outputs:    PNG / PDF figure
# Paper ref:  Figure 5
# Software:   ggtree, ggtreeExtra, tidyverse, ape, ggrepel
# ============================================================================
# ==============================================================================
# STEP 1: SETUP
# ==============================================================================
rm(list = ls())
graphics.off()

library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(ape)
library(ggnewscale)
library(tidytree)
library(ggrepel)

# 1. Load Data
tree_file <- "rename_root_pmsf_60_tips.tree"
csv_file <- "Coleoptera.csv"

tree <- read.tree(tree_file)
meta_data <- read.csv(csv_file)

# 2. Match Data
tree_data <- data.frame(label = tree$tip.label) %>%
  mutate(join_name = sub(".*-", "", label)) %>%
  left_join(meta_data, by = c("join_name" = "name"))

# ==============================================================================
# STEP 2: DEFINE GROUPS (Including Outgroup)
# ==============================================================================
# We need a clean lookup table for every tip
tip_group_map <- tree_data %>%
  select(label, Order, Superfamily) %>%
  mutate(
    group_id = case_when(
      Order != "Coleoptera" ~ "Outgroup",
      Order == "Coleoptera" & !is.na(Superfamily) & Superfamily != "" ~ Superfamily,
      TRUE ~ NA_character_ # Unassigned Coleoptera (stay black)
    )
  )

# Convert to a named vector for fast lookup
tip_lookup <- setNames(tip_group_map$group_id, tip_group_map$label)

# ==============================================================================
# STEP 3: TOP-DOWN RECURSIVE COLORING
# ==============================================================================
message("--- Applying Top-Down Coloring Logic ---")

tree_tbl <- as_tibble(tree)

# Initialize all nodes to NA (Black)
node_data <- tree_tbl %>%
  select(node) %>%
  mutate(SF_Group = NA_character_)

# FUNCTION: Top-Down Search
# node: The current node ID to inspect
paint_top_down <- function(node, tree_tbl, current_node_data, tip_lookup) {

  # 1. Identify all tips under this node
  #    (Using tidytree::offspring is easiest, though recursion can be faster for huge trees)
  node_offspring <- offspring(tree_tbl, node)

  # Filter only valid tips (exclude internal nodes)
  descendant_tips <- node_offspring %>%
    filter(!is.na(label)) %>%
    pull(label)

  # Handle Case: The node itself IS a tip (offspring returns 0 rows usually)
  if (length(descendant_tips) == 0) {
    # Check if 'node' is a tip index
    if (node <= length(tree$tip.label)) {
      tip_label <- tree$tip.label[node]
      descendant_tips <- c(tip_label)
    }
  }

  # 2. Check the Groups of these tips
  #    Get the group_id for all descendant tips
  descendant_groups <- tip_lookup[descendant_tips]

  #    Remove NAs (tips we don't want to color anyway)
  valid_groups <- na.omit(unique(descendant_groups))

  # 3. DECISION LOGIC
  if (length(valid_groups) == 1) {
    # --- PURE CLADE ---
    # All valid descendants belong to the SAME group.
    # Color this node and ALL its descendants.
    group_name <- valid_groups[1]

    nodes_to_color <- c(node, node_offspring$node)
    current_node_data$SF_Group[current_node_data$node %in% nodes_to_color] <- group_name

    # STOP recursing here (we colored the whole subtree)
    return(current_node_data)

  } else {
    # --- MIXED CLADE ---
    # Contains multiple groups (e.g., Cucujoidea + Tenebrionoidea).
    # Do NOT color this node. Leave it black.
    # RECURSE: Go down to its children and check them.
    children <- child(tree_tbl, node)$node

    for (kid in children) {
      current_node_data <- paint_top_down(kid, tree_tbl, current_node_data, tip_lookup)
    }

    return(current_node_data)
  }
}

# Start the recursion from the Root Node
root_node <- rootnode(tree_tbl)$node
node_data <- paint_top_down(root_node, tree_tbl, node_data, tip_lookup)

# ==============================================================================
# STEP 4: COLOR PALETTE (High Contrast)
# ==============================================================================
found_groups <- sort(setdiff(unique(node_data$SF_Group), c(NA, "Outgroup")))
n_groups <- length(found_groups)

set.seed(999)
raw_colors <- grDevices::rainbow(n_groups, s = 0.8, v = 0.9)
shuffled_colors <- sample(raw_colors)
names(shuffled_colors) <- found_groups

final_colors <- c("Outgroup" = "black", shuffled_colors)

# ==============================================================================
# STEP 5: BAR POSITIONS (Closest 0.1% Fit)
# ==============================================================================
message("--- Calculating Bar Positions ---")

plot_data_final <- tree_tbl %>%
  left_join(node_data, by = "node")

p_base <- ggtree(tree, layout="rectangular") %<+% plot_data_final
dt <- p_base$data
tree_width <- max(dt$x, na.rm=TRUE)

# Per-character width estimate used to push the superfamily bar past the longest
# tip label in each cluster. Tuned to leave a clear gap at the current tip size.
char_width_factor <- tree_width * 0.008

bar_df <- data.frame()

groups_to_label <- found_groups
# Add Outgroup to label list if present
if ("Outgroup" %in% unique(node_data$SF_Group)) {
  groups_to_label <- c("Outgroup", groups_to_label)
}

for (grp in groups_to_label) {
  # Get all tips belonging to this group
  tips <- names(tip_lookup)[tip_lookup == grp & !is.na(tip_lookup)]

  if (length(tips) > 0) {
    tip_data <- dt %>% filter(label %in% tips)

    if (nrow(tip_data) > 0) {
      # Cluster Detection (Handles separated chunks)
      tip_data <- tip_data %>% arrange(y)
      y_diff <- c(1, diff(tip_data$y))
      tip_data$cluster <- cumsum(y_diff > 1)

      for (clus in unique(tip_data$cluster)) {
        clus_data <- tip_data %>% filter(cluster == clus)

        ymin <- min(clus_data$y)
        ymax <- max(clus_data$y)

        clus_data$name_width <- nchar(as.character(clus_data$label)) * char_width_factor
        clus_data$total_width <- clus_data$x + clus_data$name_width

        # Small gap beyond the longest tip label
        x_bar_start <- max(clus_data$total_width, na.rm=TRUE) + (tree_width * 0.008)

        bar_df <- rbind(bar_df, data.frame(
          label = grp,
          xmin = x_bar_start,
          xmax = x_bar_start,
          ymin = ymin,
          ymax = ymax,
          color = final_colors[grp]
        ))
      }
    }
  }
}

# ==============================================================================
# STEP 6: PARSE NODE LABELS (bootstrap/sCF/gCF)
# ==============================================================================
message("--- Parsing Node Support Values ---")

# IQ-TREE .cf.tree internal node labels are: "SH-aLRT/bootstrap/gCF/sCF"
#   1st value = SH-aLRT support
#   2nd value = UFBoot (bootstrap)
#   3rd value = gCF (gene concordance factor)
#   4th value = sCF (site concordance factor)
#
# We extract bootstrap (2nd) and gCF (3rd).
#   - If bootstrap == 100 -> draw a filled circle on the node; text shows only gCF
#   - Otherwise           -> text shows "bootstrap/gCF", no circle

node_support <- dt %>%
  filter(isTip == FALSE) %>%
  mutate(
    parts = strsplit(as.character(label), "/"),
    boot  = as.numeric(sapply(parts, `[`, 2)),
    gCF   = as.numeric(sapply(parts, `[`, 3)),
    node_label = case_when(
      is.na(boot) | is.na(gCF) ~ NA_character_,
      boot == 100 ~ as.character(gCF),
      TRUE ~ paste0(boot, "/", gCF)
    )
  )

# Subset of nodes with bootstrap == 100 — these get a filled circle marker
node_full_support <- node_support %>%
  filter(!is.na(boot) & boot == 100)

# ==============================================================================
# STEP 7: FINAL PLOT
# ==============================================================================
message("--- Generating Final Plot ---")

p <- ggtree(tree, layout="rectangular") %<+% plot_data_final
p$layers <- list()  # Remove ALL default black tree layers (horizontal + vertical)
p <- p +
  geom_tree(aes(color = SF_Group), size = 0.1) +
  geom_tiplab(size = 1.3, align = FALSE, color = "black", offset = 0.002) +
  # Filled circle marker for nodes with bootstrap == 100
  geom_point(data = node_full_support,
             aes(x = x, y = y),
             shape = 16, size = 0.3, color = "black") +
  # Node support labels on internal nodes: gCF only (boot=100) or boot/gCF (otherwise)
  geom_text(data = node_support %>% filter(!is.na(node_label)),
            aes(x = x, y = y, label = node_label),
            hjust = -0.15, vjust = 0.5, size = 1.15, color = "black") +
  scale_color_manual(values = final_colors, name = "Superfamily", na.value = "black") +
  theme_tree2() +
  guides(color = "none")

p_final <- p +
  geom_segment(data = bar_df,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax),
               color = bar_df$color,
               size = 0.8) +

  geom_text(data = bar_df,
            aes(x = xmin + (tree_width * 0.008),
                y = (ymin + ymax)/2,
                label = label),
            hjust = 0,
            size = 1.9,
            fontface = "bold",
            color = "black") +

  xlim(0, max(bar_df$xmin) + (tree_width * 0.4)) +

  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

print(p_final)
ggsave("Coleoptera_Final_TopDown_Support.png", p_final,
       width = 6.5, height = 15, dpi = 1200, limitsize = FALSE)
ggsave("Coleoptera_Final_TopDown_Support.pdf", p_final,
       width = 6.5, height = 15, limitsize = FALSE)
