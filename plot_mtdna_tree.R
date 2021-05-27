library(tidyverse)
library(ape)
library(treeio)
library(ggtree)

# set working directory
setwd("")

# import tree
mt.tree <- read.iqtree("partitions.nex.treefile")

# remove outgroup
mt.tree.pruned <- drop.tip(mt.tree, "JN632674")

# check the node labels
#plot.phylo(as.phylo(mt.tree.pruned), no.margin = TRUE, cex = 0.75)
#nodelabels(cex = 0.75, frame = "circle")

# create the basic plot
mt.t1 <- ggtree(mt.tree.pruned,
                layout = "roundrect",
                ladderize = TRUE,
                right = TRUE) +
  # add scale bar
  geom_treescale(x = 0, y = 0, fontsize = 2.3) +
  # add tip labels
  geom_tiplab(size = 2.3, align = TRUE, linesize = 0.1) +
  # add node labels (support values)
  geom_nodepoint(aes(color = UFboot, subset = !is.na(as.numeric(UFboot))),
                 size = 1.5) +
  # annotate clades
  # note: Some clade annotations are knowngly wrong due to the placement of a
  #       few individuals outside its assigned (sub)species. To get the tree
  #       annotations as shown in the paper requires image editing. This was
  #       done just to facilitate the image editing process later on.
  geom_cladelabel(node = 54, label = "Northern",
                  color = c("#000000", "#000000"),
                  offset = 0.02, offset.text = 0.003,
                  barsize = 0, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 61, label = "West African",
                  color = c("#F0E442", "#000000"),
                  offset = 0.0125, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 55, label = "Kordofan",
                  color = c("#E69F00", "#000000"),
                  offset = 0.0125, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 65, label = "Nubian",
                  color = c("#D55E00", "#000000"),
                  offset = 0.0125, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 72, label = "Reticulated",
                  color = c("#CC79A7", "#000000"),
                  offset = 0.02, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 79,
                  label = expression(paste("Masai ", italic("s. l."))),
                  color = c("#000000", "#000000"),
                  offset = 0.02, offset.text = 0.003,
                  barsize = 0, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 80, label = "Luangwa",
                  color = c("#006046", "#000000"),
                  offset = 0.0125, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 87,
                  label = expression(paste("Masai ", italic("s. str."))),
                  color = c("#009E73", "#000000"),
                  offset = 0.0125, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 96, label = "Southern",
                  color = c("#000000", "#000000"),
                  offset = 0.02, offset.text = 0.003,
                  barsize = 0, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 91, label = "South African",
                  color = c("#56B4E9", "#000000"),
                  offset = 0.0125, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 96, label = "Angolan",
                  color = c("#0072B2", "#000000"),
                  offset = 0.0125, offset.text = 0.003,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  # add legend
  scale_colour_distiller("UFBoot support",
                         palette = "Spectral",
                         direction = 1) +
  theme(legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.75, "line"),
        legend.position = "bottom")

# flip clades
mt.t2 <- flip(mt.t1, 52, 78) %>%
  flip(79, 96) %>%
  flip(34, 81) %>%
  flip(27, 53) %>%
  flip(54, 72) %>%
  flip(69, 71)

# check node numbers interactively
#identify(mt.t2)

# show plot
mt.t2

# save plot in '.tif' format
ggsave("mtdna_tree.tiff", width = 87, height = 170, units = "mm", dpi = 300)
