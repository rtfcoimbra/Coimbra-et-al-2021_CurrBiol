library(tidyverse)
library(ape)
library(treeio)
library(ggtree)

# set working directory
setwd("")

# import tree
nj.tree <- read.newick("boot.tree.suptree", node.label = "support")

# remove outgroup
nj.tree.pruned <- drop.tip(nj.tree, "WOAK")

# check the node labels
#plot.phylo(as.phylo(nj.tree.pruned), no.margin = TRUE, cex = 0.75)
#nodelabels(cex = 0.75, frame = "circle")

# create the basic plot
nj.t1 <- ggtree(nj.tree.pruned,
                layout = "roundrect",
                ladderize = TRUE,
                right = TRUE) +
  # add scale bar
  geom_treescale(x = 0, y = 0, fontsize = 2.3) +
  # add tip labels
  geom_tiplab(size = 2.3, align = TRUE, linesize = 0.1) +
  # add node labels (support values)
  geom_nodepoint(aes(color = support, subset = !is.na(as.numeric(support))),
                 size = 1.5) +
  # annotate monophyletic clades
  geom_cladelabel(node = 53, label = "Northern",
                  color = c("#000000", "#000000"),
                  offset = 0.095, offset.text = 0.015,
                  barsize = 0, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 64, label = "West African",
                  color = c("#F0E442", "#000000"),
                  offset = 0.06, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 55, label = "Kordofan",
                  color = c("#E69F00", "#000000"),
                  offset = 0.06, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 59, label = "Nubian",
                  color = c("#D55E00", "#000000"),
                  offset = 0.06, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 68, label = "Reticulated",
                  color = c("#CC79A7", "#000000"),
                  offset = 0.095, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 89,
                  label = expression(paste("Masai ", italic("s. l."))),
                  color = c("#000000", "#000000"),
                  offset = 0.095, offset.text = 0.015,
                  barsize = 0, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 90,
                  label = expression(paste("Masai ", italic("s. str."))),
                  color = c("#009E73", "#000000"),
                  offset = 0.06, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 95, label = "Luangwa",
                  color = c("#006046", "#000000"),
                  offset = 0.06, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  geom_cladelabel(node = 78, label = "Southern",
                  color = c("#000000", "#000000"),
                  offset = 0.095, offset.text = 0.015,
                  barsize = 0, angle = -90, hjust = 0.5,
                  fontsize = 3.2, align = TRUE) +
  geom_cladelabel(node = 80, label = "Angolan",
                  color = c("#0072B2", "#000000"),
                  offset = 0.06, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  # annotate paraphyletic clades
  # note: This is a workaround the issue introduced by 'flip()' that requires
  #       image editing. The correct method is 'geom_strip()'. Here, I
  #       intentionaly placed the annotation of the paraphyletic South African
  #       giraffe in the wrong node to facilitate image editing later.
  geom_cladelabel(node = 72, label = "South African",
                  color = c("#56B4E9", "#000000"),
                  offset = 0.06, offset.text = 0.015,
                  barsize = 1, angle = -90, hjust = 0.5,
                  fontsize = 2.3, align = TRUE) +
  # add legend
  scale_colour_distiller("Bootstrap support",
                         palette = "Spectral",
                         direction = 1) +
  theme(legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.75, "line"),
        legend.position = "bottom")

# flip clades
nj.t2 <- flip(nj.t1, 52, 77) %>%
  flip(53, 68) %>%
  flip(90, 95)

# check node numbers interactively
#identify(nj.t2)

# show plot
nj.t2

# save plot in '.tif' format
ggsave("bionj_tree.tiff", width = 87, height = 170, units = "mm", dpi = 300)
