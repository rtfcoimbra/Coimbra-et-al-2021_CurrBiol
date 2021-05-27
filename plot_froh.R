library(tidyverse)
library(patchwork)

# set working directory
setwd("")

# read csv of realized inbreeding coefficients (F_ROH)
f.roh <- tibble(read_csv("roh_outputs_giraffe.csv"))

# set factor levels
fct.lvls <- c("West African", "Kordofan", "Nubian",
              "Reticulated",
              "Masai s. str.", "Luangwa",
              "South African", "Angolan")

# set color palette
cbPalette <- c("West African"  = "#F0E442",
               "Kordofan"      = "#E69F00",
               "Nubian"        = "#D55E00",
               "Reticulated"   = "#CC79A7",
               "Masai s. str." = "#009E73",
               "Luangwa"       = "#006046",
               "South African" = "#56B4E9",
               "Angolan"       = "#0072B2")

# convert 'subspecies' to factor and set the levels' order
f.roh <- f.roh %>%
  mutate(class = case_when(class == "A" ~ "5120",
                           class == "B" ~ "1280",
                           class == "C" ~ "320",
                           class == "D" ~ "64",
                           class == "E" ~ "16",
                           class == "F" ~ "4"),
         class = fct_relevel(class, "5120", "1280", "320", "64", "16", "4"),
         subspecies = fct_relevel(subspecies, fct.lvls))

# create basic plot for F_ROH
p1 <- ggplot(data = f.roh,
             mapping = aes(x = sample,
                           y = f_roh,
                           fill = class)) +
  # add barplot
  geom_bar(stat = "identity",
           width = 0.8) +
  # add matrix of panels defined by two column faceting variables
  facet_grid(cols = vars(subspecies),
             scales = "free_x",
             space = "free_x") +
  # change color palette
  scale_fill_manual(name = "HBD classes",
                    values = c("#003f5c", "#444e86", "#955196",
                               "#dd5182", "#ff6e54", "#ffa600")) +
  # change y axis limit
  scale_y_continuous(limits = c(0, 0.8)) +
  # rename Y label
  ylab(expression("F"[ROH])) +
  # adjust plot appearance
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 7),
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.8, "line"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, reverse = TRUE))

###############################################################################

# read csv ROH segments
hbd.segs <- tibble(read_csv("giaffe_roh_sgmt.csv"))

# convert 'subspecies' to factor and set the levels' order
hbd.segs <- hbd.segs %>%
  mutate(subspecies = fct_relevel(subspecies, fct.lvls))

p2 <- ggplot(data = hbd.segs,
             mapping = aes(x = s_roh_mbp2,
                           y = n_roh2)) +
  geom_point(mapping = aes(color = subspecies,
                           shape = subspecies),
             size = 2) +
  scale_color_manual(name = "Subspecies",
                     labels = fct.lvls,
                     values = cbPalette) +
  scale_shape_manual(name = "Subspecies",
                     labels = fct.lvls,
                     values = c(18, 18, 18, 15, 17, 17, 16, 16)) +
  labs(x = "Sum of ROH (Mbp)",
       y = "Number of ROH") +
  theme(axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5),
        legend.title = element_text(face = "bold", size = 7),
        legend.text = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(1, "line"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

###############################################################################

# generate composed plot
p1 / p2 + plot_layout(nrow = 2, heights = c(1/3, 2/3))

# save plot in '.pdf' format
ggsave("froh.pdf", width = 174, height = 154, units = "mm")
