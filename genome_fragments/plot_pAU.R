library(tidyverse)
library(ggtree)
library(patchwork)

# set working directory
setwd("")

# read tree topologies
trees <- read.tree("topologies.tree")

# plot topologies
trees.plot <- ggtree(trees, ladderize = TRUE, right = TRUE) +
  facet_wrap(~.id, ncol = 3) +
  geom_tiplab(size = 2.5) +
  xlim(0, 6.2) +
  theme(strip.text = element_text(size = 8))

# read p-AU table
data <- read_tsv("combined.au") %>%
  rename(frag_size = Fragment,
         topo_id = Topology,
         p_au = pAU)

# calculates mean, standard deviation (SD), standard error (SE),
# and confidence interval (CI)
data.sum <- data %>%
  group_by(frag_size, topo_id) %>%
  summarize(n = n(),
            mean = mean(p_au),
            median = median(p_au),
            sd = sd(p_au)
  ) %>%
  mutate(se = sd/sqrt(n)) %>%
  mutate(ci = se*qt((1-0.05)/2+0.5, n-1))

# fix order of topology IDs
data.sum$topo_id <- factor(data.sum$topo_id,
                           levels = c("Top1", "Top2", "Top3",
                                      "Top4", "Top5", "Top6",
                                      "Top7", "Top8", "Top9",
                                      "Top10", "Top11", "Top12",
                                      "Top13", "Top14", "Top15"))

# check if significant rejection was reached (change "TopN")
#head(subset(data.sum, topo_id == "Top1" & mean + ci <= 0.05))

# plot p-AU over fragment sizes
pAU.plot <- ggplot(data.sum, aes(x = frag_size/1000,
                                 y = mean,
                                 color = factor(topo_id))) +
  geom_point(size = 1.2) +
  geom_line(size = 0.5) +
  geom_errorbar(aes(ymin = mean - ci,
                    ymax = mean + ci),
                width = 0.5) +
  geom_hline(aes(yintercept = 0.95),
             color = "green",
             lty = 2) +
  geom_hline(aes(yintercept = 0.05),
             color = "red",
             lty = 2) +
  labs(x = "Fragment size (kbp)",
       y = expression(paste(italic("p"), "-AU"))) +
  scale_colour_viridis_d(direction = -1) +
  theme_minimal() +
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key = element_blank(),
        legend.key.size = unit(0.7, "lines"),
        legend.background = element_blank())

# plot composition
pAU.plot / trees.plot +
  plot_layout(heights = c(1,3)) +
  plot_annotation(tag_levels = "A")

# save plot in '.pdf' format
ggsave("au_test.pdf", width = 210, height = 297, units = "mm")
