library(tidyverse)

# set working directory
setwd("")

# find input files
files <- list.files(pattern = "*.sfs")

# create a tibble with the number of homozygous, heterozygous, and total sites
# for all samples and calculate heterozygosity
tbl <- files %>%
  set_names(str_remove(files, ".sfs")) %>%
  map_dfr(read_table, .id = "sample", col_names = FALSE) %>%
  rename(hom_sites = X1,
         het_sites = X2) %>%
  mutate(total_sites = rowSums(across(where(is.numeric))),
         heterozygosity = het_sites/total_sites,
         subspecies = case_when(
           str_detect(sample, regex("WA"))                   ~ "West African",
           str_detect(sample, regex("GNP|SNR|ZNP"))          ~ "Kordofan",
           str_detect(sample, regex("ETH|MF"))               ~ "Nubian",
           str_detect(sample, regex("ISC|RET"))              ~ "Reticulated",
           str_detect(sample, regex("MA|SGR"))               ~ "Masai s. str.",
           str_detect(sample, regex("LVNP"))                 ~ "Luangwa",
           str_detect(sample, regex("BNP|KKR|MTNP|SUN|V23")) ~ "South African",
           str_detect(sample, regex("ENP|HNB"))              ~ "Angolan",
           TRUE                                              ~ NA_character_
         )
  )

# calculates mean, standard deviation (SD), standard error (SE),
# and confidence interval (CI)
tbl.sum <- tbl %>%
  group_by(sample, subspecies) %>%
  summarise(n = n(),
            mean = mean(heterozygosity),
            sd = sd(heterozygosity)
  ) %>%
  mutate(se = sd/sqrt(n)) %>%
  mutate(ci = se*qt((1-0.05)/2+0.5, n-1))

# set factor levels
fct.lvls <- c("West African", "Kordofan", "Nubian",
              "Reticulated",
              "Masai s. str.", "Luangwa",
              "South African", "Angolan")

# convert 'subspecies' to factor and set levels' order
tbl.sum$subspecies <- factor(tbl.sum$subspecies, levels = fct.lvls)

# set color palette
cbPalette <- c("West African"  = "#F0E442",
               "Kordofan"      = "#E69F00",
               "Nubian"        = "#D55E00",
               "Reticulated"   = "#CC79A7",
               "Masai s. str." = "#009E73",
               "Luangwa"       = "#006046",
               "South African" = "#56B4E9",
               "Angolan"       = "#0072B2")

# create basic plot for heterozygosity
ggplot(data = tbl.sum,
       mapping = aes(x = sample, y = mean)) +
  # add barplot
  geom_col(mapping = aes(fill = subspecies),
           width = 0.8) +
  # add error bars
  geom_errorbar(mapping = aes(ymin = mean-sd, ymax = mean+sd),
                position = "dodge",
                width = 0.4,
                size = 0.25) +
  # add matrix of panels defined by two column faceting variables
  facet_grid(cols = vars(subspecies),
             scales = "free_x",
             space = "free_x") +
  # change color palette
  scale_fill_manual(values = cbPalette) +
  # change y axis limit
  scale_y_continuous(limits = c(0, 6e-04)) +
  # rename Y label
  ylab("Heterozygosity") +
  # adjust plot appearance
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 7, angle = 90, hjust = 0.5),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 7),
        legend.position = "none")

# save plot in '.pdf' format
ggsave("heterozygosity.pdf", width = 174, height = 77, units = "mm")
