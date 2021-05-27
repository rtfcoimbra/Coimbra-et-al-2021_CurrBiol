library(ggplot2)
library(patchwork)

# set working directory
setwd("~/Sync/rcoimbra_phd/project_speciation/results/ngsadmix")
# read likelihoods list
fin <- read.table("likelihoods.list")
# add column with value of K
data <- cbind(c(rep("1", 100), rep("2", 100), rep("3", 100), rep("4", 100), rep("5", 100),
                rep("6", 100), rep("7", 100), rep("8", 100), rep("9", 100), rep("10", 100)),
              fin)
# add headers
colnames(data) <- c("K", "Likelihoods")

# write file formatted for CLUMPAK
write.table(data[, c(1, 2)],
            "clumpak.file",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# convert the variable K to a factor variable
data$K <- as.factor(data$K)

# create basic plot
p1 <- ggplot(data, aes(x = reorder(K, Likelihoods, FUN = median), y = Likelihoods)) +
  # add boxplot
  geom_boxplot(outlier.shape = NA, lwd = 0.4) +
  # add stripchart
  geom_jitter(aes(colour = "red", alpha = 0.5), size = 0.6, width = 0.15) +
  # change title and axis font size and remove legend
  theme(axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        legend.position = "none")

# get Delta K values from clumpak output
# bash: grep -Po 'Delta\(K=\d+\) = \K\d+.\d+' output.log > delta_k

# read list of delta K values
deltak <- read.table("delta_k")
# create data frame
df <- data.frame(K = c(seq(2, 9)),
                 DeltaK = deltak$V1)
# add rows for K=1 and K=10 with NA
df <- rbind(c(1, NA), df, c(10, NA))
# convert the variable K to a factor variable
df$K <- as.factor(df$K)

# create basic plot
p2 <- ggplot(df, aes(x = K, y = DeltaK, group = 1)) +
  # add lineplot
  geom_line(color = "red", alpha = 0.5) +
  # add scatterplot
  geom_point() +
  # change scale
  scale_y_continuous(breaks = c(0,2,4)) +
  # change title and axis font size and remove legend
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(hjust = 0.5),
        legend.position = "none")

# set plot composition
p1 / p2 + plot_layout(nrow = 2, heights = c(2, 1))

# save plot in '.pdf' format
ggsave("likelihoods_K.pdf", width = 170, height = 105, units = "mm")
