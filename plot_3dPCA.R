library(ggplot2)
library(plot3D)

# set working directory
setwd("")

# read bamlist used with ANGSD
bams <- read.table("bamlist")[,1]

# extract sample names from bamlist
samples <- sub("/.*/", "", bams)
samples <- sub(".clean.bam", "", bams)

# read covariance matrix generated with PCAngsd
giraffe_cov <- as.matrix(read.table("snps.ld_pruned.hwe_filter.cov"))

# append sample names as row and column names to covariance matrix
dimnames(giraffe_cov) <- list(samples, samples)

# perform PCA
pca <- prcomp(giraffe_cov, scale = TRUE)

# scree plot
eigenval <- pca$sdev^2
explained_var <- 100*(eigenval/sum(eigenval))
df1 <- data.frame(prin_comp = c(seq(1, length(eigenval))),
                  explained_var)
ggplot(df1, aes(prin_comp, explained_var)) +
  geom_col(fill = c(rep("steelblue", 3),
                    rep("grey40", length(eigenval)-3))) +
  xlab("Principal components") +
  ylab("Explained variance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
# save scree plot in PDF
ggsave("scree_plot.pdf", width = 170, height = 92, units = "mm")

# create data frame
df2 <- as.data.frame(pca$x)

# add column for species
df2$sp <- ""
df2$sp[c(grep("WA|GNP|SNR|ZNP|ETH|MF", row.names(df2)))] <- "Northern"
df2$sp[c(grep("ISC|RET", row.names(df2)))] <- "Reticulated"
df2$sp[c(grep("LVNP|MA|SGR", row.names(df2)))] <- "Masai s. l."
df2$sp[c(grep("BNP|KKR|MTNP|SUN|V23|ENP|HNB", row.names(df2)))] <- "Southern"
# change the order of groups
df2$sp <- factor(df2$sp,
                 levels = c("Northern", "Reticulated", "Masai s. l.", "Southern"))

# add column for subspecies
df2$subsp <- ""
df2$subsp[c(grep("WA", row.names(df2)))] <- "West African"
df2$subsp[c(grep("GNP|SNR|ZNP", row.names(df2)))] <- "Kordofan"
df2$subsp[c(grep("ETH|MF", row.names(df2)))] <- "Nubian"
df2$subsp[c(grep("ISC|RET", row.names(df2)))] <- "Reticulated"
df2$subsp[c(grep("MA|SGR", row.names(df2)))] <- "Masai s. str."
df2$subsp[c(grep("LVNP", row.names(df2)))] <- "Luangwa"
df2$subsp[c(grep("BNP|KKR|MTNP|SUN|V23", row.names(df2)))] <- "South African"
df2$subsp[c(grep("ENP|HNB", row.names(df2)))] <- "Angolan"
# change the order of groups
df2$subsp <- factor(df2$subsp,
                    levels = c("West African", "Kordofan", "Nubian",
                               "Reticulated",
                               "Masai s. str.", "Luangwa",
                               "South African", "Angolan"))

# set color blind friendly palette
cbPalette <- c("West African" = "#F0E442",
               "Kordofan" = "#E69F00",
               "Nubian" = "#D55E00",
               "Reticulated" = "#CC79A7",
               "Masai s. str." = "#009E73",
               "Luangwa" = "#006046",
               "South African" = "#56B4E9",
               "Angolan" = "#0072B2")

# set point shapes
shapes <- c(18, 15, 17, 16)

# save PCA plot in PDF
pdf("3dpca_plot.pdf", width = 8, height = 6)

# 3D PCA plot
layout(matrix(c(1, 1, 1, 0,
                1, 1, 1, 0,
                1, 1, 1, 2,
                1, 1, 1, 2),
              nrow = 4, ncol = 4, byrow = TRUE))
par(mar = c(0, 0, 0, 0))
scatter3D(df2$PC1, df2$PC2, df2$PC3, bty = "g", theta = 30, phi = 45,
          colvar = NULL, colkey = FALSE, col = cbPalette[as.factor(df2$subsp)],
          pch = shapes[as.factor(df2$sp)], cex = 2.5, cex.lab = 2,
          xlab = paste("PC1 (", round(explained_var[1], 2), "%)", sep = ""),
          ylab = paste("PC2 (", round(explained_var[2], 2), "%)", sep = ""),
          zlab = paste("PC3 (", round(explained_var[3], 2), "%)", sep = ""))
plot.new()
legend("left",
       c("West African", "Kordofan", "Nubian",
         "Reticulated",
         "Masai s. str.", "Luangwa",
         "South African", "Angolan"),
       col = cbPalette,
       pch = c(rep(18, 3), 15, rep(17, 2), rep(16, 2)),
       pt.cex = 2.5, cex = 1.5, y.intersp = 1.5, bty = "n")

invisible(dev.off())
