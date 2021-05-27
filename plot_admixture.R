# Multiple plot function
#
# Source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# load libraries
library(ggplot2)
library(reshape2)

# set working directory
setwd("")

# read bamlist used with ANGSD
bams <- read.table("bamlist")[,1]

# extract sample names from bamlist
samples <- sub("/.*/", "", bams)
samples <- sub(".clean.bam", "", samples)

# create a rank to keep the sample order of the bamlist
rank <- c(seq(1, length(samples)))

# get input file names
file.names <- dir("./", pattern = ".qopt")

# generate an admixture barplot for each K
for(i in 1:length(file.names)){
  # read input Q-matrix
  q.matrix <- read.table(paste0("./", file.names[i]))
  # create a data frame with rank, sample names, and Q-matrix
  data <- cbind(rank, samples, q.matrix)
  # reformat the data frame for plotting
  data.m <- melt(data, id = c("rank", "samples"),
                 value.name = "proportion",
                 variable.name = "ancestry")

  # assign plot to a variable
  assign(paste0("p", i),
         # read data frame and map aesthetics for plotting
         ggplot(data.m, aes(x = rank, y = proportion, fill = ancestry)) +
           # add barplot
           geom_bar(stat = "identity") +
           # add solid vertical lines between species
           geom_vline(xintercept = c(16.5, 26.5, 38.5),
                      color = "black",
                      lwd = 1.25) +
           # add dashed vertical lines between species
           geom_vline(xintercept = c(5.5, 10.5, 44.5),
                      color = "black",
                      lty = 2) +
           # add y label
           ylab(paste("K =", i)) +
           # change axis elements, remove legend and background panel
           theme(axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 16),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 legend.position = "none",
                 panel.background = element_blank()))
}

# render multiple plots in a single page (colors set by trial and error)
multiplot(p2 + scale_fill_manual(values = c("#56B4E9", "#F0E442")),
          p3 + scale_fill_manual(values = c("#009E73", "#56B4E9", "#F0E442")),
          p4 + scale_fill_manual(values = c("#009E73", "#56B4E9", "#F0E442", "#CC79A7")),
          p5 + scale_fill_manual(values = c("#E69F00", "#F0E442", "#009E73", "#56B4E9", "#CC79A7")),
          p6 + scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "#006046", "#CC79A7", "#F0E442")),
          p7 + scale_fill_manual(values = c("#D55E00", "#006046", "#E69F00", "#CC79A7", "#F0E442", "#009E73", "#56B4E9")),
          p8 + scale_fill_manual(values = c("#009E73", "#E69F00", "#F0E442", "#D55E00", "#56B4E9", "#CC79A7", "#006046", "#0072B2")),
          p9 + scale_fill_manual(values = c("#666666", "#F0E442", "#0072B2", "#006046", "#E69F00", "#CC79A7", "#009E73", "#D55E00", "#56B4E9")),
          p1 + scale_fill_manual(values = c("#F0E442", "#CC79A7", "#873b00", "#E69F00", "#0072B2", "#56B4E9", "#006046", "#D55E00", "#009E73", "#666666")) +
          # fix y label
          ylab("K = 10")
         )
