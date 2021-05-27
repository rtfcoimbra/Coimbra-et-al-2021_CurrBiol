library(RZooROH)

# set working directory
setwd("")

# read input file of genotype probabilities
data <- zoodata(genofile = "snps.filtered.gen",
                zformat = "gp",
                samplefile = "snps.filtered.samples")

# define a model with pre-defined rates for 7 classes (6 HBD and 1 non-HBD)
model <- zoomodel(predefined = TRUE,
                  K = 7,
                  krates = c(4, 16, 64, 320, 1280, 5120, 5120))

# estimate the parameters of the model, the global and local realized autozygosity,
# partition it in the different HBD classes, and identify the HBD segments
results <- zoorun(model, data)
