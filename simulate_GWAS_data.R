# simulate GWAS data in R
# https://web.stanford.edu/group/candes/knockoffs/tutorials/gwas_tutorial.html#creating_knockoffs

## make sure you have BioConductor:
# if (!requireNamespace("BiocManager", quietly = TRUE)) 
#   install.packages("BiocManager")
# BiocManager::install(version = "3.12")

# install the snpStats package from the Bioconductor repository
instal
source("https://bioconductor.org/biocLite.R")
biocLite("snpStats")
library(snpStats)

dat.cases = read.impute("data_sim/sim.out.cases.gen")
