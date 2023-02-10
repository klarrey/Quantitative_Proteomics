#Install bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

#Install DEP Package
BiocManager::install("DEP")

#load library
library("DEP")


setwd("C:/Users/3077080/OneDrive - University of Arkansas for Medical Sciences/Qproteomics/dep")

head(DiUbi)
dim(DiUbi)

#############################################################
#add_rejections Mark significant proteins
#Examples

# Load example
data <- UbiLength
data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Make SummarizedExperiment
columns <- grep("LFQ.", colnames(data_unique))
exp_design <- UbiLength_ExpDesign
se <- make_se(data_unique, columns, exp_design)
head(assays(se)[[1]])

# Filter, normalize and impute missing values
filt <- filter_missval(se, thr = 0)
norm <- normalize_vsn(filt)
norm.2 <- as.matrix(assays(filt)[[1]])        #Gives better normalization for this data
imputed <- impute(norm, fun = "MinProb", q = 0.01)

# Test for differentially expressed proteins
diff <- test_diff(imputed, "control", "Ctrl")
dep <- add_rejections(diff, alpha = 0.05, lfc = 1)
