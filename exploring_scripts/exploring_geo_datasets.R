#   Data plots for selected GEO samples

install.packages("BiocManager")
BiocManager::install("GEOquery")
install.packages("preprocessCore")
install.packages("limma")

library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

library(GEOquery)
gds <- getGEO("GSE50760")
eset <- gds[[1]]
metadata <- pData(eset) #gets the covariates
exprs <- exprs(eset)
featureNames(eset)

getGEOSuppFiles("GSE50760", makeDirectory = FALSE, baseDir = "sample_datasets/GEOSuppDownloadsRAW")
gunzip("sample_datasets/GEOSuppDownloads/GSE50760_RAW.tar", destname = "sample_datasets/GSE50760")

data <- read.table("sample_datasets/GSE62944_01_27_15_TCGA_20_420_Clinical_Variables_7706_Samples.txt", sep='\t',col.names=NA, quote=F)

clinicals<-t(read.delim("sample_datasets/GSE62944_01_27_15_TCGA_20_420_Clinical_Variables_7706_Samples.txt",sep='\t',header=1, row.names=1,check.names=F))

gse <- getGEO("GSE50760", GSEMatrix=TRUE)
exprs <- exprs(gse[[1]])
