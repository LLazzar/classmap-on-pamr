install.packages("GEOquery")
library(GEOquery)

# set the GSE ID of the dataset you want to download
gse_id <- "GSE1234" # replace with the actual GSE ID

# download the dataset from GEO and parse it
gse <- getGEO(gse_id, GSEMatrix = TRUE)
