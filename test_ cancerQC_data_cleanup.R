# Test data cleanup for CancerQC_v2 shiny app

library(dplyr)
library(data.table)
library(ggplot2)

#### Load and prepare objects

# Tumor sample groups
tumor <- c("FF", "FFPE")

# Load and clean QC metrics data (NEEDS MANUAL UPDATE)
QC <- read.csv("./Data/ready_to_plot.05-02-17.csv")
QC <- QC[!duplicated(QC),]  # remove exact duplicates
QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),] # WARNING: this table has also WELL_ID duplicates where second entries are empty
QC$COLLECTING_DATE <- as.Date(QC$COLLECTING_DATE, format = "%Y-%m-%d")
QC_tumor <- QC %>% filter(GROUP %in% tumor)
QC_tumor$GROUP <- as.character(QC_tumor$GROUP)
QC_tumor$TUMOUR_TYPE <- as.character(QC_tumor$TUMOUR_TYPE)
QC_tumor[QC_tumor$TUMOUR_TYPE == "",]$TUMOUR_TYPE <- "Unknown"
QC_tumor[QC_tumor$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- "Unknown"
QC_tumor[QC_tumor$TUMOUR_TYPE == "Testicular Germ Cell Tumours",]$TUMOUR_TYPE <- "Testicular"
QC_tumor[QC_tumor$TUMOUR_TYPE == "Malignant Melanoma",]$TUMOUR_TYPE <- "Melanoma"
QC_tumor[QC_tumor$TUMOUR_TYPE == "Endometrial Carcinoma",]$TUMOUR_TYPE <- "Endometrial"
QC_tumor[QC_tumor$TUMOUR_TYPE == "Adult Glioma",]$TUMOUR_TYPE <- "Glioma"

# Load GMC codes
GMCs <- read.csv("./Data/GMCs.csv")
GMCs <- GMCs %>% select(CODE, GMC)
GMCs$CODE <- as.character(GMCs$CODE)

# Save test data
save(GMCs, QC, QC_tumor, file = "test_QC_data.RData")

