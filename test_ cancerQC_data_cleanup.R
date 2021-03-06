# Test data cleanup for CancerQC_v2 shiny app

library(dplyr)

#### Load and prepare objects

# Tumor sample groups
tumor <- c("FF", "FFPE")

# Load and clean QC metrics data
QC <- read.csv("./Data/ready_to_plot.05-02-17.csv")
QC <- QC[!duplicated(QC),]  # remove exact duplicates
QC <- QC[!duplicated(QC$WELL_ID, fromLast = T),]
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


# Create annonymized GMC variable
QC_tumor$CENTER_annon <- QC_tumor$CENTER
levels(QC_tumor$CENTER_annon)
levels(QC_tumor$CENTER_annon) <- paste0("GMC", 1:length(levels(QC_tumor$CENTER)))
levels(QC_tumor$CENTER_annon)
table(QC_tumor$CENTER, QC_tumor$CENTER_annon)
# Remove unnecessary fields
QC_tumor <- QC_tumor %>% select(-(WELL_ID), -(CENTER), -(DIVERSITY), -(GbQ30NoDupsNoClip), -(perc_bases_ge_15x_mapQ_ge11))
# Save test data
save(QC_tumor, file = "test_QC_data.RData")

# Prepare a dummy QC data table
load("test_QC_data.RData")
# Remove unnecessary fields
QC_table_dummy <- QC_tumor %>% select(-(WELL_ID), -(CENTER), -(DIVERSITY), -(GbQ30NoDupsNoClip), -(perc_bases_ge_15x_mapQ_ge11))
# Save new table as test data
save(QC_table_dummy, file = "test_QC_data.RData")
# Remove some GMCs
QC_table_dummy <- QC_table_dummy %>% filter(!CENTER_annon %in% c("GMC9", "GMC10"))
# Write the dummy table
write.csv(QC_table_dummy, file = "./QC_dummy_table.csv", row.names = F)


