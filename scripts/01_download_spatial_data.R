#Human DLPFC Spatial Transcriptomics Analysis
#Script: 01_download_spatial_data.R
#Purpose: Download and inspect the spatial transcriptomics data

library(spatialLIBD)
library(SpatialExperiment)
library(SummarizedExperiment)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

spe <- fetch_data(type = "spe")

#Inspect the SpatialExperiment Object
print(class(spe))
print(dim(spe))
print(assayNames(spe))
print(colnames(rowData(spe)))
print(colnames(colData(spe)))
print(head(spatialCoords(spe)))
print(imgData(spe))

#Select the representative DLPFC section:151673
stopifnot("151673" %in% spe$sample_id)
spe151673 <- spe[, spe$sample_id == "151673"]


#Verify the subset
print(dim(spe151673))
print(unique(spe151673$sample_id))
print(table(spe151673$sample_id))
print(class(spe151673))
print(imgData(spe151673))


#Save
saveRDS(
  spe151673,"data/processed/spe151673_raw.rds")
stopifnot(file.exists("data/processed/spe151673_raw.rds"))

writeLines(
  capture.output(sessionInfo()),
  con = "results/tables/01_sessionInfo.txt"
)
cat("\nSpatial data downloaded and section 151673 saved successfully.\n")