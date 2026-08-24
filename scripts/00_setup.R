#Project:
#Integrated snRNA-seq and Spatial Transcriptomics Analysis 
#of the Human Dorsolateral Prefrontal Cortex

#1.Set project directory 

stopifnot(basename(getwd()) == "DLPFC_scRNA_spatial_project")

for(d in c(
  "data/raw",
  "data/processed",
  "results/figures",
  "results/tables"
)){
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


#2.Install CRAN packages

cran_packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "patchwork",
  "tidyr",
  "readr",
  "Matrix",
  "remotes",
  "mclust",
  "aricode",
  "pheatmap",
  "scales"
)

for (pkg in cran_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

#3.Install BiocManager

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

#4.Install Bioconductor packages

bioc_packages <- c(
  "spatialLIBD",
  "SpatialExperiment",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scuttle",
  "scran",
  "scater",
  "BayesSpace",
  "HDF5Array"
)
for (pkg in bioc_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

#4b. GitHub packages
if(!requireNamespace("pak", quietly = TRUE)) install.packages("pak", repos = "https://cloud.r-project.org")
if(!requireNamespace("spacexr", quietly = TRUE)) pak::pak("dmcable/spacexr")

#5. Load packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(spatialLIBD)
library(SingleCellExperiment)
library(scater)
library(scran)
library(scuttle)

#6.Create a package-version table
package_names <- c(cran_packages, bioc_packages, "spacexr")

package_versions <- data.frame(
  package = package_names,
  version = sapply(package_names, function(pkg)
    tryCatch(as.character(packageVersion(pkg)),
             error = function(e) NA_character_))
)
print(package_versions)

#7.Save package versions
write.csv(
  package_versions,
  file = "results/tables/package_versions.csv",
  row.names = FALSE
)

#8.Save the complete R environment information
writeLines(
  capture.output(sessionInfo()),
  con = "results/tables/session_info.txt"
)
cat("\nSetup completed successfully.\n")
