# Spatial Transcriptomics Project
# 03: Seurat Preprocessing and Clustering  (Plan C)
#
# Design notes:
#   - spe151673_qc is the SINGLE SOURCE OF TRUTH. The Seurat object is a
#     derived view. Filtering happens only on the SPE; results flow back
#     to the SPE, never the other way around.
#   - Ground-truth layer labels are STRIPPED from the Seurat object here
#     and stored separately. They are not opened until script 05.
#   - Resolution is selected from a pre-specified resolution scan without
#     using ground-truth labels or external validation scores. The selected
#     resolution (0.47) is retained as the expression-based clustering baseline.

library(SpatialExperiment)
library(SingleCellExperiment)
library(scuttle)
library(spatialLIBD)
library(Seurat)
library(ggplot2)
library(scales)

set.seed(42)

spe151673_qc <- readRDS("data/processed/spe151673_qc.rds")
dim(spe151673_qc)



# 1. Quarantine the ground truth
# The manual annotations ship inside colData(). If they ride along into
# the Seurat object they are one table() away from contaminating every
# parameter choice made below. Extract, save, remove.

truth <- data.frame(
  barcode = colnames(spe151673_qc),
  layer_full = spe151673_qc$layer_guess_reordered,
  layer_short = spe151673_qc$layer_guess_reordered_short,
  stringsAsFactors = FALSE
)
saveRDS(truth, "data/processed/ground_truth_layers.rds")

spot_metadata <- as.data.frame(colData(spe151673_qc))

truth_cols <- grepl(
  "layer_guess|^spatialLIBD$|^Maynard$",
  colnames(spot_metadata),
  ignore.case = TRUE
)

cat("Removing ground-truth columns:",
    colnames(spot_metadata)[truth_cols], "\n")
spot_metadata <- spot_metadata[, !truth_cols, drop = FALSE]

colnames(spot_metadata)   # confirm nothing layer-related survived
stopifnot(!any(grepl("layer|spatialLIBD|Maynard|guess",
                     colnames(spot_metadata), ignore.case = TRUE)))


# 2. Ensembl IDs -> gene symbols
counts_matrix <- assay(spe151673_qc, "counts")

rownames(counts_matrix) <- scuttle::uniquifyFeatureNames(
  ID    = rownames(spe151673_qc),
  names = rowData(spe151673_qc)$gene_name
)

# canonical laminar markers must be findable for script 04
markers_needed <- c("SNAP25","MBP","PCP4","RORB","KRT17","MOBP")
missing_markers <- setdiff(markers_needed, rownames(counts_matrix))
if (length(missing_markers)) cat("MISSING:", missing_markers, "\n")
stopifnot(length(missing_markers) == 0)



# 3. Build the Seurat object (derived view)
dLPFC_seurat <- CreateSeuratObject(
  counts       = counts_matrix,
  assay        = "Spatial",
  meta.data    = spot_metadata,
  project      = "Human_DLPFC_151673",
  min.cells    = 0,
  min.features = 0
)

stopifnot(identical(colnames(dLPFC_seurat), colnames(spe151673_qc)))
stopifnot(!any(grepl("layer|spatialLIBD|guess",
                     colnames(dLPFC_seurat@meta.data), ignore.case = TRUE)))

dLPFC_seurat


# 4. SCTransform

dLPFC_seurat <- SCTransform(
  object              = dLPFC_seurat,
  assay               = "Spatial",
  new.assay.name      = "SCT",
  vst.flavor          = "v2",
  variable.features.n = 3000,
  verbose             = TRUE
)

DefaultAssay(dLPFC_seurat)
length(VariableFeatures(dLPFC_seurat))



# 5. PCA

dLPFC_seurat <- RunPCA(dLPFC_seurat, npcs = 50, verbose = FALSE)

p_elbow <- ElbowPlot(dLPFC_seurat, ndims = 50)
p_elbow
ggsave("results/figures/03_elbow_plot.png", p_elbow,
       width = 6, height = 4, dpi = 300)

# Set this AFTER looking at the elbow plot. Record the reason in README.
n_pcs <- 20



# 6. Resolution scan
# Scan a predefined range of Louvain resolutions without consulting
# ground-truth labels or external validation metrics (e.g., ARI/NMI).
# Resolution 0.47 is retained as the expression-based baseline used
# throughout downstream analyses.

dLPFC_seurat <- FindNeighbors(dLPFC_seurat, dims = 1:n_pcs, verbose = FALSE)

res_grid <- c(seq(0.1, 1.0, by = 0.1), seq(0.40, 0.50, by = 0.01))
res_grid <- sort(unique(res_grid))
res_scan <- data.frame(resolution = res_grid, n_clusters = NA_integer_)

for (i in seq_along(res_grid)) {
  tmp <- FindClusters(dLPFC_seurat,
                      resolution = res_grid[i],
                      verbose    = FALSE)
  res_scan$n_clusters[i] <- length(unique(Idents(tmp)))
}

print(res_scan)
write.csv(res_scan, "results/tables/03_resolution_scan.csv", row.names = FALSE)

chosen_res <- 0.47
cat("Chosen resolution:", chosen_res, "\n")



# 7. Final clustering + UMAP
dLPFC_seurat <- FindClusters(dLPFC_seurat,
                             resolution = chosen_res,
                             verbose    = FALSE)
dLPFC_seurat <- RunUMAP(dLPFC_seurat, dims = 1:n_pcs, verbose = FALSE)

table(dLPFC_seurat$seurat_clusters)

p_umap <- DimPlot(dLPFC_seurat, reduction = "umap", label = TRUE) +
  ggtitle(paste0("Louvain clusters (res = ", chosen_res, ")"))
p_umap
ggsave("results/figures/03_umap_clusters.png", p_umap,
       width = 6, height = 5, dpi = 300)


# 8. Write results back to the SPE  (Plan C)
# No VisiumV1 needed: the SPE already holds correct coordinates and the
# H&E image. Match by barcode, never by position.

stopifnot(identical(colnames(dLPFC_seurat), colnames(spe151673_qc)))

spe151673_qc$seurat_cluster <- factor(dLPFC_seurat$seurat_clusters)

n_clus <- nlevels(spe151673_qc$seurat_cluster)

p_spatial <- vis_clus(
  spe        = spe151673_qc,
  clustervar = "seurat_cluster",
  sampleid   = "151673",
  colors     = setNames(hue_pal()(n_clus),
                        levels(spe151673_qc$seurat_cluster)),
  ... = " Louvain clusters"
)
p_spatial
ggsave("results/figures/03_spatial_clusters.png", p_spatial,
       width = 6, height = 6, dpi = 300)


# 9. Save
saveRDS(dLPFC_seurat, "data/processed/dLPFC_seurat_processed.rds")
saveRDS(spe151673_qc, "data/processed/spe151673_clustered.rds")

cat("\n03 completed. resolution =", chosen_res,
    "| n_clusters =", n_clus,
    "| PCs =", n_pcs, "\n")
cat("Ground truth remains sealed until script 05.\n")

stopifnot(
  file.exists("data/processed/dLPFC_seurat_processed.rds"),
  file.exists("data/processed/spe151673_clustered.rds"),
  file.exists("data/processed/ground_truth_layers.rds")
)

writeLines(capture.output(sessionInfo()), "results/tables/03_sessionInfo.txt")

out <- c(list.files("data/processed", pattern = "^dLPFC_seurat_processed|^spe151673_clustered|^ground_truth", full.names = TRUE),
         list.files("results/tables",  pattern = "^03", full.names = TRUE),
         list.files("results/figures", pattern = "^03", full.names = TRUE))
print(data.frame(file = basename(out), bytes = file.size(out)))
stopifnot(all(file.size(out) > 0))