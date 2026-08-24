# 04_marker_genes

library(Seurat)
library(SpatialExperiment)
library(spatialLIBD)
library(dplyr)
library(pheatmap)

set.seed(42)

dLPFC_seurat <- readRDS(
  "data/processed/dLPFC_seurat_processed.rds"
)

spe151673 <- readRDS(
  "data/processed/spe151673_clustered.rds"
)

stopifnot(
  identical(
    colnames(dLPFC_seurat),
    colnames(spe151673)
  )
)

table(dLPFC_seurat$seurat_clusters)

#01_Find each cluster feature gene
DefaultAssay(dLPFC_seurat) <- "SCT"
dLPFC_seurat <- PrepSCTFindMarkers(dLPFC_seurat, verbose = TRUE)

clusters <- levels(Idents(dLPFC_seurat))

markers_list <- lapply(clusters, function(cl) {
  df <- FindMarkers(
    dLPFC_seurat,
    ident.1         = cl,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    recorrect_umi   = FALSE
  )
  df$gene    <- rownames(df)
  df$cluster <- cl
  df
})

markers <- do.call(rbind, markers_list)
markers <- subset(markers, avg_log2FC > 0)
rownames(markers) <- NULL
write.csv(markers, "results/tables/04_cluster_markers.csv", row.names = FALSE)



# Cluster labels are assigned from marker-gene profiles and layer-signature
# scores only. Ground-truth anatomical annotations remain sealed until script 05.
#02_marker panel + marker existence check
layer_markers <- list(
  L1 = c("AQP4", "RELN", "FABP7"),
  L2 = c("HPCAL1", "LAMP5"),
  L3 = c("ADCYAP1", "CARTPT", "FREM3"),
  L4 = c("RORB", "PVALB"),
  L5 = c("PCP4", "TRABD2A"),
  L6 = c("KRT17", "NTNG2"),
  WM = c("MOBP", "MBP", "PLP1")
)
present <- lapply(layer_markers, function(g) g[g %in% rownames(dLPFC_seurat)])
missing <- setdiff(unlist(layer_markers), rownames(dLPFC_seurat))

sapply(present, length)
missing

#03_Calculate each spot layer signature score

dLPFC_seurat <- AddModuleScore(
  dLPFC_seurat,
  features = present,
  name     = "sig_",
  ctrl     = 50,
  seed.    = 42
)

score_cols <- paste0("sig_", seq_along(present))
colnames(dLPFC_seurat@meta.data)[
  match(score_cols, colnames(dLPFC_seurat@meta.data))
] <- names(present)

sig_mat <- sapply(names(present), function(layer){
  tapply(dLPFC_seurat@meta.data[[layer]],
         dLPFC_seurat$seurat_clusters,
         mean)
})
  
round(sig_mat, 3)


#04_Standardization + Preliminary Judgment /Signature-based annotation
sig_scaled <- scale(sig_mat)
auto_call <- colnames(sig_scaled)[apply(sig_scaled, 1, which.max)]
names(auto_call) <- rownames(sig_scaled)

auto_call

margin <- apply(sig_scaled, 1, function(x){
  s <- sort(x, decreasing = TRUE); s[1] -s[2]
})
round(margin, 3)


# Save annotation evidence
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  round(sig_mat, 4),
  "results/tables/04_layer_signature_scores.csv"
)

write.csv(
  data.frame(
    cluster   = names(auto_call),
    auto_call = unname(auto_call),
    margin    = unname(margin)
  ),
  "results/tables/04_signature_annotation_summary.csv",
  row.names = FALSE
)

#05_Marker_informed cluster annotation
cluster_annotation <- c(
  "0" = "L2/L3",
  "1" = "WM",
  "2" = "L6",
  "3" = "L5",
  "4" = "Non-laminar",
  "5" = "L1"
)
stopifnot(all(nzchar(cluster_annotation)))
stopifnot(setequal(names(cluster_annotation),
                   levels(dLPFC_seurat$seurat_clusters)))

lp <- cluster_annotation[as.character(dLPFC_seurat$seurat_clusters)]
names(lp) <- colnames(dLPFC_seurat)

dLPFC_seurat$layer_pred <- factor(
  lp,
  levels = c("L1", "L2/L3", "L5", "L6", "WM", "Non-laminar")
)

table(dLPFC_seurat$layer_pred)

stopifnot(!any(is.na(dLPFC_seurat$layer_pred)))
table(dLPFC_seurat$layer_pred)

stopifnot(identical(colnames(spe151673), colnames(dLPFC_seurat)))
spe151673$layer_pred <- dLPFC_seurat$layer_pred


#06_Spatial visualization
#6.1 Overall: The annotated stratified diagram of the space
lc <- spatialLIBD::libd_layer_colors

pred_cols <- c(
  "L1"          = unname(lc["Layer1"]),   # #F0027F 洋红
  "L2/L3"       = unname(lc["Layer3"]),   # #4DAF4A 绿
  "L5"          = unname(lc["Layer5"]),   # #FFD700 金黄
  "L6"          = unname(lc["Layer6"]),   # #FF7F00 橙
  "WM"          = unname(lc["WM"]),       # #1A1A1A 近黑
  "Non-laminar" = "#BDBDBD"               # 灰
)
stopifnot(!any(is.na(pred_cols)),
          setequal(names(pred_cols), levels(spe151673$layer_pred)))

p_clus <- vis_clus(spe = spe151673, sampleid = "151673",
                   clustervar = "layer_pred",
                   colors = pred_cols, point_size = 1.5) +
  ggplot2::ggtitle("Predicted — expression-only clustering(marker annotated")

png("results/figures/04_spatial_layer_pred.png",
    width = 1600, height = 1600, res = 200)
print(p_clus); dev.off()
p_clus

#6.2 Single-gene marker spatial map
marker_genes <- c(
  "AQP4",
  "CUX2",
  "PCP4",
  "KRT17",
  "MBP",
  "COL1A1"
)

stopifnot(all(marker_genes %in% rownames(dLPFC_seurat)))

for(g in marker_genes){
  colData(spe151673)[[g]] <- as.numeric(
    GetAssayData(dLPFC_seurat, assay = "SCT", layer = "data")[g, ]
  )
  p <- vis_gene(spe = spe151673, sampleid = "151673",
                geneid = g, point_size = 1.5)
  png(sprintf("results/figures/04_marker_%s.png",g), width = 1600, height = 1600, res = 200)
  print(p)
  dev.off()
}

#07: save

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

names(markers_list) <- levels(dLPFC_seurat$seurat_clusters)
stopifnot(!is.null(names(markers_list)),
          length(markers_list) == nlevels(dLPFC_seurat$seurat_clusters))

cat("cluster 1 (WM) top genes:\n");
print(head(rownames(markers_list[["1"]]), 8))
cat("cluster 4 (Non-laminar) top genes:\n")
print(head(rownames(markers_list[["4"]]), 8))



all_markers <- do.call(rbind, lapply(names(markers_list), function(cl) {
  df <- markers_list[[cl]]
  df$cluster <- cl
  df$gene    <- rownames(df)
  df
}))

stopifnot(nrow(all_markers) > 0)
print(table(all_markers$cluster))
cat("total marker rows:", nrow(all_markers), "\n")   # 预期 9067

write.csv(all_markers, "results/tables/04_all_markers.csv", row.names = FALSE)
saveRDS(markers_list, "data/processed/04_markers_list.rds")

gene_cols <- intersect(
  c("AQP4", "CUX2", "PCP4", "KRT17", "MBP", "COL1A1",
    "CLDN5", "FLT1", "PECAM1"),
  colnames(colData(spe151673))
)
if (length(gene_cols)) colData(spe151673)[gene_cols] <- NULL


stopifnot(!any(is.na(dLPFC_seurat$layer_pred)),
          identical(colnames(spe151673), colnames(dLPFC_seurat)))

saveRDS(dLPFC_seurat, "data/processed/dLPFC_seurat_04_annotated.rds")
saveRDS(spe151673,    "data/processed/spe151673_04_annotated.rds")


write.csv(
  data.frame(
    barcode    = colnames(dLPFC_seurat),
    cluster    = as.character(dLPFC_seurat$seurat_clusters),
    layer_pred = as.character(dLPFC_seurat$layer_pred)
  ),
  "results/tables/04_spot_annotations.csv",
  row.names = FALSE
)


write.csv(
  data.frame(cluster    = names(cluster_annotation),
             prediction = unname(cluster_annotation)),
  "results/tables/04_layer_predictions.csv",
  row.names = FALSE
)

# 7.6 sessionInfo
writeLines(capture.output(sessionInfo()), "results/tables/04_sessionInfo.txt")

# 7.7 
out <- c(list.files("data/processed", pattern = "^04|_04_", full.names = TRUE),
         list.files("results/tables",  pattern = "^04",      full.names = TRUE))
print(data.frame(file = basename(out), bytes = file.size(out)))