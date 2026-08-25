# 06: Spatial-aware clustering benchmark
# Three Arm design：
#   A  SCT     + PCA + Louvain      
#   B  logNorm + PCA + Louvain      
#   C  logNorm + PCA + BayesSpace  

library(SpatialExperiment)
library(Seurat)
library(BayesSpace)
library(mclust)
library(aricode)
library(ggplot2)

#load+ Reconstruction evaluation vector
spe151673 <- readRDS("data/processed/spe151673_04_annotated.rds")
dLPFC_seurat <- readRDS("data/processed/dLPFC_seurat_04_annotated.rds")

truth <- spe151673$layer_guess_reordered_short
keep <- !is.na(truth)
truth_eval <- droplevels(factor(truth[keep]))
stopifnot(
  length(truth_eval) ==sum(keep))

#Arm A: 
metrics_A <- read.csv("results/tables/05_baseline_metrics.csv")
print(metrics_A)

# Arm B: logNorm + PCA + Louvain, res = 0.60 (7 clusters)
set.seed(42)
DefaultAssay(dLPFC_seurat) <- "Spatial"
Assays(dLPFC_seurat)

dLPFC_logn <- NormalizeData(dLPFC_seurat, verbose = FALSE)
dLPFC_logn <- FindVariableFeatures(dLPFC_logn, nfeatures = 3000, verbose = FALSE)
dLPFC_logn <- ScaleData(dLPFC_logn, verbose = FALSE)
dLPFC_logn <- RunPCA(dLPFC_logn, npcs = 50, verbose = FALSE)
dLPFC_logn <- FindNeighbors(dLPFC_logn, dims = 1:20, verbose = FALSE)


length(VariableFeatures(dLPFC_logn))   


scan_B <- do.call(rbind, lapply(seq(0.3, 1.2, by = 0.1), function(r) {
  so <- FindClusters(dLPFC_logn, resolution = r, verbose = FALSE)
  data.frame(resolution = r,
             n_clusters = nlevels(droplevels(so$seurat_clusters)))
}))
print(scan_B)
write.csv(scan_B, "results/tables/06_resolution_scan_lognorm.csv", row.names = FALSE)
stopifnot(!6 %in% scan_B$n_clusters)  

#Arm B : logNorm + PCA + Louvain, res = 0.6
dLPFC_logn <- FindClusters(dLPFC_logn, resolution = 0.6, verbose = FALSE)

common_B <- intersect(
  colnames(dLPFC_logn),
  colnames(spe151673)[keep]
)

cl_B <- droplevels(factor(
  dLPFC_logn$seurat_clusters[common_B]
))

truth_B <- droplevels(factor(
  truth[match(common_B, colnames(spe151673))]
))

stopifnot(
  length(cl_B) == length(truth_B),
  !anyNA(cl_B),
  !anyNA(truth_B)
)

ari_B <- mclust::adjustedRandIndex(cl_B, truth_B)
nmi_B <- aricode::NMI(cl_B, truth_B)

cat(sprintf(
  "Arm B - logNorm res 0.6: %d clusters | n = %d | ARI = %.4f | NMI = %.4f\n",
  nlevels(cl_B), length(cl_B), ari_B, nmi_B
))

#ArmC : logNorm + PCA + BayesSpace
set.seed(42)

head(colData(spe151673)[, c("array_row", "array_col")])

spe_bs <- spe151673
colData(spe_bs)$row <- colData(spe_bs)$array_row
colData(spe_bs)$col <- colData(spe_bs)$array_col

spe_bs <- spatialPreprocess(spe_bs, platform = "Visium",
                            n.PCs = 20, n.HVGs = 3000, log.normalize = TRUE)
#spatial clustering, q = 7
spe_bs <- spatialCluster(spe_bs, q = 7, platform = "Visium",
                         d = 20, init.method = "mclust",
                         model = "t", gamma = 3,
                         nrep = 50000, burn.in = 1000,
                         save.chain = FALSE)
cl_C <- droplevels(factor(spe_bs$spatial.cluster[keep]))
ari_C <- mclust::adjustedRandIndex(cl_C, truth_eval)
nmi_C <- aricode::NMI(cl_C, truth_eval)

cat(sprintf("Arm C - BayesSpace q = 7:%d clusters | ARI = %.4f | NMI = %.4f\n",
            nlevels(cl_C), ari_C, nmi_C))

#Arm C confusion matrix
cm_C <- table(Predicted = cl_C, Truth = truth_eval)
print(cm_C)
round(prop.table(cm_C, margin = 2), 3)

#Arm C : Cluster -> layer tag
bs_label <- c("1" = "L6", "2" = "L6/WM boundary", "3" = "L3",
              "4" = "L2", "5" = "L5", "6" = "L1", "7" = "WM")
lab_C <- factor(bs_label[as.character(cl_C)],
                levels = c("L1", "L2", "L3", "L5", "L6", "L6/WM boundary", "WM"))
stopifnot(!any(is.na(lab_C)))
cm_C2 <- table(Predicted = lab_C, Truth = truth_eval)
print(cm_C2)

#recall
recall_C <- data.frame(
  truth_layer = c("L1", "L2", "L3", "L4", "L5", "L6", "WM"),
  n_truth = as.integer(colSums(cm_C2)),
  n_correct = as.integer(c(cm_C2["L1", "L1"], cm_C2["L2", "L2"], cm_C2["L3", "L3"],
                           0, cm_C2["L5", "L5"], cm_C2["L6", "L6"], cm_C2["WM", "WM"])),
  resolved = c("Yes","Yes","Yes", "No(not recovered", "Yes",
              "Yes (plus boundary domain)","Yes(split with boundary domain)")
  
)
recall_C$recall <- round(recall_C$n_correct / recall_C$n_truth, 3)
print(recall_C)

write.csv(cm_C2, "results/tables/06_confusion_matrix_armC.csv")
write.csv(recall_C, "results/tables/06_per_layer_recall_armC.csv", row.names = FALSE)

summary_arms <- data.frame(
  arm        = c("A", "B", "C"),
  norm       = c("SCTransform", "logNormalize", "logNormalize"),
  method     = c("Louvain", "Louvain", "BayesSpace"),
  param      = c("res 0.47", "res 0.60", "q = 7"),
  n_clusters = c(
    metrics_A$n_clusters[1],
    nlevels(cl_B),
    nlevels(cl_C)
  ),
  ARI = c(
    metrics_A$ARI[1],
    round(ari_B, 4),
    round(ari_C, 4)
  ),
  NMI = c(
    metrics_A$NMI[1],
    round(nmi_B, 4),
    round(nmi_C, 4)
  )
)
print(summary_arms)
write.csv(summary_arms, "results/tables/06_arm_summary.csv", row.names = FALSE)


cm_df_C <- as.data.frame(cm_C2)
cm_df_C$Prop <- as.vector(prop.table(cm_C2, margin = 2))

p_cmC <- ggplot(cm_df_C, aes(Truth, Predicted, fill = Prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#B2182B",
                      name = "Column\nproportion") +
  scale_y_discrete(limits = rev(levels(cm_df_C$Predicted))) +
  labs(title = "151673 — BayesSpace (spatial prior) vs manual annotation",
       subtitle = sprintf("ARI = %.3f | NMI = %.3f | n =  %d | q = 7",
                          ari_C, nmi_C, length(truth_eval)),
       x = "Manual annotation (ground truth)",
       y = "Predicted (BayesSpace, q = 7)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

png("results/figures/06_confusion_matrix_armC.png",
    width = 1800, height = 1400, res = 200)
print(p_cmC); dev.off()

p_cmC

#Arm C spatial graph 
lab_all <- factor(bs_label[as.character(spe_bs$spatial.cluster)],
                  levels = c("L1","L2","L3","L5","L6","L6/WM boundary","WM"))
spe_bs$bayesspace_label <- lab_all
stopifnot(!any(is.na(lab_all)))

bs_cols <- c("L1"="#F7F7A0", "L2"="#F1A340", "L3"="#E7298A",
             "L5"="#66C2A5", "L6"="#3288BD",
             "L6/WM boundary"="#B15928", "WM"="#525252")

p_spC <- vis_clus(spe = spe_bs, sampleid = "151673",
                  clustervar = "bayesspace_label",
                  colors = bs_cols, point_size = 1.5)

png("results/figures/06_spatial_armC.png",
    width = 1600, height = 1600, res = 200)
print(p_spC); dev.off()

p_spC

saveRDS(spe_bs, "data/processed/spe151673_06_bayesspace.rds")

write.csv(
  data.frame(barcode    = colnames(spe_bs),
             bs_cluster = spe_bs$spatial.cluster,
             bs_label   = as.character(spe_bs$bayesspace_label)),
  "results/tables/06_spot_labels_armC.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "results/tables/06_sessionInfo.txt")

out <- c(list.files("results/tables",  pattern = "^06", full.names = TRUE),
         list.files("results/figures", pattern = "^06", full.names = TRUE))
print(data.frame(file = basename(out), bytes = file.size(out)))