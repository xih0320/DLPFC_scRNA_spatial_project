#05：Baseline Benchmark - ARI /NMI vs Manual Annotation
#ground truth

library(SpatialExperiment)
library(Seurat)
library(mclust)
library(aricode)
library(tidyr)
library(spatialLIBD)
library(ggplot2)

spe151673 <- readRDS("data/processed/spe151673_04_annotated.rds")
dLPFC_seurat <-readRDS("data/processed/dLPFC_seurat_04_annotated.rds")

dim(spe151673)
table(spe151673$layer_pred)

stopifnot(identical(colnames(spe151673), colnames(dLPFC_seurat)))

#1.Locate the "ground truth"column
grep("layer|guess|ground", colnames(colData(spe151673)),
     value = TRUE, ignore.case = TRUE)

#2.Check three ground truth
table(spe151673$layer_guess, useNA = "ifany")
table(spe151673$layer_guess_reordered, useNA = "ifany")
table(spe151673$layer_guess_reordered_short, useNA = "ifany")

table(spe151673$layer_guess, spe151673$layer_guess_reordered)

#3.Construct paired vectors for evaluation
truth <- spe151673$layer_guess_reordered_short
pred <- spe151673$layer_pred

keep <- !is.na(truth)
sum(keep)
stopifnot(sum(keep) == 3580 - 27)

truth_eval <- droplevels(factor(truth[keep]))
pred_eval <- droplevels(factor(pred[keep]))

stopifnot(length(truth_eval) == length(pred_eval))
table(truth_eval)
table(pred_eval)
     
#04. ARI / NMI 
ari <- mclust::adjustedRandIndex(pred_eval, truth_eval)
nmi <- aricode::NMI(pred_eval, truth_eval)

cat(sprintf("ARI = %.4f\nNMI = %.4f\nn = %d\n",
            ari, nmi, length(truth_eval)))
     
     
#05.confusion matrix
cm <- table(Predicted = pred_eval, Truth = truth_eval)
print(cm)

round(prop.table(cm, margin = 2), 3)

#06. Save confusion matrix + ggplot
write.csv(as.data.frame.matrix(cm),
           "results/tables/05_confusion_matrix_counts.csv")
write.csv(round(as.data.frame.matrix(prop.table(cm, margin = 2)), 4),
          "results/tables/05_confusion_matrix_colprop.csv")
cm_df <- as.data.frame(cm)
cm_df$Prop <- as.vector(prop.table(cm,margin = 2))

p_cm <- ggplot(cm_df, aes(x = Truth, y = Predicted, fill = Prop)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 3.5) +
  scale_fill_gradient(low = "white", high = "#2166AC",
                      name = "Column\nproportion") +
  scale_y_discrete(limits = rev(levels(cm_df$Predicted))) +
  labs(title = "151673 - Expression - only clustering vs manual annotation",
       subtitle = sprintf("ARI = %.3f·NMI = %.3f·n = %d", ari, nmi, length(truth_eval)),
       x = "Manual annotation (ground truth)",
       y = "Predicted(res0.47, 6 clusters") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

png("results/figures/05_confusion_matrix.png",
    width = 1800, height = 1400, res = 200)
print(p_cm)
dev.off()

#07. Sensitivity check: Resolution VS ARI
names(dLPFC_seurat@graphs)
set.seed(42)
DefaultAssay(dLPFC_seurat) <- "SCT"
res_grid <- c(0.30, 0.40, 0.45, 0.47, 0.50, 0.60, 0.80)

sens <- do.call(rbind, lapply(res_grid, function(r){
  so <- FindClusters(dLPFC_seurat, graph.name = "SCT_snn",resolution = r, verbose = FALSE)
  cl <- droplevels(so$seurat_clusters[keep])
  data.frame(
    resolution = r,
    n_cluster = nlevels(cl),
    ARI = mclust::adjustedRandIndex(cl, truth_eval),
    NMI = aricode::NMI(cl, truth_eval)
  )
}))
print(sens)
write.csv(sens, "results/tables/05_resolution_sensitivity.csv", row.names = FALSE)


sens_long <- pivot_longer(sens, c(ARI, NMI), names_to = "metric", values_to = "value")

p_sens <- ggplot(sens_long, aes(resolution, value, colour = metric)) +
  geom_line() + geom_point(size = 2) +
  geom_vline(xintercept = 0.47, linetype = "dashed", colour = "grey40") +
  geom_text(data = sens, aes(resolution, 0.29, label = n_cluster),
            inherit.aes = FALSE, size = 3, colour = "grey30") +
  labs(title = "Resolution sensitivity (labels below = cluster count)",
       subtitle = "Dashed line: res 0.47, chosen in script 03 before unsealing annotations",
       x = "Louvain resolution", y = "Score") +
  theme_minimal(base_size = 11)
png("results/figures/05_resolution_sensitivity.png",
    width = 1600, height = 1000, res = 200)
print(p_sens); dev.off()

p_sens

#08.Each layer recall
recall_tbl <- data.frame(
  truth_layer = c("L1","L2","L3","L4","L5","L6","WM"),
  n_truth     = as.integer(colSums(cm)),
  n_correct   = as.integer(c(cm["L1","L1"], cm["L2/L3","L2"], cm["L2/L3","L3"], 0,
                             cm["L5","L5"], cm["L6","L6"], cm["WM","WM"])),
  resolved    = c("Yes","No (merged with L3)","No (merged with L2)",
                  "No (not recovered)","Yes","Yes","Yes")
)
recall_tbl$recall <- round(recall_tbl$n_correct / recall_tbl$n_truth, 3)
print(recall_tbl)

write.csv(recall_tbl, "results/tables/05_per_layer_recall.csv", row.names = FALSE)

keep_excl <- keep & pred != "Non-laminar"
ari_excl <- mclust::adjustedRandIndex(
              droplevels(factor(pred[keep_excl])),
              droplevels(factor(truth[keep_excl])))
cat(sprintf("ARI excl Non-laminar = %.4f (n = %d)\n", ari_excl, sum(keep_excl)))
stopifnot(sum(keep_excl) == 3192)


#09. Summary Indicators
metrics <- data.frame(
  method              = "Louvain (expression-only)",
  resolution          = 0.47,
  n_clusters          = nlevels(pred_eval),
  n_spots             = length(truth_eval),
  ARI                 = round(ari, 4),
  NMI                 = round(nmi, 4),
  ARI_excl_nonlaminar = round(ari_excl, 4),
  n_spots_excl        = sum(keep_excl)
)
write.csv(metrics, "results/tables/05_baseline_metrics.csv", row.names = FALSE)


lc <- spatialLIBD::libd_layer_colors
stopifnot(all(c("Layer1","Layer2","Layer3","Layer4",
                "Layer5","Layer6","WM") %in% names(lc)))
# Truth table
truth_cols <- c(
  "L1" = unname(lc["Layer1"]), "L2" = unname(lc["Layer2"]),
  "L3" = unname(lc["Layer3"]), "L4" = unname(lc["Layer4"]),
  "L5" = unname(lc["Layer5"]), "L6" = unname(lc["Layer6"]),
  "WM" = unname(lc["WM"])
)
stopifnot(!any(is.na(truth_cols)))

p_truth <- vis_clus(spe = spe151673, sampleid = "151673",
                    clustervar = "layer_guess_reordered_short",
                    colors = truth_cols, point_size = 1.5) +
  ggplot2::ggtitle("Ground truth — expert manual annotation")

png("results/figures/05_ground_truth_layers.png",
    width = 1600, height = 1600, res = 200)
print(p_truth); dev.off()
p_truth

# ---------------------------------------------------------------
# 09b. Retrospective QC check: was the mito >= 0.20 filter layer-biased?
#
# This section uses ground truth, so it cannot live in script 02.
# The mito criterion was dropped in 02 on truth-free evidence only:
# unimodal ratio distribution (median 0.167, max 0.330), no damaged-spot
# subpopulation, and removal tracking sequencing depth. What follows is a
# retrospective confirmation made AFTER unsealing -- it is reported as
# such, and was never the reason the filter was dropped.
#
# Denominator is the PRE-QC object: these are per-layer fractions of all
# spots before any filtering. Running this on the post-QC object would
# silently produce different numbers.
# ---------------------------------------------------------------

pre_qc_path <- "data/processed/spe151673_raw.rds"   # <- confirm the name written by 01
stopifnot(file.exists(pre_qc_path))
spe_pre <- readRDS(pre_qc_path)

cd_pre <- as.data.frame(colData(spe_pre))
stopifnot(all(c("expr_chrM_ratio", "layer_guess_reordered_short")
              %in% colnames(cd_pre)))

mito_thr <- 0.20

# Anchor assertion: catches a wrong object or a changed threshold before
# any number reaches the README.
n_flagged_all <- sum(cd_pre$expr_chrM_ratio >= mito_thr, na.rm = TRUE)
cat(
  "Spots with mitochondrial ratio >=", mito_thr,
  "that would have been removed:", n_flagged_all, "\n"
)

# mito alone flags 741. The 792 in earlier notes was the union of all three
# conventional criteria (mito 741, UMI<=1000 59, gene<=500 11 subset of UMI;
# overlap 8). 733 spots are mito-only.
stopifnot(n_flagged_all == 741)

ok        <- !is.na(cd_pre$layer_guess_reordered_short)
layer_pre <- droplevels(factor(cd_pre$layer_guess_reordered_short[ok]))
flagged   <- cd_pre$expr_chrM_ratio[ok] >= mito_thr

bias_tbl <- data.frame(
  layer     = levels(layer_pre),
  n_spots   = as.integer(table(layer_pre)),
  n_flagged = as.integer(table(layer_pre[flagged]))
)
bias_tbl$frac_flagged <- round(bias_tbl$n_flagged / bias_tbl$n_spots, 3)
print(bias_tbl)



# The finding itself: removal is strongly layer-dependent, not uniform.
stopifnot(bias_tbl$frac_flagged[bias_tbl$layer == "L1"] >
            bias_tbl$frac_flagged[bias_tbl$layer == "WM"])

write.csv(bias_tbl, "results/tables/05_qc_mito_bias_by_layer.csv",
          row.names = FALSE)

p_mito_bias <- ggplot(bias_tbl, aes(x = layer, y = frac_flagged)) +
  geom_col(fill = "grey35") +
  geom_text(aes(label = n_flagged), vjust = -0.4, size = 3, colour = "grey20") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = sprintf(
      "Spots with mitochondrial ratio >= %.2f that would have been removed",
      mito_thr
    ),
    subtitle = "Retrospective — annotations were sealed when this filter was dropped in 02",
    x = NULL,
    y = "Proportion of spots in layer"
  ) +
  theme_classic(base_size = 11)

ggsave("results/figures/05_qc_threshold_bias_retrospective.png", p_mito_bias,
       width = 6, height = 4, dpi = 300)





#10. sessionInfo
writeLines(capture.output(sessionInfo()), "results/tables/05_sessionInfo.txt")

out <- c(list.files("results/tables",  pattern = "^05", full.names = TRUE),
         list.files("results/figures", pattern = "^05", full.names = TRUE))
print(data.frame(file = basename(out), bytes = file.size(out)))




