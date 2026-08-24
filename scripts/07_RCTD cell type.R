#Spatial Transcriptomics Project
#07: RCTD cell-type devonvolution

library(SpatialExperiment)
library(spatialLIBD)
library(spacexr)
library(SingleCellExperiment)
library(HDF5Array)
library(ggplot2)
library(Matrix)

#01 :Download spatial data 
spe_bs <- readRDS("data/processed/spe151673_06_bayesspace.rds")
dim(spe_bs)
table(spe_bs$bayesspace_label)
table(spe_bs$layer_pred)

#02 : Reference data : LIBD DLPFC snRNA-seq
sce_path_zip <- fetch_data("spatialDLPFC_snRNAseq")
sce_path <- unzip(sce_path_zip, exdir = tempdir())
sce <- HDF5Array::loadHDF5SummarizedExperiment(
  file.path(tempdir(), "sce_DLPFC_annotated"))
dim(sce)
table(sce$cellType_broad_hc, useNA = "ifany")

#03 : Get rid of Ambiguous
sce <- sce[, sce$cellType_broad_hc != "Ambiguous"]
sce$cellType_broad_hc <- droplevels(sce$cellType_broad_hc)
table(sce$cellType_broad_hc)

#04 : Check the genes name align
# Use Ensembl gene IDs for alignment with spatial data
stopifnot("gene_id" %in% colnames(rowData(sce)))

gene_id <- rowData(sce)$gene_id

stopifnot(
  length(gene_id) == nrow(sce),
  !anyNA(gene_id)
)

sum(duplicated(rowData(sce)$gene_id))
rownames(sce) <- gene_id

cat("spatial gene IDs:\n"); print(head(rownames(spe_bs), 3))
cat("reference gene IDs:\n"); print(head(rownames(sce), 3))

n_shared <- length(intersect(rownames(spe_bs), rownames(sce)))
cat(sprintf("shared genes: %d(spatial %d, reference %d)\n",n_shared,nrow(spe_bs), nrow(sce)))
stopifnot(n_shared > 0)
                   
#05: Build the Reference Object(snRNA-seq side)
set.seed(42)
cells_keep <- unlist(lapply(split(colnames(sce), sce$cellType_broad_hc), function(x) {
  sample(x, min(length(x), 2000))
}))
sce_sub <- sce[, cells_keep]
table(sce_sub$cellType_broad_hc)

ref_counts <- as(as.matrix(assay(sce_sub, "counts")), "dgCMatrix")
ref_celltype <- sce_sub$cellType_broad_hc
names(ref_celltype) <- colnames(sce_sub)
ref_nUMI <- colSums(ref_counts)

reference <- Reference(
  counts     = ref_counts,
  cell_types = ref_celltype,
  nUMI       = ref_nUMI
)

print(reference)
table(reference@cell_types)


# 06: Build the SpatialRNA object(Visium side)
sp_counts <- as(counts(spe_bs), "dgCMatrix")
class(sp_counts)
format(object.size(sp_counts), "MB")

sp_coords <- as.data.frame(spatialCoords(spe_bs))
colnames(sp_coords) <- c("x","y")
rownames(sp_coords) <- colnames(spe_bs)

sp_nUMI <- colSums(sp_counts)

query <- SpatialRNA(
  coords = sp_coords,
  counts = sp_counts,
  nUMI  = sp_nUMI
)

dim(query@counts)
head(query@coords)
identical(colnames(query@counts), rownames(query@coords))

#07: Create and run RCTD
myRCTD <- create.RCTD(
  spatialRNA = query,
  reference = reference,
  max_cores = 4,
  CELL_MIN_INSTANCE =25
)
myRCTD <- run.RCTD(myRCTD, doublet_mode = "full")

#08: Extract cell-type proportions
results <- myRCTD@results
prop_matrix <- normalize_weights(results$weights)

dim(prop_matrix)
head(prop_matrix)

stopifnot(all(abs(rowSums(prop_matrix) - 1) < 1e-6))

#Alignment check
stopifnot(identical(rownames(prop_matrix), colnames(spe_bs)))
pm <- as.matrix(prop_matrix)

print(round(apply(pm, 2, summary), 4))

#Aggregation
layer <- colData(spe_bs)$layer_guess_reordered_short
by_layer <- aggregate(pm, by = list(layer = layer), FUN = mean)
stopifnot(nrow(by_layer) == 7)

global_summary <- as.data.frame(t(apply(pm, 2, summary)))
global_summary$cell_type <- rownames(global_summary)
global_summary <- global_summary[, c("cell_type", setdiff(names(global_summary), "cell_type"))]
global_summary[, -1] <- round(global_summary[, -1], 4)
print(global_summary)

write.csv(global_summary, "results/tables/07_global_proportion_summary.csv",
          row.names = FALSE)
stopifnot(file.size("results/tables/07_global_proportion_summary.csv") > 0)
          

#09： Sanity check - Oligo should be high in WM, low elsewhere
spe_bs$prop_Oligo <- pm[, "Oligo"]
spe_bs$prop_Excit <- pm[, "Excit"]

print(by_layer[, c("layer", "Oligo", "Excit")], digits = 3)
wm <- by_layer$layer == "WM"
stopifnot(by_layer$Oligo[wm] == max(by_layer$Oligo),
          by_layer$Excit[wm] == min(by_layer$Excit))


dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
plots <- list()
for (ct in c("Oligo", "Excit")) {
  p <- vis_gene(spe = spe_bs, sampleid = "151673",
                geneid = paste0("prop_", ct),
                point_size = 1.5, spatial = TRUE)
  png(sprintf("results/figures/08_prop_%s.png", ct),
      width = 1400, height = 1400, res = 200)
  print(p); dev.off()
  plots[[ct]] <- p
}
stopifnot(all(file.size(sprintf("results/figures/08_prop_%s.png",
                                c("Oligo", "Excit"))) > 0))

plots$Oligo
plots$Excit
print(by_layer[, c("layer", "Oligo", "Excit")], digits = 3)

#10: Q1: What are the 388 Non-laminar spots?
prop_df <- as.data.frame(pm)             
prop_df$barcode    <- rownames(prop_df)
prop_df$layer_pred <- spe_bs$layer_pred   

nl   <- which(prop_df$layer_pred == "Non-laminar")
rest <- which(!prop_df$layer_pred %in% c("Non-laminar", "WM"))
stopifnot(length(nl) == 388,
          length(rest) > 0,
          !any(prop_df$layer_pred[rest] %in% c("Non-laminar", "WM")))

celltypes <- colnames(pm)
q1_table <- do.call(rbind, lapply(celltypes, function(ct){
  x   <- prop_df[[ct]]
  q75 <- quantile(x[c(nl, rest)], 0.75)
  data.frame(
    cell_type    = ct,
    mean_NL      = mean(x[nl]),
    mean_rest    = mean(x[rest]),
    median_NL    = median(x[nl]),
    median_rest  = median(x[rest]),
    frac_hi_NL   = mean(x[nl]   > q75),
    frac_hi_rest = mean(x[rest] > q75),
    ratio        = mean(x[nl]) / mean(x[rest])
  )
}))
q1_table[, -1] <- round(q1_table[, -1], 4)
q1_table <- q1_table[order(-q1_table$ratio), ]
print(q1_table)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(q1_table, "results/tables/07_nonlaminar_celltype_comparison_ctxonly.csv",
          row.names = FALSE)
stopifnot(file.size("results/tables/07_nonlaminar_celltype_comparison_ctxonly.csv") > 0)




#11. Q2: Is L4's cell-type composition distinct from L3/L5?
truth <- spe_bs$layer_guess_reordered_short
idx <- function(L) which(truth == L)          

q2_table <- do.call(rbind, lapply(celltypes, function(ct){
  x <- prop_df[[ct]]
  data.frame(cell_type = ct,
             L3 = mean(x[idx("L3")]),
             L4 = mean(x[idx("L4")]),
             L5 = mean(x[idx("L5")]))
}))
q2_table[, -1] <- round(q2_table[, -1], 4)
q2_table$min_gap <- round(pmin(abs(q2_table$L4 - q2_table$L3),
                               abs(q2_table$L4 - q2_table$L5)), 4)
stopifnot(length(idx("L4")) == 218)           
print(q2_table)

write.csv(q2_table, "results/tables/07_L4_celltype_composition.csv", row.names = FALSE)
stopifnot(file.size("results/tables/07_L4_celltype_composition.csv") > 0)

#12. Q3: Is the L6/WM boundary domain real?
dom <- spe_bs$bayesspace_label
stopifnot(!is.null(dom), length(dom) == nrow(prop_df))

q3_table <- data.frame(
  domain     = names(table(dom)),
  n          = as.vector(table(dom)),
  mean_Oligo = round(as.vector(tapply(prop_df$Oligo, dom, mean)), 3),
  mean_Excit = round(as.vector(tapply(prop_df$Excit, dom, mean)), 3)
)
q3_table <- q3_table[order(q3_table$mean_Oligo), ]
print(q3_table)

write.csv(q3_table, "results/tables/07_bayesspace_domain_composition.csv",
          row.names = FALSE)
stopifnot(file.size("results/tables/07_bayesspace_domain_composition.csv") > 0)

#13. Save everything

saveRDS(myRCTD, "data/processed/07_rctd_full.rds")
write.csv(as.matrix(prop_matrix), "results/tables/07_rctd_proportions.csv")
write.csv(by_layer, "results/tables/07_prop_by_true_layer.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/tables/07_sessionInfo.txt")

for (f in list.files("results/tables", pattern = "^07_", full.names = TRUE))
  stopifnot(file.size(f) > 0)
