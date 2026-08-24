# Spatial Transcriptomics Project
#02: Quality Control(QC)

library(SpatialExperiment)
library(ggplot2)
library(spatialLIBD)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

spe151673 <- readRDS(
  "data/processed/spe151673_raw.rds"
)
imgData(spe151673)
head(spatialCoords(spe151673))
head(colData(spe151673))
head(rowData(spe151673))

ggplot(
  as.data.frame(spatialCoords(spe151673)),
  aes(
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres
  )
) +
  geom_point(size = 0.6) +
  coord_fixed() +
  scale_y_reverse() +
  theme_classic()

# QC summary 
summary(colData(spe151673)$sum_umi)
summary(colData(spe151673)$sum_gene)
summary(colData(spe151673)$expr_chrM_ratio)

#UMI distribution
ggplot(as.data.frame(colData(spe151673)),
       aes(sum_umi)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "UMI per Spot")

#Gene distribution
ggplot(as.data.frame(colData(spe151673)),
       aes(sum_gene)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Genes per Spot")

#Mito ratio
ggplot(as.data.frame(colData(spe151673)),
       aes(expr_chrM_ratio)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Mitochondrial Ratio")

# UMI

plot(
  spatialCoords(spe151673),
  col = colorRampPalette(c("blue", "yellow", "red"))(100)[
    cut(
      colData(spe151673)$sum_umi,
      breaks  = 100
    )
  ],
  pch = 16,
  cex = 0.7,
  asp = 1
)

title("UMI Counts")

# Gene

plot(
  spatialCoords(spe151673),
  col = colorRampPalette(c("blue", "yellow", "red"))(100)[
    cut(
      colData(spe151673)$sum_gene,
      breaks  = 100
    )
  ],
  pch = 16,
  cex = 0.7,
  asp = 1
)

title("Detected genes")

#Mito
plot(
  spatialCoords(spe151673),
  col = colorRampPalette(c("blue","yellow","red"))(100)[
    cut(
      colData(spe151673)$expr_chrM_ratio,
      breaks = 100
    )
  ],
  pch = 16,
  cex = 0.7,
  asp = 1
)

title("Mito ratio")

#QC Filtering
par(mfrow = c(1, 3))

boxplot(colData(spe151673)$sum_umi,
        main = "UMI")

boxplot(colData(spe151673)$sum_gene,
        main = "Genes")

boxplot(colData(spe151673)$expr_chrM_ratio,
        main = "Mito ratio")
par(mfrow = c(1, 1))
# UMI vs Gene

plot(
  colData(spe151673)$sum_umi,
  colData(spe151673)$sum_gene,
  pch = 16,
  cex = 0.5,
  col = rgb(0,0,0,0.3),
  xlab = "Total UMI",
  ylab = "Detected Genes",
  main = "QC: UMI vs Gene"
)

plot(
  colData(spe151673)$sum_umi,
  colData(spe151673)$expr_chrM_ratio,
  pch = 16,
  cex = 0.5,
  col = rgb(1,0,0,0.3),
  xlab = "Total UMI",
  ylab = "Mito Ratio",
  main = "QC: UMI vs Mito Ratio"
)
# Save: UMI vs mito ratio (evidence for skipping the mito filter)
png("results/figures/02_qc_umi_vs_mito.png",
    width = 6, height = 5, units = "in", res = 300)
plot(
  colData(spe151673)$sum_umi,
  colData(spe151673)$expr_chrM_ratio,
  pch = 16, cex = 0.5, col = rgb(1, 0, 0, 0.3),
  xlab = "Total UMI", ylab = "Mito Ratio",
  main = "QC: UMI vs Mito Ratio"
)
abline(h = 0.20, lty = 2, col = "blue")
dev.off()






cd <- colData(spe151673)
cat("Spots failing each criterion: \n")
cat("UMI <= 1000 :", sum(cd$sum_umi <= 1000), "\n")
cat("Gene <= 500 :", sum(cd$sum_gene <= 500), "\n")
cat("Mito >=0.2 :", sum(cd$expr_chrM_ratio >= 0.2), "\n")
summary(cd$expr_chrM_ratio)


cd_df <- as.data.frame(colData(spe151673))

p1 <- ggplot(cd_df, aes(expr_chrM_ratio)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = 0.20, linetype = "dashed", color = "red") +
  labs(title = "Mitochondrial ratio: unimodal, no damaged-spot population", subtitle = "Dashed line = conventional 0.20 threshold; not applied", x = NULL) +
  theme_classic()

ggsave("results/figures/02_mito_ratio_distribution.png", p1, width = 6, height = 4, dpi = 300)

# mito ratio is unimodal (median 0.167, max 0.330) with no separate
# damaged-spot population, so the conventional 0.20 threshold would cut
# the tail of a single distribution rather than remove damaged spots.
# No mito filter applied. Retrospective layer-bias check in 05.

keep <- colData(spe151673)$sum_umi > 1000 &
  colData(spe151673)$sum_gene > 500 
  

table(keep)
spe151673_qc <- spe151673[, keep]
spe151673_qc <- spe151673_qc[rowSums(counts(spe151673_qc)) > 0,]
dim(spe151673_qc)
saveRDS(
  spe151673_qc,
  file = "data/processed/spe151673_qc.rds"
)

stopifnot(file.exists("data/processed/spe151673_qc.rds"))

writeLines(
  capture.output(sessionInfo()),
  con = "results/tables/02_sessionInfo.txt"
)
cat("\nQC completed. Filtered object saved.\n")