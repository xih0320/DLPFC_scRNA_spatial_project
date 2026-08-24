#1.Ground truth VS Expression-Only
# Final comparison figure using raw spatial coordinates
library(SpatialExperiment)
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)

spe151673 <- readRDS(
  "data/processed/spe151673_04_annotated.rds"
)
spe_bs <- readRDS(
  "data/processed/spe151673_07_bayesspace.rds"
)
# Build plotting data
plot_df <- as.data.frame(spatialCoords(spe151673))

plot_df$truth <- factor(
  spe151673$layer_guess_reordered_short,
  levels = c("L1", "L2", "L3", "L4", "L5", "L6", "WM")
)

plot_df$pred <- factor(
  spe151673$layer_pred,
  levels = c("L1", "L2/L3", "L5", "L6", "WM", "Non-laminar")
)


# Unified biological palette
truth_cols <- c(
  "L1" = "#F0027F",
  "L2" = "#A6D854",
  "L3" = "#4DAF4A",
  "L4" = "#377EB8",
  "L5" = "#FFD700",
  "L6" = "#FF7F00",
  "WM" = "#1A1A1A"
)

pred_cols <- c(
  "L1" = "#F0027F",
  "L2/L3" = "#4DAF4A",
  "L5" = "#FFD700",
  "L6" = "#FF7F00",
  "WM" = "#1A1A1A",
  "Non-laminar" = "#BDBDBD"
)


# ------------------------------------------------------------
# Ground truth
# ------------------------------------------------------------

p_truth <- ggplot(
  plot_df,
  aes(
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres,
    color = truth
  )
) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = truth_cols,
    na.value = "#E5E5E5",
    drop = FALSE
  ) +
  scale_y_reverse() +
  coord_fixed() +
  ggtitle("Ground truth — expert manual annotation") +
  labs(color = NULL) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    legend.text = element_text(size = 10)
  )

# Expression-only
p_pred <- ggplot(
  plot_df,
  aes(
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres,
    color = pred
  )
) +
  geom_point(size = 1.5) +
  scale_color_manual(
    values = pred_cols,
    drop = FALSE
  ) +
  scale_y_reverse() +
  coord_fixed() +
  ggtitle("Expression-only clustering — marker annotated") +
  labs(color = NULL) +
  theme_void() +
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    legend.text = element_text(size = 10)
  )

# Combine

p_compare <- p_truth + p_pred +
  plot_layout(ncol = 2)

p_compare


ggsave(
  "results/figures/09_truth_vs_expression_only.png",
  p_compare,
  width = 14,
  height = 7,
  dpi = 300
)

#2. Ground truth VS BayesSpace spatial-aware clustering
# Load BayesSpace result

spe_bs <- readRDS(
  "data/processed/spe151673_07_bayesspace.rds"
)

# Check labels
levels(factor(spe_bs$layer_guess_reordered_short))
levels(factor(spe_bs$bayesspace_label))

table(spe_bs$bayesspace_label, useNA = "ifany")

# Build plotting dataframe

plot_bs <- as.data.frame(spatialCoords(spe_bs))

# Ground truth:
# Keep the original 7-layer annotation
plot_bs$truth <- factor(
  spe_bs$layer_guess_reordered_short,
  levels = c(
    "L1", "L2", "L3", "L4",
    "L5", "L6", "WM"
  )
)

# BayesSpace:
# Keep the labels generated in script 07
plot_bs$bayes <- factor(
  spe_bs$bayesspace_label,
  levels = c(
    "L1",
    "L2",
    "L3",
    "L5",
    "L6",
    "L6/WM boundary",
    "WM"
  )
)

# Unified color palette
#
# IMPORTANT:
# Same biological layer = same color as Figure 1

truth_cols <- c(
  "L1" = "#F0027F",   # pink
  "L2" = "#A6D854",   # light green
  "L3" = "#4DAF4A",   # green
  "L4" = "#377EB8",   # blue
  "L5" = "#FFD700",   # yellow
  "L6" = "#FF7F00",   # orange
  "WM" = "#1A1A1A"    # black
)

bayes_cols <- c(
  "L1" = "#F0027F",
  "L2" = "#A6D854",
  "L3" = "#4DAF4A",
  
  # BayesSpace did not recover L4 as a distinct domain
  
  "L5" = "#FFD700",
  "L6" = "#FF7F00",
  
  # Unique spatial transition domain
  "L6/WM boundary" = "#984EA3",
  
  "WM" = "#1A1A1A"
)

# Ground truth

p_truth_bs <- ggplot(
  plot_bs,
  aes(
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres,
    color = truth
  )
) +
  geom_point(size = 1.5) +
  
  scale_color_manual(
    values = truth_cols,
    na.value = "#D9D9D9",
    drop = FALSE,
    na.translate = TRUE
  ) +
  
  scale_y_reverse() +
  coord_fixed() +
  
  ggtitle(
    "Ground truth — expert manual annotation"
  ) +
  
  labs(color = NULL) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text = element_text(
      size = 10
    )
  )


# BayesSpace

p_bayes <- ggplot(
  plot_bs,
  aes(
    x = pxl_col_in_fullres,
    y = pxl_row_in_fullres,
    color = bayes
  )
) +
  
  geom_point(size = 1.5) +
  
  scale_color_manual(
    values = bayes_cols,
    drop = FALSE
  ) +
  
  scale_y_reverse() +
  coord_fixed() +
  
  ggtitle(
    "BayesSpace — spatial-aware clustering"
  ) +
  
  labs(color = NULL) +
  
  theme_void() +
  
  theme(
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    legend.text = element_text(
      size = 10
    )
  )

# Combine

p_truth_vs_bayes <- p_truth_bs + p_bayes +
  plot_layout(
    ncol = 2,
    widths = c(1, 1)
  )

p_truth_vs_bayes

# Save

ggsave(
  "results/figures/09_truth_vs_bayesspace.png",
  p_truth_vs_bayes,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat(
  "\nGround truth vs BayesSpace figure completed.\n"
)

#3. A/B/C ARI +Benchmark
bench <- read.csv(
  "results/tables/07_arm_summary.csv",
  stringsAsFactors = FALSE
)
print(bench)

#Create informative method labels
bench$label <- c(
  "A\nSCT + Louvain",
  "B\nlogNorm + Louvain",
  "C\nlogNorm + BayesSpace"
)
bench$label <- factor(
  bench$label,
  levels = c(
    "A\nSCT + Louvain",
    "B\nlogNorm + Louvain",
    "C\nlogNorm + BayesSpace"
  )
)
# Convert ARI/NMI to long format
bench_long <- bench %>%
  select(label, ARI, NMI) %>%
  pivot_longer(
    cols = c(ARI, NMI),
    names_to = "Metric",
    values_to = "Score"
    
  )
print(bench_long)

# Grouped bar plot
p_benchmark <- ggplot(
  bench_long,
  aes(
    x = label, 
    y = Score,
    fill = Metric
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65
  ) +
# Add exact values above bars
  geom_text(
    aes(label = sprintf("%.3f", Score)),
    position = position_dodge(width = 0.75),
    vjust = -0.45,
    size = 4.2
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    title = "Quantitative benchmark of cortical-layer recovery",
    subtitle = "Agreement with expert manual annotation",
    x = NULL,
    y = "Agreement score",
    fill = NULL
  ) +
  theme_classic(base_size = 13) +
  
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      size = 12,
      hjust = 0.5
    ),
    
    axis.text.x = element_text(
      size = 11
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    axis.title.y = element_text(
      size = 12
    ),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 11
    )
  )
p_benchmark

ggsave(
  "results/figures/09_ABC_ARI_NMI_benchmark.png",
  p_benchmark,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)
p_arms +
  annotate("segment", x = 2.19, xend = 2.81, y = 0.60, yend = 0.60) +
  annotate("segment", x = 2.19, xend = 2.19, y = 0.58, yend = 0.60) +
  annotate("segment", x = 2.81, xend = 2.81, y = 0.58, yend = 0.60) +
  annotate("text", x = 2.5, y = 0.625, size = 5, fontface = "bold",
           label = "+0.178 ARI\nspatial prior only")

cat("\nA/B/C ARI-NMI benchmark figure completed.\n")