library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(irlba)
library(patchwork)
library(grid)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/PCA/gen0_diff_PCA"))


################################################################################
# PCA of allele-frequency differences relative to Gen0
################################################################################
# - calculate ALT-frequency differences between each Gen1-Gen3 population and
#   independent Gen0 reference
# - run PCA on these difference profiles
################################################################################

# read merged allele-frequency table
dt <- fread(
  paste0(git_path, "/general_analysis/vcftools/freq/LD_all_pop_existing_sites_miss0.9.csv.gz"), 
  header = TRUE)

# calculate ALT-frequency differences relative to Gen0
dt <- dt %>%
  mutate(
    ALT_FREQ_DIFF_gen1_ND = gen1_ND - gen0,
    ALT_FREQ_DIFF_gen1_Control = gen1_Control - gen0,
    ALT_FREQ_DIFF_gen2_ND = gen2_ND - gen0,
    ALT_FREQ_DIFF_gen2_Control = gen2_Control - gen0,
    ALT_FREQ_DIFF_gen3_ND = gen3_ND - gen0,
    ALT_FREQ_DIFF_gen3_Control = gen3_Control - gen0
  ) %>% 
  select(
    ALT_FREQ_DIFF_gen1_ND, ALT_FREQ_DIFF_gen1_Control,
    ALT_FREQ_DIFF_gen2_ND, ALT_FREQ_DIFF_gen2_Control,
    ALT_FREQ_DIFF_gen3_ND, ALT_FREQ_DIFF_gen3_Control,
  )

# transpose table so generation-treatment groups are rows + variants are columns
pca_dt <- t(dt)

# run PCA without scaling (alues are already allele-frequency differences)
p1 <- prcomp(pca_dt, center = TRUE, scale. = FALSE)
summary(p1)
# Importance of components:
#   PC1    PC2    PC3    PC4     PC5       PC6
# Standard deviation     2.4663 1.9568 1.3677 1.3231 1.05328 1.016e-14
# Proportion of Variance 0.4154 0.2615 0.1278 0.1196 0.07577 0.000e+00
# Cumulative Proportion  0.4154 0.6769 0.8047 0.9242 1.00000 1.000e+00

# export PCA scores for each generation-treatment comparison
pca_scores <- as.data.frame(p1$x)
pca_scores$sample <- rownames(pca_scores)
pca_scores <- pca_scores[, c("sample", setdiff(names(pca_scores), "sample"))]
write.csv(pca_scores, "pca_gen0_diffs_scores.csv", row.names = FALSE)

# export variance proportions
pca_var <- data.frame(
  PC = paste0("PC", seq_along(p1$sdev)),
  Standard_deviation = p1$sdev,
  Proportion_of_variance = (p1$sdev^2) / sum(p1$sdev^2),
  Cumulative_proportion = cumsum((p1$sdev^2) / sum(p1$sdev^2))
)
write.csv(pca_var, "pca_gen0_diffs_variance_summary.csv", row.names = FALSE)


################################################################################
# Save screeplot
################################################################################

# calculate variance explained per PC
prop_var <- (p1$sdev^2) / sum(p1$sdev^2)
scree_df <- data.frame(
  PC_num = seq_along(prop_var),
  PC = paste0("PC", seq_along(prop_var)),
  Proportion = prop_var
)

# plot variance explained by each PC
scree_plot <- ggplot(scree_df, aes(x = PC_num, y = Proportion)) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_x_continuous(
    breaks = scree_df$PC_num,
    labels = scree_df$PC
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "", y = "Explained Variance") +
  theme_classic()
ggsave("pca_gen0_diffs_screeplot.pdf", plot = scree_plot, width = 3, height = 3)


################################################################################
# Save PCA scatter plot
################################################################################

# add generation and treatment labels from row names
plot_df <- as.data.frame(p1$x)
plot_df$sample <- rownames(plot_df)
plot_df$generation <- sub(".*_(gen[0-9]+)_.*", "\\1", plot_df$sample)
plot_df$treatment  <- sub(".*_(ND|Control)$", "\\1", plot_df$sample)

# set plotting order and display labels
plot_df$generation <- factor(
  plot_df$generation,
  levels = c("gen1", "gen2", "gen3"),
  labels = c("Gen1", "Gen2", "Gen3")
)
plot_df$treatment <- factor(
  plot_df$treatment,
  levels = c("ND", "Control"),
  labels = c("ND", "Control")
)

# create PC axis labels with explained variance
var_expl <- (p1$sdev^2) / sum(p1$sdev^2) * 100
xlab_pc1 <- paste0("PC1 (", round(var_expl[1], 1), "%)")
ylab_pc2 <- paste0("PC2 (", round(var_expl[2], 1), "%)")
ylab_pc3 <- paste0("PC3 (", round(var_expl[3], 1), "%)")

# shared colors, shapes, and theme settings for both PCA panels
common_scales <- list(
  scale_color_manual(
    name = "Generation",
    values = c("Gen1" = "darkgreen", "Gen2" = "darkred", "Gen3" = "darkblue")
  ),
  scale_shape_manual(
    name = "Treatment",
    values = c("ND" = 16, "Control" = 17),
    labels = c("ND" = "ND", "Control" = "control")
  ),
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      direction = "horizontal",
      title.position = "top",
      theme = theme(
        legend.background = element_rect(color = "black", fill = "white"),
        legend.margin = margin(4, 4, 4, 4)
      )
    ),
    shape = guide_legend(
      nrow = 1,
      byrow = TRUE,
      direction = "horizontal",
      title.position = "top",
      theme = theme(
        legend.background = element_rect(color = "black", fill = "white"),
        legend.margin = margin(4, 4, 4, 4)
      )
    )
  ),
  theme_classic(),
  theme(aspect.ratio = 1)
)

# plot PC1/PC2 and PC1/PC3
p_pc12 <- ggplot(plot_df, aes(PC1, PC2, color = generation, shape = treatment)) +
  geom_point(size = 4) +
  labs(x = xlab_pc1, y = ylab_pc2) +
  common_scales
p_pc13 <- ggplot(plot_df, aes(PC1, PC3, color = generation, shape = treatment)) +
  geom_point(size = 4) +
  labs(x = xlab_pc1, y = ylab_pc3) +
  common_scales

# combine both panels and collect one shared legend
p <- (p_pc12 + p_pc13) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.spacing.x = unit(1.5, "cm")
    )
  )
ggsave("pca_gen0_diffs_plot.pdf", plot = p, width = 5, height = 5)
