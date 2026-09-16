library(vcfR)
library(adegenet)
library(StAMPP)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(tidyr)
library(tidyverse)
library(grid)
library(rstatix)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/Fst_StAMPP"))

# Path to file that can be requested from the authors:
# - LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar.vcf.gz
vcf_path <- readline(
  prompt = "Please enter the path to LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar.vcf.gz: "
)
if (!file.exists(vcf_path)) {
  stop("File not found: ", vcf_path)
}
vcf_file <- normalizePath(vcf_path)


################################################################################
# Pairwise FST analysis with StAMPP
################################################################################
# - calculate pairwise FST values with bootstrap and confidence intervals
# - test diffrence significance with Wilcoxon test
# - plot pairwise comparison of FST differences and ratios
################################################################################


################################################################################
# Calculate pairwise FST with bootstrap confidence intervals
################################################################################

# read VCF and convert to genlight format
vcf <- read.vcfR(vcf_path)
gl <- vcfR2genlight(vcf)

# match sample IDs from the VCF to population ID files
id_dir <- paste0(git_path,"/additional_data/id_data")
id_files <- list.files(id_dir, pattern = "_IDs\\.txt$", full.names = TRUE)
sample_ids <- trimws(indNames(gl))
pop_vector <- rep(NA, length(sample_ids))
for (f in id_files) {
  ids <- readLines(f, warn = FALSE)
  ids <- trimws(ids)
  ids <- ids[ids != ""]
  pop_name <- basename(f)
  pop_name <- sub("_IDs\\.txt$", "", pop_name)
  pop_vector[sample_ids %in% ids] <- pop_name
}

# assign populations to genlight object and check unmatched individuals
pop(gl) <- as.factor(pop_vector)
table(pop(gl), useNA = "ifany")
indNames(gl)[is.na(pop(gl))]
# remove individuals without population assignment
gl <- gl[!is.na(pop(gl))]

# run pairwise FST analysis
fst_res <- stamppFst(gl, nboots = 1000, percent = 95, nclusters = 1)


# save estimates, p-values, and bootstrap confidence intervals
write.csv(fst_res$Fsts,
          file = file.path("pairwise_FST_matrix_StAMPP.csv"),
          row.names = TRUE)
write.csv(fst_res$Pvalues,
          file = file.path("pairwise_Pvalues_matrix_StAMPP.csv"),
          row.names = TRUE)
write.csv(fst_res$Bootstraps,
          file = file.path("pairwise_FST_bootstraps_CI_StAMPP.csv"),
          row.names = FALSE)
write.csv(fst_res$Bootstraps_unsorted,
          file = "pairwise_FST_bootstraps_unsorted_StAMPP.csv",
          row.names = FALSE)


################################################################################
# Plot selected pairwise FST values with confidence intervals
################################################################################

# extract FST estimates, confidence intervals, and p-values
boot <- read.csv("pairwise_FST_bootstraps_CI_StAMPP.csv",
                 header = TRUE, 
                 check.names = FALSE)
fst_summary <- boot[, c(
  "Population1",
  "Population2",
  "Fst",
  "Lower bound CI limit",
  "Upper bound CI limit",
  "p-value"
)]
names(fst_summary) <- c(
  "pop1", "pop2", "Fst", "CI_low", "CI_high", "p_value"
)

# define selected population contrasts and their plotting order
plot_order <- tibble(
  group_label = c(
    rep("Gen0 vs control treatment", 3),
    rep("Gen0 vs ND", 3),
    rep("control treatment over generations", 3),
    rep("ND over generations", 3),
    rep("ND vs control treatment", 3)
  ),
  group_color = c(
    rep("steelblue", 3),
    rep("coral", 3),
    rep("steelblue", 3),
    rep("coral", 3),
    rep("#7B3294", 3)
  ),
  pop1 = c(
    "gen0", "gen0", "gen0",
    "gen0", "gen0", "gen0",
    "gen0", "gen1_Control", "gen2_Control",
    "gen0", "gen1_ND", "gen2_ND",
    "gen1_ND", "gen2_ND", "gen3_Control"
  ),
  pop2 = c(
    "gen1_Control", "gen2_Control", "gen3_Control",
    "gen1_ND", "gen2_ND", "gen3_ND",
    "gen1_Control", "gen2_Control", "gen3_Control",
    "gen1_ND", "gen2_ND", "gen3_ND",
    "gen1_Control", "gen2_Control", "gen3_ND"
  ),
  short_label = c(
    # gen0 vs Control
    "Gen1", "Gen2", "Gen3",
    # gen0 vs ND
    "Gen1", "Gen2", "Gen3",
    # Control over generations
    "Gen0 vs Gen1", "Gen1 vs Gen2", "Gen2 vs Gen3",
    # ND over generations
    "Gen0 vs Gen1", "Gen1 vs Gen2", "Gen2 vs Gen3",
    # ND vs Control
    "Gen1", "Gen2", "Gen3"
  ),
  order_in_group = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)
)

# add FST estimates and confidence intervals to selected contrasts
plot_df <- plot_order %>%
  left_join(fst_summary, by = c("pop1", "pop2")) %>%
  mutate(
    group_label = factor(
      group_label,
      levels = c(
        "Gen0 vs control treatment",
        "Gen0 vs ND",
        "control treatment over generations",
        "ND over generations",
        "ND vs control treatment"
      )
    )
  )

# create plotting labels and colors for grouped FST contrasts
plot_df <- plot_df %>%
  mutate(panel_x = paste(group_label, order_in_group, sep = "__"))
plot_df$panel_x <- factor(plot_df$panel_x, levels = rev(plot_df$panel_x))
label_map <- setNames(plot_df$short_label, plot_df$panel_x)
group_colors <- c(
  "Gen0 vs control treatment" = "coral",
  "Gen0 vs ND" = "steelblue",
  "control treatment over generations" = "coral",
  "ND over generations" = "steelblue",
  "ND vs control treatment" = "#7B3294"
)

# plot FST estimates with bootstrap confidence intervals
p1 <- ggplot(plot_df, aes(x = panel_x, y = Fst, color = group_label)) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high), size = 0.5) +
  coord_flip() +
  facet_wrap(~ group_label, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_discrete(labels = label_map) +
  scale_color_manual(values = group_colors, guide = "none") +
  theme_bw() +
  labs(x = NULL, y = expression(F[ST])) +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.8, "lines")
  )
ggsave("Fst_vals_with_CIs.pdf", p1, width = 5, height = 6)


################################################################################
# Test differences between FST bootstrap distributions
################################################################################

boot <- read.csv("pairwise_FST_bootstraps_CI_StAMPP.csv",
                 header = TRUE, 
                 check.names = FALSE)
# bootstrap replicate columns are named 1 to 1000 in the StAMPP output
boot_cols <- as.character(1:1000)

# define FST contrasts to compare
fst_pairs <- tibble(
  pop1 = c(
    "gen0", "gen0", "gen0",
    "gen0", "gen0", "gen0",
    "gen1_Control", "gen2_Control",
    "gen1_ND", "gen2_ND",
    "gen1_ND", "gen2_ND", "gen3_Control"
  ),
  pop2 = c(
    "gen1_Control", "gen2_Control", "gen3_Control",
    "gen1_ND", "gen2_ND", "gen3_ND",
    "gen2_Control", "gen3_Control",
    "gen2_ND", "gen3_ND",
    "gen1_Control", "gen2_Control", "gen3_ND"
  ),
  label = c(
    "Gen0 vs Gen1 control",
    "Gen0 vs Gen2 control",
    "Gen0 vs Gen3 control",
    "Gen0 vs Gen1 ND",
    "Gen0 vs Gen2 ND",
    "Gen0 vs Gen3 ND",
    "Gen1 vs Gen2 control",
    "Gen2 vs Gen3 control",
    "Gen1 vs Gen2 ND",
    "Gen2 vs Gen3 ND",
    "Gen1 ND vs control",
    "Gen2 ND vs control",
    "Gen3 ND vs control"
  )
)

# extract bootstrap replicate values for one population pair
get_boot_values <- function(pop1, pop2) {
  row <- boot %>%
    filter(Population1 == pop1, Population2 == pop2)
  as.numeric(row[, boot_cols])
}

# convert bootstrap replicate values to long format for pairwise testing
fst_long <- fst_pairs %>%
  mutate(
    Fst_boot = map2(pop1, pop2, get_boot_values)
  ) %>%
  select(label, Fst_boot) %>%
  unnest_longer(
    Fst_boot,
    values_to = "Fst",
    indices_to = "bootstrap_id"
  )

# test whether selected FST contrasts differ across bootstrap replicates
wilcox_results <- fst_long %>%
  pairwise_wilcox_test(
    Fst ~ label,
    p.adjust.method = "BH"
  )
wilcox_results <- wilcox_results %>%
  mutate(significant = p.adj < 0.05)


# plot which FST contrasts differ significantly
ggplot(wilcox_results, aes(x = group2, y = group1, fill = significant)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "darkred")) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Different?") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
# Result: Diffrences are significant between every pair





################################################################################
# Plot pairwise FST differences and log2 ratios
################################################################################


boot <- read.csv(
  "pairwise_FST_bootstraps_CI_StAMPP.csv",
  header = TRUE,
  check.names = FALSE
)

# define contrasts to compare
fst_pairs <- tibble(
  pop1 = c(
    "gen0", "gen0", "gen0",
    "gen0", "gen0", "gen0",
    "gen1_Control", "gen2_Control",
    "gen1_ND", "gen2_ND",
    "gen1_ND", "gen2_ND", "gen3_Control"
  ),
  pop2 = c(
    "gen1_Control", "gen2_Control", "gen3_Control",
    "gen1_ND", "gen2_ND", "gen3_ND",
    "gen2_Control", "gen3_Control",
    "gen2_ND", "gen3_ND",
    "gen1_Control", "gen2_Control", "gen3_ND"
  ),
  label = c(
    "Gen0 vs Gen1 control",
    "Gen0 vs Gen2 control",
    "Gen0 vs Gen3 control",
    "Gen0 vs Gen1 ND",
    "Gen0 vs Gen2 ND",
    "Gen0 vs Gen3 ND",
    "Gen1 vs Gen2 control",
    "Gen2 vs Gen3 control",
    "Gen1 vs Gen2 ND",
    "Gen2 vs Gen3 ND",
    "Gen1 ND vs control",
    "Gen2 ND vs control",
    "Gen3 ND vs control"
  )
)

# add mean FST estimates to selected contrasts
fst_values <- fst_pairs %>%
  left_join(
    boot %>% select(Population1, Population2, Fst),
    by = c("pop1" = "Population1", "pop2" = "Population2")
  )

# create all pairwise combinations of FST contrasts
labels <- fst_values$label
n <- length(labels)
grid_df <- expand.grid(i = seq_len(n), j = seq_len(n)) %>%
  mutate(
    Comparison_A = labels[i],
    Comparison_B = labels[j],
    Fst_A = fst_values$Fst[i],
    Fst_B = fst_values$Fst[j],
    mean_diff = Fst_A - Fst_B,
    log2_ratio = log2(Fst_A / Fst_B)
  )

# split matrix into upper and lower triangles for separate color scales
upper_df <- grid_df %>%
  filter(i < j)
lower_df <- grid_df %>%
  filter(i > j)

# plot absolute differences in upper triangle and log2 ratios in lower triangle
# upper triangle: FST_A - FST_B
# lower triangle: log2(FST_A / FST_B)
p_combined <- ggplot() +
  geom_tile(
    data = upper_df,
    aes(x = j, y = i, fill = mean_diff),
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_gradient2(
    low = "forestgreen",
    mid = "white",
    high = "red2",
    midpoint = 0,
    name = "A - B"
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = lower_df,
    aes(x = j, y = i, fill = log2_ratio),
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkgoldenrod1",
    midpoint = 0,
    name = "log2(A / B)"
  ) +
  annotate(
    "segment",
    x = 0.5, y = 0.5,
    xend = n + 0.5, yend = n + 0.5,
    colour = "black",
    linewidth = 1.5,
    lineend = "round"
  ) +
  annotate(
    "segment",
    x = 0.5, y = 0.5,
    xend = n + 0.5, yend = n + 0.5,
    colour = "white",
    linewidth = 0.75,
    lineend = "round"
  ) +
  coord_fixed() +
  scale_x_continuous(
    breaks = seq_len(n),
    labels = labels,
    expand = c(0, 0)
  ) +
  scale_y_reverse(
    breaks = seq_len(n),
    labels = labels,
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.4, "cm")
  )
ggsave("Fst_differences_and_log2_ratio.pdf", p_combined, width = 5, height = 5)
