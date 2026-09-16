library(pcadapt)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(grid)
library(ggtext)
library(concaveman)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/PCA/indiv_pcadapt_PCA"))


################################################################################
# Plotting helper functions
################################################################################


# create a screeplot from a pcadapt result object
pcadapt_screeplot <- function(result) {
  prop_var <- (result$singular.values^2)
  df <- data.frame(
    PC = seq_along(prop_var),
    proportion = prop_var
  )
  p <- ggplot(df, aes(x = PC, y = proportion)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue") +
    labs(
      x = "Principal component",
      y = "Variance proportion",
      title = ""
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.grid        = element_blank()
    )
  return(p)
}


# calculate concave hull polygons for each group in PCA space
concave_hull_df <- function(df, x = "PC1", y = "PC2", group = "GEN", concavity = 4) {
  out <- lapply(split(df, df[[group]]), function(d) {
    # only individuals with complete coordinates for the PCs
    d <- d[complete.cases(d[, c(x, y)]), ]
    
    # calculate concave hull around the points of this group
    hull_mat <- concaveman(as.matrix(d[, c(x, y)]),concavity = concavity)
    # convert hull coordinates back to a data frame for ggplot
    hull_df <- as.data.frame(hull_mat)
    colnames(hull_df) <- c(x, y)
    
    # keep group label so hulls can be colored/faceted by group
    hull_df[[group]] <- d[[group]][1]
    hull_df
  })
  do.call(rbind, out)
}


# create PC1/PC2 + PC1/PC3 scatter plots from pcadapt result
# - individuals are colored by generation and shaped by treatment
# - concave hulls can be drawn around generation groups
pcadapt_pc_scatter <- function(result, legend_pos = c(0.99, 0.9), fam_file, id_map_file, hull = TRUE) {
  
  # squared singular values for PC axis labels
  prop_var <- result$singular.values^2
  lab_pc1 <- paste0("PC1 (", round(prop_var[1], 3), ")")
  lab_pc2 <- paste0("PC2 (", round(prop_var[2], 3), ")")
  lab_pc3 <- paste0("PC3 (", round(prop_var[3], 3), ")")
  
  # read PLINK fam file to match scores to individuals
  id_map <- read.csv(id_map_file, header = TRUE,stringsAsFactors = FALSE)
  fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
  ids <- fam$IID
  # convert pcadapt scores to tibble and add individual IDs from fam file
  colnames(result$scores) <- paste0("PC", seq_len(ncol(result$scores)))
  scores <- as_tibble(result$scores) %>% mutate(ID = ids, .before = 1)
  
  # add generation and treatment metadata + remove incomplete records
  scores <- scores %>%
    left_join(id_map %>% select(ID, GEN, TREAT2), by = "ID") %>%
    filter(!is.na(GEN), !is.na(TREAT2), !is.na(PC1), !is.na(PC2), !is.na(PC3))
  # set plotting order and labels for generation and treatment
  scores <- scores %>%
    mutate(
      GEN = factor(
        GEN,
        levels = c("gen0", "gen1", "gen2", "gen3"),
        labels = c("Gen0", "Gen1", "Gen2", "Gen3")
      ),
      TREAT = factor(
        TREAT2,
        levels = c("gen0", "ND", "Control"),
        labels = c("Gen0", "ND", "Control")
      )
    )
  
  # define colors for generations and shapes for treatments
  gen_colours <- c(
    "Gen0" = "grey42",
    "Gen1" = "darkgreen",
    "Gen2" = "darkred",
    "Gen3" = "darkblue"
  )
  treat_shapes <- c(
    "Gen0" = 15,
    "ND"   = 16,
    "Control"  = 17
  )
  
  # common ggplot settings for both PCA panels
  base_opts <- list(
    scale_color_manual(values = gen_colours, drop = TRUE, na.translate = FALSE),
    scale_shape_manual(
      values = treat_shapes,
      labels = c(
        "Gen0" = "Gen0",
        "ND" = "ND",
        "Control" = "control"
      ),
      drop = TRUE,
      na.translate = FALSE
    ),
    theme_classic(base_size = 12),
    theme(
      legend.title      = element_text(face = "bold"),
      plot.margin       = margin(4, 2, 2, 2),
      panel.spacing     = unit(1, "mm"),
      axis.line         = element_line(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(2, "pt"),
      panel.grid        = element_blank()
    )
  )
  
  # calculate hulls around generations for both plotted PC combinations
  hull_pc12 <- concave_hull_df(scores, x = "PC1", y = "PC2", group = "GEN", concavity = 4)
  hull_pc13 <- concave_hull_df(scores, x = "PC1", y = "PC3", group = "GEN", concavity = 4)
  
  # PC1 vs PC2 panel without legend
  p12 <- ggplot(scores, aes(PC1, PC2, color = GEN, shape = TREAT)) +
    geom_point(alpha = 0.5, size = 1.9) +
    labs(
      x = lab_pc1,
      y = lab_pc2,
      color = "Generation",
      shape = "Treatment"
    ) +
    base_opts +
    theme(
      legend.position = "none",
      plot.title = element_markdown(size = 14)
    )
  # PC1 vs PC3 panel with legend
  p13 <- ggplot(scores, aes(PC1, PC3, color = GEN, shape = TREAT)) +
    geom_point(alpha = 0.5, size = 1.9) +
    labs(
      x = lab_pc1,
      y = lab_pc3,
      color = "Generation",
      shape = "Treatment"
    ) +
    base_opts +
    theme(
      legend.position      = legend_pos,
      legend.justification = c(1, 0.5),
      legend.background    = element_rect(fill = scales::alpha("white", 0.7), colour = NA),
      legend.key.size      = unit(4, "mm"),
      plot.title           = element_markdown(size = 14)
    ) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 3)),
      shape = guide_legend(override.aes = list(alpha = 1, size = 3))
    )
  
  # Add concave hull outlines if requested
  if (hull == TRUE) {
    p12 <- p12 +
      geom_polygon(
        data = hull_pc12,
        aes(PC1, PC2, group = GEN, color = GEN),
        fill = NA,
        linewidth = 0.6,
        inherit.aes = FALSE
      )
    p13 <- p13 +
      geom_polygon(
        data = hull_pc13,
        aes(PC1, PC3, group = GEN, color = GEN),
        fill = NA,
        linewidth = 0.6,
        inherit.aes = FALSE
      )
  }
  
  # Combine both PCA panels into one figure
  p <- plot_grid(p12, p13, ncol = 2, align = "hv")
  return(p)
}




################################################################################
# Individual-level PCA with pcadapt
################################################################################
# PCA is run with pcadapt on the LD-pruned genotype dataset.
#  - first run uses K = 20 to inspect the screeplot
#  - second run with K = 3 -> export PCA scores and create the final plots.
################################################################################


pcadapt_data <- read.pcadapt("LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar.bed", 
                             type = "bed")

# run PCA with many PCs to inspect the screeplot
K <- 20
MAF <- 0
method <- "mahalanobis"
result <- pcadapt(pcadapt_data, K = K, min.maf = MAF, method = method)
p <- pcadapt_screeplot(result)
pdf("pca_indv_screeplot_pcadapt.pdf", 
    width = 4, height = 4)
print(p)
dev.off()

# re-run PCA with selected number of PCs for final output
K <- 3
result <- pcadapt(pcadapt_data, K = K, min.maf = MAF ,method = method)
saveRDS(result, file = "pcadapt_result.rds")

# Add individual IDs from the PLINK fam file to the PCA scores 
# + export PCA scores for downstream use
fam <- read.table(
  "LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar.fam",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(fam) <- c("FID","IID","PAT","MAT","SEX","PHENO")
pca_scores <- as.data.frame(result$scores)
stopifnot(nrow(pca_scores) == nrow(fam))
pca_scores$ID <- fam$IID
pca_scores <- pca_scores[, c("ID", setdiff(names(pca_scores), "ID"))]
write.csv(pca_scores, "pcadapt_scores_outlier_pcadapt.csv", row.names = FALSE)

# create PCA scatter plot with generation and treatment metadata
fam_file <- "LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar.fam"
id_map_file <- paste0(git_path,"/additional_data/id_data/samples_id_gen.csv")
p <- pcadapt_pc_scatter(result,legend_pos = c(0.99, 0.7), fam_file = fam_file, id_map_file = id_map_file)
pdf(paste0("pca_pcadapt_outlier_plot.pdf"), 
    width = 6.5, height = 3.5)
print(p)
dev.off()




################################################################################
# Calculate PCA variance within generations
################################################################################
# calculate how much individuals vary within each generation in PC1/PC2 space
# calculating: 
# - within-generation means
# - variances and covariance
# - total PC1/PC2 variance.
################################################################################


result <- readRDS("pcadapt_result.rds")
fam <- read.table(
  "LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar.fam",
  header = FALSE,
  stringsAsFactors = FALSE
)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO")
id_map <- read.csv(
  paste0(git_path, "/additional_data/id_data/samples_id_gen.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# add PC column names and combine scores with sample metadata
colnames(result$scores) <- paste0("PC", seq_len(ncol(result$scores)))
scores <- as_tibble(result$scores) %>%
  mutate(ID = fam$IID, .before = 1) %>%
  left_join(id_map %>% select(ID, GEN, TREAT2), by = "ID") %>%
  filter(!is.na(GEN), !is.na(PC1), !is.na(PC2)) %>%
  mutate(
    GEN = factor(
      GEN,
      levels = c("gen0", "gen1", "gen2", "gen3"),
      labels = c("Gen0", "Gen1", "Gen2", "Gen3")
    )
  )

# calculate within-generation spread in PC1/PC2 space
generation_variance <- scores %>%
  group_by(GEN) %>%
  summarise(
    n = n(),
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2),
    var_PC1 = var(PC1),
    var_PC2 = var(PC2),
    cov_PC1_PC2 = cov(PC1, PC2),
    total_PC1_PC2_variance = var_PC1 + var_PC2,
    .groups = "drop"
  )
write.csv(
  generation_variance,
  "variance_within_generation_PC1_PC2.csv",
  row.names = FALSE
)
