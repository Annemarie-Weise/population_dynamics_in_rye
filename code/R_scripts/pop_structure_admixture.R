library(tidyverse)
library(ggplot2)
library(dplyr)
library(clue)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path,"/general_analysis/Admixture"))


################################################################################
# Plot ADMIXTURE cross-validation errors and ancestry proportions
################################################################################
# read ADMIXTURE CV errors + Q-matrix output files
#  - CV-error plot across K values
#  - plot ancestry proportions for selected K values
################################################################################


# read ADMIXTURE cross-validation errors
cv <- read.csv("CV_errors.csv", check.names = FALSE)
colnames(cv) <- c("K", "CV_error")

# plot CV error across tested K values
p <- ggplot(cv, aes(x = K, y = CV_error)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2) +
  labs(
    x = "K",
    y = "CV error",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line         = element_line(color = "black"),
    axis.ticks        = element_line(color = "black"),
    axis.ticks.length = unit(2, "pt"),
    panel.grid        = element_blank()
  )
pdf("cv_errors_admixture.pdf", 
    width = 4, height = 4)
print(p)
dev.off()


# define ADMIXTURE output prefix and metadata input files
prefix <- "output/LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar_chr."
fam_file <- "LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar_chr.fam"
group_file <- paste0(git_path, "/additional_data/id_data/samples_id_gen.csv")

# read PLINK sample IDs and group metadata
fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
groups <- read.csv(group_file, header = TRUE, stringsAsFactors = FALSE)
names(groups) <- c("Old_ID", "group1", "group", "treat")

cols <- c(
  "steelblue", "coral", "lightgreen", "#984EA3", 
  "darkblue", "darkred", "darkgreen", "#F781BF", "#999999"
)


# read Q matrix
read_Q <- function(K, prefix, fam) {
  q <- as.matrix(read.table(paste0(prefix, K, ".Q"), header = FALSE))
  rownames(q) <- fam[, 1]
  return(q)
}

# align components between consecutive K values
align_Q <- function(previous_Q, current_Q) {
  sim <- cor(previous_Q, current_Q)
  sim[is.na(sim)] <- -1
  sim_scaled <- (sim + 1) / 2
  matched <- as.integer(
    solve_LSAP(sim_scaled, maximum = TRUE)
  )
  unmatched <- setdiff(
    seq_len(ncol(current_Q)),
    matched
  )
  component_order <- c(matched, unmatched)
  list(
    Q = current_Q[, component_order, drop = FALSE],
    component_order = component_order
  )
}

# align K = 2 to K = 9
Q_aligned <- list()
component_order <- list()
Q_aligned[["2"]] <- read_Q(2, prefix, fam)
component_order[["2"]] <- seq_len(2)
for (K in 3:9) {
  current_Q <- read_Q(K, prefix, fam)
  aligned <- align_Q(
    Q_aligned[[as.character(K - 1)]],
    current_Q
  )
  Q_aligned[[as.character(K)]] <- aligned$Q
  component_order[[as.character(K)]] <- aligned$component_order
}

component_colors <- bind_rows(
  lapply(2:9, function(K) {
    data.frame(
      K = K,
      aligned_component = seq_len(K),
      original_Q_column = component_order[[as.character(K)]],
      color = cols[seq_len(K)]
    )
  })
)
write.csv(component_colors, "admixture_component_colors.csv",row.names = FALSE)


# fixed individual order based on aligned K = 2
q2 <- as.data.frame(Q_aligned[["2"]])
q2$Old_ID <- rownames(Q_aligned[["2"]])
names(q2)[1:2] <- paste0("V", 1:2)
sample_order <- q2 %>%
  inner_join(groups[, c("Old_ID", "group")], by = "Old_ID") %>%
  mutate(group = factor(group, levels = unique(groups$group))) %>%
  arrange(group, V1) %>%
  pull(Old_ID)


# read and prepare the ADMIXTURE Q matrix for one K value
#  - Q matrix is combined with sample metadata and ordered by group
#  - individuals appear in the same order in all ancestry plots
prep_tbl <- function(K, Q_aligned, groups, sample_order) {
  q <- as.data.frame(Q_aligned[[as.character(K)]])
  q$Old_ID <- rownames(Q_aligned[[as.character(K)]])
  names(q)[seq_len(K)] <- paste0("V", seq_len(K))
  q <- q %>%
    inner_join(
      groups[, c("Old_ID", "group")],
      by = "Old_ID"
    ) %>%
    mutate(
      group = factor(group, levels = unique(groups$group)),
      order_id = match(Old_ID, sample_order)
    ) %>%
    arrange(order_id)
  return(q)
}

# plot ancestry proportions for one K value
plot_K <- function(K, Q_aligned, groups, sample_order,
                   draw_x_labels = FALSE,
                   shared_vpos = NULL) {
  
  # Prepare Q matrix
  tbl <- prep_tbl(K, Q_aligned, groups, sample_order)
  mat <- t(as.matrix(tbl[, paste0("V", seq_len(K))]))
  
  # Draw stacked ancestry barplot
  bp <- barplot(
    mat,
    col = cols[seq_len(K)],
    names.arg = rep("", ncol(mat)),
    ylab = paste0("Anc. Prop., K = ", K),
    border = NA,
    space = 0,
    axes = TRUE
  )
  
  # Find group centers and group boundaries along the x-axis
  runs <- rle(as.character(tbl$group))
  len <- runs$lengths
  grp <- runs$values
  centers_idx <- cumsum(len) - (len - 1) / 2
  centers_x <- bp[centers_idx]
  if (length(len) > 1) {
    boundaries_idx <- cumsum(len)
    vpos_local <- (
      bp[boundaries_idx[-length(boundaries_idx)]] +
        bp[boundaries_idx[-length(boundaries_idx)] + 1]
    ) / 2
  } else {
    vpos_local <- numeric(0)
  }
  
  # use shared group boundaries if supplied, otherwise use local boundaries
  vpos <- if (!is.null(shared_vpos)) shared_vpos else vpos_local
  
  # draw vertical lines between metadata groups
  if (length(vpos) > 0) {
    usr <- par("usr")
    segments(
      x0 = vpos,
      y0 = usr[3] + 0.01,
      x1 = vpos,
      y1 = usr[4],
      col = "black",
      lwd = 1,
      lend = 1
    )
  }
  
  # add group labels only to the bottom panel
  grp <- grp %>%
    str_replace("^gen", "Gen") %>%
    str_replace("_", " ") %>%
    str_replace("Control", "control")
  if (draw_x_labels) {
    par(xpd = NA)
    text(
      centers_x,
      par("usr")[3] - 0.03 * diff(par("usr")[3:4]),
      labels = grp,
      srt = 90,
      cex = 0.9,
      font = 3,
      adj = 1
    )
    par(xpd = FALSE)
  }
  invisible(vpos_local)
}


# plot ADMIXTURE results for K = 6 to K = 9
# - group boundaries are taken from the first panel and reused
pdf("admixture_K6_to_K9.pdf", width = 7, height = 5) #width = 6, height = 5)
par(
  mfcol = c(4, 1),
  mar = c(0.5, 3.0, 0.5, 0.2),
  oma = c(9, 1, 1, 1),
  mgp = c(1.4, 0.35, 0),
  cex.axis = 0.8,
  las = 1,
  xaxs = "i",
  tcl = -0.2
)
vpos_ref <- plot_K(9, Q_aligned, groups, sample_order, draw_x_labels = FALSE)
plot_K(8, Q_aligned, groups, sample_order, draw_x_labels = FALSE, shared_vpos = vpos_ref)
plot_K(7, Q_aligned, groups, sample_order, draw_x_labels = FALSE, shared_vpos = vpos_ref)
plot_K(6, Q_aligned, groups, sample_order, draw_x_labels = TRUE, shared_vpos = vpos_ref)
dev.off()


# Plot ADMIXTURE results for K = 2 to K = 5
pdf("admixture_K2_to_K5.pdf", width = 7, height = 5)# width = 6, height = 5)
par(
  mfcol = c(4, 1),
  mar = c(0.5, 3.0, 0.5, 0.2),
  oma = c(9, 1, 1, 1),
  mgp = c(1.4, 0.35, 0),
  cex.axis = 0.8,
  las = 1,
  xaxs = "i",
  tcl = -0.2
)
vpos_ref <- plot_K(5, Q_aligned, groups, sample_order, draw_x_labels = FALSE)
plot_K(4, Q_aligned, groups, sample_order, draw_x_labels = FALSE, shared_vpos = vpos_ref)
plot_K(3, Q_aligned, groups, sample_order, draw_x_labels = FALSE, shared_vpos = vpos_ref)
plot_K(2, Q_aligned, groups, sample_order, draw_x_labels = TRUE, shared_vpos = vpos_ref)
dev.off()



################################################################################
# Calculate individual admixture statistics from an ADMIXTURE Q files
################################################################################

calc_admixture <- function(q_file, fam_file, id_col = 1) {
  q <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
  
  # determine K from number of ancestry components
  K <- ncol(q)
  names(q) <- paste0("V", seq_len(K))
  
  # read FAM file
  fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  q$Old_ID <- fam[, id_col]
  
  # matrix containing only ancestry proportions
  q_mat <- as.matrix(q[, paste0("V", seq_len(K))])
  
  # calculate individual ancestry statistics
  q$dominant_component <- paste0("V", max.col(q_mat, ties.method = "first"))
  q$dominant_proportion <- apply(q_mat, 1, max)
  q$clusteredness <- sqrt((K / (K - 1)) * rowSums((q_mat - 1 / K)^2))
  
  q$K <- K
  q <- q[, c(
    "Old_ID",
    "K",
    paste0("V", seq_len(K)),
    "dominant_component",
    "dominant_proportion",
    "clusteredness"
  )]
  return(q)
}


# calculate individual ancestry statistics for K values
fam_file <- "LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar_chr.fam"
q_prefix <- "output/LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar_chr."

admixture_stats <- list()
for (K in 2:9) {
  q_file <- paste0(q_prefix, K, ".Q")
  admixture_stats[[as.character(K)]] <- calc_admixture(
    q_file = q_file,
    fam_file = fam_file
  )
}

# combine all K values into one table
admixture_stats <- bind_rows(admixture_stats)
write.csv(admixture_stats, "admixture_individual_stats_K2_K9.csv", row.names = FALSE)

# add group metadata
group_file <- paste0(git_path, "/additional_data/id_data/samples_id_gen.csv")
groups <- read.csv(group_file, header = TRUE, stringsAsFactors = FALSE)
names(groups) <- c("Old_ID", "group1", "group", "treat")
admixture_stats <- admixture_stats %>% left_join(groups, by = "Old_ID")

# calculate summary statistics per K and generation
admixture_summary <- admixture_stats %>%
  group_by(K, group1) %>%
  summarise(
    n = n(),
    mean_clusteredness = mean(clusteredness, na.rm = TRUE),
    median_clusteredness = median(clusteredness, na.rm = TRUE),
    mean_dominant_proportion = mean(dominant_proportion, na.rm = TRUE),
    median_dominant_proportion = median(dominant_proportion, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(admixture_summary, "admixture_generation_summary_K2_K9.csv", row.names = FALSE)



################################################################################
# Quantify dark-green ancestry component across generations for K = 7
################################################################################

# read component-color assignment
component_colors <- read.csv("admixture_component_colors.csv", stringsAsFactors = FALSE)

# identify dark-green component for K = 7
green_col <- component_colors %>%
  filter(K == 7, color == "darkgreen") %>%
  pull(original_Q_column)

q7 <- read.table(paste0(q_prefix, "7.Q"), header = FALSE, stringsAsFactors = FALSE)
names(q7) <- paste0("V", seq_len(ncol(q7)))
fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
q7$Old_ID <- fam[, 1]

# add group metadata and dark-green ancestry proportion
q7 <- q7 %>%
  left_join(groups, by = "Old_ID") %>%
  mutate(darkgreen = .data[[paste0("V", green_col)]])

# summarize dark-green ancestry component per generation
darkgreen_summary <- q7 %>%
  group_by(group1) %>%
  summarise(
    n = n(),
    mean_darkgreen = mean(darkgreen),
    median_darkgreen = median(darkgreen),
    Q1_darkgreen = quantile(darkgreen, 0.25),
    Q3_darkgreen = quantile(darkgreen, 0.75),
    percent_ge_5 = mean(darkgreen >= 0.05) * 100,
    percent_ge_10 = mean(darkgreen >= 0.10) * 100,
    .groups = "drop"
  )
write.csv(darkgreen_summary, "K7_darkgreen_generation_summary.csv", row.names = FALSE)









