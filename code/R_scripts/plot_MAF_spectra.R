library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(purrr)
library(broom)
library(gtools)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/vcftools/freq"))


################################################################################
# Plot minor allele frequency spectra across generations
################################################################################
# - combine MAF values across generations and treatments
# - remove zero-frequency sites
# - create dodged MAF histograms at two zoom levels
################################################################################


################################################################################
# Plotting helper functions
################################################################################

# create dodged MAF histogram for all generations and treatments
# -> generations shown by color and treatments by transparency
plot_maf_hist <- function(df, x_max, bin_width) {
  # define minor and major MAF bin positions for visual guides
  bin_edges <- seq(0, x_max, by = bin_width)
  major_edges <- seq(0, x_max, by = 0.1)
  
  ggplot(df, aes(x = MAF)) +
    geom_vline(
      xintercept = bin_edges,
      linetype = "dashed",
      colour = "grey85",
      linewidth = 0.25
    ) +
    geom_vline(
      xintercept = major_edges,
      linetype = "dashed",
      colour = "grey65",
      linewidth = 0.5
    ) +
    geom_histogram(
      aes(fill = gen, alpha = trt2, group = grp),
      binwidth = bin_width,
      boundary = 0,
      position = position_dodge2(preserve = "single", padding = 0.05),
      color = "black",
      linewidth = 0.2
    ) +
    scale_x_continuous(
      limits = c(0, x_max),
      breaks = major_edges,
      minor_breaks = bin_edges
    ) +
    scale_fill_manual(
      values = base_cols,
      labels = c(
        gen0 = "Gen0",
        gen1 = "Gen1",
        gen2 = "Gen2",
        gen3 = "Gen3"
      ),
      name = NULL
    ) +
    scale_alpha_manual(
      values = c(Control = 1.0, ND = 0.5),
      labels = c(Control = "control", ND = "ND"),
      name = NULL,
      na.translate = FALSE
    ) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
      alpha = guide_legend(
        order = 2,
        override.aes = list(fill = "grey50", colour = "black")
      )
    ) +
    labs(x = "MAF", y = "Variant count") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(3, "pt"),
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.box = "vertical",
      legend.background = element_rect(fill = "white", colour = "grey70"),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.text = element_text(size = 11)
    )
}


################################################################################
# Create MAF histograms at different zoom levels
################################################################################


pops <- c(
  "gen1_ND", "gen1_Control",
  "gen2_ND", "gen2_Control",
  "gen3_ND", "gen3_Control",
  "gen0"
)
base_cols <- c(
  gen0 = "grey42",
  gen1 = "darkgreen",
  gen2 = "darkred",
  gen3 = "darkblue"
)

# read all population-specific MAF tables and add plotting groups
maf_all <- rbindlist(lapply(pops, function(pop) {
  dt <- fread(paste0("LD_freq_", pop, "_MAF_0.9miss.csv.gz"))
  dt[, pop := pop]
  dt
})) %>%
  # keep only polymorphic sites with finite MAF values
  # -> needed to not overestimate the number of rare alleles because of zeros
  filter(!is.na(MAF), MAF > 0) %>%
  mutate(
    gen = if_else(pop == "gen0", "gen0", str_extract(pop, "^gen\\d+")),
    trt = case_when(
      pop == "gen0" ~ "gen0",
      str_detect(pop, "_Control$") ~ "Control",
      str_detect(pop, "_ND$") ~ "ND"
    ),
    trt2 = factor(if_else(trt %in% c("ND", "Control"), trt, NA_character_),
                  levels = c("Control", "ND")),
    # set plotting order for generations, treatments, and groups
    grp = if_else(pop == "gen0", "gen0", paste0(gen, "_", trt)),
    gen = factor(gen, levels = c("gen0", "gen1", "gen2", "gen3")),
    trt = factor(trt, levels = c("ND", "Control", "gen0")),
    pop = factor(pop, levels = pops),
    grp = factor(
      grp,
      levels = c(
        "gen0",
        "gen1_Control", "gen1_ND",
        "gen2_Control", "gen2_ND",
        "gen3_Control", "gen3_ND"
      )
    )
  )

# full MAF range up to 0.5
ggsave(
  "LD_MAF_hist_dodged_no-zero_miss0.9.pdf",
  plot_maf_hist(maf_all, x_max = 0.5, bin_width = 0.02),
  width = 8,
  height = 5.5,
  units = "in"
)

# zoom into rare and low-frequency variants
ggsave(
  "LD_MAF_hist_zoomed_in_no-zero_miss0.9.pdf",
  plot_maf_hist(maf_all, x_max = 0.2, bin_width = 0.01),
  width = 8,
  height = 5.5,
  units = "in"
)



################################################################################
# Summarize variant counts in low-MAF bins
################################################################################

maf_low_counts <- maf_all %>%
  filter(gen %in% c("gen1", "gen2", "gen3"), MAF <= 0.03) %>%
  mutate(
    MAF_bin = case_when(
      MAF <= 0.01 ~ "0.00-0.01",
      MAF <= 0.02 ~ "0.01-0.02",
      MAF <= 0.03 ~ "0.02-0.03"
    )
  ) %>%
  count(gen, trt, MAF_bin, name = "n_variants") %>%
  arrange(MAF_bin, trt, gen)
write.csv(maf_low_counts, "MAF_low_frequency_bin_counts.csv", row.names = FALSE)



################################################################################
# Combine temporal and treatment differences in one table
################################################################################

maf_low_differences <- bind_rows(
  # temporal differences
  maf_low_counts %>%
    select(gen, trt, MAF_bin, n_variants) %>%
    pivot_wider(names_from = gen, values_from = n_variants) %>%
    mutate(
      Gen1_to_Gen2 = gen2 - gen1,
      Gen2_to_Gen3 = gen3 - gen2,
      Gen1_to_Gen3 = gen3 - gen1
    ) %>%
    select(MAF_bin, trt, Gen1_to_Gen2, Gen2_to_Gen3, Gen1_to_Gen3) %>%
    pivot_longer(
      cols = starts_with("Gen"),
      names_to = "comparison",
      values_to = "difference"
    ) %>%
    mutate(comparison_type = "temporal", group = as.character(trt)) %>%
    select(MAF_bin, comparison_type, group, comparison, difference),
  # treatment differences
  maf_low_counts %>%
    select(gen, trt, MAF_bin, n_variants) %>%
    pivot_wider(
      names_from = trt,
      values_from = n_variants
    ) %>%
    transmute(
      MAF_bin,
      comparison_type = "treatment",
      group = as.character(gen),
      comparison = "ND_minus_control",
      difference = ND - Control
    )
)
write.csv(maf_low_differences, "MAF_low_frequency_bin_differences.csv", row.names = FALSE)
