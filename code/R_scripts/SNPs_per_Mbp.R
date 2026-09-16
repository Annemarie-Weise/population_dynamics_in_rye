library(tidyverse)
library(readr)
library(ggplot2)
library(gtools)
library(patchwork)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/vcftools/freq"))

################################################################################
# Plot SNP counts per 25 Mbp window
################################################################################
# - count all observed and significant SNPs in 25 Mbp genomic windows
# - plot total SNP counts as generation-specific lines 
# - plot fraction of significant SNPs compared to observed SNPs
################################################################################

# window size used for genomic binning
win_bp <- 25e6

# read merged allele-frequency table and assign SNPs to 25 Mbp windows
all_sites <- read_csv(
  "LD_all_pop_existing_sites_miss0.9.csv.gz",
  show_col_types = FALSE
) %>%
  select(CHROM, POS, ALT_AC_gen0,
         ALT_AC_gen1_ND,  ALT_AC_gen2_ND,  ALT_AC_gen3_ND,
         ALT_AC_gen1_Control, ALT_AC_gen2_Control, ALT_AC_gen3_Control
  ) %>%
  mutate(
    CHROM = factor(CHROM, levels = mixedsort(unique(CHROM))),
    POS = as.numeric(POS),
    BIN_START = floor(POS / win_bp) * win_bp,
    BIN_END = BIN_START + win_bp,
    midBP = BIN_START + win_bp / 2
  )

# convert allele-count columns to long format + add generation/treatment labels
counts_all <- all_sites %>%
  pivot_longer(
    cols = starts_with("ALT_AC_"),
    names_to = "variable",
    values_to = "ALT_AC"
  ) %>%
  mutate(
    Generation = case_when(
      variable == "ALT_AC_gen0" ~ "Gen0",
      str_detect(variable, "gen1") ~ "Gen1",
      str_detect(variable, "gen2") ~ "Gen2",
      str_detect(variable, "gen3") ~ "Gen3"
    ),
    Treatment = case_when(
      variable == "ALT_AC_gen0" ~ "Gen0",
      str_detect(variable, "ND") ~ "ND",
      str_detect(variable, "Control") ~ "Control"
    )
  ) %>%
  # count only SNPs where the ALT allele is present in population
  filter(ALT_AC > 0) %>%
  mutate(
    Generation = factor(Generation, levels = c("Gen0", "Gen1", "Gen2", "Gen3"))
  )

# duplicate Gen0 into both treatment panels
counts_all <- bind_rows(
  counts_all %>%
    filter(Treatment != "Gen0"),
  counts_all %>%
    filter(Treatment == "Gen0") %>%
    mutate(Treatment = "ND"),
  counts_all %>%
    filter(Treatment == "Gen0") %>%
    mutate(Treatment = "Control")
)

# count all ALT-present SNPs per chromosome, treatment, generation, and window
df_counts <- counts_all %>%
  group_by(CHROM, Treatment, Generation, BIN_START, BIN_END, midBP) %>%
  summarise(
    n_all = n(),
    .groups = "drop"
  )

# read significant SNPs and assign them to the same 25 Mbp windows
special_snps <- read_csv(
  "logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz",
  show_col_types = FALSE
) %>%
  select(CHROM, POS, SIG_ND, SIG_Control
  ) %>%
  mutate(
    CHROM = factor(CHROM, levels = levels(all_sites$CHROM)),
    POS = as.numeric(POS),
    BIN_START = floor(POS / win_bp) * win_bp,
    BIN_END = BIN_START + win_bp,
    midBP = BIN_START + win_bp / 2
  )

# count significant SNPs per treatment and window
df_special <- bind_rows(
  special_snps %>%
    filter(SIG_ND) %>%
    mutate(Treatment = "ND"),
  special_snps %>%
    filter(SIG_Control) %>%
    mutate(Treatment = "Control")
) %>%
  group_by(CHROM, Treatment, BIN_START, BIN_END, midBP) %>%
  summarise(
    n_special = n(),
    width_bp = max(BIN_END - BIN_START),
    .groups = "drop"
  )

# scale significant-SNP bars to fit on the primary y-axis
# -> secondary axis converts scaled values back to raw counts
ymax_all <- max(df_counts$n_all, na.rm = TRUE)
ymax_special <- max(df_special$n_special, na.rm = TRUE)
if (ymax_special == 0) ymax_special <- 1
scale_fac <- (ymax_all * 0.8) / ymax_special
df_special <- df_special %>%
  mutate(n_special_scaled = n_special * scale_fac)

# plot significant SNP counts as grey bars and all SNP counts as lines
p <- ggplot() +
  geom_col(
    data = df_special,
    aes(
      x = midBP,
      y = n_special_scaled,
      width = width_bp
    ),
    fill = "grey70",
    alpha = 0.7
  ) +
  geom_line(
    data = df_counts %>% arrange(Treatment, CHROM, Generation, midBP),
    aes(
      x = midBP,
      y = n_all,
      color = Generation,
      group = interaction(Treatment, CHROM, Generation)
    ),
    linewidth = 0.2
  ) +
  facet_grid(
    rows = vars(Treatment),
    cols = vars(CHROM),
    scales = "free_x",
    labeller = labeller(
      RowGroup = c(
        "Control" = "control",
        "ND" = "ND"
      )
    )
  ) +
  scale_x_continuous(labels = function(x) x / 1e6) +
  scale_y_continuous(
    name = "Number of SNPs per window",
    sec.axis = sec_axis(
      ~ . / scale_fac,
      name = "Number of significant SNPs per window"
    )
  ) +
  scale_color_manual(
    values = c(
      "Gen0" = "grey42",
      "Gen1" = "darkgreen",
      "Gen2" = "darkred",
      "Gen3" = "darkblue"
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    strip.text.x = element_text(size = 9),
    strip.text.y = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    x = "Genome position (Mbp, 25 Mbp window midpoint)",
    color = "Generation"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(linewidth = 1.2)
    )
  )


# Plot fraction of significant SNPs per 25 Mbp window
# - separate significant SNPs into control only, ND only, and both

# count all tested SNPs per chromosome and window
df_tested <- all_sites %>%
  distinct(CHROM, POS, BIN_START, BIN_END, midBP) %>%
  group_by(CHROM, BIN_START, BIN_END, midBP) %>%
  summarise(
    n_tested = n(),
    .groups = "drop"
  )

# classify significant SNPs by treatment-specific significance
sig_groups <- special_snps %>%
  mutate(
    sig_group = case_when(
      SIG_Control & !SIG_ND ~ "Control only",
      SIG_ND & !SIG_Control ~ "ND only",
      SIG_Control & SIG_ND ~ "Both",
      TRUE ~ NA
    ),
    sig_group = factor(
      sig_group,
      levels = c("Control only", "ND only", "Both")
    )
  ) %>%
  filter(!is.na(sig_group))

# count significant SNPs per significance group and window
df_sig_groups <- sig_groups %>%
  group_by(CHROM, BIN_START, BIN_END, midBP, sig_group) %>%
  summarise(
    n_significant = n(),
    .groups = "drop"
  )

# add all significance groups to every tested window
# -> windows without significant SNPs receive a count of zero
df_sig_fraction <- df_tested %>%
  crossing(
    sig_group = factor(
      c("Control only", "ND only", "Both"),
      levels = c("Control only", "ND only", "Both")
    )
  ) %>%
  left_join(
    df_sig_groups,
    by = c(
      "CHROM",
      "BIN_START",
      "BIN_END",
      "midBP",
      "sig_group"
    )
  ) %>%
  mutate(
    n_significant = replace_na(n_significant, 0),
    significant_fraction = n_significant / n_tested,
    width_bp = BIN_END - BIN_START
  )

# plot fraction of significant SNPs as stacked bars
p_sig_fraction <- ggplot(
  df_sig_fraction,
  aes(
    x = midBP,
    y = significant_fraction,
    fill = sig_group
  )
) +
  geom_col(
    aes(width = width_bp),
    position = "stack"
  ) +
  facet_grid(
    cols = vars(CHROM),
    scales = "free_x"
  ) +
  scale_x_continuous(
    labels = function(x) x / 1e6
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c(
      "Control only" = "steelblue",
      "ND only" = "coral",
      "Both" = "#7B3294"
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    x = "Genome position (Mbp, 25 Mbp window midpoint)",
    y = "Significant / tested SNPs",
    fill = "Significance"
  )



################################################################################
# Combine SNP-count and significant-SNP fraction plots
################################################################################

# remove x-axis information from upper plot because it is shown below
p_top <- p +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# combine both plots into one figure
p_combined <- p_top /
  p_sig_fraction +
  plot_layout(
    heights = c(2, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )

# save combined plot
ggsave(
  "logistic_trend/SNPs_per_25Mb_window_ND_Control.pdf",
  plot = p_combined,
  width = 14,
  height = 6.5
)
