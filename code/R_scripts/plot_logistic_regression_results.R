library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(purrr)
library(broom)
library(gtools)
library(VennDiagram)
library(scales)
library(cowplot)
library(UpSetR)
library(grid)
library(gridExtra)
library(ggnewscale)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path,"/general_analysis/vcftools/freq"))


################################################################################
# Plot and summarize significant SNPs from logistic modelling
################################################################################
# - read empirical p-values from the logistic/Wright-Fisher analysis
# - create QQ plots
# - identify SNPs significant in ND and/or Control
# - plot their allele-frequency trajectories
# - export summary tables + overview plots
################################################################################


################################################################################
# Helper function
################################################################################

# classify Gen1-to-Gen3 allele-frequency change direction per treatment
summarise_snp_set <- function(df, set_name) {
  d <- df %>%
    mutate(
      dir_Control = case_when(
        is.na(Delta_gen1_gen3_Control) ~ NA,
        Delta_gen1_gen3_Control > 0 ~ "positive",
        Delta_gen1_gen3_Control < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      dir_ND = case_when(
        is.na(Delta_gen1_gen3_ND) ~ NA,
        Delta_gen1_gen3_ND > 0 ~ "positive",
        Delta_gen1_gen3_ND < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      direction_relation = case_when(
        is.na(dir_Control) | is.na(dir_ND) ~ NA,
        dir_Control == "zero" | dir_ND == "zero" ~ "zero_in_at_least_one",
        dir_Control == dir_ND ~ "same_direction",
        dir_Control != dir_ND ~ "different_direction"
      )
    )
  
  # count direction classes and summarize delta distributions
  tibble(
    set = set_name,
    n_snps = nrow(d),
    n_complete_both_deltas = sum(!is.na(d$Delta_gen1_gen3_Control) & !is.na(d$Delta_gen1_gen3_ND)),
    n_zero_in_at_least_one = sum(d$direction_relation == "zero_in_at_least_one", na.rm = TRUE),
    # direction comparison between Control and ND
    n_same_direction = sum(d$direction_relation == "same_direction", na.rm = TRUE),
    n_different_direction = sum(d$direction_relation == "different_direction", na.rm = TRUE),
    n_pos_pos = sum(d$dir_Control == "positive" & d$dir_ND == "positive", na.rm = TRUE),
    n_neg_neg = sum(d$dir_Control == "negative" & d$dir_ND == "negative", na.rm = TRUE),
    n_pos_neg = sum(d$dir_Control == "positive" & d$dir_ND == "negative", na.rm = TRUE),
    n_neg_pos = sum(d$dir_Control == "negative" & d$dir_ND == "positive", na.rm = TRUE),
    # positive/negative allele-frequency changes per treatment
    n_positive_Control = sum(d$Delta_gen1_gen3_Control > 0, na.rm = TRUE),
    n_negative_Control = sum(d$Delta_gen1_gen3_Control < 0, na.rm = TRUE),
    n_positive_ND = sum(d$Delta_gen1_gen3_ND > 0, na.rm = TRUE),
    n_negative_ND = sum(d$Delta_gen1_gen3_ND < 0, na.rm = TRUE),
    # combined across both treatments 
    # -> each treatment-specific delta is counted separately
    n_positive_total_treatments = 
      sum(d$Delta_gen1_gen3_Control > 0, na.rm = TRUE) +
      sum(d$Delta_gen1_gen3_ND  > 0, na.rm = TRUE),
    n_negative_total_treatments = 
      sum(d$Delta_gen1_gen3_Control < 0, na.rm = TRUE) +
      sum(d$Delta_gen1_gen3_ND  < 0, na.rm = TRUE),
    # delta distribution in Control
    mean_delta_Control   = mean(d$Delta_gen1_gen3_Control, na.rm = TRUE),
    median_delta_Control = median(d$Delta_gen1_gen3_Control, na.rm = TRUE),
    sd_delta_Control     = sd(d$Delta_gen1_gen3_Control, na.rm = TRUE),
    min_delta_Control    = min(d$Delta_gen1_gen3_Control, na.rm = TRUE),
    q25_delta_Control    = quantile(d$Delta_gen1_gen3_Control, 0.25, na.rm = TRUE),
    q75_delta_Control    = quantile(d$Delta_gen1_gen3_Control, 0.75, na.rm = TRUE),
    max_delta_Control    = max(d$Delta_gen1_gen3_Control, na.rm = TRUE),
    # delta distribution in ND
    mean_delta_ND   = mean(d$Delta_gen1_gen3_ND, na.rm = TRUE),
    median_delta_ND = median(d$Delta_gen1_gen3_ND, na.rm = TRUE),
    sd_delta_ND     = sd(d$Delta_gen1_gen3_ND, na.rm = TRUE),
    min_delta_ND    = min(d$Delta_gen1_gen3_ND, na.rm = TRUE),
    q25_delta_ND    = quantile(d$Delta_gen1_gen3_ND, 0.25, na.rm = TRUE),
    q75_delta_ND    = quantile(d$Delta_gen1_gen3_ND, 0.75, na.rm = TRUE),
    max_delta_ND    = max(d$Delta_gen1_gen3_ND, na.rm = TRUE)
  )
}


################################################################################
# Plotting for sig. SNPs in logistic setting
################################################################################

# read refined empirical p-value results for both treatments
emp_control <- fread(
  "logistic_trend/LD_empirical_pvalues_Control_refined_full_0.9miss.csv.gz"
) %>% as.data.frame()
emp_nd <- fread(
  "logistic_trend/LD_empirical_pvalues_ND_refined_full_0.9miss.csv.gz"
) %>% as.data.frame()


################################################################################
# QQ plot of empirical p-values
################################################################################

# combine ND and Control empirical p-values + calculate expected/observed -log10(p)
qq_nd <- emp_nd %>%
  as.data.frame() %>%
  mutate(track = "ND")
qq_control <- emp_control %>%
  as.data.frame() %>%
  mutate(track = "Control")
qq_all <- bind_rows(qq_nd, qq_control) %>%
  filter(!is.na(p_empirical_final)) %>%
  filter(p_empirical_final > 0, p_empirical_final <= 1) %>%
  group_by(track) %>%
  arrange(p_empirical_final, .by_group = TRUE) %>%
  mutate(
    exp_p    = (row_number() - 0.5) / n(),
    obs_logp = -log10(p_empirical_final),
    exp_logp = -log10(exp_p)
  ) %>%
  ungroup()

# plot p-value distributions against uniform expectation
p <- ggplot(qq_all, aes(x = exp_logp, y = obs_logp, color = track)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "ND" = "coral",
      "Control" = "steelblue"
    ),
    labels = c(
      "ND" = "ND",
      "Control" = "control"
    )
  ) +
  theme_bw() +
  labs(
    x = "Expected -log10(p)",
    y = "Observed -log10(p)",
    title = "",
    color = NULL
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    )
  ) +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "grey70"),
    legend.key = element_rect(fill = "white")
  )
ggsave("logistic_trend/LD_qq_empirical_both_0.9miss.pdf", p,width = 4,height = 4,units = "in")


################################################################################
# Save significant SNPs with frequencies and deltas
################################################################################

key_cols <- c("CHROM", "POS", "ID", "REF_BASE", "ALT_BASE")
alpha <- 0.05

# keep SNPs significant after BH correction in each treatment line
sig_control_key <- emp_control %>%
  filter(!is.na(p_adj_BH), p_adj_BH < alpha) %>%
  select(all_of(key_cols)) %>%
  distinct() %>%
  mutate(SIG_Control = TRUE)
sig_nd_key <- emp_nd %>%
  filter(!is.na(p_adj_BH), p_adj_BH < alpha) %>%
  select(all_of(key_cols)) %>%
  distinct() %>%
  mutate(SIG_ND = TRUE)

# read merged allele-frequency table
global <- fread("LD_all_pop_existing_sites_miss0.9.csv.gz") %>%
  as.data.frame()

# combine significance status with allele counts + calculate exact frequencies/deltas
snps_sig_table <- global %>%
  select(
    all_of(key_cols),
    N_CHROM_gen0, ALT_AC_gen0,
    N_CHROM_gen1_Control, ALT_AC_gen1_Control,
    N_CHROM_gen2_Control, ALT_AC_gen2_Control,
    N_CHROM_gen3_Control, ALT_AC_gen3_Control,
    N_CHROM_gen1_ND,  ALT_AC_gen1_ND,
    N_CHROM_gen2_ND,  ALT_AC_gen2_ND,
    N_CHROM_gen3_ND,  ALT_AC_gen3_ND
  ) %>%
  left_join(sig_control_key, by = key_cols) %>%
  left_join(sig_nd_key, by = key_cols) %>%
  # add observed slopes for both treatments
  left_join(
    emp_control %>%
      select(all_of(key_cols), slope_Control = slope_obs),
    by = key_cols
  ) %>%
  left_join(
    emp_nd %>%
      select(all_of(key_cols), slope_ND = slope_obs),
    by = key_cols
  ) %>%
  # mark whether each SNP is significant in Control, ND, or both
  mutate(
    SIG_Control = replace_na(SIG_Control, FALSE),
    SIG_ND  = replace_na(SIG_ND, FALSE),
    SIG_BOTH = SIG_Control & SIG_ND
  ) %>%
  filter(SIG_Control | SIG_ND) %>%
  # add slope diffrences
  mutate(
    Delta_slope_ND_Control = slope_ND - slope_Control,
    Abs_Delta_slope_ND_Control = abs(Delta_slope_ND_Control)
  ) %>%
  # calculate exact ALT allele frequencies from allele counts (not the vcftools rounded ones)
  mutate(
    ALT_FREQ_gen0     = ALT_AC_gen0 / N_CHROM_gen0,
    ALT_FREQ_gen1_Control = ALT_AC_gen1_Control / N_CHROM_gen1_Control,
    ALT_FREQ_gen2_Control = ALT_AC_gen2_Control / N_CHROM_gen2_Control,
    ALT_FREQ_gen3_Control = ALT_AC_gen3_Control / N_CHROM_gen3_Control,
    ALT_FREQ_gen1_ND  = ALT_AC_gen1_ND / N_CHROM_gen1_ND,
    ALT_FREQ_gen2_ND  = ALT_AC_gen2_ND / N_CHROM_gen2_ND,
    ALT_FREQ_gen3_ND  = ALT_AC_gen3_ND / N_CHROM_gen3_ND
  ) %>%
  # calculate allele-frequency changes between generations
  mutate(
    Delta_gen0_gen3_Control = ALT_FREQ_gen3_Control - ALT_FREQ_gen0,
    Delta_gen0_gen3_ND  = ALT_FREQ_gen3_ND  - ALT_FREQ_gen0,
    Delta_gen1_gen3_Control = ALT_FREQ_gen3_Control - ALT_FREQ_gen1_Control,
    Delta_gen1_gen3_ND  = ALT_FREQ_gen3_ND  - ALT_FREQ_gen1_ND,
    Delta_gen0_gen1_Control = ALT_FREQ_gen1_Control - ALT_FREQ_gen0,
    Delta_gen0_gen1_ND  = ALT_FREQ_gen1_ND  - ALT_FREQ_gen0,
    Delta_gen1_gen2_Control = ALT_FREQ_gen2_Control - ALT_FREQ_gen1_Control,
    Delta_gen1_gen2_ND  = ALT_FREQ_gen2_ND  - ALT_FREQ_gen1_ND,
    Delta_gen2_gen3_Control = ALT_FREQ_gen3_Control - ALT_FREQ_gen2_Control,
    Delta_gen2_gen3_ND  = ALT_FREQ_gen3_ND  - ALT_FREQ_gen2_ND
  ) %>%
  # calculate treatment differences within generations
  mutate(
    Delta_gen1_Control_ND = ALT_FREQ_gen1_Control - ALT_FREQ_gen1_ND,
    Delta_gen2_Control_ND = ALT_FREQ_gen2_Control - ALT_FREQ_gen2_ND,
    Delta_gen3_Control_ND = ALT_FREQ_gen3_Control - ALT_FREQ_gen3_ND
  ) %>%
  # compare generation-to-generation changes between Control and ND
  mutate(
    Delta_diff_Control_vs_ND_gen0_gen3 = Delta_gen0_gen3_Control - Delta_gen0_gen3_ND,
    Delta_diff_Control_vs_ND_gen1_gen3 = Delta_gen1_gen3_Control - Delta_gen1_gen3_ND,
    Delta_diff_Control_vs_ND_gen0_gen1 = Delta_gen0_gen1_Control - Delta_gen0_gen1_ND,
    Delta_diff_Control_vs_ND_gen1_gen2 = Delta_gen1_gen2_Control - Delta_gen1_gen2_ND,
    Delta_diff_Control_vs_ND_gen2_gen3 = Delta_gen2_gen3_Control - Delta_gen2_gen3_ND
  )
fwrite(
  snps_sig_table,
  file = "logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz",
  compress = "gzip"
)


################################################################################
# Plot allele-frequency trajectories of significant SNPs
################################################################################

# read significant SNP table
snps_sig_table <- fread(
  "logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz"
) %>%
  as.data.frame()

# prepare long-format allele-frequency trajectories for Control and ND
plot_freq <- bind_rows(
  snps_sig_table %>%
    filter(SIG_Control) %>%
    select(
      all_of(key_cols),
      SIG_Control, SIG_ND, SIG_BOTH,
      ALT_FREQ_gen0,
      ALT_FREQ_gen1_Control, ALT_FREQ_gen2_Control, ALT_FREQ_gen3_Control
    ) %>%
    pivot_longer(
      cols = starts_with("ALT_FREQ_"),
      names_to = "generation_raw",
      values_to = "ALT_FREQ"
    ) %>%
    mutate(
      treatment = "Control",
      generation = case_when(
        generation_raw == "ALT_FREQ_gen0"     ~ "Gen0",
        generation_raw == "ALT_FREQ_gen1_Control" ~ "Gen1",
        generation_raw == "ALT_FREQ_gen2_Control" ~ "Gen2",
        generation_raw == "ALT_FREQ_gen3_Control" ~ "Gen3",
        TRUE ~ NA_character_
      )
    ),
  snps_sig_table %>%
    filter(SIG_ND) %>%
    select(
      all_of(key_cols),
      SIG_Control, SIG_ND, SIG_BOTH,
      ALT_FREQ_gen0, ALT_FREQ_gen1_ND, ALT_FREQ_gen2_ND, ALT_FREQ_gen3_ND
    ) %>%
    pivot_longer(
      cols = starts_with("ALT_FREQ_"),
      names_to = "generation_raw",
      values_to = "ALT_FREQ"
    ) %>%
    mutate(
      treatment = "ND",
      generation = case_when(
        generation_raw == "ALT_FREQ_gen0"    ~ "Gen0",
        generation_raw == "ALT_FREQ_gen1_ND" ~ "Gen1",
        generation_raw == "ALT_FREQ_gen2_ND" ~ "Gen2",
        generation_raw == "ALT_FREQ_gen3_ND" ~ "Gen3",
        TRUE ~ NA_character_
      )
    )
) %>%
  filter(!is.na(ALT_FREQ), !is.na(POS), !is.na(CHROM),!is.na(generation)) %>%
  mutate(
    sig_class = case_when(
      SIG_ND & !SIG_Control ~ "ND only",
      SIG_Control & !SIG_ND ~ "Control only",
      SIG_ND & SIG_Control ~ "Both",
      TRUE ~ NA_character_
    ),
    generation = factor(generation, levels = c("Gen0", "Gen1", "Gen2", "Gen3")),
    treatment = factor(treatment, levels = c("ND", "Control")),
    CHROM = factor(CHROM, levels = mixedsort(unique(CHROM)))
  )

# SNPs significant in only one treatment are highlighted
plot_freq_only_one <- plot_freq %>%
  filter(sig_class %in% c("ND only", "Control only"))

# Plot allele-frequency trajectories by genomic position
p_freq <- ggplot(
  plot_freq,
  aes(x = POS, y = ALT_FREQ)
) +
  # highlight trajectories significant only in ND/Control
  geom_line(
    data = plot_freq_only_one,
    aes(
      group = interaction(CHROM, POS, ID, treatment),
      color = sig_class
    ),
    alpha = 0.7,
    linewidth = 2
  ) +
  scale_color_manual(
    name = "Significance:",
    values = c(
      "ND only" = "coral",
      "Control only" = "steelblue"
    ),
    labels = c(
      "ND only" = "ND only",
      "Control only" = "control only"
    )
  ) +
  ggnewscale::new_scale_color() +
  
  # add all significant trajectories as thin black background lines
  geom_line(
    aes(group = interaction(CHROM, POS, ID, treatment)),
    color = "black",
    linewidth = 0.2,
    alpha = 0.45
  ) +
  geom_point(
    aes(color = generation, shape = treatment),
    size = 1.8,
    alpha = 0.9
  ) +
  facet_wrap(~ CHROM, scales = "free_x", ncol = 4) +
  scale_x_continuous(labels = function(x) x / 1e6) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    name = "Generation:",
    values = c(
      "Gen0" = "grey42",
      "Gen1" = "darkgreen",
      "Gen2" = "darkred",
      "Gen3" = "darkblue"
    )
  ) +
  scale_shape_manual(
    name = "Treatment:",
    values = c(
      "ND" = 16,
      "Control" = 17
    ),
    labels = c(
      "ND" = "ND",
      "Control" = "control"
    )
  ) +
  labs(
    x = "Genomic position (Mbp)",
    y = "ALT allele frequency",
    title = "Allele-frequency trajectories of significant SNPs"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    
    legend.background = element_rect(
      fill = alpha("white", 0.85),
      color = "grey70"
    ),
    legend.box.background = element_rect(
      fill = alpha("white", 0.0),
      color = NA
    )
  ) +
  guides(
    color = guide_legend(
      order = 2,
      override.aes = list(size = 3, alpha = 1)
    ),
    shape = guide_legend(
      order = 3,
      override.aes = list(size = 3, alpha = 1)
    ),
    color_ggnewscale_1 = guide_legend(
      order = 1,
      override.aes = list(linewidth = 2, alpha = 0.7)
    )
  )

ggsave(
  filename = "logistic_trend/LD_significant_ALTfreq_trajectories_by_chr.pdf",
  plot = p_freq,
  width = 15,
  height = 9
)


################################################################################
# Summarize significant SNP sets
################################################################################


snps_sig_table <- fread(
  "logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz"
)
summary_table <- bind_rows(
  summarise_snp_set(snps_sig_table, "all_snps"),
  summarise_snp_set(filter(snps_sig_table, SIG_BOTH), "sig_BOTH"),
  summarise_snp_set(filter(snps_sig_table, SIG_Control & !SIG_ND), "only_sig_Control"),
  summarise_snp_set(filter(snps_sig_table, SIG_ND & !SIG_Control), "only_sig_ND")
)
fwrite(
  summary_table,
  file = "logistic_trend/significant_SNPs_summary_overview.csv"
)


################################################################################
# UpSet plot of significant SNP direction classes
################################################################################

# convert direction classes to binary membership columns for UpSetR
upset_df_upsetr <- snps_sig_table %>%
  mutate(
    `control positive` = as.integer(SIG_Control & !is.na(Delta_gen1_gen3_Control) & Delta_gen1_gen3_Control > 0),
    `control negative` = as.integer(SIG_Control & !is.na(Delta_gen1_gen3_Control) & Delta_gen1_gen3_Control < 0),
    `ND positive` = as.integer(SIG_ND & !is.na(Delta_gen1_gen3_ND) & Delta_gen1_gen3_ND > 0),
    `ND negative` = as.integer(SIG_ND & !is.na(Delta_gen1_gen3_ND) & Delta_gen1_gen3_ND < 0)
  ) %>%
  distinct(CHROM, POS, ID, REF_BASE, ALT_BASE, .keep_all = TRUE) %>%
  select(
    `control positive`,
    `control negative`,
    `ND positive`,
    `ND negative`
  ) %>%
  as.data.frame()

# plot upset plot with colored side bars
p_upset <- upset(
  upset_df_upsetr,
  sets = c(
    "control positive",
    "control negative",
    "ND positive",
    "ND negative"
  ),
  sets.bar.color = c(
    "steelblue",
    "steelblue",
    "coral",
    "coral"
  ),
  order.by = "freq",
  keep.order = TRUE,
  mainbar.y.label = "Number of SNPs",
  sets.x.label = "Number of SNPs",
  set_size.show = TRUE,
  set_size.scale_max =
    max(colSums(upset_df_upsetr)) * 1.20,
  text.scale = c(
    1.2,
    1.2,
    1.2,
    1.2,
    1.4,
    1.2
  )
)
pdf(
  "logistic_trend/LD_significant_SNPs_upset.pdf",
  width = 5,
  height = 4
)
print(p_upset)
dev.off()


################################################################################
# Timeline plots of significant SNP allele-frequency trajectories
################################################################################


snps_sig_table <- fread(
  "logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz"
) %>%
  mutate(
    snp = paste(CHROM, POS, sep = "_")
  )

# prepare long-format allele-frequency table for timeline plots
# -> Gen0 is duplicated for Control and ND as an independent visual reference
plot_df <- bind_rows(
  snps_sig_table %>%
    transmute(
      CHROM, POS, ID, snp,
      generation = "Gen0",
      Control = ALT_FREQ_gen0,
      ND = ALT_FREQ_gen0
    ) %>%
    pivot_longer(
      cols = c(Control, ND),
      names_to = "treatment",
      values_to = "allele_freq"
    ),
  snps_sig_table %>%
    transmute(
      CHROM, POS, ID, snp,
      Gen1_Control = ALT_FREQ_gen1_Control, 
      Gen2_Control = ALT_FREQ_gen2_Control,
      Gen3_Control = ALT_FREQ_gen3_Control,
      Gen1_ND = ALT_FREQ_gen1_ND,
      Gen2_ND = ALT_FREQ_gen2_ND,
      Gen3_ND = ALT_FREQ_gen3_ND
    ) %>%
    pivot_longer(
      cols = starts_with("Gen"),
      names_to = c("generation", "treatment"),
      names_sep = "_",
      values_to = "allele_freq"
    )
) %>%
  mutate(
    generation = factor(generation, levels = c("Gen0", "Gen1", "Gen2", "Gen3")),
    treatment = factor(treatment, levels = c("Control", "ND"))
  )

# plot allele-frequency trajectories across generations
make_freq_plot <- function(df, title = NULL) {
  ggplot(
    df,
    aes(
      generation,
      allele_freq,
      group = interaction(snp, treatment),
      color = treatment
    )
  ) +
    geom_line(alpha = 0.3, linewidth = 0.75) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
    scale_color_manual(
      values = c("Control" = "steelblue", "ND" = "coral"),
      labels = c("ND" = "ND", "Control" = "control")
      ) +
    scale_x_discrete(
      expand = c(0, 0),
      labels = c("Gen0" = "0", "Gen1" = "1", "Gen2" = "2", "Gen3" = "3")
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = title,
      x = "Generation",
      y = "ALT allele frequency",
      color = "Treatment"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 15, hjust = 0.7),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
}

# SNPs significant only in Control
p_control_only <- make_freq_plot(
  plot_df %>%
    filter(
      treatment == "Control",
      snp %in% snps_sig_table$snp[
        snps_sig_table$SIG_Control & !snps_sig_table$SIG_ND
      ]
    ),
  "Significant in control only"
) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# SNPs significant only in ND
p_nd_only <- make_freq_plot(
  plot_df %>%
    filter(
      treatment == "ND",
      snp %in% snps_sig_table$snp[
        snps_sig_table$SIG_ND & !snps_sig_table$SIG_Control
      ]
    ),
  "Significant in ND only"
) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


# SNPs significant in both treatments - Control trajectories
p_both_control <- make_freq_plot(
  plot_df %>%
    filter(
      treatment == "Control",
      snp %in% snps_sig_table$snp[snps_sig_table$SIG_BOTH]
    ),
  "Significant in both: control"
) +
  theme(
    legend.position = "none"
  )


# SNPs significant in both treatments - ND trajectories
p_both_nd <- make_freq_plot(
  plot_df %>%
    filter(
      treatment == "ND",
      snp %in% snps_sig_table$snp[snps_sig_table$SIG_BOTH]
    ),
  "Significant in both: ND"
) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


# Combine timeline panels
combined_plot <- (
  (p_control_only + p_nd_only) /
    (p_both_control + p_both_nd)
)

ggsave(
  "logistic_trend/LD_SNPs_sig_altfreqs_per_timeline.pdf",
  plot = combined_plot,
  width = 7.5,
  height = 8
)




################################################################################
# Summary of treatment-specific allele-frequency changes in significant SNPs
################################################################################

# load table containing SNP significance classes and allele-frequency changes
snps_sig_table <- fread(
  "logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz"
) %>%
  mutate(
    snp = paste(CHROM, POS, sep = "_")
  )

# summarize Gen1--Gen3 allele-frequency changes for significant SNP groups
#
#SNPs significant in only one treatment:
# - compare absolute allele-frequency change in significant treatment
#   with change of the same SNP in other treatment
summary_table <- snps_sig_table %>%
  mutate(
    group = case_when(
      SIG_Control & !SIG_ND ~ "control only",
      SIG_ND & !SIG_Control ~ "ND only",
      SIG_BOTH ~ "both"
    ),
    
    # absolute change in the treatment in which the SNP was significant
    abs_delta_significant = case_when(
      group == "ND only" ~ abs(Delta_gen1_gen3_ND),
      group == "control only" ~ abs(Delta_gen1_gen3_Control),
      TRUE ~ NA
    ),
    
    # absolute change of the same SNP in the other treatment
    abs_delta_other = case_when(
      group == "ND only" ~ abs(Delta_gen1_gen3_Control),
      group == "control only" ~ abs(Delta_gen1_gen3_ND),
      TRUE ~ NA
    ),
    
    # relative reduction in the other treatment
    # calculated only when the allele-frequency change is larger
    # in the treatment in which the SNP was significant
    relative_reduction_in_other_treatment = if_else(
      !is.na(abs_delta_significant) & !is.na(abs_delta_other) &
        abs_delta_significant > 0 & abs_delta_significant > abs_delta_other,
      (abs_delta_significant - abs_delta_other) / abs_delta_significant,
      NA
    )
  ) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    median_abs_delta_control =
      median(abs(Delta_gen1_gen3_Control), na.rm = TRUE),
    q25_abs_delta_control =
      quantile(abs(Delta_gen1_gen3_Control), 0.25, na.rm = TRUE),
    q75_abs_delta_control =
      quantile(abs(Delta_gen1_gen3_Control), 0.75, na.rm = TRUE),
    median_abs_delta_ND =
      median(abs(Delta_gen1_gen3_ND), na.rm = TRUE),
    q25_abs_delta_ND =
      quantile(abs(Delta_gen1_gen3_ND), 0.25, na.rm = TRUE),
    q75_abs_delta_ND =
      quantile(abs(Delta_gen1_gen3_ND), 0.75, na.rm = TRUE),
    median_abs_treatment_diff =
      median(abs(Delta_diff_Control_vs_ND_gen1_gen3), na.rm = TRUE),
    q25_abs_treatment_diff =
      quantile(abs(Delta_diff_Control_vs_ND_gen1_gen3), 0.25, na.rm = TRUE),
    q75_abs_treatment_diff =
      quantile(abs(Delta_diff_Control_vs_ND_gen1_gen3), 0.75, na.rm = TRUE),
    median_relative_reduction_in_other_treatment =
      median(relative_reduction_in_other_treatment, na.rm = TRUE),
    n_relative_reduction =
      sum(!is.na(relative_reduction_in_other_treatment)),
    proportion_stronger_in_significant_treatment =
      mean(abs_delta_significant > abs_delta_other, na.rm = TRUE),
    q25_relative_reduction =
      quantile(relative_reduction_in_other_treatment, 0.25, na.rm = TRUE),
    q75_relative_reduction =
      quantile(relative_reduction_in_other_treatment, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

fwrite(
  summary_table,
  "logistic_trend/significant_SNPs_treatment_delta_summary.csv"
)



# Summary of significant SNP groups per chromosome
# classify significant SNPs into treatment-specific significance groups
chromosome_summary <- snps_sig_table %>%
  mutate(
    group = case_when(
      SIG_Control & !SIG_ND ~ "control only",
      SIG_ND & !SIG_Control ~ "ND only",
      SIG_BOTH ~ "both"
    )
  ) %>%
  # count SNPs of each significance group per chromosome
  count(CHROM, group, name = "n") %>%
  group_by(group) %>%
  mutate(
    percent = 100 * n / sum(n)
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from = group,
    values_from = c(n, percent),
    values_fill = 0,
    names_glue = "{group}_{.value}"
  ) %>%
  arrange(mixedorder(CHROM))

fwrite(
  chromosome_summary,
  "logistic_trend/significant_SNPs_per_chromosome.csv"
)
