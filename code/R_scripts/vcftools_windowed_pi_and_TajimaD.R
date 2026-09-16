library(data.table)
library(tidyverse)
library(rstatix)
library(rprojroot)

git_path <- find_root(has_dir(".git"))

################################################################################
# Plotting and pairwise Wilcoxon tests for vcftools Tajima's D and pi
################################################################################
#  - create violin/boxplots 
#  - test pairwise group differences
#  - plot chromosome-wise profiles
################################################################################


################################################################################
# Plotting helper functions
################################################################################

# add generation, treatment, and plotting-group labels for violin/boxplots
add_distribution_groups <- function(df) {
  df %>%
    mutate(
      Generation = case_when(
        Population == "gen0" ~ "Gen0",
        str_detect(Population, "^gen1") ~ "Gen1",
        str_detect(Population, "^gen2") ~ "Gen2",
        str_detect(Population, "^gen3") ~ "Gen3",
        TRUE ~ NA_character_
      ),
      Treatment = case_when(
        Population == "gen0" ~ "",
        str_detect(Population, "_ND$") ~ "ND",
        str_detect(Population, "_Control$") ~ "Control",
        TRUE ~ NA_character_
      ),
      Group = if_else(Population == "gen0", "Gen0", paste(Generation, Treatment, sep = "_")),
      Group = factor(Group, levels = group_levels),
      Treatment_fill = if_else(Generation == "Gen0", "Gen0", Treatment),
      Treatment_fill = factor(Treatment_fill, levels = c("Gen0", "ND", "Control"))
    )
}

# add generation, treatment, chromosome, and midpoint information for line plots
add_lineplot_groups <- function(df) {
  df <- df %>%
    mutate(
      Generation = case_when(
        Population == "gen0" ~ "gen0",
        str_detect(Population, "^gen1") ~ "gen1",
        str_detect(Population, "^gen2") ~ "gen2",
        str_detect(Population, "^gen3") ~ "gen3",
        TRUE ~ NA_character_
      ),
      Treatment = case_when(
        Population == "gen0" ~ "gen0",
        str_detect(Population, "_ND$") ~ "ND",
        str_detect(Population, "_Control$") ~ "Control",
        TRUE ~ NA_character_
      ),
      Window_bp = case_when(
        Window == "1Mbp" ~ 1e6,
        Window == "10Mbp" ~ 10e6,
        Window == "25Mbp" ~ 25e6,
        TRUE ~ NA_real_
      ),
      CHR = as.factor(CHROM)
    )
  if ("BIN_END" %in% names(df)) {
    df <- df %>% mutate(midBP = BIN_START + (BIN_END - BIN_START) / 2)
  } else {
    df <- df %>% mutate(midBP = BIN_START + Window_bp / 2)
  }
  df
}

# prepare chromosome-wise line-plot data for one window size.
# - Gen0 is duplicated into both treatment rows 
prepare_lineplot_data <- function(df, window_label, value_col) {
  df_gen <- df %>% add_lineplot_groups()
  df_gen %>%
    filter(Window == window_label) %>%
    bind_rows(
      df_gen %>% filter(Window == window_label, Treatment == "gen0") %>%
        mutate(Treatment = "Control", Population = "gen0_in_Control"),
      df_gen %>% filter(Window == window_label, Treatment == "gen0") %>%
        mutate(Treatment = "ND", Population = "gen0_in_ND")
    ) %>%
    filter(Treatment != "gen0", !is.na(.data[[value_col]])) %>%
    mutate(
      PlotRow = factor(Treatment, levels = c("Control", "ND")),
      Generation = factor(Generation, levels = c("gen0", "gen1", "gen2", "gen3"))
    )
}

# create and save violin/boxplots for one statistic and window size
plot_violin_box <- function(df, window_label, value_col, y_label, output_prefix) {
  # keep selected window size and remove missing values
  plot_df <- df %>%
    add_distribution_groups() %>%
    filter(Window == window_label, !is.na(.data[[value_col]]), !is.na(Group))
  
  # calculate generation-specific median lines using Gen0 and ND groups
  median_lines <- plot_df %>%
    filter(Group == "Gen0" | Treatment == "ND") %>%
    group_by(Generation) %>%
    summarise(median_value = median(.data[[value_col]], na.rm = TRUE), .groups = "drop")

  # positions for generation labels below the x-axis
  generation_labels <- data.frame(
    x = c(1, 2.5, 4.5, 6.5),
    y = -Inf,
    Generation = c("Gen0", "Gen1", "Gen2", "Gen3"),
    label = c("Gen0", "Gen1", "Gen2", "Gen3")
  )

  p <- ggplot(plot_df, aes(x = Group, y = .data[[value_col]], fill = Treatment_fill)) +
    geom_violin(trim = FALSE, alpha = 0.7, linewidth = 0.4) +
    geom_hline(
      data = median_lines,
      aes(yintercept = median_value, color = Generation),
      linewidth = 0.7,
      inherit.aes = FALSE
    ) +
    geom_boxplot(width = 0.12, outlier.alpha = 0.25, outlier.size = 0.8, linewidth = 0.35, fill = "white") +
    scale_x_discrete(labels = c(
      "Gen0" = "",
      "Gen1_ND" = "ND", "Gen1_Control" = "control",
      "Gen2_ND" = "ND", "Gen2_Control" = "control",
      "Gen3_ND" = "ND", "Gen3_Control" = "control"
    )) +
    scale_fill_manual(values = c("Gen0" = "grey42", 
                                 "ND" = "coral", 
                                 "Control" = "steelblue"), 
                      guide = "none") +
    scale_color_manual(values = c("Gen0" = "grey42", 
                                  "Gen1" = "darkgreen", 
                                  "Gen2" = "darkred", 
                                  "Gen3" = "darkblue"), 
                      guide = "none") +
    geom_text(data = generation_labels, 
              aes(x = x, y = y, label = label, color = Generation), 
              vjust = 3.2, size = 5, inherit.aes = FALSE) +
    labs(x = NULL, y = y_label, title = paste("Window size:", window_label)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 13),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = grid::unit(2, "pt"),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.margin = margin(5.5, 5.5, 35, 5.5)
    )

  ggsave(
    paste0(output_prefix, "_", window_label, "_violin_boxplot_NDvsControl.pdf"), 
    p, width = 5, height = 7, units = "in")
  p
}


# plot chromosome-wise profiles separated by treatment and colored by generation
plot_by_treatment <- function(df, value_col, y_label, title = "") {
  ggplot(df, aes(
    x = midBP,
    y = .data[[value_col]],
    colour = Generation,
    group = interaction(PlotRow, Generation)
  )) +
    geom_line(linewidth = 0.35, na.rm = TRUE) +
    facet_grid(
      rows = vars(PlotRow),
      cols = vars(CHR),
      scales = "free_x",
      labeller = labeller(
        PlotRow = c(
          "Control" = "control",
          "ND" = "ND"
        )
      )
    ) +
    scale_x_continuous(labels = function(x) x / 1e6) +
    scale_colour_manual(
      values = c(
        gen0 = "grey42",
        gen1 = "darkgreen",
        gen2 = "darkred",
        gen3 = "darkblue"
      ),
      labels = c(
        gen0 = "Gen0",
        gen1 = "Gen1",
        gen2 = "Gen2",
        gen3 = "Gen3"
      ),
      drop = FALSE
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(linewidth = 1.5),
        ncol = 1
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.direction = "vertical",
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Genome position (Mbp, window midpoint)",
      y = y_label,
      colour = "Generation",
      title = title
    )
}


# calculate group summaries and run pairwise Wilcoxon tests with BH correction
run_pairwise_wilcox <- function(df, window_label, value_col, output_prefix) {
  test_df <- df %>%
    add_distribution_groups() %>%
    filter(Window == window_label, !is.na(.data[[value_col]]), !is.na(Group))
  summary_df <- test_df %>%
    group_by(Group) %>%
    summarise(
      n_windows = n(),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
  pairwise_df <- test_df %>%
    pairwise_wilcox_test(as.formula(paste(value_col, "~ Group")), 
                         p.adjust.method = "BH")

  write.csv(summary_df, 
            paste0(output_prefix, "_", window_label, "_summary_by_group.csv"), 
            row.names = FALSE)
  write.csv(pairwise_df, 
            paste0(output_prefix, "_", window_label, "_pairwise_wilcox_tests.csv"), 
            row.names = FALSE)
  list(summary = summary_df, pairwise = pairwise_df)
}



################################################################################
# General settings
################################################################################


tajima_dir <- paste0(git_path,"/general_analysis/vcftools/tajimaD")
pi_dir     <- paste0(git_path,"/general_analysis/vcftools/windowed_pi")
window_label <- "25Mbp"

# order of groups in violin/boxplots
group_levels <- c(
  "Gen0",
  "Gen1_ND", "Gen1_Control",
  "Gen2_ND", "Gen2_Control",
  "Gen3_ND", "Gen3_Control"
)



################################################################################
# Tajima's D: read, plot, test, chromosome profile
################################################################################

setwd(tajima_dir)

# read all Tajima's D files and extract population/window size from file names
tajima_files <- list.files(pattern = "\\.Tajima\\.D$", full.names = TRUE)
tajima_df <- map_dfr(tajima_files, function(f) {
  file_name <- basename(f)
  info <- str_match(file_name, "tajimaD_(.*)_(1Mbp|10Mbp|25Mbp)\\.Tajima\\.D$")
  fread(f) %>% mutate(Population = info[, 2], Window = info[, 3])
})

# create violin/boxplot for Tajima's D
tajima_violin <- plot_violin_box(
  df = tajima_df,
  window_label = window_label,
  value_col = "TajimaD",
  y_label = "Tajima's D",
  output_prefix = "tajimaD"
)

# run pairwise Wilcoxon tests between generation-treatment groups
tajima_wilcox <- run_pairwise_wilcox(
  df = tajima_df,
  window_label = window_label,
  value_col = "TajimaD",
  output_prefix = "tajimaD"
)

# prepare and save chromosome-wise Tajima's D profile
plot_df_tajima <- prepare_lineplot_data(
  df = tajima_df,
  window_label = window_label,
  value_col = "TajimaD"
)
p_tajima_25mbp <- plot_by_treatment(
  df = plot_df_tajima,
  value_col = "TajimaD",
  y_label = "Tajima's D",
  title = "Tajima's D per 25 Mbp window"
)
ggsave("tajimaD_25Mbp_by_generation_treatment.pdf", 
       p_tajima_25mbp, 
       width = 14, 
       height = 4.4, 
       units = "in")

################################################################################
# Nucleotide diversity: read, plot, test, chromosome profile
################################################################################

setwd(pi_dir)

# read all windowed pi files and extract population/window size from file names
pi_files <- list.files(pattern = "\\.windowed\\.pi$", full.names = TRUE)
pi_df <- map_dfr(pi_files, function(f) {
  file_name <- basename(f)
  pop <- str_remove(file_name, "_window_pi_.*$")
  window <- str_extract(file_name, "(1Mbp|10Mbp|25Mbp)")

  fread(f) %>% mutate(Population = pop, Window = window)
})

# create violin/boxplot for nucleotide diversity
pi_violin <- plot_violin_box(
  df = pi_df,
  window_label = window_label,
  value_col = "PI",
  y_label = expression(Nucleotide~diversity~pi),
  output_prefix = "pi"
)

# run pairwise Wilcoxon tests between generation-treatment groups
pi_wilcox <- run_pairwise_wilcox(
  df = pi_df,
  window_label = window_label,
  value_col = "PI",
  output_prefix = "pi"
)

# prepare and save chromosome-wise nucleotide diversity profile
plot_df_pi <- prepare_lineplot_data(
  df = pi_df,
  window_label = window_label,
  value_col = "PI"
)
p_pi_25mbp <- plot_by_treatment(
  df = plot_df_pi,
  value_col = "PI",
  y_label = expression(Nucleotide~diversity~pi),
  title = "Nucleotide diversity per 25 Mbp window"
)
ggsave("windowed_pi_25Mbp_by_generation_treatment.pdf", 
       p_pi_25mbp, 
       width = 14, 
       height = 4.4, 
       units = "in")
