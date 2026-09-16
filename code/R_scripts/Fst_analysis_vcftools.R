library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/vcftools/Fst"))

################################################################################
# Plot windowed FST across generations and treatments
################################################################################
# read vcftools windowed Weir-FST files for 25 Mbp windows
# create chromosome-wise FST line plots for:
#   - each generation compared to Gen0
#   - consecutive generation comparisons
#   - ND versus Control within generations
################################################################################


################################################################################
# Plotting helper functions
################################################################################

# convert chromosome labels to integer codes for sorting and plotting
normalize_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("(?i)^chr([1-7])r$", "\\1", x, perl = TRUE)
  x <- sub("(?i)^chrun$", "8", x, perl = TRUE)
  x <- sub("(?i)^chrb$", "9", x, perl = TRUE)
  as.integer(x)
}

# convert integer chromosome codes back to ordered chromosome labels
chr_factor <- function(chr_int) {
  factor(
    chr_int,
    levels = 1:9,
    labels = c(paste0("chr", 1:7, "R"), "chrUn", "chrB")
  )
}

# read windowed file and add comparison metadata
read_fst_win <- function(path, series, row_group = NA) {
  read_tsv(path, show_col_types = FALSE) %>%
    transmute(
      Series = series,
      RowGroup = row_group,
      CHR = normalize_chr(CHROM),
      BIN_START = BIN_START,
      BIN_END = BIN_END,
      midBP = (BIN_START + BIN_END) / 2,
      FST = pmax(0, as.numeric(WEIGHTED_FST))
    ) %>%
    filter(!is.na(CHR), is.finite(FST)) %>%
    mutate(
      CHR_f = chr_factor(CHR)
    ) %>%
    filter(!CHR_f %in% c("chrUn", "chrB"))
}

# read + combine many FST files listed in a metadata table
read_many_fst_win <- function(meta) {
  bind_rows(lapply(seq_len(nrow(meta)), function(i) {
    read_fst_win(
      path = meta$file[i],
      series = meta$Series[i],
      row_group = meta$RowGroup[i]
    )
  }))
}

# plot chromosome-wise FST profiles for one set of comparisons
plot_fst_lines <- function(fst_df, cols, title = "", one_row = FALSE) {
  
  # set facet order: Control first, then ND
  if ("RowGroup" %in% names(fst_df)) {
    fst_df <- fst_df %>%
      mutate(
        RowGroup = factor(RowGroup, levels = c("Control", "ND"))
      )
  }
  p <- ggplot(
    fst_df,
    aes(midBP, FST, colour = Series, group = Series)
  ) +
    geom_line(linewidth = 0.45) +
    scale_x_continuous(labels = function(x) x / 1e6) +
    scale_colour_manual(values = cols) +
    guides(
      colour = guide_legend(
        override.aes = list(linewidth = 1.5),
        ncol = 1
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9),
      strip.text.x = element_text(size = 9),
      
      # reverse direction of vertical Control/ND labels
      strip.text.y.right = element_text(
        angle = 270,
        size = 9
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 8),
      #axis.title.x = element_text(size = 9),
      #axis.title.y = element_text(size = 9)
    ) +
    labs(
      x = "Genome position (Mbp)",
      y = expression("Weighted " ~ F[ST]),
      colour = NULL,
      title = title
    )
  
  if (one_row) {
    p + facet_grid(
      cols = vars(CHR_f),
      scales = "free_x"
    )
  } else {
    p + facet_grid(
      rows = vars(RowGroup),
      cols = vars(CHR_f),
      scales = "free_x",
      labeller = labeller(
        RowGroup = c(
          "Control" = "control",
          "ND" = "ND"
        )
      )
    )
  }
}


################################################################################
# General settings
################################################################################

# colors for generation-wise comparisons
gen_cols <- c(
  "Gen1" = "darkgreen",
  "Gen2" = "darkred",
  "Gen3" = "darkblue"
)

# colors for consecutive-generation comparisons
step_cols <- c(
  "Gen0 vs Gen1" = "darkgreen",
  "Gen1 vs Gen2" = "darkred",
  "Gen2 vs Gen3" = "darkblue"
)


################################################################################
# FST of each generation versus Gen0
################################################################################

# input files and plotting labels for Gen0 comparisons
meta_vs_gen0 <- data.frame(
  file = c(
    "Fst_gen0_vs_gen1_Control_25000kb.windowed.weir.fst",
    "Fst_gen0_vs_gen2_Control_25000kb.windowed.weir.fst",
    "Fst_gen0_vs_gen3_Control_25000kb.windowed.weir.fst",
    "Fst_gen0_vs_gen1_ND_25000kb.windowed.weir.fst",
    "Fst_gen0_vs_gen2_ND_25000kb.windowed.weir.fst",
    "Fst_gen0_vs_gen3_ND_25000kb.windowed.weir.fst"
  ),
  Series = c(
    "Gen1", "Gen2", "Gen3",
    "Gen1", "Gen2", "Gen3"
  ),
  RowGroup = c(
    "Control", "Control", "Control",
    "ND", "ND", "ND"
  )
)

# read all Gen0 comparison files
fst_vs_gen0 <- read_many_fst_win(meta_vs_gen0)

# plot chromosome-wise FST profiles
p1 <- plot_fst_lines(
  fst_vs_gen0,
  cols = gen_cols,
  title = expression("F"["ST"] ~ "versus Gen0")
)
ggsave("fst_25000kb_vs_gen0.pdf", p1, width = 14, height = 4.4)


################################################################################
# FST between consecutive generations
################################################################################

# input files and plotting labels for consecutive-generation comparisons
meta_steps <- data.frame(
  file = c(
    "Fst_gen0_vs_gen1_Control_25000kb.windowed.weir.fst",
    "Fst_gen1_Control_vs_gen2_Control_25000kb.windowed.weir.fst",
    "Fst_gen2_Control_vs_gen3_Control_25000kb.windowed.weir.fst",
    "Fst_gen0_vs_gen1_ND_25000kb.windowed.weir.fst",
    "Fst_gen1_ND_vs_gen2_ND_25000kb.windowed.weir.fst",
    "Fst_gen2_ND_vs_gen3_ND_25000kb.windowed.weir.fst"
  ),
  Series = c(
    "Gen0 vs Gen1", "Gen1 vs Gen2", "Gen2 vs Gen3",
    "Gen0 vs Gen1", "Gen1 vs Gen2", "Gen2 vs Gen3"
  ),
  RowGroup = c(
    "Control", "Control", "Control",
    "ND", "ND", "ND"
  )
)

# read all consecutive-generation comparison files
fst_steps <- read_many_fst_win(meta_steps)

# plot chromosome-wise FST profiles
p2 <- plot_fst_lines(
  fst_steps,
  cols = step_cols,
  title = expression("F"["ST"] ~ "between consecutive generations")
)
ggsave("fst_25000kb_consecutive_generations.pdf", p2, width = 15, height = 4.4)


################################################################################
# FST between ND and Control within generations
################################################################################

# input files and plotting labels for treatment comparisons
meta_control_vs_nd <- data.frame(
  file = c(
    "Fst_gen1_ND_vs_gen1_Control_25000kb.windowed.weir.fst",
    "Fst_gen2_ND_vs_gen2_Control_25000kb.windowed.weir.fst",
    "Fst_gen3_ND_vs_gen3_Control_25000kb.windowed.weir.fst"
  ),
  Series = c(
    "Gen1",
    "Gen2",
    "Gen3"
  ),
  RowGroup = NA
)

# read all ND versus Control comparison files
fst_control_vs_nd <- read_many_fst_win(meta_control_vs_nd)

# plot chromosome-wise FST profiles -> one row
p3 <- plot_fst_lines(
  fst_control_vs_nd,
  cols = gen_cols,
  title = expression("F"["ST"] ~ "control versus ND within generations"),
  one_row = TRUE
)
ggsave("fst_25000kb_Control_vs_ND.pdf", p3, width = 14, height = 2.8)


