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
source(paste0(git_path, "/code/R_scripts/wright_fisher_simulation_functions.R"))
setwd(paste0(git_path,"/general_analysis/vcftools/freq"))

################################################################################
# Logistic allele-count trends and Wright-Fisher drift test
################################################################################
# - fit SNP-wise logistic regressions of ALT allele counts across generations
# - use Wright-Fisher simulations to test obs allele-frequency trends against 
#   neutral drift
# calculatest:
#   - empirical GLM slope and p-value per SNP
#   - simulation-based empirical p-values
#   - refined empirical p-values for initially significant SNPs
################################################################################


################################################################################
# Helper functions
################################################################################

# fit SNP-wise binomial GLMs for one treatment line
run_empirical_glm <- function(global, treatment) {
  pops <- c("gen1", "gen2", "gen3")
  count_cols_alt <- paste0("ALT_AC_", pops, "_", treatment)
  count_cols_n   <- paste0("N_CHROM_", pops, "_", treatment)
  
  # select ALT allele counts and chromosome counts for Gen1-Gen3
  glm_input <- global %>%
    select(
      CHROM, POS, ID, REF_BASE, ALT_BASE,
      all_of(count_cols_alt),
      all_of(count_cols_n)
    )
  
  # convert wide count table to long format for GLM fitting
  maf_long <- glm_input %>%
    pivot_longer(
      cols = starts_with("ALT_AC_") | starts_with("N_CHROM_"),
      names_to = c(".value", "gen"),
      names_pattern = "(ALT_AC|N_CHROM)_(gen\\d+)"
    ) %>%
    mutate(time = as.numeric(str_extract(gen, "\\d+")))
  
  # fit one binomial GLM per SNP: ALT counts ~ generation time
  maf_long %>%
    group_by(CHROM, POS, ID, REF_BASE, ALT_BASE) %>%
    group_modify(~ {
      fit <- glm(
        cbind(ALT_AC, N_CHROM - ALT_AC) ~ time,
        family = binomial,
        data = .x
      )
      null <- glm(
        cbind(ALT_AC, N_CHROM - ALT_AC) ~ 1,
        family = binomial,
        data = .x
      )
      s <- summary(fit)$coefficients
      data.frame(
        intercept = s[1, 1],
        slope = s[2, 1],
        p_value = s[2, 4],
        pseudo_r2 = 1 - as.numeric(logLik(fit) / logLik(null))
      )
    }) %>%
    ungroup()
}


# prepare allele counts and starting frequencies for simulations
prepare_sim_data <- function(global, treatment) {
  key_cols <- c("CHROM", "POS", "ID", "REF_BASE", "ALT_BASE")
  global %>%
    select(
      all_of(key_cols),
      all_of(paste0("ALT_AC_gen1_", treatment)),
      all_of(paste0("ALT_AC_gen2_", treatment)),
      all_of(paste0("ALT_AC_gen3_", treatment)),
      all_of(paste0("N_CHROM_gen1_", treatment)),
      all_of(paste0("N_CHROM_gen2_", treatment)),
      all_of(paste0("N_CHROM_gen3_", treatment))
    ) %>%
    mutate(
      p0 = .data[[paste0("ALT_AC_gen1_", treatment)]] /
        .data[[paste0("N_CHROM_gen1_", treatment)]]
    ) %>%
    filter(
      is.finite(p0),
      .data[[paste0("N_CHROM_gen1_", treatment)]] > 0,
      .data[[paste0("N_CHROM_gen2_", treatment)]] > 0,
      .data[[paste0("N_CHROM_gen3_", treatment)]] > 0
    )
}


# calculate observed logistic slopes from empirical allele counts
get_obs_emp_slopes <- function(dat, treatment) {
  alt_cols <- c(
    paste0("ALT_AC_gen1_", treatment),
    paste0("ALT_AC_gen2_", treatment),
    paste0("ALT_AC_gen3_", treatment)
  )
  n_cols <- c(
    paste0("N_CHROM_gen1_", treatment),
    paste0("N_CHROM_gen2_", treatment),
    paste0("N_CHROM_gen3_", treatment)
  )
  K_obs <- as.matrix(dat[, alt_cols])
  N_obs <- as.matrix(dat[, n_cols])
  fit_all_snps_fast(
    K_mat = K_obs,
    N_mat = N_obs,
    time = 1:3
  )
}


# simulate one Wright-Fisher trajectory and fit logistic slopes
simulate_one_run_glm <- function(dat, Ne_vec, treatment = c("Control", "ND")) {
  treatment <- match.arg(treatment)
  
  # Select observed chromosome counts for the simulated generations
  if (treatment == "Control") {
    n_mat <- cbind(dat$N_CHROM_gen2_Control, dat$N_CHROM_gen3_Control)
    K_gen1 <- dat$ALT_AC_gen1_Control
    N_gen1 <- dat$N_CHROM_gen1_Control
  } else {
    n_mat <- cbind(dat$N_CHROM_gen2_ND, dat$N_CHROM_gen3_ND)
    K_gen1 <- dat$ALT_AC_gen1_ND
    N_gen1 <- dat$N_CHROM_gen1_ND
  }
  
  # simulate ALT allele counts from Gen1 through Gen3 under drift
  sim_list <- simulate_wf_trajectory(
    p0 = dat$p0,
    Ne_vec = Ne_vec,
    n_chrom_mat = n_mat,
    return = "count"
  )
  
  # combine observed Gen1 counts with simulated Gen2-Gen3 counts
  sim_mat <- do.call(cbind, sim_list)
  K_mat <- cbind(K_gen1, sim_mat)
  N_mat <- cbind(N_gen1, n_mat)
  
  # fit logistic slopes to the simulated count trajectories
  fit_all_snps_fast(K_mat = K_mat,N_mat = N_mat, time = 1:3)
}


# run repeated Wright-Fisher simulations and calculate empirical p-values
run_wf_empirical_p <- function(dat, obs, treatment, Ne_vec, n_sims) {
  n_snps <- nrow(dat)
  abs_obs <- abs(obs)
  
  # count how often simulated slopes are at least as extreme as observed slopes
  extreme_counts <- integer(n_snps)
  valid_counts <- integer(n_snps)
  
  for (s in seq_len(n_sims)) {
    cat("Simulation", s, "of", n_sims, treatment, "\n")
    
    # simulate one drift replicate and fit logistic slopes
    beta_sim_s <- simulate_one_run_glm(
      dat = dat,
      Ne_vec = Ne_vec,
      treatment = treatment
    )
    
    # compare absolute simulated and observed slopes
    valid <- !is.na(beta_sim_s) & !is.na(obs)
    extreme <- valid & abs(beta_sim_s) >= abs_obs
    valid_counts[valid] <- valid_counts[valid] + 1L
    extreme_counts[extreme] <- extreme_counts[extreme] + 1L
  }
  
  # empirical p-value with +1 correction to avoid zero p-values
  p_empirical <- rep(NA, n_snps)
  has_valid <- valid_counts > 0
  p_empirical[has_valid] <- (extreme_counts[has_valid] + 1) / (valid_counts[has_valid] + 1)
  
  data.frame(
    extreme_count = extreme_counts,
    valid_count = valid_counts,
    p_empirical = p_empirical
  )
}


# correlation test between glm values and custom IRLS values
cor_glm_wf <- function(treatment) {
  glm <- fread(
    paste0("logistic_trend/LD_ALT_AC_logistic_trends_", treatment, "_0.9miss.csv.gz")
  ) %>%
    select(all_of(key_cols), glm_slope = slope)
  wf <- fread(
    paste0("logistic_trend/LD_empirical_pvalues_", treatment, "_refined_full_0.9miss.csv.gz")
  ) %>%
    select(all_of(key_cols), wf_slope = slope_obs)
  dat <- inner_join(glm, wf, by = key_cols)
  data.frame(
    treatment = treatment,
    pearson = cor(dat$glm_slope, dat$wf_slope, use = "complete.obs", method = "pearson"),
    spearman = cor(dat$glm_slope, dat$wf_slope, use = "complete.obs", method = "spearman")
  )
}



################################################################################
# Empirical logistic regression with glm()
################################################################################

# read merged allele frequency table
global <- fread("LD_all_pop_existing_sites_miss0.9.csv.gz") %>%
  as.data.frame()

# run empirical logistic regressions separately for both treatment lines
glm_control <- run_empirical_glm(global, "Control")
glm_nd  <- run_empirical_glm(global, "ND")
fwrite(
  glm_control,
  file = "logistic_trend/LD_ALT_AC_logistic_trends_Control_0.9miss.csv.gz",
  compress = "gzip"
)
fwrite(
  glm_nd,
  file = "logistic_trend/LD_ALT_AC_logistic_trends_ND_0.9miss.csv.gz",
  compress = "gzip"
)

# QQ plot of glm p-values
#########################
# combine ND and Control glm p-values + calculate expected/observed -log10(p)
qq_nd_glm <- glm_nd %>%
  as.data.frame() %>%
  mutate(track = "ND")
qq_control_glm <- glm_control %>%
  as.data.frame() %>%
  mutate(track = "Control")
qq_glm_all <- bind_rows(qq_nd_glm, qq_control_glm) %>%
  filter(!is.na(p_value)) %>%
  filter(p_value > 0, p_value <= 1) %>%
  group_by(track) %>%
  arrange(p_value, .by_group = TRUE) %>%
  mutate(
    exp_p    = (row_number() - 0.5) / n(),
    obs_logp = -log10(p_value),
    exp_logp = -log10(exp_p)
  ) %>%
  ungroup()

# plot glm p-value distributions against uniform expectation
p_glm <- ggplot(qq_glm_all, aes(x = exp_logp, y = obs_logp, color = track)) +
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
ggsave("logistic_trend/LD_qq_glm_pvalues_both_0.9miss.pdf", p_glm, width = 4, 
       height = 4, units = "in")


################################################################################
# Wright-Fisher simulation setup
################################################################################
# - prepare simulation input tables and define Ne trajectories for ND and Control
# - simulations start from Gen1 allele frequencies 
################################################################################

global <- fread("LD_all_pop_existing_sites_miss0.9.csv.gz") %>%
  as.data.frame()

# prepare simulation data for both treatment lines
dat_control <- prepare_sim_data(global, "Control")
dat_nd  <- prepare_sim_data(global, "ND")

# effective population size trajectories used for drift simulations
Ne_nd <- round(c(159, 159) * 2.5)
Ne_control <- round(c(134, 134) * 2.5)


################################################################################
# Initial Wright-Fisher simulations
################################################################################

# Number of simulations for the genome-wide initial screen
initial_sims <- 2000
set.seed(123)

# observed logistic slopes for Control
obs_control <- get_obs_emp_slopes(dat_control, "Control")
# simulated null distribution and empirical p-values for Control
sim_counts_control <- run_wf_empirical_p(
  dat = dat_control,
  obs = obs_control,
  treatment = "Control",
  Ne_vec = Ne_control,
  n_sims = initial_sims
)
# combine SNP metadata, observed slopes, empirical p-values, and BH correction
results_control <- dat_control %>%
  select(CHROM, POS, ID, REF_BASE, ALT_BASE) %>%
  mutate(
    slope_obs = obs_control
  ) %>%
  bind_cols(sim_counts_control) %>%
  mutate(
    p_adj = p.adjust(p_empirical, method = "BH")
  )
fwrite(
  results_control,
  file = "logistic_trend/LD_empirical_pvalues_Control_0.9miss.csv.gz",
  compress = "gzip"
)


# observed logistic slopes for ND
obs_nd <- get_obs_emp_slopes(dat_nd, "ND")
# simulated null distribution and empirical p-values for ND
sim_counts_nd <- run_wf_empirical_p(
  dat = dat_nd,
  obs = obs_nd,
  treatment = "ND",
  Ne_vec = Ne_nd,
  n_sims = initial_sims
)
# combine SNP metadata, observed slopes, empirical p-values, and BH correction
results_nd <- dat_nd %>%
  select(CHROM, POS, ID, REF_BASE, ALT_BASE) %>%
  mutate(
    slope_obs = obs_nd
  ) %>%
  bind_cols(sim_counts_nd) %>%
  mutate(
    p_adj = p.adjust(p_empirical, method = "BH")
  )
fwrite(
  results_nd,
  file = "logistic_trend/LD_empirical_pvalues_ND_0.9miss.csv.gz",
  compress = "gzip"
)


################################################################################
# Refine empirical p-values for initially significant SNPs
################################################################################

key_cols <- c("CHROM", "POS", "ID", "REF_BASE", "ALT_BASE")
set.seed(123)

# keep initially significant Control SNPs for refined simulation
dat_control_sig <- dat_control %>%
  semi_join(
    results_control %>%
      filter(p_adj < 0.05) %>%
      select(all_of(key_cols)),
    by = key_cols
  )
# recalculate observed slopes and empirical p-values with more simulations
obs_control_sig <- get_obs_emp_slopes(dat_control_sig, "Control")
sim_counts_control_refined <- run_wf_empirical_p(
  dat = dat_control_sig,
  obs = obs_control_sig,
  treatment = "Control",
  Ne_vec = Ne_control,
  n_sims = 100000
)
refined_table_control <- dat_control_sig %>%
  select(all_of(key_cols)) %>%
  bind_cols(sim_counts_control_refined) %>%
  rename(
    extreme_count_refined = extreme_count,
    valid_count_refined = valid_count,
    p_empirical_refined = p_empirical
  )
# merge refined p-values back into the full Control result table
results_control_refined <- results_control %>%
  left_join(refined_table_control, by = key_cols) %>%
  mutate(
    p_empirical_final = coalesce(p_empirical_refined, p_empirical),
    p_adj_BH = p.adjust(p_empirical_final, method = "BH"),
    p_adj_BY = p.adjust(p_empirical_final, method = "BY")
  )
fwrite(
  results_control_refined,
  file = "logistic_trend/LD_empirical_pvalues_Control_refined_full_0.9miss.csv.gz",
  compress = "gzip"
)


# keep initially significant ND SNPs for refined simulation
dat_nd_sig <- dat_nd %>%
  semi_join(
    results_nd %>%
      filter(p_adj < 0.05) %>%
      select(all_of(key_cols)),
    by = key_cols
  )
# recalculate observed slopes and empirical p-values with more simulations
obs_nd_sig <- get_obs_emp_slopes(dat_nd_sig, "ND")
sim_counts_nd_refined <- run_wf_empirical_p(
  dat = dat_nd_sig,
  obs = obs_nd_sig,
  treatment = "ND",
  Ne_vec = Ne_nd,
  n_sims = 100000
)
refined_table_nd <- dat_nd_sig %>%
  select(all_of(key_cols)) %>%
  bind_cols(sim_counts_nd_refined) %>%
  rename(
    extreme_count_refined = extreme_count,
    valid_count_refined = valid_count,
    p_empirical_refined = p_empirical
  )
# merge refined p-values back into the full ND result table
results_nd_refined <- results_nd %>%
  left_join(refined_table_nd, by = key_cols) %>%
  mutate(
    p_empirical_final = coalesce(p_empirical_refined, p_empirical),
    p_adj_BH = p.adjust(p_empirical_final, method = "BH"),
    p_adj_BY = p.adjust(p_empirical_final, method = "BY")
  )
fwrite(
  results_nd_refined,
  file = "logistic_trend/LD_empirical_pvalues_ND_refined_full_0.9miss.csv.gz",
  compress = "gzip"
)



################################################################################
# Correlation: glm slope vs fast IRLS slope
################################################################################

key_cols <- c("CHROM", "POS", "ID", "REF_BASE", "ALT_BASE")
cor_summary <- bind_rows(
  cor_glm_wf("Control"),
  cor_glm_wf("ND")
)
fwrite(
  cor_summary,
  "logistic_trend/glm_vs_fastIRLS_slope_correlations.csv"
)

