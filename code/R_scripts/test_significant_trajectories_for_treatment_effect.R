library(data.table)
library(dplyr)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
source(paste0(git_path, "/code/R_scripts/wright_fisher_simulation_functions.R"))
setwd(paste0(git_path, "/general_analysis/vcftools/freq"))


################################################################################
# Wright-Fisher test of treatment differences in significant SNP trajectories
################################################################################
# - simulate neutral allele-frequency trajectories separately for ND and Control
#   for the sig. SNPs
# - fit logistic slopes to both simulated trajectories
# - compare simulated treatment differences in slopes with the observed
#   treatment difference
# calculates:
#   - observed difference between ND and Control logistic slopes
#   - simulation-based empirical p-value for treatment divergence
#   - BH-adjusted empirical p-value
################################################################################


################################################################################
# Helper functions
################################################################################

# simulate one Wright-Fisher replicate for both treatments
# and calculate the difference between simulated logistic slopes
simulate_one_run_treatment_diff <- function(dat, Ne_control, Ne_nd) {
  # Control
  #--------
  # observed Gen1 starting frequencies
  p0_control <- dat$ALT_AC_gen1_Control / dat$N_CHROM_gen1_Control
  
  # observed chromosome counts for simulated generations
  n_mat_control <- cbind(dat$N_CHROM_gen2_Control, dat$N_CHROM_gen3_Control)
  
  # simulate Gen2 and Gen3 under neutral drift
  sim_control <- simulate_wf_trajectory(
    p0 = p0_control,
    Ne_vec = Ne_control,
    n_chrom_mat = n_mat_control,
    return = "count"
  )
  
  # combine observed Gen1 with simulated Gen2 and Gen3
  K_control <- cbind(dat$ALT_AC_gen1_Control, do.call(cbind, sim_control))
  N_control <- cbind(dat$N_CHROM_gen1_Control,n_mat_control)
  
  # fit logistic slope to simulated Control trajectory
  beta_control <- fit_all_snps_fast(
    K_mat = K_control,
    N_mat = N_control,
    time = 1:3
  )
  
  # ND
  #---
  # observed Gen1 starting frequencies
  p0_nd <- dat$ALT_AC_gen1_ND / dat$N_CHROM_gen1_ND
  
  # observed chromosome counts for simulated generations
  n_mat_nd <- cbind(dat$N_CHROM_gen2_ND, dat$N_CHROM_gen3_ND)
  
  # simulate Gen2 and Gen3 under neutral drift
  sim_nd <- simulate_wf_trajectory(
    p0 = p0_nd,
    Ne_vec = Ne_nd,
    n_chrom_mat = n_mat_nd,
    return = "count"
  )
  
  # combine observed Gen1 with simulated Gen2 and Gen3
  K_nd <- cbind(dat$ALT_AC_gen1_ND, do.call(cbind, sim_nd))
  N_nd <- cbind(dat$N_CHROM_gen1_ND, n_mat_nd)
  
  # fit logistic slope to simulated ND trajectory
  beta_nd <- fit_all_snps_fast(K_mat = K_nd, N_mat = N_nd, time = 1:3)

  
  # Treatment difference
  # signed difference in simulated logistic slopes
  beta_nd - beta_control
}


# run repeated Wright-Fisher simulations and calculate empirical p-values
# for the difference between ND and Control trajectories
run_wf_treatment_diff <- function(dat ,abs_obs_diff ,Ne_control ,Ne_nd ,n_sims) {
  n_snps <- nrow(dat)
  
  # count how often simulated treatment differences are at least
  # as extreme as the observed treatment difference
  extreme_counts <- integer(n_snps)
  valid_counts <- integer(n_snps)
  
  for (s in seq_len(n_sims)) {
    # print progress every 1000 simulations
    if (s %% 1000 == 0) {cat("Simulation", s, "of", n_sims, "\n")}
    
    # simulate one independent drift replicate for ND and Control
    # and calculate difference between simulated slopes
    beta_diff_sim <- simulate_one_run_treatment_diff(
      dat = dat,
      Ne_control = Ne_control,
      Ne_nd = Ne_nd
    )
    
    # compare absolute simulated and observed slope differences
    valid <- !is.na(beta_diff_sim) & !is.na(abs_obs_diff)
    extreme <- valid & abs(beta_diff_sim) >= abs_obs_diff
    valid_counts[valid] <- valid_counts[valid] + 1L
    extreme_counts[extreme] <- extreme_counts[extreme] + 1L
  }
  
  # empirical p-value with +1 correction to avoid zero p-values
  p_empirical <- rep(NA_real_, n_snps)
  has_valid <- valid_counts > 0
  p_empirical[has_valid] <- (extreme_counts[has_valid] + 1) / (valid_counts[has_valid] + 1)
  
  data.frame(
    extreme_count = extreme_counts,
    valid_count = valid_counts,
    p_empirical = p_empirical
  )
}


################################################################################
# Wright-Fisher treatment-difference simulation setup
################################################################################

# read table containing all 703 SNPs significant in at least one treatment
sig_snps <- fread(
  "logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz"
) %>%
  as.data.frame()

key_cols <- c("CHROM", "POS", "ID", "REF_BASE", "ALT_BASE")

# effective population sizes used in the original drift simulations
Ne_nd <- round(c(159, 159) * 2.5)
Ne_control <- round(c(134, 134) * 2.5)


################################################################################
# Observed treatment differences
################################################################################

# signed observed difference between treatment-specific slopes
abs_obs_diff <- sig_snps$Abs_Delta_slope_ND_Control

# Wright-Fisher simulations of treatment differences
n_sims <- 100000
set.seed(123)
sim_treatment_diff <- run_wf_treatment_diff(
  dat = sig_snps,
  abs_obs_diff = abs_obs_diff,
  Ne_control = Ne_control,
  Ne_nd = Ne_nd,
  n_sims = n_sims
)

# Combine observed and simulated results
results_treatment_diff <- sig_snps %>%
  select(
    all_of(key_cols),
    SIG_Control,
    SIG_ND,
    SIG_BOTH,
    slope_Control,
    slope_ND,
    Delta_slope_ND_Control,
    Abs_Delta_slope_ND_Control
  ) %>%
  bind_cols(sim_treatment_diff) %>%
  mutate(
    p_adj_BH = p.adjust(
      p_empirical,
      method = "BH"
    )
  )
fwrite(
  results_treatment_diff,
  file =
    "logistic_trend/LD_WF_treatment_difference_703_SNPs.csv.gz",
  compress = "gzip"
)