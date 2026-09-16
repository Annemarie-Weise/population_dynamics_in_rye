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
source(paste0(git_path, "/code/R_scripts/wright_fisher_simulation_functions.R"))
setwd(paste0(git_path, "/general_analysis/vcftools/freq"))


################################################################################
# Compare empirical MAF shifts with Wright-Fisher drift simulations
################################################################################
# - calculate pairwise Kolmogorov-Smirnov statistics between empirical MAF dist
# - estimate Ne from allele-frequency changes (Gen0->Gen1) 
# - simulate neutral Wright-Fisher drift under different Ne scenarios
# - compare simulated MAF shifts (Gen1 ->...) with the observed KS statistics
################################################################################


################################################################################
# Helper functions
################################################################################

# calculate the KS statistic between MAF distributions of two populations
run_pairwise_ks <- function(data, pop1, pop2) {
  x <- data$MAF[data$pop == pop1]
  y <- data$MAF[data$pop == pop2]
  test <- suppressWarnings(ks.test(x, y))
  data.frame(
    pop1 = pop1,
    pop2 = pop2,
    statistic = unname(test$statistic),
    stringsAsFactors = FALSE
  )
}


# return only the KS statistic for two numeric vectors
ks_stat_only <- function(x, y) {
  suppressWarnings(as.numeric(ks.test(x, y)$statistic))
}


# avoid exact allele frequencies of 0 or 1 for Ne estimation
clip_freq <- function(p, eps = 1e-4) {
  pmin(pmax(p, eps), 1 - eps)
}


# calculate empirical one-sided p-value from simulated and observed statistics
emp_p <- function(sim, obs) {
  (sum(sim >= obs) + 1) / (length(sim) + 1)
}


# estimate Ne from allele-frequency change using Jonas Plan I correction
# - adapted for individual-level sequencing
estimate_Ne_plan1 <- function(p0, pt, n_chrom0, n_chromt, N, t = 1,
                              min_freq = 0.05,
                              max_freq = 0.95) {
  # avoid exact 0 or 1 allele frequencies
  p0 <- clip_freq(p0)
  pt <- clip_freq(pt)
  
  # allele-frequency change between time points
  delta <- pt - p0
  
  # z = mean allele frequency between the two time points
  z <- (p0 + pt) / 2
  
  # keep only SNPs with intermediate frequency for Ne calculation
  keep <- z >= min_freq & z <= max_freq
  p0 <- p0[keep]
  pt <- pt[keep]
  z  <- z[keep]
  n_chrom0 <- n_chrom0[keep]
  n_chromt <- n_chromt[keep]
  delta <- delta[keep]
  
  # Jorde-Ryman/Jonas denominator for biallelic SNPs
  D <- z - p0 * pt
  
  # Uncorrected standardized allele-frequency change
  Fc = sum((p0 - pt)^2) / sum(z - p0*pt)
  
  # sampling correction for individual-level sequencing:
  # Jonas uses Ct = 1/(2*St) for individual sampling
  # Since n_chrom = 2*St, this becomes Ct = 1/n_chrom_t
  C0 <- 1 / n_chrom0
  Ct <- 1 / n_chromt
  
  # Plan I correction from Jonas et al., without Pool-seq read correction
  # F_corrected = sum( D * [ Fc_i * (1 - 1/(2N)) - C0 - Ct + 1/N ] )
  #               --------------------------------------------------
  #                              sum( D * [1 - Ct] )
  # Since Fc_i = delta^2 / D, this becomes
  numerator <- sum(
    delta^2 * (1 - 1 / (2 * N)) -
      D * (C0 + Ct - 1 / N)
  )
  denominator <- sum(
    D * (1 - Ct)
  )
  F_corrected <- numerator / denominator
  
  # convert corrected drift variance to Ne
  Ne <- -t / (2 * log(1 - F_corrected))
  data.frame(
    F_corrected = F_corrected,
    Ne = round(Ne),
    F_uncorrected = Fc
  )
}


# simulate Gen0 -> Gen1 under one Ne multiplier
# and compare the observed Gen0-Gen1 KS statistic with the simulated null
run_calibration_scenario <- function(multiplier, n_sims = 1000) {
  
  Ne_nd <- round(base_Ne_nd * multiplier)
  Ne_control <- round(base_Ne_control * multiplier)

  sim_res <- replicate(n_sims, {
    # ND: simulate one drift step from Gen0 to Gen1
    sim_nd <- simulate_wf_trajectory(
      p0 = dat$ALT_FREQ_gen0,
      Ne_vec = Ne_nd,
      n_chrom_mat = cbind(dat$N_CHROM_gen1_ND),
      return = "freq",
      maf = TRUE
    )
    
    # Control: simulate one drift step from Gen0 to Gen1
    sim_control <- simulate_wf_trajectory(
      p0 = dat$ALT_FREQ_gen0,
      Ne_vec = Ne_control,
      n_chrom_mat = cbind(dat$N_CHROM_gen1_Control),
      return = "freq",
      maf = TRUE
    )
  
    c(
      ND = ks_stat_only(maf_gen0, sim_nd[[1]]),
      Control = ks_stat_only(maf_gen0, sim_control[[1]])
    )
  }) |> t() |> as.data.frame()
  
  data.frame(
    multiplier = multiplier,
    treatment = c("ND", "Control"),
    Ne = c(Ne_nd, Ne_control),
    observed_ks = c(obs_nd_g1, obs_control_g1),
    empirical_p = c(
      emp_p(sim_res$ND, obs_nd_g1),
      emp_p(sim_res$Control, obs_control_g1)
    ),
    sim_median_ks = c(
      median(sim_res$ND),
      median(sim_res$Control)
    )
  )
}

get_obs_ks <- function(p1, p2) {
  obs_res[obs_res$pop1 == p1 & obs_res$pop2 == p2, ]$statistic
}



################################################################################
# Calculate empirical KS statistics against Gen0
################################################################################

pops <- c(
  "gen0",
  "gen1_ND", "gen1_Control",
  "gen2_ND", "gen2_Control",
  "gen3_ND", "gen3_Control"
)

# read all population-specific MAF tables
maf_list <- lapply(pops, function(pop) {
  dt <- fread(paste0("LD_freq_", pop, "_MAF_0.9miss.csv.gz"))
  dt[, pop := pop]
  dt
})
maf_all <- rbindlist(maf_list)

# add generation and treatment labels
maf_all <- maf_all %>%
  mutate(
    gen = case_when(
      pop == "gen0" ~ "gen0",
      TRUE ~ str_extract(pop, "^gen\\d+")
    ),
    trt = case_when(
      pop == "gen0" ~ "gen0",
      str_detect(pop, "_Control$") ~ "Control",
      str_detect(pop, "_ND$")  ~ "ND"
    )
  )

# Observed Gen0-Gen1 KS statistics used for calibration
obs_nd_g1 <- run_pairwise_ks(maf_all,"gen0","gen1_ND")$statistic
obs_control_g1 <- run_pairwise_ks(maf_all,"gen0","gen1_Control")$statistic


################################################################################
# Estimate initial drift-calibration values
################################################################################

# read ALT frequencies and chromosome counts for all shared sites
dat <- fread("LD_all_pop_existing_sites_miss0.9.csv.gz") %>%
  select(
    ID,
    ALT_FREQ_gen0 = gen0,
    ALT_FREQ_gen1_ND = gen1_ND,
    ALT_FREQ_gen2_ND = gen2_ND,
    ALT_FREQ_gen3_ND = gen3_ND,
    ALT_FREQ_gen1_Control = gen1_Control,
    ALT_FREQ_gen2_Control = gen2_Control,
    ALT_FREQ_gen3_Control = gen3_Control,
    N_CHROM_gen0,
    N_CHROM_gen1_ND,
    N_CHROM_gen2_ND,
    N_CHROM_gen3_ND,
    N_CHROM_gen1_Control,
    N_CHROM_gen2_Control,
    N_CHROM_gen3_Control
  ) %>%
  filter(
    if_all(
      starts_with("ALT_FREQ_"),
      is.finite
    )
  )

# initial calibration values from Gen0-Gen1 allele-frequency differences
Ne_nd_01 <- estimate_Ne_plan1(
  p0 = dat$ALT_FREQ_gen0,
  pt = dat$ALT_FREQ_gen1_ND,
  n_chrom0 = dat$N_CHROM_gen0,
  n_chromt = dat$N_CHROM_gen1_ND,
  N = 1500,
  t = 1
)

Ne_control_01 <- estimate_Ne_plan1(
  p0 = dat$ALT_FREQ_gen0,
  pt = dat$ALT_FREQ_gen1_Control,
  n_chrom0 = dat$N_CHROM_gen0,
  n_chromt = dat$N_CHROM_gen1_Control,
  N = 1500,
  t = 1
)

base_Ne_nd <- Ne_nd_01$Ne
base_Ne_control <- Ne_control_01$Ne

################################################################################
# Gen0-Gen1 sensitivity analysis
################################################################################

Ne_multipliers <- seq(1, 4, by = 0.25)
set.seed(123)
n_sims <- 1000

# observed Gen0 MAF distribution used as calibration reference
maf_gen0 <- pmin(
  dat$ALT_FREQ_gen0,
  1 - dat$ALT_FREQ_gen0
)

# run Gen0-Gen1 calibration for each multiplier
calibration_results <- bind_rows(
  lapply(
    Ne_multipliers,
    run_calibration_scenario,
    n_sims = n_sims
  )
)

fwrite(calibration_results,"LD_drift_calibration_Gen0_Gen1.csv")




################################################################################
# Final Wright-Fisher drift test starting from Gen1
################################################################################

# selected multiplier from Gen0-Gen1 sensitivity analysis
final_multiplier <- 2.5

# final calibrated Ne values
# two drift steps: Gen1 -> Gen2 -> Gen3
Ne_final <- list(
  ND = rep(round(base_Ne_nd * final_multiplier), 2),             # 398
  Control = rep(round(base_Ne_control * final_multiplier), 2)    # 335
)


# Observed KS statistics: Gen1 versus Gen2 and Gen3
obs_nd_g2 <- run_pairwise_ks(maf_all, "gen1_ND", "gen2_ND")$statistic
obs_nd_g3 <- run_pairwise_ks(maf_all, "gen1_ND", "gen3_ND")$statistic

obs_control_g2 <- run_pairwise_ks(maf_all, "gen1_Control", "gen2_Control")$statistic
obs_control_g3 <- run_pairwise_ks(maf_all, "gen1_Control", "gen3_Control")$statistic

# Gen1 MAF distributions used as simulation reference
maf_gen1_nd <- maf_all$MAF[maf_all$pop == "gen1_ND"]
maf_gen1_control <- maf_all$MAF[maf_all$pop == "gen1_Control"]

# Wright-Fisher simulations starting from observed Gen1 allele frequencies
set.seed(123)
n_sims <- 1000
sim_res_gen1 <- replicate(n_sims, {
  # ND: Gen1 -> Gen2 -> Gen3
  sim_nd <- simulate_wf_trajectory(
    p0 = dat$ALT_FREQ_gen1_ND,
    Ne_vec = Ne_final$ND,
    n_chrom_mat = cbind(
      dat$N_CHROM_gen2_ND,
      dat$N_CHROM_gen3_ND
    ),
    return = "freq",
    maf = TRUE
  )
  
  # Control: Gen1 -> Gen2 -> Gen3
  sim_control <- simulate_wf_trajectory(
    p0 = dat$ALT_FREQ_gen1_Control,
    Ne_vec = Ne_final$Control,
    n_chrom_mat = cbind(
      dat$N_CHROM_gen2_Control,
      dat$N_CHROM_gen3_Control
    ),
    return = "freq",
    maf = TRUE
  )
  
  c(
    nd_g2 = ks_stat_only(maf_gen1_nd, sim_nd[[1]]),
    nd_g3 = ks_stat_only(maf_gen1_nd, sim_nd[[2]]),
    control_g2 = ks_stat_only(maf_gen1_control, sim_control[[1]]),
    control_g3 = ks_stat_only(maf_gen1_control, sim_control[[2]])
  )
}) |> t() |> as.data.frame()



# Empirical p-values
obs_gen1 <- c(
  nd_g2 = obs_nd_g2,
  nd_g3 = obs_nd_g3,
  control_g2 = obs_control_g2,
  control_g3 = obs_control_g3
)

final_drift_results <- data.frame(
  comparison = names(obs_gen1),
  observed_ks = as.numeric(obs_gen1),
  empirical_p = mapply(
    function(sim_col, obs_val) {
      emp_p(sim_res_gen1[[sim_col]], obs_val)
    },
    names(obs_gen1),
    obs_gen1
  ),
  sim_median_ks = sapply(
    names(obs_gen1),
    function(x) median(sim_res_gen1[[x]])
  )
)
fwrite(final_drift_results,"LD_drift_test_final_vsGen1.csv")