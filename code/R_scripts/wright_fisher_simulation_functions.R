################################################################################
# Wright-Fisher simulation helper functions
################################################################################
# - simulating neutral allele-frequency drift under Wright-Fisher model
# - sampling observed allele counts/allele frequencies from simulated frequency
################################################################################

# Simulate one generation of Wright-Fisher drift
# p = current allele frequencies
# Ne = effective population size used for binomial sampling of 2 * Ne chrs
drift_step <- function(p, Ne) {
  rbinom(length(p), size = 2 * Ne, prob = p) / (2 * Ne)
}


# sample observed ALT allele counts from true allele frequencies
# n_chrom = number of sampled chromosomes per SNP
sample_counts <- function(p, n_chrom) {
  rbinom(length(p), size = n_chrom, prob = p)
}


# Sample observed ALT allele frequencies from true allele frequencies
# - maf = TRUE: return MAFs instead of ALT frequencies
sample_freq <- function(p, n_chrom, maf = FALSE) {
  out <- sample_counts(p, n_chrom) / n_chrom
  if (maf) {
    out <- pmin(out, 1 - out)
  }
  out
}


# Simulate a Wright-Fisher trajectory over multiple generations
simulate_wf_trajectory <- function(p0, Ne_vec, n_chrom_mat = NULL,
                                   return = c("freq", "count"),
                                   maf = FALSE) {
  return <- match.arg(return)
  
  # initial allele frequencies
  p <- p0
  # one output vector per generation step
  out <- vector("list", length(Ne_vec))
  
  for (i in seq_along(Ne_vec)) {
    # apply one generation of drift
    p <- drift_step(p, Ne_vec[i])
    if (is.null(n_chrom_mat)) {
      # return true population frequencies if no sampling scheme is given
      out[[i]] <- p
      
    } else if (return == "count") {
      # sample observed allele counts using generation-specific chromosome counts
      out[[i]] <- sample_counts(p, n_chrom_mat[, i])
      
    } else {
      # sample observed allele frequencies, optionally converted to MAF
      out[[i]] <- sample_freq(p, n_chrom_mat[, i], maf = maf)
    }
  }
  out
}


# fast SNP-wise logistic regression using custom IRLS
fit_all_snps_fast <- function(K_mat, N_mat, time = 0:3, maxit = 25, tol = 1e-8) {
  n_snps <- nrow(K_mat)
  
  # time matrix in the same shape as K_mat/N_mat
  t_mat  <- matrix(rep(time, each = n_snps), nrow = n_snps)
  t2_mat <- t_mat^2
  
  # starting values from corrected allele proportions
  p <- (K_mat + 0.5) / (N_mat + 1)
  eta <- qlogis(p)
  
  # initial beta0/beta1 values from linear regression on logit frequencies
  t_bar <- mean(time)
  beta1 <- rowSums((eta - rowMeans(eta)) * (t_mat - t_bar)) / sum((time - t_bar)^2)
  beta0 <- rowMeans(eta) - beta1 * t_bar
  
  # IRLS loop for fitting logistic regression per SNP
  for (it in seq_len(maxit)) {
    # current linear predictor and fitted allele frequency
    eta <- beta0 + beta1 * t_mat
    eta <- pmin(pmax(eta, -30), 30)  # avoid numerical overflow
    mu <- plogis(eta)
    
    # binomial weights for each SNP/time point
    w  <- N_mat * mu * (1 - mu)
    w[w < 1e-12] <- 1e-12  # avoid division by zero
    
    # working response for weighted least squares update
    z <- eta + (K_mat - N_mat * mu) / w
    
    # weighted sums needed for closed-form update of beta0 and beta1
    sum_w   <- rowSums(w)
    sum_wx  <- rowSums(w * t_mat)
    sum_wx2 <- rowSums(w * t2_mat)
    sum_wz  <- rowSums(w * z)
    sum_wxz <- rowSums(w * z * t_mat)
    det <- sum_w * sum_wx2 - sum_wx^2
    
    beta0_new <- beta0
    beta1_new <- beta1
    
    # update only SNPs with a non-singular weighted regression system
    good <- abs(det) > 1e-12
    beta0_new[good] <- (sum_wz[good] * sum_wx2[good] - sum_wxz[good] * sum_wx[good]) / det[good]
    beta1_new[good] <- (sum_w[good] * sum_wxz[good] - sum_wx[good] * sum_wz[good]) / det[good]
    
    # stop when parameter updates are smaller than the tolerance
    delta <- pmax(abs(beta0_new - beta0), abs(beta1_new - beta1))
    beta0 <- beta0_new
    beta1 <- beta1_new
    
    if (max(delta[is.finite(delta)], na.rm = TRUE) < tol) break
  }
  
  # set uninformative cases to NA
  p_obs <- K_mat / N_mat
  bad <- apply(p_obs, 1, function(x) length(unique(x)) < 2)
  # Safety cap - prevent rare numerical explosions in IRLS fit
  beta1 <- pmin(pmax(beta1, -30), 30)
  beta1[bad | !is.finite(beta1)] <- NA
  beta1
}


