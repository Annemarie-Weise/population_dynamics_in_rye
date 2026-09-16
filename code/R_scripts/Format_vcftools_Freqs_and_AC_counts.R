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
# Parse vcftools allele-frequency files
################################################################################
# MAF is reported in a seperate column
# Sites with non-finite allele frequencies, e.g. NA, are removed.

read_maf_vcftools_frq <- function(frq_gz) {
  dt <- fread(frq_gz, header = FALSE)
  names(dt) <- c("CHROM","POS","N_ALLELES","N_CHROM","REF_AF","ALT_AF")
  
  # vcftools reports "BASE:FREQ" -> split into base and frequency
  dt$REF_BASE <- sub(":.*", "", dt$REF_AF)
  dt$REF_FREQ <- as.numeric(sub(".*:", "", dt$REF_AF))
  dt$ALT_BASE <- sub(":.*", "", dt$ALT_AF)
  dt$ALT_FREQ <- as.numeric(sub(".*:", "", dt$ALT_AF))
  # drop combined allele-frequency columns
  dt <- dt[, -(5:6)]
  gc()
  
  # remove nan/non-finite sites
  keep_non_nan <- is.finite(dt$REF_FREQ) & is.finite(dt$ALT_FREQ)
  dt <- dt[keep_non_nan]
  
  # add MAF and site ID
  dt$MAF <- pmin(dt$REF_FREQ, dt$ALT_FREQ)
  dt$ID <- paste(dt$CHROM, dt$POS, sep="_")
  return(dt)
}

# Convert each population-specific vcftools frequency file to a cleaned table
pops <- c("gen1_ND","gen1_Control","gen2_ND","gen2_Control","gen3_ND","gen3_Control","gen0")
for (pop in pops) {
  cur_pop <- read_maf_vcftools_frq(paste0("LD_freq_",pop,"_miss0.9_vcftools.frq"))
  fwrite(cur_pop, file = paste0("LD_freq_",pop,"_MAF_0.9miss.csv.gz"),compress = "gzip")
}



################################################################################
# Add allele counts from the corresponding *.frq.count files
################################################################################
# Source VCF was already filtered to gen0-gen3 individuals and
# filtered with --mac 1 -> fixed-REF sites are already removed
# Check if all sites:
#   - are present in all gen0-gen3 population frequency files
#   - are not fixed for the ALT allele across all populations

pops <- c(
  "gen0",
  "gen1_ND", "gen1_Control",
  "gen2_ND", "gen2_Control",
  "gen3_ND", "gen3_Control"
)

# read ALT freqs for each population and add population label
dt_all <- rbindlist(lapply(pops, function(pop) {
  fread(
    paste0("LD_freq_", pop, "_MAF_0.9miss.csv.gz"),
    select = c("CHROM", "POS", "ID", "REF_BASE", "ALT_BASE", "ALT_FREQ")
  )[, POP := pop]
}))
# convert to one row per variant, with ALT_FREQ columns for each population
dt_wide <- dcast(
  dt_all,
  CHROM + POS + ID + REF_BASE + ALT_BASE ~ POP,
  value.var = "ALT_FREQ"
)
# keep sites present in all populations and remove sites fixed for ALT
alt_min <- do.call(pmin, c(dt_wide[, ..pops], na.rm = FALSE))
dt_wide <- dt_wide[!is.na(alt_min) &alt_min < 1]
rm(dt_all, alt_min)
gc()

# Read vcftools *.frq.count files and add N_CHROM and ALT_AC
read_ac_vcftools <- function(pop) {
  dt <- fread(
    paste0("AC_per_site_LD_miss0.9_", pop, "_vcftools.frq.count"),
    header = FALSE
  )
  # count output stores alleles as "BASE:COUNT" -> extract numeric ALT_AC
  setnames(dt, c("CHROM", "POS", "N_ALLELES", "N_CHROM", "REF_AC", "ALT_AC"))
  dt[, ALT_AC := as.numeric(sub(".*:", "", ALT_AC))]
  # remove sites with invalid chromosome or allele counts
  dt <- dt[is.finite(N_CHROM) & is.finite(ALT_AC)]
  
  # match count table to frequency table by CHROM_POS ID
  dt[, ID := paste(CHROM, POS, sep = "_")]
  dt <- dt[, .(ID, N_CHROM, ALT_AC)]
  
  # make count columns population-specific before merging
  setnames(
    dt,
    c("N_CHROM", "ALT_AC"),
    paste0(c("N_CHROM_", "ALT_AC_"), pop)
  )
  dt
}

# read count files for all populations and merge them by site ID
ac_list <- lapply(pops, read_ac_vcftools)
dt_ac <- Reduce(
  function(x, y) merge(x, y, by = "ID", all = TRUE),
  ac_list
)
# final table of shared gen0-gen3 variants with frequencies and counts
dt_wide <- merge(dt_wide, dt_ac, by = "ID", all.x = TRUE)
fwrite(
  dt_wide,
  file = "LD_all_pop_existing_sites_miss0.9.csv.gz",
  compress = "gzip"
)
