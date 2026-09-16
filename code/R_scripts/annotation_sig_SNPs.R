library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/annotation_and_GO/"))

# Path to feature table from: https://zenodo.org/records/15801361
# - SECALE.CEREALE.Lo7V3.pgsb.r2.Feb2025.featuretable.csv
featuretable_path <- readline(
  prompt = "Please enter the path to SECALE.CEREALE.Lo7V3.pgsb.r2.Feb2025.featuretable.csv: "
)
if (!file.exists(featuretable_path)) {
  stop("File not found: ", featuretable_path)
}
featuretable_file <- normalizePath(featuretable_path)


################################################################################
# Annotate significant SNPs and prepare gene tables for GO analysis
################################################################################
# - convert significant SNPs into a minimal VCF file
# - parse the annotated VCF and join gene feature information
# - create grouped gene tables for downstream GO-term analysis
################################################################################


################################################################################
# Create minimal VCF with significant SNP positions
################################################################################

# read significant SNPs from logistic trend analysis
snps <- read_csv(paste0(git_path, "/general_analysis/vcftools/freq/logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz")
                 , show_col_types = FALSE)

# convert SNP table to minimal VCF columns
sites_vcf <- snps %>%
  transmute(
    CHROM = as.character(CHROM),
    POS = as.integer(POS),
    ID = ID,
    REF = REF_BASE,
    ALT = ALT_BASE,
    QUAL = ".",
    FILTER = ".",
    INFO = "."
  ) %>%
  arrange(CHROM, POS, REF, ALT)

# write VCF header
con <- file("significant_SNPs.sites.vcf", open = "wt")
writeLines("##fileformat=VCFv4.2", con)
writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", con)
close(con)
write_tsv(sites_vcf,"significant_SNPs.sites.vcf",append = TRUE,col_names = FALSE)


################################################################################
# Parse annotated VCF and join gene feature information
################################################################################
# - after external annotation of significant_SNPs.sites.vcf -> extract ANN fields
# - split annotation entries into separate columns
# - add gene-level information from the rye feature table

# read annotated VCF while skipping metadata header lines
vcf <- read_tsv(
  "significant_SNPs.sites.annotated.vcf",
  comment = "##",
  col_names = FALSE,
  show_col_types = FALSE
)

# recover VCF column names from the #CHROM header line
header_line <- read_lines("significant_SNPs.sites.annotated.vcf") %>%
  .[str_starts(., "#CHROM")]
vcf_names <- str_split(header_line, "\t")[[1]]
vcf_names[1] <- "CHROM"
names(vcf) <- vcf_names

# extract raw ANN annotation string from the INFO field
ann_raw <- vcf %>%
  mutate(
    ANN = str_extract(INFO, "(?<=ANN=)[^;]+")
  ) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, ANN)

# split multiple annotations per SNP + separate SnpEff ANN fields
ann_long <- ann_raw %>%
  separate_rows(ANN, sep = ",") %>%
  separate(
    ANN,
    into = c(
      "ANN_ALLELE",
      "EFFECT",
      "IMPACT",
      "GENE_NAME",
      "GENE_ID",
      "FEATURE_TYPE",
      "FEATURE_ID",
      "BIOTYPE",
      "RANK",
      "HGVS_C",
      "HGVS_P",
      "CDNA_POS",
      "CDS_POS",
      "AA_POS",
      "DISTANCE",
      "ERRORS_WARNINGS"
    ),
    sep = "\\|",
    fill = "right",
    extra = "merge"
  )

# read rye feature table with gene descriptions, InterPro, and GO terms
feature <- read_csv(featuretable_file, show_col_types = FALSE, na = c("", "NA", "NaN", "nan"))

# collapse transcript-level feature information to one row per gene
# collapse transcript-level feature information to one row per gene
feature_gene <- feature %>%
  group_by(gene_id) %>%
  summarise(
    confidence = paste(sort(unique(na.omit(confidence))), collapse = "|"),
    description = paste(sort(unique(na.omit(description))), collapse = "|"),
    accession = paste(sort(unique(na.omit(accession))), collapse = "|"),
    family = paste(sort(unique(na.omit(Famliy))), collapse = "|"),
    InterPro = paste(sort(unique(na.omit(InterPro))), collapse = "|"),
    GO = paste(sort(unique(na.omit(GO))), collapse = "|"),
    n_transcripts_featuretable = n_distinct(transcript_id),
    .groups = "drop"
  )

# add gene-level feature information to SNP annotations
ann_long <- ann_long %>%
  left_join(
    feature_gene,
    by = c("GENE_ID" = "gene_id")
  )
write_csv( ann_long, "significant_SNPs.annotation.csv")


################################################################################
# Split annotated SNPs into gene tables for GO-term analysis
################################################################################


sig <- read_csv(
  paste0(git_path,"/general_analysis/vcftools/freq/logistic_trend/significant_SNPs_ND_or_Control_with_freqs_and_deltas.csv.gz"),
  show_col_types = FALSE
  ) %>% mutate(key = paste(CHROM, POS, REF_BASE, ALT_BASE, sep = "_"))

# new direct ND-Control trajectory-difference test
treatment_diff <- read_csv(
  paste0(git_path,"/general_analysis/vcftools/freq/logistic_trend/LD_WF_treatment_difference_sig_SNPs.csv.gz"),
  show_col_types = FALSE
) %>%
  mutate(key = paste(CHROM, POS, REF_BASE, ALT_BASE, sep = "_")) %>%
  select(
    key,
    extreme_count_treatment_diff = extreme_count,
    valid_count_treatment_diff = valid_count,
    p_treatment_diff_empirical = p_empirical,
    p_treatment_diff_BH = p_adj_BH
  )

# functional annotation
ann <- read_csv("significant_SNPs.annotation.csv",show_col_types = FALSE) %>%
  mutate(key = paste(CHROM, POS, REF, ALT, sep = "_"))

# add treatment-difference results to the large SNP table
sig_full <- sig %>%
  left_join(treatment_diff, by = "key") %>%
  mutate(
    SIG_treatment_diff =
      !is.na(p_treatment_diff_BH) &
      p_treatment_diff_BH < 0.05
  )

# combine annotation with significance and frequency-change information
ann_full <- ann %>%
  left_join(sig_full, by = "key", suffix = c("_ANN", ""))

# define SNP groups for downstream GO-term analysis
outdir <- "grouped_annotations"
groups <- list(
  only_in_ND_significant_positives = ann_full %>%
    filter(SIG_ND, !SIG_Control, Delta_gen1_gen3_ND > 0),
  only_in_Control_significant_positives = ann_full %>%
    filter(SIG_Control, !SIG_ND, Delta_gen1_gen3_Control > 0),
  only_in_ND_significant_negatives = ann_full %>%
    filter(SIG_ND, !SIG_Control, Delta_gen1_gen3_ND < 0),
  only_in_Control_significant_negatives = ann_full %>%
    filter(SIG_Control, !SIG_ND, Delta_gen1_gen3_Control < 0),
  in_both_TREATs_positive = ann_full %>%
    filter(SIG_ND, SIG_Control, Delta_gen1_gen3_ND > 0, Delta_gen1_gen3_Control > 0),
  in_both_TREATs_negative = ann_full %>%
    filter(SIG_ND, SIG_Control, Delta_gen1_gen3_ND < 0, Delta_gen1_gen3_Control < 0)
)

# summarise SNP annotations to one row per gene for each SNP group
gene_col <- "GENE_ID"
effect_col <- "EFFECT"
for (group_name in names(groups)) {
  # save SNP-level information for each group
  snp_annotations <- groups[[group_name]] %>%
    filter(
      !is.na(.data[[gene_col]]),
      .data[[gene_col]] != "",
      .data[[effect_col]] != "intergenic_region"
    ) %>%
    group_by(.data[[gene_col]], key) %>%
    summarise(
      across(
        c(EFFECT, description, confidence, IMPACT, AA_POS, HGVS_P),
        ~ paste(unique(na.omit(.x)), collapse = "|")
      ),
      across(
        c(
          starts_with("ALT_FREQ_gen"),
          starts_with("Delta_gen1_gen3"),
          Delta_diff_Control_vs_ND_gen1_gen3,
          p_treatment_diff_empirical,
          p_treatment_diff_BH,
          SIG_treatment_diff
        ),
        dplyr::first
      ),
      .groups = "drop"
    ) %>%
    arrange(.data[[gene_col]], key)
  write_tsv(
    snp_annotations,
    file.path(outdir, paste0(group_name, ".snps.tsv"))
  )
  
  gene_annotations <- snp_annotations %>%
    # collapse all SNP annotations per gene
    group_by(.data[[gene_col]]) %>%
    summarise(
      effects = paste(EFFECT, collapse = ";"),
      description = paste(description, collapse = ";"),
      confidence = paste(confidence, collapse = ";"),
      AA_pos = paste(AA_POS, collapse = ";"),
      HGVS_P = paste(HGVS_P, collapse = ";"),
      IMPACT = paste(IMPACT, collapse = ";"),
      n_snps = n_distinct(key),
      snps = paste(key, collapse = ";"),
      .groups = "drop"
    ) %>%
    arrange(.data[[gene_col]])
  # save one gene table per SNP group
  write_tsv(
    gene_annotations,
    file.path(outdir, paste0(group_name, ".genes.tsv"))
  )
  cat(group_name, ":", nrow(gene_annotations), "genes\n")
}

#only_in_ND_significant_positives : 32 genes                                                                                                                          
#only_in_Control_significant_positives : 14 genes                                                                                                                     
#only_in_ND_significant_negatives : 178 genes                                                                                                                         
#only_in_Control_significant_negatives : 113 genes                                                                                                                    
#in_both_TREATs_positive : 26 genes                                                                                                                                   
#in_both_TREATs_negative : 537 genes      
