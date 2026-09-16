library(readr)
library(dplyr)
library(tidyr)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/annotation_and_GO"))

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
# GO-term enrichment for significant SNP gene groups
################################################################################
# - build a GO TERM2GENE table from the rye feature annotation
# - run Fisher's exact test GO enrichment for all groups
################################################################################


################################################################################
# Helper function
################################################################################

# Go-Term Enrichment -> Fisher's exact test
run_go_hypergeom <- function(target_genes, go_terme, background_genes, alpha = 0.05) {
  # Keep unique target/background genes and restrict targets to the background
  target_genes <- unique(target_genes)
  background_genes <- unique(background_genes)
  target_genes <- intersect(target_genes, background_genes)
  
  # all genes
  N <- length(background_genes)  
  # all traget genes
  n <- length(target_genes)  
  # unique Go terms
  go_terms <- unique(go_terme$GO) 
  
  result_list <- list()
  for (go in go_terms) {
    genes_with_go <- unique(go_terme$gene_id[go_terme$GO == go])
    
    # background genes with this GO term
    K <- length(genes_with_go)
    # target genes with this GO term
    x <- length(intersect(target_genes, genes_with_go))
    
    # Test only GO terms represented in the target set
    if (x > 0) {
      p_value <- phyper(x - 1, K, N - K, n,lower.tail = FALSE)
      result_list[[go]] <- data.frame(
        GO = go,
        x = x,
        K = K,
        GeneRatio = paste0(x, "/", n),
        BgRatio = paste0(K, "/", N),
        p_value = p_value,
        genes_in_target = paste(sort(intersect(target_genes, genes_with_go)), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
  results <- bind_rows(result_list)
  
  # Adjust p-values and sort by significance
  results <- results %>%
    mutate(
      p_adjust = p.adjust(p_value, method = "BH"),
      significant = p_adjust < alpha
    ) %>%
    arrange(p_adjust, p_value)
  return(results)
}


################################################################################
# Build TERM2GENE table
################################################################################

# read rye feature annotation table
feature <- read_csv(featuretable_file, show_col_types = FALSE,na = c("", "NA", "NaN", "nan"))

# split multiple GO terms per gene into separate rows
gene2go <- feature %>%
  filter(!is.na(GO), GO != "") %>%
  select(gene_id, GO) %>%
  separate_rows(GO, sep = "\\|") %>%
  filter(!is.na(GO), GO != "") %>%
  distinct(gene_id, GO)

# save GO-to-gene mapping
term2gene <- gene2go %>%
  select(GO, gene_id)
write_tsv(term2gene, "Lo7V3_TERM2GENE.tsv")


################################################################################
# Run GO enrichment for grouped gene lists
################################################################################

# read TERM2GENE table and define background genes
go_terme <- read_tsv("Lo7V3_TERM2GENE.tsv")
background_genes <- go_terme %>%
  distinct(gene_id) %>%
  pull(gene_id)
n_background_genes <- n_distinct(go_terme$gene_id)
n_background_go_terms <- n_distinct(go_terme$GO)

# find grouped gene annotation files
gene_files <- list.files("grouped_annotations", pattern = "\\.genes\\.tsv$", full.names = TRUE)

# run enrichment separately for each gene group
for (gene_file in gene_files) {
  group_name <- basename(gene_file) %>%
    sub("\\.genes\\.tsv$", "", .)
  
  # read target genes for this group
  target_genes <- read_tsv(
    gene_file,
    show_col_types = FALSE
  ) %>%
    pull(GENE_ID) %>%
    unique()
  
  # test GO-term enrichment
  res <- run_go_hypergeom(target_genes, go_terme, background_genes)
  write_tsv(
    res,
    file.path(paste0(group_name, "_GO_res.tsv"))
  )
  target_genes_used <- intersect(target_genes, background_genes)
  
  n_genes_used <- length(target_genes_used)
  
  n_go_terms_hit <- go_terme %>%
    filter(gene_id %in% target_genes_used) %>%
    summarise(n = n_distinct(GO)) %>%
    pull(n)
  cat(
    group_name, ":",
    length(target_genes), "genes total,",
    n_genes_used, "GO-annotated genes used,",
    n_go_terms_hit, "GO terms with hits\n"
  )
}
cat(
  "Background:",
  n_background_genes, "genes,",
  n_background_go_terms, "unique GO terms\n"
)

#in_both_TREATs_negative : 537 genes total, 287 GO-annotated genes used, 208 GO terms with hits                                                                                                            
#in_both_TREATs_positive : 26 genes total, 8 GO-annotated genes used, 14 GO terms with hits                                                                                                                
#only_in_Control_significant_negatives : 113 genes total, 57 GO-annotated genes used, 71 GO terms with hits                                                                                                
#only_in_Control_significant_positives : 14 genes total, 7 GO-annotated genes used, 10 GO terms with hits                                                                                                  
#only_in_ND_significant_negatives : 178 genes total, 90 GO-annotated genes used, 89 GO terms with hits                                                                                                     
#only_in_ND_significant_positives : 32 genes total, 17 GO-annotated genes used, 25 GO terms with hits        
#Background: 27058 genes, 1563 unique GO terms


