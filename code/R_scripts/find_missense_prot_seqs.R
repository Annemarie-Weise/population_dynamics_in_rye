library(Biostrings)
library(data.table)
library(rprojroot)

git_path <- find_root(has_dir(".git"))
setwd(paste0(git_path, "/general_analysis/annotation_and_GO"))

# Path to protein FASTA from: https://zenodo.org/records/15801361
# - SECALE.CEREALE.Lo7V3.pgsb.r2.Feb2025_all.aa.fa
protein_fa_path <- readline(
  prompt = "Please enter the path to SECALE.CEREALE.Lo7V3.pgsb.r2.Feb2025_all.aa.fa: "
)
if (!file.exists(protein_fa_path)) {
  stop("File not found: ", protein_fa_path)
}
protein_fa_file <- normalizePath(protein_fa_path)


################################################################################
# Create alternative protein sequences for missense SNPs
################################################################################
# - read missense SNP annotation tables
# - extract the corresponding reference protein sequences from the rye FASTA
# - create alternative protein sequences carrying missense variant
################################################################################


################################################################################
# Helper functions
################################################################################

# convert three-letter amino acid codes to one-letter codes
convert_aa <- function(x) {
  if (is.na(x)) return(NA)
  if (nchar(x) == 1) return(x)
  aa3_to_aa1[x]
}

# parse HGVS protein notation into REF_AA, POS, and ALT_AA
parse_hgvs_p <- function(hgvs) {
  x <- hgvs
  x <- gsub("^p\\.", "", x)
  x <- gsub("[()]", "", x)
  m <- regexec("^([A-Za-z*]{1,3})([0-9]+)([A-Za-z*]{1,3})$", x)
  r <- regmatches(x, m)[[1]]
  if (length(r) != 4) {
    return(data.table(ref = NA, pos = NA, alt = NA))
  }
  data.table(
    ref = convert_aa(r[2]),
    pos = as.integer(r[3]),
    alt = convert_aa(r[4])
  )
}

# replace the REF_AA with ALT_AA and return the variant protein sequence
make_alt_sequence <- function(seq, pos, ref, alt, gene_id) {
  seq <- as.character(seq)
  if (is.na(pos) || is.na(ref) || is.na(alt)) {
    warning(paste("Could not parse HGVS_P for gene", gene_id))
    return(NA)
  }
  
  # check that the HGVS reference amino acid matches the FASTA sequence
  observed_ref <- substr(seq, pos, pos)
  if (observed_ref != ref) {
    warning(paste(
      "Reference mismatch for", gene_id,"at position", pos,
      ": HGVS says", ref,"but FASTA has", observed_ref
    ))
    return(NA)
  }
  paste0(
    substr(seq, 1, pos - 1),
    alt,
    substr(seq, pos + 1, nchar(seq))
  )
}

# add reference and alternative protein sequences to one missense SNP table
process_missense_file <- function(input_file, output_file) {
  tab <- fread(input_file)
  
  tab[, GENE_ID := paste0(GENE_ID, ".1")]
  tab[, AA_seq_REF := NA_character_]
  tab[, AA_seq_ALT := NA_character_]
  
  for (i in seq_len(nrow(tab))) {
    gene <- tab$GENE_ID[i]
    hgvs <- tab$HGVS_P[i]
    
    ref_seq <- sub("\\*$", "", as.character(faa[[gene]]))
    
    parsed <- parse_hgvs_p(hgvs)
    tab$AA_seq_REF[i] <- ref_seq
    tab$AA_seq_ALT[i] <- make_alt_sequence(
      seq = ref_seq,
      pos = parsed$pos,
      ref = parsed$ref,
      alt = parsed$alt,
      gene_id = gene
    )
  }
  
  fwrite(tab, output_file, sep = "\t")
}


################################################################################
# Read protein FASTA and define amino acid code conversion
################################################################################

# reference protein sequences
faa <- readAAStringSet(protein_fa_file)

# three-letter to one-letter amino acid code map
aa3_to_aa1 <- c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
  Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
  Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
  Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
  Ter = "*", Stop = "*"
)

# create alternative protein sequences for missense SNP groups
process_missense_file(
  input_file = "ND_only_missense/ND_only_sig_neg_prot_missense.csv",
  output_file = "ND_only_missense/ND_only_neg_missense_vars.tsv"
)
process_missense_file(
  input_file = "ND_only_missense/ND_only_sig_pos_prot_missense.csv",
  output_file = "ND_only_missense/ND_only_pos_missense_vars.tsv"
)