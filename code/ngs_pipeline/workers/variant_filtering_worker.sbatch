#!/bin/bash

# Worker script for variant filtering using vcftools.
# Must be launched by the module script, which provides all variables.
#
# Required variables:
#   VCF_IN              – Input VCF (.vcf.gz)
#   OUT_VCF             - Full name (with path) for output
#   MIN_DP              – Minimum depth per site
#   MAX_DP              – Maximum depth per site
#   MAF                 – Minor allele frequency
#   MAX_MISSING         – Maximum missingness
#   MIN_Q               – Minimum QUAL
#   RM_INDELS           - 1 = remove indels, 0 = keep indels
#   RM_MULTIALLELIC     - 1 = remove multiallelic sites, 0 = keep multiallelic sites

set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load vcftools

OPTS=(
  --gzvcf "$VCF_IN"
  --minDP "$MIN_DP"
  --maxDP "$MAX_DP"
  --maf "$MAF"
  --max-missing "$MAX_MISSING"
  --minQ "$MIN_Q"
  --recode
  --recode-INFO-all
)

if [[ "${RM_INDELS:-0}" -eq 1 ]]; then
  OPTS+=(--remove-indels)
fi
if [[ "${RM_MULTIALLELIC:-0}" -eq 1 ]]; then
  OPTS+=(--min-alleles 2 --max-alleles 2)
fi

vcftools "${OPTS[@]}" --stdout | gzip -c > "${OUT_VCF}"
