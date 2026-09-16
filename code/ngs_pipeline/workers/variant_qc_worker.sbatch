#!/bin/bash
#
# Worker: VCF statistics using vcftools.
#
# Required variables:
#   VCF_OUT          – Prefix for output
#   VCF_IN           – VCF file path (vcf.gz)
#   TAJIMAD_WINDOW   – Window size for --TajimaD
#   PI_WINDOW        – Window size for --window-pi
#
# Output:
#   Produces vcftools output files:
#     *.ldepth, *.sitemean, *.idepth, *.lmiss, *.imiss, *.frq, *.hwe, *.Tajima.D, *.windowed.pi
#
# Note:
#   Not meant to be run standalone without the pipeline environment — it expects
#   variables and paths prepared by the main pipeline wrapper.

set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load vcftools

vcftools --gzvcf "$VCF_IN" --site-depth --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --site-mean-depth --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --depth --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --missing-site --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --missing-indv --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --freq --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --hardy --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --TajimaD "$TAJIMAD_WINDOW" --out "$VCF_OUT"
vcftools --gzvcf "$VCF_IN" --window-pi "$PI_WINDOW" --out "$VCF_OUT"

gzip "${VCF_OUT}".*

