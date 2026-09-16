#!/bin/bash

# Worker script for concatenating region-split VCF files.
# Must be launched by the main SLURM submission script, which provides
# all required variables and resource settings.
#
# Required variables:
#   MERGE_VCF_CORES – Number of threads to use for bcftools concat
#   CHROM_LIST     – Path to a sorted list of chromosomes/regions
#                     (one entry per line, order defines concat order)
#   VCF_PRE        – Path and filename prefix of the split VCF files
#                     (expects ${VCF_PRE}.<region>.vcf.gz files)
#
# Function:
#   Builds an ordered list of region-specific VCF files based on CHROM_LIST
#   and concatenates them using `bcftools concat` into a single compressed
#   VCF file.
#
# Output:
#   Concatenated VCF written to ${VCF_PRE}.vcf.gz
#
# Note:
#   Requires SLURM resources (cpus-per-task, mem, time, array range)
#   set by the main submission script.


set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load bcftools

bcftools concat --threads $MERGE_VCF_CORES -n -Oz -o "${VCF_PRE}.vcf.gz" \
  $(sed "s|^|${VCF_PRE}.|; s|$|.vcf.gz|" "$CHROM_LIST")

