#!/bin/bash

# Worker script for generating bcftools statistics for a single VCF file.
# Must be launched by the main SLURM submission script, which provides
# all required variables and resource settings.
#
# Required variables:
#   VCF_PRE – Path and filename prefix of the input VCF
#             (expects ${VCF_PRE}.vcf.gz to exist)
#
# Function:
#   Runs bcftools stats on the specified compressed VCF file
#   and writes the resulting statistics to a .stats file with
#   the same prefix.
#
# Output:
#   bcftools statistics written to ${VCF_PRE}.stats
#
# Note:
#   Requires SLURM resources (cpus-per-task, mem, time, array range)
#   set by the main submission script.


set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load bcftools

bcftools stats "${VCF_PRE}.vcf.gz" > "${VCF_PRE}.stats"
