#!/bin/bash

# Worker script for variant calling using bcftools.
# Must be launched by the main SLURM submission script, which provides
# all required variables and resource settings.
#
# Required variables:
#   BAMLIST               – File containing BAM paths (one per line)
#   REFERENCE_GENOME      – Reference genome FASTA
#   OUT                   – Output directory for VCF files
#   VCF_NAME              – Base name for output VCF
#   VARIANT_CALLING_CORES – Threads for bcftools mpileup
#   PRE_FILTERING_CORES   – Threads for bcftools call/filter
#   MAX_READ_DEPTH        – Max read depth per position
#   MAP_QUAL              – Minimum mapping quality
#   BASE_QUAL             – Minimum base quality
#   ANNOT_TAGS            – mpileup annotation tags (comma-separated, optional)
#   VARIANT_QUAL          – Minimum variant QUAL
#   MIN_BASE_DEPTH        – Minimum total depth (DP)
#   MAX_BASE_DEPTH        – Maximum total depth (DP)
#   MIN_ALT_DEPTH         – Minimum alternate allele depth
#
# Optional variables:
#   REGION                – If set to 1, process one region per array task
#   CHROM_LIST            – List of chromosomes/regions (used if REGION=1)
#
# Function:
#   Runs a bcftools mpileup -> call -> filter pipeline, optionally restricted
#   to a single chromosome/region when used in a SLURM array job.
#
# Output:
#   Compressed VCF written to: ${OUT}/${VCF_NAME}.vcf.gz
#   or (if REGION=1): ${OUT}/${VCF_NAME}.<REGION>.vcf.gz
#
# Note:
#   SLURM array settings (cpus-per-task, mem, time, array range)
#   must be configured in the main submission script.


set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load bcftools

ANNOT_TAGS="${ANNOT_TAGS:-AD,DP}"

ANNOT_TAGS="${ANNOT_TAGS// /}"
TAGS=",$ANNOT_TAGS,"
FILTER_EXPR="(QUAL < ${VARIANT_QUAL})"
if [[ "$TAGS" == *",DP,"* || "$TAGS" == *"FORMAT/DP"* ]]; then
  FILTER_EXPR+=" || (INFO/DP < ${MIN_BASE_DEPTH}) || (INFO/DP > ${MAX_BASE_DEPTH})"
fi
if [[ "$TAGS" == *",AD,"* || "$TAGS" == *"FORMAT/AD"* ]]; then
  FILTER_EXPR+=" || (MAX(FORMAT/AD[*:1]) < ${MIN_ALT_DEPTH})"
fi

# chromosome/region handling
REGION_FLAG="${REGION:-0}"
REGION_OPT=()
OUT_SUFFIX=""
if [[ "$REGION_FLAG" -eq 1 ]]; then
  CHR="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$CHROM_LIST")"
  if [[ -z "$CHR" ]]; then
    echo "ERROR: No chromosome found for task ${SLURM_ARRAY_TASK_ID} in ${CHROM_LIST}." >&2
    exit 1
  fi
  REGION_OPT=(-r "$CHR")
  OUT_SUFFIX=".${CHR}"
fi

bcftools mpileup "${REGION_OPT[@]}" \
  --threads "$VARIANT_CALLING_CORES" -d "$MAX_READ_DEPTH" -q "$MAP_QUAL" -Q "$BASE_QUAL" -a "$ANNOT_TAGS" -b "$BAMLIST" -f "$REFERENCE_GENOME" -Ou \
  | bcftools call --threads "$PRE_FILTERING_CORES" -mv -Ou \
  | bcftools filter --threads "$PRE_FILTERING_CORES" -e "$FILTER_EXPR" -Oz -o "${OUT}/${VCF_NAME}${OUT_SUFFIX}.vcf.gz"

