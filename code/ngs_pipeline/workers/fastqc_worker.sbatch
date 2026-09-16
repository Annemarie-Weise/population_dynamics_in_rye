#!/bin/bash

# Worker script for per-file FastQC analysis.
# Not meant to be run directly — it must be launched by the main SLURM
# submission script, which provides all required variables and resources.
#
# Required variables:
#   GLOB_PATTERN     – Pattern matching the FASTQ files to process
#   OUT              – Output directory where FastQC reports are written
#
# Function:
#   Processes one FASTQ file per SLURM array task using FastQC
#   with $SLURM_CPUS_PER_TASK threads.
#
# Output:
#   $OUT/<sample>_fastqc.html
#   $OUT/<sample>_fastqc.zip
#
# Note:
#   SLURM array settings (cpus-per-task, mem, time, array range)
#   must be configured in the main submission script.

set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load fastqc

FILES=( $GLOB_PATTERN )

idx=${SLURM_ARRAY_TASK_ID:?}
(( idx < ${#FILES[@]} )) || { echo "Index $idx out of ${#FILES[@]}"; exit 1; }

FILE="${FILES[$idx]}"

fastqc -o "$OUT" --threads "$SLURM_CPUS_PER_TASK" --noextract "$FILE"
