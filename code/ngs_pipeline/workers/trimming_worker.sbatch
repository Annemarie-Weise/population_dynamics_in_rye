#!/bin/bash

# Worker script for per-file read trimming with TrimGalore.
# Launched only by the main SLURM script, which provides all variables.
#
# Required variables:
#   GLOB_PATTERN – Pattern to locate input FASTQ files
#   OUT          – Output directory for trimmed reads + FastQC
#   CORES        – Threads for TrimGalore
#   QUALITY      - PhredScore to trimm low Quality ends
#   MIN_LENGTH   - Discard reads that are shorter than this length
#   CLIP3        – Bases to clip from 3' end of the reads
#   CLIP5        – Bases to clip from 5' end of the reads
#
# Optional variables:
#   ADAPTER      - Adapter sequence for TrimGalore, if empty => autodetect
#
# Function:
#   Selects one FASTQ file based on $SLURM_ARRAY_TASK_ID and runs TrimGalore
#   with clipping, adapter/quality trimming, multithreading, and FastQC.
#
# Output:
#   Trimmed FASTQ + FastQC reports in $OUT.
#
# Note:
#   Requires SLURM resources (cpus-per-task, mem, time) set by the main script.


set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load TrimGalore
module load fastqc


FILES=( $GLOB_PATTERN )

idx=${SLURM_ARRAY_TASK_ID:?}
(( idx < ${#FILES[@]} )) || { echo "Index $idx out of ${#FILES[@]}"; exit 1; }

FILE="${FILES[$idx]}"
BASENAME=$(basename "$FILE")
BASENAME="${BASENAME%.fastq.gz}"
BASENAME="${BASENAME%.fq.gz}"
BASENAME="${BASENAME%.fastq}"
BASENAME="${BASENAME%.fq}"

CMD=( trim_galore
  --output_dir "$OUT"
  --cores "$CORES"
  -q "$QUALITY"
  --trim-n
  --length "$MIN_LENGTH"
  --fastqc
  --basename "$BASENAME"
)

if (( CLIP5 > 0 )); then
  CMD+=( --clip_R1 "$CLIP5" )
fi

if (( CLIP3 > 0 )); then
  CMD+=( --three_prime_clip_R1 "$CLIP3" )
fi

if [[ -n "${ADAPTER:-}" ]]; then
  CMD+=( -a "$ADAPTER" )
fi

CMD+=( "$FILE" )

"${CMD[@]}"

mv "${OUT}/${BASENAME}_trimmed.fq.gz" "${OUT}/${BASENAME}_trim.fq.gz"
