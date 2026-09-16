#!/bin/bash

# Pipeline module script for read trimming orchestration (SLURM).
# Submits and chains:
#   1) TrimGalore trimming array job
#   2) MultiQC summary report
#
# Required variables (provided by the pipeline/main script):
#   JOB_DIR, JOB_NAME – Directory/name prefix for SLURM logs and jobid tracking
#   CODE_DIR          – Base path containing workers/ (sbatch scripts)
#   GLOB_PATTERN      – Glob expression selecting input FASTQ files
#   OUT               – Output directory for trimmed reads and reports
#   ARRAY_RANGE       – SLURM array range for trimming jobs
#
# Trimming parameters (passed through to the worker):
#   TRIMMING_CORES, TRIMMING_MEM, TRIMMING_TIME_LIMIT
#   QUALITY, MIN_LENGTH, CLIP5, CLIP3, ADAPTER
#
# Output:
#   - Trimmed FASTQ files written to $OUT
#   - MultiQC report written to $OUT/multiqc
#   - SLURM stdout/stderr logs written to $JOB_DIR
#   - Trimming array jobid written to:
#       $JOB_DIR/${JOB_NAME}_trimming_array.jobid
#
# Note:
#   Not meant to be run standalone without the pipeline environment — it expects
#   variables and paths prepared by the main pipeline wrapper.


#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --time=00:05:00
#SBATCH --auks=yes

########################################################
# General script settings
########################################################

set -euo pipefail

########################################################

########################################################
# This Block calls the trimming worker
########################################################

mkdir -p "$JOB_DIR" "$OUT"

CPUS_PER_TASK=$(( 2 * TRIMMING_CORES + 3 ))               # peak threads = (2×TRIMMING_CORES)+3 for SE

TRIMMING_JOBID=$(sbatch \
    --job-name="${JOB_NAME}_trimming" \
    --auks=yes \
    --partition=cpu \
    --ntasks=1 \
    --array="$ARRAY_RANGE" \
    --cpus-per-task="$CPUS_PER_TASK" \
    --mem="$TRIMMING_MEM" \
    --nodes=1 \
    --output="$JOB_DIR/${JOB_NAME}_trimming_%A_%a.out" \
    --error="$JOB_DIR/${JOB_NAME}_trimming_%A_%a.err" \
    --mail-type=ALL \
    --mail-user="$MAIL_USER" \
    --time="$TRIMMING_TIME_LIMIT" \
    --export=ALL,GLOB_PATTERN="$GLOB_PATTERN",OUT="$OUT",CORES="$TRIMMING_CORES",CLIP5="$CLIP5",CLIP3="$CLIP3",QUALITY="$QUALITY",MIN_LENGTH="$MIN_LENGTH",ADAPTER="$ADAPTER" \
    "$CODE_DIR"/workers/trimming_worker.sbatch | awk '{print $NF}')

echo "$TRIMMING_JOBID" > "$JOB_DIR/${JOB_NAME}_trimming_array.jobid"

########################################################

########################################################
# This block calls the MultiQC worker
########################################################

MULTIQC_OUT="${OUT}/multiqc"
mkdir -p "$MULTIQC_OUT"

MULTIQC_JOBID=$(sbatch \
    --job-name="${JOB_NAME}_trimming_multiqc" \
    --dependency=afterok:${TRIMMING_JOBID} \
    --ntasks=1 \
    --nodes=1 \
    --partition=cpu \
    --auks=yes \
    --cpus-per-task=1 \
    --mem="$MULTIQC_MEM" \
    --output="${JOB_DIR}/${JOB_NAME}_trimming_multiqc_%j.out" \
    --error="${JOB_DIR}/${JOB_NAME}_trimming_multiqc_%j.err" \
    --mail-type=FAIL \
    --mail-user="$MAIL_USER" \
    --time="$MULTIQC_TIME_LIMIT" \
    --export=ALL,INPUT="$OUT",OUT="$MULTIQC_OUT",GLOB_PATTERN="",NAME="${JOB_NAME}_trimming_multiqc",CONFIG="$MULTIQC_CONFIG" \
    "$CODE_DIR"/workers/multiqc_worker.sbatch | awk '{print $NF}')

########################################################
