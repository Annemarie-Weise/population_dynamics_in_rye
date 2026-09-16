#!/bin/bash

# Pipeline module script for FastQC quality control (SLURM).
# Submits and chains:
#   1) FastQC array job (per input FASTQ)
#   2) MultiQC summary report
#
# Required variables (provided by the pipeline/main script):
#   JOB_DIR, JOB_NAME – Directory/name prefix for SLURM logs
#   CODE_DIR          – Base path containing workers/ (sbatch scripts)
#   GLOB_PATTERN      – Glob expression selecting input FASTQ files
#   OUT               – Output directory for FastQC and MultiQC results
#   ARRAY_RANGE       – SLURM array range for FastQC jobs
#
# Resource/config variables:
#   FASTQC_CORES, FASTQC_MEM, FASTQC_TIME_LIMIT
#   MULTIQC_MEM, MULTIQC_TIME_LIMIT, MULTIQC_CONFIG
#   MAIL_USER
#
# Output:
#   - FastQC reports written to $OUT
#   - MultiQC report written to $OUT/multiqc
#   - SLURM stdout/stderr logs written to $JOB_DIR
#
# Note:
#   Not meant to be run standalone without the pipeline environment — it expects
#   variables and paths prepared by the main pipeline wrapper.


#SBATCH --partition=cpu
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --auks=yes

########################################################
# General script settings
########################################################

set -euo pipefail

########################################################

########################################################
# This Block calls the FastQC worker
########################################################

mkdir -p "$JOB_DIR" "$OUT"

FASTQC_JOBID=$(sbatch \
  --job-name="${JOB_NAME}_fastqc" \
  --ntasks=1 \
  --nodes=1 \
  --partition=cpu \
  --auks=yes \
  --array="$ARRAY_RANGE" \
  --cpus-per-task="$FASTQC_CORES" \
  --mem="$FASTQC_MEM" \
  --output="${JOB_DIR}/${JOB_NAME}_fastqc_%A_%a.out" \
  --error="${JOB_DIR}/${JOB_NAME}_fastqc_%A_%a.err" \
  --mail-type=ALL \
  --mail-user="$MAIL_USER" \
  --time="$FASTQC_TIME_LIMIT" \
  --export=ALL,GLOB_PATTERN="$GLOB_PATTERN",OUT="$OUT" \
  "$CODE_DIR"/workers/fastqc_worker.sbatch | awk '{print $NF}')

########################################################

########################################################
# This block calls the MultiQC worker
########################################################

MULTIQC_OUT="${OUT}/multiqc"
mkdir -p "$MULTIQC_OUT"

MULTIQC_JOBID=$(sbatch \
    --job-name="${JOB_NAME}_fastqc_multiqc" \
    --dependency=afterok:${FASTQC_JOBID} \
    --ntasks=1 \
    --nodes=1 \
    --partition=cpu \
    --auks=yes \
    --cpus-per-task=1 \
    --mem="$MULTIQC_MEM" \
    --output="${JOB_DIR}/${JOB_NAME}_fastqc_multiqc_%j.out" \
    --error="${JOB_DIR}/${JOB_NAME}_fastqc_multiqc_%j.err" \
    --mail-type=FAIL \
    --mail-user="$MAIL_USER" \
    --time="$MULTIQC_TIME_LIMIT" \
    --export=ALL,INPUT="$OUT",OUT="$MULTIQC_OUT",GLOB_PATTERN="",NAME="${JOB_NAME}_fastqc_multiqc",CONFIG="$MULTIQC_CONFIG" \
    "$CODE_DIR"/workers/multiqc_worker.sbatch | awk '{print $NF}')

########################################################
