#!/bin/bash

# Pipeline module script for read alignment orchestration (SLURM).
# Submits and chains:
#   1) Alignment array job (per input FASTQ)
#   2) MultiQC summary report
#   3) Optional BAM merging and indexing
#
# Required variables:
#   JOB_DIR, JOB_NAME – Directory/name prefix for SLURM logs and jobid tracking
#   CODE_DIR          – Base path containing workers/ (sbatch scripts)
#   OUT               – Output directory for BAMs and reports
#   GLOB_PATTERN      – Glob expression selecting input FASTQ files (if TRIMMING_OUT not used)
#   ARRAY_RANGE       – SLURM array range for alignment jobs
#   REFERENCE_GENOME  – Reference FASTA (used to build index if INDEX is empty)
#
# Optional variables:
#   TRIMMING_OUT      – If set to an existing directory, use FASTQs from there and
#                       attach dependency on $JOB_DIR/${JOB_NAME}_trimming_array.jobid
#   INDEX             – Path to minimap2 index (.mmi); if empty, index is created in $OUT
#   MERGING_FLAG      – If set to 1, submit a merge/index job after alignment
#
# Resource/config variables:
#   ALIGNMENT_CORES, SAMTOOLS_CORES, ALIGNMENT_MEM, ALIGNMENT_TIME_LIMIT
#   MAP_QUAL
#   MULTIQC_MEM, MULTIQC_TIME_LIMIT, MULTIQC_CONFIG
#   MERGING_CORES, MERGING_MEM, MERGING_TIME_LIMIT, MERGING_BATCH_SIZE   [if MERGING_FLAG=1]
#   INDEX_SIZE, MAIL_USER
#
# Output:
#   - Alignment array jobid written to:
#       $JOB_DIR/${JOB_NAME}_alignment_array.jobid
#   - BAMs and alignment QC files written to $OUT
#   - MultiQC report written to $OUT/multiqc
#   - If MERGING_FLAG=1: merged BAM written under $OUT/merged_bam
#   - SLURM stdout/stderr logs written to $JOB_DIR
#
# Note:
#   Not meant to be run standalone without the pipeline environment — it expects
#   variables and paths prepared by the main pipeline wrapper.


#SBATCH --partition=cpu
#SBATCH --cpus-per-task=3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL
#SBATCH --auks=yes

########################################################
# General script settings
########################################################

set -euo pipefail
# init module environment
. /etc/profile.d/modules.sh

########################################################

########################################################
# This Block calls the alignment worker
########################################################

mkdir -p "$JOB_DIR" "$OUT"

if [[ -n "${TRIMMING_OUT:-}" && -d "$TRIMMING_OUT" ]]; then
    echo "Using trimmed FASTQs from: $TRIMMING_OUT" >&2
    GLOB_PATTERN=( "${TRIMMING_OUT}/*.fq" "${TRIMMING_OUT}/*.fq.gz" )
    TRIMMING_JOBID=$(<"${JOB_DIR}/${JOB_NAME}_trimming_array.jobid")
else
    echo "Using FASTQs matching pattern: $GLOB_PATTERN" >&2
    TRIMMING_JOBID=""
fi

CPUS_PER_TASK=$(( ALIGNMENT_CORES + 2 + SAMTOOLS_CORES * 4 ))               # peak threads ≈ ALIGNMENT_CORES + 2 for I/O + SAMTOOLS_CORES for 4 samtools operations

if [ -z "$INDEX" ]; then #change later if needed
        module load minimap2
        INDEX="${OUT}/$(basename "$REFERENCE_GENOME").mmi"
        minimap2 -t 3 -I "$INDEX_SIZE" -d "$INDEX" "$REFERENCE_GENOME"
fi

DEP_TRIM=""
if [[ -n "$TRIMMING_JOBID" ]]; then
    TRIMMING_STATE="$(sacct -j "$TRIMMING_JOBID" --format=State --noheader | tr -d ' ')"
    if [[ "$TRIMMING_STATE" == *"FAILED"* ]]; then
        echo "ERROR: Something with trimming went wrong." >&2
        exit 1
    fi
    if [[ "$TRIMMING_STATE" != *"COMPLETED"* ]]; then
        DEP_TRIM="--dependency=afterok:${TRIMMING_JOBID}"
    fi
fi


ALIGN_JOBID=$(sbatch \
  --job-name="${JOB_NAME}_alignment" \
  $DEP_TRIM \
  --ntasks=1 \
  --nodes=1 \
  --partition=cpu \
  --auks=yes \
  --array="$ARRAY_RANGE" \
  --cpus-per-task="$CPUS_PER_TASK" \
  --mem="$ALIGNMENT_MEM" \
  --output="${JOB_DIR}/${JOB_NAME}_alignment_%A_%a.out" \
  --error="${JOB_DIR}/${JOB_NAME}_alignment_%A_%a.err" \
  --mail-type=ALL \
  --mail-user="$MAIL_USER" \
  --time="$ALIGNMENT_TIME_LIMIT" \
  --export=ALL,GLOB_PATTERN="${GLOB_PATTERN[*]}",OUT="$OUT",REF="$REFERENCE_GENOME",ALIGNMENT_CORES="$ALIGNMENT_CORES",MAP_QUAL="$MAP_QUAL",SAMTOOLS_CORES="$SAMTOOLS_CORES" \
  "$CODE_DIR"/single_use_scripts/bwa_test.sbatch | awk '{print $NF}')

echo "$ALIGN_JOBID" > "$JOB_DIR/${JOB_NAME}_alignment_array.jobid"

########################################################

########################################################
# This block calls the MultiQC worker
########################################################

MULTIQC_OUT="${OUT}/multiqc"
mkdir -p "$MULTIQC_OUT"

MULTIQC_JOBID=$(sbatch \
    --job-name="${JOB_NAME}_alignment_multiqc" \
    --dependency=afterok:${ALIGN_JOBID} \
    --ntasks=1 \
    --nodes=1 \
    --partition=cpu \
    --auks=yes \
    --cpus-per-task=1 \
    --mem="$MULTIQC_MEM" \
    --output="${JOB_DIR}/${JOB_NAME}_alignment_multiqc_%j.out" \
    --error="${JOB_DIR}/${JOB_NAME}_alignment_multiqc_%j.err" \
    --mail-type=FAIL \
    --mail-user="$MAIL_USER" \
    --time="$MULTIQC_TIME_LIMIT" \
    --export=ALL,INPUT="$OUT",OUT="$MULTIQC_OUT",GLOB_PATTERN="",NAME="${JOB_NAME}_alignment_multiqc",CONFIG="$MULTIQC_CONFIG" \
    "$CODE_DIR"/workers/multiqc_worker.sbatch | awk '{print $NF}')

########################################################

########################################################
# This block calls the merging and indexing worker
########################################################

if [[ "$MERGING_FLAG" -eq 1 ]]; then

  MERGED="${OUT}/merged_bam"
  mkdir -p "$MERGED"

  CPUS_PER_TASK=$(( MERGING_CORES + 2 ))   # CPU calc for the cluster (peak threads = cores + I/O threads)

  sbatch \
    --job-name="${JOB_NAME}_merge" \
    --dependency=afterok:${ALIGN_JOBID} \
    --partition=cpu \
    --ntasks=1 \
    --nodes=1 \
    --auks=yes \
    --cpus-per-task="$CPUS_PER_TASK" \
    --mem="$MERGING_MEM" \
    --mail-type=FAIL \
    --mail-user="$MAIL_USER" \
    --time="$MERGING_TIME_LIMIT" \
    --output="${JOB_DIR}/${JOB_NAME}_merge_%j.out" \
    --error="${JOB_DIR}/${JOB_NAME}_merge_%j.err" \
    --export=ALL,INPUT="$OUT",MERGED="$MERGED",NAME_MERGED_BAM="$JOB_NAME",MERGING_CORES="$MERGING_CORES",MERGING_BATCH_SIZE="$MERGING_BATCH_SIZE" \
    "$CODE_DIR"/workers/merge_bams_worker.sbatch
fi

########################################################
