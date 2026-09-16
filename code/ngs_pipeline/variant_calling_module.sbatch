#!/bin/bash

# Pipeline module script for variant calling orchestration (SLURM).
# Submits and chains multiple worker jobs:
#   1) Variant calling (bcftools mpileup/call/filter) [optionally per-chromosome array]
#   2) Optional VCF concatenation (if REGION=1)
#   3) bcftools stats on the final VCF
#   4) MultiQC summary report
#
# Required variables (provided by the pipeline/main script):
#   JOB_DIR, JOB_NAME   – Directory/name prefix for log files and jobid tracking
#   OUT                 – Output directory for VCFs and reports
#   CODE_DIR            – Base path containing workers/ (sbatch scripts)
#   REFERENCE_GENOME    – Reference FASTA (creates .fai if missing)
#   GLOB_PATTERN        – Glob expression selecting input BAMs (if ALIGNMENT_OUT not used)
#   VCF_NAME            – Base name/prefix for VCF-related outputs
#
# Optional variables:
#   ALIGNMENT_OUT       – If set to an existing directory, use BAMs from there and
#                         attach dependency on the recorded alignment jobid
#   REGION              – If set to 1, run variant calling as chromosome array and
#                         concatenate per-chromosome VCFs afterwards
#
# Resource/config variables:
#   VARIANT_CALLING_* (CORES/MEM/TIME_LIMIT), PRE_FILTERING_CORES
#   MERGE_VCF_* (CORES/MEM/TIME_LIMIT)        [only if REGION=1]
#   STATS_* (MEM/TIME_LIMIT)
#   MULTIQC_* (MEM/TIME_LIMIT/CONFIG)
#   MAIL_USER
#
# Output:
#   - Creates a BAM list file: ${OUT}/${VCF_NAME}.bamlist.txt
#   - Produces VCF(s) and derived stats in ${OUT}
#   - Writes MultiQC report to: ${OUT}/multiqc
#   - Writes SLURM stdout/stderr logs to: ${JOB_DIR}/
#
# Note:
#   Not meant to be run standalone without the pipeline environment — it expects
#   variables and paths prepared by the main pipeline wrapper.


#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
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
# This Block calls the variant calling worker
########################################################

mkdir -p "$JOB_DIR" "$OUT"

if [[ -n "${ALIGNMENT_OUT:-}" && -d "$ALIGNMENT_OUT" ]]; then
    echo "Using aligned BAMs from: $ALIGNMENT_OUT" >&2
    GLOB_PATTERN=( "${ALIGNMENT_OUT}/*.bam" )
    ALIGNMENT_JOBID=$(<"${JOB_DIR}/${JOB_NAME}_alignment_array.jobid")
else
    echo "Using BAMs matching pattern: $GLOB_PATTERN" >&2
    ALIGNMENT_JOBID=""
fi

FILES=( $GLOB_PATTERN )
BAMLIST="${OUT}/${VCF_NAME}.bamlist.txt"
: > "$BAMLIST"
for f in "${FILES[@]}"; do
  echo "$f" >> "$BAMLIST"
done

CPUS_PER_TASK=$(( VARIANT_CALLING_CORES + PRE_FILTERING_CORES * 2 ))

if [[ ! -f "${REFERENCE_GENOME}.fai" ]]; then
    module load samtools

    samtools faidx "$REFERENCE_GENOME"
fi

SBATCH_OPTS=()
if [[ "${REGION:-0}" -eq 1 ]]; then
  # Build chromosome list
  CHROM_LIST="${OUT}/${VCF_NAME}.chromosomes.txt"
  cut -f1 "${REFERENCE_GENOME}.fai" > "$CHROM_LIST"
  NCHR=$(wc -l < "$CHROM_LIST")
  ARRAY_RANGE_CHR="0-$((NCHR-1))"
  SBATCH_OPTS+=(--array="$ARRAY_RANGE_CHR")
  export CHROM_LIST REGION
fi

if [[ -n "$ALIGNMENT_JOBID" ]]; then
    ALIGNMENT_STATE="$(sacct -j "$ALIGNMENT_JOBID" --format=State --noheader | tr -d ' ')"
    if [[ "$ALIGNMENT_STATE" == *"FAILED"* ]]; then
        echo "ERROR: Something with alignment went wrong." >&2
        exit 1
    fi
    if [[ "$ALIGNMENT_STATE" != *"COMPLETED"* ]]; then
        SBATCH_OPTS+=(--dependency="afterok:${ALIGNMENT_JOBID}")
    fi
fi

export OUT REFERENCE_GENOME VARIANT_CALLING_CORES PRE_FILTERING_CORES \
    BAMLIST MAP_QUAL MAX_READ_DEPTH BASE_QUAL ANNOT_TAGS\
    VARIANT_QUAL MIN_BASE_DEPTH MAX_BASE_DEPTH MIN_ALT_DEPTH \
    VCF_NAME

VARCAL_JOBID=$(sbatch \
  --job-name="${JOB_NAME}_variant_calling" \
  --ntasks=1 \
  --nodes=1 \
  "${SBATCH_OPTS[@]}" \
  --partition=cpu \
  --auks=yes \
  --cpus-per-task="$CPUS_PER_TASK" \
  --mem="$VARIANT_CALLING_MEM" \
  --output="${JOB_DIR}/${JOB_NAME}_variant_calling_%A_%a.out" \
  --error="${JOB_DIR}/${JOB_NAME}_variant_calling_%A_%a.err" \
  --mail-type=ALL \
  --mail-user="$MAIL_USER" \
  --time="$VARIANT_CALLING_TIME_LIMIT" \
  --export=ALL \
  "$CODE_DIR"/workers/variant_calling_worker.sbatch | awk '{print $NF}')

########################################################

########################################################
# This Block calls the concat vcfs worker
########################################################

if [[ "${REGION:-0}" -eq 1 ]]; then
    VAR_MERGE_JOBID=$(sbatch \
    --job-name="${JOB_NAME}_concat_vcfs" \
    --ntasks=1 \
    --nodes=1 \
    --dependency=afterok:"$VARCAL_JOBID" \
    --partition=cpu \
    --auks=yes \
    --cpus-per-task="$MERGE_VCF_CORES" \
    --mem="$MERGE_VCF_MEM" \
    --output="${JOB_DIR}/${JOB_NAME}_concat_vcfs_%A.out" \
    --error="${JOB_DIR}/${JOB_NAME}_concat_vcfs_%A.err" \
    --mail-type=FAIL \
    --mail-user="$MAIL_USER" \
    --time="$MERGE_VCF_TIME_LIMIT" \
    --export=ALL,CHROM_LIST="$CHROM_LIST",VCF_PRE="${OUT}/${VCF_NAME}" \
    "$CODE_DIR"/workers/concat_vcfs_worker.sbatch | awk '{print $NF}')
fi

########################################################

########################################################
# This block calls the bcftools stats worker
########################################################

DEP_OPT=""
if [[ "${REGION:-0}" -eq 1 ]]; then
    DEP_OPT="--dependency=afterok:${VAR_MERGE_JOBID}"
else
    DEP_OPT="--dependency=afterok:${VARCAL_JOBID}"
fi

STATS_JOBID=$(sbatch \
    --job-name="${JOB_NAME}_bcftools_stats" \
    $DEP_OPT \
    --ntasks=1 \
    --nodes=1 \
    --partition=cpu \
    --auks=yes \
    --cpus-per-task=1 \
    --mem="$STATS_MEM" \
    --output="${JOB_DIR}/${JOB_NAME}_bcftools_stats_%j.out" \
    --error="${JOB_DIR}/${JOB_NAME}_bcftools_stats_%j.err" \
    --mail-type=FAIL \
    --mail-user="$MAIL_USER" \
    --time="$STATS_TIME_LIMIT" \
    --export=ALL,VCF_PRE="${OUT}/${VCF_NAME}" \
    "$CODE_DIR"/workers/bcftools_stats_worker.sbatch | awk '{print $NF}')

########################################################

########################################################
# This block calls the MultiQC worker
########################################################

MULTIQC_OUT="${OUT}/multiqc"
mkdir -p "$MULTIQC_OUT"

MULTIQC_JOBID=$(sbatch \
    --job-name="${JOB_NAME}_variant_calling_multiqc" \
    --dependency=afterok:${STATS_JOBID} \
    --ntasks=1 \
    --nodes=1 \
    --partition=cpu \
    --auks=yes \
    --cpus-per-task=1 \
    --mem="$MULTIQC_MEM" \
    --output="${JOB_DIR}/${JOB_NAME}_variant_calling_multiqc_%j.out" \
    --error="${JOB_DIR}/${JOB_NAME}_variant_calling_multiqc_%j.err" \
    --mail-type=FAIL \
    --mail-user="$MAIL_USER" \
    --time="$MULTIQC_TIME_LIMIT" \
    --export=ALL,INPUT="${OUT}/${VCF_NAME}.stats",OUT="$MULTIQC_OUT",GLOB_PATTERN="",NAME="${JOB_NAME}_variant_calling_multiqc",CONFIG="$MULTIQC_CONFIG" \
    "$CODE_DIR"/workers/multiqc_worker.sbatch | awk '{print $NF}')

########################################################
