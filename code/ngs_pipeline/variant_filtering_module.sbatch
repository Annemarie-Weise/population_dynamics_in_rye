#!/bin/bash

# Pipeline module script for variant filtering (+QC) orchestration (SLURM).
# Submits and chains multiple worker jobs:
#   1) Variant filtering (vcftools recode) -> writes filtered VCF to the SAME directory as input VCF
#   2) Variant QC statistics (vcftools)    -> writes QC output files to ${OUT}
#
# Required variables:
#   JOB_DIR, JOB_NAME   – Directory/name prefix for log files and jobid tracking
#   OUT                 – Output directory for QC outputs and MultiQC report
#   CODE_DIR            – Base path containing workers/ (sbatch scripts)
#   MAIL_USER           – Email for SLURM notifications
#   VCF_IN              – Input VCF path (must be .vcf.gz)
#   MIN_DP              – Minimum depth per site
#   MAX_DP              – Maximum depth per site
#   MAF                 – Minor allele frequency threshold (use 0 for no MAF filtering)
#   MAX_MISSING         – Maximum missing data proportion allowed (use 0 to allow fully-missing sites)
#   MIN_Q               – Minimum site QUAL
#   RM_INDELS           – If 1, remove indels
#   RM_MULTIALLELIC     – If 1, keep only biallelic sites
#   QC_STATS_FLAG       - Should QC statistics such as TajimaD,... be calculated? (0 = no, 1 = yes)
#   TAJIMAD_WINDOW      – Window size for vcftools --TajimaD
#   PI_WINDOW           – Window size for vcftools --window-pi
#
# Resource/config variables:
#   VARIANT_FILTERING_MEM, VARIANT_FILTERING_TIME_LIMIT
#   VARIANT_QC_MEM, VARIANT_QC_TIME_LIMIT
#
# Output:
#   - Filtered VCF: same directory as $VCF_IN
#       ${VCF_DIR}/${JOB_NAME}_minDP${MIN_DP}_maxDP${MAX_DP}_maf${MAF}_maxMiss${MAX_MISSING}_minQ${MIN_Q}[...].vcf.gz
#   - vcftools QC outputs in ${OUT} with prefix:
#       ${OUT}/${PREFIX}.*  (e.g. .idepth, .imiss, .frq, .hwe, .Tajima.D, .windowed.pi, ...)
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

set -euo pipefail

mkdir -p "$JOB_DIR" "$OUT"

########################################################
# This Block calls the variant filtering worker
########################################################

VCF_DIR="$(dirname "$VCF_IN")"
PREFIX="${JOB_NAME}_minDP${MIN_DP}_maxDP${MAX_DP}_maf${MAF}_maxMiss${MAX_MISSING}_minQ${MIN_Q}"
if [[ "${RM_INDELS:-0}" -eq 1 ]]; then
  PREFIX="${PREFIX}_noIndels"
fi
if [[ "${RM_MULTIALLELIC:-0}" -eq 1 ]]; then
  PREFIX="${PREFIX}_biallelic"
fi

export VCF_IN MIN_DP MAX_DP MAF MAX_MISSING \
    MIN_Q RM_INDELS RM_MULTIALLELIC \
    OUT_VCF="${VCF_DIR}/${PREFIX}.vcf.gz"

VARFIL_JOBID=$(sbatch \
  --job-name="${JOB_NAME}_variant_filtering" \
  --ntasks=1 \
  --nodes=1 \
  --partition=cpu \
  --auks=yes \
  --cpus-per-task=1 \
  --mem="$VARIANT_FILTERING_MEM" \
  --output="${JOB_DIR}/${JOB_NAME}_variant_filtering_%A.out" \
  --error="${JOB_DIR}/${JOB_NAME}_variant_filtering_%A.err" \
  --mail-type=ALL \
  --mail-user="$MAIL_USER" \
  --time="$VARIANT_FILTERING_TIME_LIMIT" \
  --export=ALL \
  "$CODE_DIR"/workers/variant_filtering_worker.sbatch | awk '{print $NF}')

########################################################

########################################################
# This Block calls the variant qc worker
########################################################
if (( QC_STATS_FLAG == 1 )); then

    export VCF_OUT="${OUT}/${PREFIX}" \
        VCF_IN="${VCF_DIR}/${PREFIX}.vcf.gz" \
        TAJIMAD_WINDOW PI_WINDOW

    VARQC_JOBID=$(sbatch \
        --job-name="${JOB_NAME}_variant_qc" \
        --ntasks=1 \
        --nodes=1 \
        --dependency=afterok:"$VARFIL_JOBID" \
        --partition=cpu \
        --auks=yes \
        --cpus-per-task=1 \
        --mem="$VARIANT_QC_MEM" \
        --output="${JOB_DIR}/${JOB_NAME}_variant_qc_%A.out" \
        --error="${JOB_DIR}/${JOB_NAME}_variant_qc_%A.err" \
        --mail-type=FAIL \
        --mail-user="$MAIL_USER" \
        --time="$VARIANT_QC_TIME_LIMIT" \
        --export=ALL \
        "$CODE_DIR"/workers/variant_qc_worker.sbatch | awk '{print $NF}')

fi

########################################################
