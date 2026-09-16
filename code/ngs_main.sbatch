#!/bin/bash

#SBATCH --job-name=submit_ngs
#SBATCH --partition=cpu
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=/path/where/job_log/should/be/stored/submit_ngs_%j.out   # change
#SBATCH --error=/path/where/job_log/should/be/stored/submit_ngs_%j.err    # change
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yourname@ipk-gatersleben.de                             # change
#SBATCH --auks=yes


########################################################
# General Information
########################################################
#
# Main SLURM submission script for an NGS mini-pipeline.
# This script submits module jobs via sbatch (it does not run the heavy
# computation itself).
#
# Main modules (dependency chain if enabled together):
#   1) TRIMMING          – TrimGalore trimming + MultiQC
#   2) ALIGNMENT         – minimap2/samtools alignment + MultiQC (+ optional merging)
#   3) VARIANT_CALLING   – bcftools calling (+ optional per-region mode), stats + MultiQC
#
# Optional independent modules (do not affect dependencies):
#   - INITIAL_FASTQC     – FastQC on input FASTQs + MultiQC
#   - INITIAL_MULTIQC    – MultiQC summary of existing QC reports
#   - VARIANT_FILTERING  - Filtering of a defined VCF file and calculating of many vcftools statistics
#
# Behavior:
#   - Enabled main modules are chained automatically (TRIMMING -> ALIGNMENT -> VARIANT_CALLING).
#   - Inputs are collected via the GLOB_PATTERN_* settings.
#   - ARRAY_RANGE is derived from the number of matched files (limited by ARRAY_CAP).
#   - Outputs are written under ALL_OUT into per-module folders:
#       ${ALL_OUT}/${JOB_NAME}_<module>/
#     SLURM logs are written to JOB_DIR.
#
# Notes:
#   - If VARIANT_CALLING is run on external BAMs, they must be coordinate-sorted and indexed (bai/csi).
#   - It is recommended to run modules independently when fine-grained control over
#     intermediate results is desired, allowing quality-based filtering between steps.
#
########################################################


########################################################
# Select which modules of the pipeline should run.
# Set each flag to:  0 = disabled   |   1 = enabled
########################################################

# Generate FastQC reports from FASTQ files and summarize them with MultiQC
INITIAL_FASTQC=0

# Summarize existing FastQC reports with MultiQC (no new FastQC run)
INITIAL_MULTIQC=0

# Trim FASTQ files and summarize trimming reports with MultiQC
TRIMMING=0

# Align reads, generate sorted BAM files, and summarize mapping reports with MultiQC
ALIGNMENT=0

# Perform variant calling
VARIANT_CALLING=0

# Filter and run vcftools statistics if wanted for a specific VCF file
VARIANT_FILTERING=0

# Base directory containing the pipeline scripts and module directory
CODE_DIR="/path/to/population_dynamics_in_rye/code/ngs_pipeline"

########################################################


########################################################
# Dataset-specific settings (must be adjusted for every run)
########################################################

# Short job name used for Slurm logs and module output folders
JOB_NAME="Testjob"

# Email address for Slurm notifications
MAIL_USER="yourname@ipk-gatersleben.de"

# Base output directory for all module output
# Each module creates a subfolder: ${JOB_NAME}_<module>
# WARNING: If JOB_NAME or this directory is not changed, old output may be overwritten
ALL_OUT="/path/to/output/directory"

# Directory for log files (will be created if missing)
JOB_DIR="${ALL_OUT}/job"

# Input path, depending on which modules are used:
# - Quality reports (for INITIAL_MULTIQC)
# - FASTQ files (for INITIAL_FASTQC, TRIMMING, ALIGNMENT)
# - BAM files (for VARIANT_CALLING)
INPUT="/path/to/input/directory"

# File pattern for INPUT:
# Modify according to file types and directory structure
# - Pattern to subfolders of INPUT which contain report files to summarize (leave empty if the files are in INPUT)
GLOB_PATTERN_REPORTS="$INPUT"/Reports
# - FASTQ pattern for FASTQC, TRIMMING, ALIGNMENT
GLOB_PATTERN_FASTQ=(
  "${INPUT}/*.fastq"
  "${INPUT}/*.fastq.gz"
  "${INPUT}/*.fq"
  "${INPUT}/*.fq.gz"
)
# - BAM pattern for VARIANT_CALLING
GLOB_PATTERN_BAM=(
    "${INPUT}/bams_1/*.bam"
    "${INPUT}/bams_2/*.bam"
)
# VCF/VCF.GZ file used for VARIANT_FILTERING module
VCF_IN="${INPUT}/path/from/inputs/to/your_vcf.vcf.gz"

# Path to the reference genome (FASTA or FASTA.GZ)
# Ignored for ALIGNMENT if a minimap2 index (.mmi) is provided below
# Is needed for VARIANT_CALLING
# -> fai file be generated in the same dir if not already present so make sure the pipline can write in the dir
REFERENCE_GENOME="/path/to/reference/genome.fasta"

# Path to an existing minimap2 index
# If empty, the index will be generated automatically
INDEX="/path/to/reference/index/genome.fasta.mmi"

########################################################


########################################################
# Optional Module settings
########################################################

# Trimming settings
#------------------
QUALITY=20                                              # Minimum quality for read end trimming (Phred33)
MIN_LENGTH=20                                           # Minimum length for reads to be kept
CLIP3=0                                                 # Clip bases from 3'-end (should be 0 if not wanted)
CLIP5=0                                                 # Clip bases from 5'-end (should be 0 if not wanted)
ADAPTER=""                                              # Adapter sequence to be detected z.B.: "AGATCGGAAGAGC", leave empty for autodetection

# Alignment
#----------
MAP_QUAL=30                                             # Minimum mapping quality for reads to be retained (30 is standard)
MERGING_FLAG=0                                          # Merging flag: 0 = no merging, 1 = merge BAM files

# Variant calling
#----------------
VAR_MAP_QUAL=30                                         # Minimum MAPQ retained after filtering
BASE_QUAL=30                                            # Minimum Quality for Bases
ANNOT_TAGS="DP,AD"                                      # Annotation-Tags for mpileup (comma separated)
VARIANT_QUAL=40                                         # Minimum variant quality
MAX_BASE_DEPTH=260000                                   # Max. base depth for variants (over all individuals)
MIN_BASE_DEPTH=10                                       # Min. base depth for variants (over all individuals)
MIN_ALT_DEPTH=0                                         # Min. depth of the first alternative allele in at least one sample
REGION=1                                                # Flag for Variant calling by region, makes calling faster (1 = yes, 0= no)

# Variant filtering
#------------------
MIN_DP=4
MAX_DP=100
MAF=0
MAX_MISSING=0
MIN_Q=40
RM_INDELS=1
RM_MULTIALLELIC=1
# Vcftools stats settings
QC_STATS_FLAG=0                                         # QC stats flag: 0 = no QC statistics will be calculated, 1 = QC statistics will be calculated
TAJIMAD_WINDOW=10000
PI_WINDOW=10000

########################################################


########################################################
# Optional technical settings — modify only if needed
########################################################

# Fastqc settings
#----------------
FASTQC_CORES=2                                          # Number of CPU threads for each fastqc job
FASTQC_MEM="10G"                                        # Memory allocated for fastqc jobs
FASTQC_TIME_LIMIT="00:20:00"                            # Time limit for each fastqc job (depends on cluster load)

# Multiqc settings
#-----------------
MULTIQC_MEM="60G"                                       # Memory allocated for multiqc
MULTIQC_TIME_LIMIT="01:00:00"                           # Time limit for the multiqc job
MULTIQC_CONFIG="${CODE_DIR}/custom_multiqc_config.yaml" # Multiqc config file, e.g. number formatting

# Trimming settings
#------------------
TRIMMING_CORES=4                                        # Number of CPU threads for Trim Galore
TRIMMING_MEM="10G"                                      # Memory allocated for trimming jobs
TRIMMING_TIME_LIMIT="00:30:00"                          # Time limit for each trimming job (depends on cluster load)

# Alignment settings
#-------------------
INDEX_SIZE=30                                           # Memory used by minimap2 when building the index in G (please give an integer for this one)
INDEX_TIME=30                                           # Time limit for creating the Index in Minutes (please give an integer for this one)
ALIGNMENT_CORES=16                                      # Number of CPU threads for minimap2
ALIGNMENT_MEM="40G"                                     # Memory allocated for alignment jobs
ALIGNMENT_TIME_LIMIT="00:30:00"                         # Time limit for each alignment job (depends on cluster load)
SAMTOOLS_CORES=4                                        # Number of CPU threads for each samtools job (4 is standard)
# Merging settings - ignore if not wanted
MERGING_CORES=8                                         # Number of CPU threads for merging
MERGING_MEM="40G"                                       # Memory allocated for merging
MERGING_TIME_LIMIT="03:00:00"                           # Time limit for merging (usually much less)
MERGING_BATCH_SIZE=100                                  # Number of BAMs to merge per batch

# Varaint calling settings
#-------------------------
MAX_READ_DEPTH=3000000                                  # limits number examined reads when the coverage is really high (performance reasons)
VARIANT_CALLING_CORES=16                                # Number of CPU threads used for bcftools mpileup
PRE_FILTERING_CORES=4                                   # Number of CPU threads used for bcftools call and filter
VARIANT_CALLING_MEM="200G"                              # Memory allocated for variant calling
VARIANT_CALLING_TIME_LIMIT="150:00:00"                  # Time limit for variant calling
VCF_NAME="${JOB_NAME}\
_q${VAR_MAP_QUAL}Q${BASE_QUAL}\
_QUAL${VARIANT_QUAL}\
_DP${MIN_BASE_DEPTH}-${MAX_BASE_DEPTH}\
_AD${MIN_ALT_DEPTH}"                                    # Präfix for the vcf file, can be changed if wanted
# Merging settings - acitive if REGION is 1
MERGE_VCF_CORES=4                                       # Number of CPU threads for concatting the vcf files
MERGE_VCF_MEM="100G"                                    # Memory allocated for concatting the vcfs
MERGE_VCF_TIME_LIMIT="1:00:00"                          # Time limit for concatting the vcfs
# bcftools stats settings
STATS_MEM="100G"                                        # Memory allocated for creating a stats file for the final vcf
STATS_TIME_LIMIT="00:30:00"                             # Time limit for creating a stats file for the final vcf

# Variant filtering settings
#---------------------------
VARIANT_FILTERING_MEM="100G"
VARIANT_FILTERING_TIME_LIMIT="48:00:00"
# vcftools stats settings
VARIANT_QC_MEM="100G"
VARIANT_QC_TIME_LIMIT="78:00:00"

# Other cluster settings
#-----------------------
ARRAY_CAP=25                                            # Maximum number of files processed simultaneously

########################################################


########################################################
# General script settings
########################################################

set -euo pipefail
shopt -s nullglob

if (( INITIAL_FASTQC == 1 || TRIMMING == 1 || ALIGNMENT == 1 )); then
    GLOB_PATTERN="${GLOB_PATTERN_FASTQ[*]}"
    FILES=( $GLOB_PATTERN )
    N=${#FILES[@]}
    if (( N == 0 )); then
        echo "No input FASTQ files found matching: $GLOB_PATTERN_FASTQ" >&2
        exit 1
    fi
    echo "${N} files found." >&2
    ARRAY_RANGE="0-$((N-1))%${ARRAY_CAP}"
elif (( VARIANT_CALLING == 1 )); then
    echo "len(GLOB_PATTERN_BAM[*])=${#GLOB_PATTERN_BAM[*]}" >&2
    GLOB_PATTERN="${GLOB_PATTERN_BAM[*]}"
    echo "AFTER  len(GLOB_PATTERN)=${#GLOB_PATTERN}" >&2
    FILES=( $GLOB_PATTERN )
    N=${#FILES[@]}
    if (( N == 0 )); then
        echo "No input BAM files found matching: $GLOB_PATTERN_BAM" >&2
        exit 1
    fi
    echo "${N} files found." >&2
    ARRAY_RANGE="0-$((N-1))%${ARRAY_CAP}"
fi

########################################################

########################################################
# This Block calls the FASTQC module
########################################################

if (( INITIAL_FASTQC == 1 )); then

    export OUT="${ALL_OUT}/${JOB_NAME}_fastqc" \
       GLOB_PATTERN="${GLOB_PATTERN_FASTQ[*]}" \
       JOB_DIR JOB_NAME MAIL_USER CODE_DIR N\
       FASTQC_CORES ARRAY_RANGE FASTQC_MEM FASTQC_TIME_LIMIT \
       MULTIQC_MEM MULTIQC_TIME_LIMIT MULTIQC_CONFIG

    FASTQC_JOBID=$(sbatch \
        --job-name="${JOB_NAME}_submit_fastqc" \
        --output="${JOB_DIR}/${JOB_NAME}_submit_fastqc_%j.out" \
        --error="${JOB_DIR}/${JOB_NAME}_submit_fastqc_%j.err" \
        --mail-user="$MAIL_USER" \
        --export=ALL \
        "$CODE_DIR"/fastqc_module.sbatch | awk '{print $NF}')
fi

########################################################

########################################################
# This Block calls the MULTIQC worker
########################################################

if (( INITIAL_MULTIQC == 1 )); then

    export OUT="${ALL_OUT}/${JOB_NAME}_multiqc" \
       NAME="${JOB_NAME}_multiqc" \
       CONFIG="${MULTIQC_CONFIG}" \
       GLOB_PATTERN=$GLOB_PATTERN_REPORTS \
       INPUT

    mkdir -p "$JOB_DIR" "$OUT"

    MULTIQC_JOBID=$(sbatch \
        --job-name="${JOB_NAME}_multiqc" \
        --ntasks=1 \
        --nodes=1 \
        --partition=cpu \
        --auks=yes \
        --cpus-per-task=1 \
        --mem="$MULTIQC_MEM" \
        --output="${JOB_DIR}/${JOB_NAME}_multiqc_%j.out" \
        --error="${JOB_DIR}/${JOB_NAME}_multiqc_%j.err" \
        --mail-type=ALL \
        --mail-user="$MAIL_USER" \
        --time="$MULTIQC_TIME_LIMIT" \
        --export=ALL \
        "$CODE_DIR"/workers/multiqc_worker.sbatch | awk '{print $NF}')
fi

########################################################

########################################################
# This Block calls the TRIMMING module
########################################################

if (( TRIMMING == 1 )); then

    TRIMMING_OUT="${ALL_OUT}/${JOB_NAME}_trimming"
    TRIMMING_JOBID_FILE="${JOB_DIR}/${JOB_NAME}_trimming_array.jobid"
    rm -f "$TRIMMING_JOBID_FILE"

    export OUT="$TRIMMING_OUT" \
       GLOB_PATTERN="${GLOB_PATTERN_FASTQ[*]}" \
       JOB_DIR JOB_NAME MAIL_USER CODE_DIR N \
       QUALITY MIN_LENGTH CLIP3 CLIP5 ADAPTER \
       TRIMMING_CORES ARRAY_RANGE TRIMMING_MEM TRIMMING_TIME_LIMIT \
       MULTIQC_MEM MULTIQC_TIME_LIMIT MULTIQC_CONFIG

    TRIMMING_SUBMIT_ID=$(sbatch \
        --job-name="${JOB_NAME}_submit_trimming" \
        --output="${JOB_DIR}/${JOB_NAME}_submit_trimming_%j.out" \
        --error="${JOB_DIR}/${JOB_NAME}_submit_trimming_%j.err" \
        --mail-user="$MAIL_USER" \
        --export=ALL \
        "$CODE_DIR"/trimming_module.sbatch | awk '{print $NF}')

fi

########################################################

########################################################
# This Block calls the ALIGNMENT module
########################################################

if (( ALIGNMENT == 1 )); then

    ALIGNMENT_JOBID_FILE="${JOB_DIR}/${JOB_NAME}_alignment_array.jobid"
    rm -f "$ALIGNMENT_JOBID_FILE"
    INDEX_T=$(printf "%02d:%02d:00" $(((INDEX_TIME + 10) / 60)) $(((INDEX_TIME + 10) % 60)))
    MEM=$(( INDEX_SIZE + 15 ))
    DEP_OPT=""
    export GLOB_PATTERN="${GLOB_PATTERN_FASTQ[*]}"

    if (( TRIMMING == 1 )); then
        export TRIMMING_OUT="${ALL_OUT}/${JOB_NAME}_trimming"
        DEP_OPT="--dependency=afterok:${TRIMMING_SUBMIT_ID}"
    fi

    export OUT="${ALL_OUT}/${JOB_NAME}_alignment" \
        JOB_DIR JOB_NAME MAIL_USER CODE_DIR N \
        INDEX_SIZE="${INDEX_SIZE}G" INDEX REFERENCE_GENOME MAP_QUAL MERGING_FLAG \
        MERGING_CORES MERGING_MEM MERGING_TIME_LIMIT MERGING_BATCH_SIZE \
        ALIGNMENT_CORES ARRAY_RANGE ALIGNMENT_MEM ALIGNMENT_TIME_LIMIT \
        SAMTOOLS_CORES MULTIQC_MEM MULTIQC_TIME_LIMIT MULTIQC_CONFIG

    ALIGNMENT_SUBMIT_JOBID=$(sbatch \
       --job-name="${JOB_NAME}_submit_alignment" \
       $DEP_OPT \
       --time="$INDEX_T" \
       --mem="${MEM}G"\
       --output="${JOB_DIR}/${JOB_NAME}_submit_alignment_%j.out" \
       --error="${JOB_DIR}/${JOB_NAME}_submit_alignment_%j.err" \
       --mail-user="$MAIL_USER" \
       --export=ALL \
       "$CODE_DIR"/alignment_module.sbatch | awk '{print $NF}')

fi

########################################################

########################################################
# This Block calls the VARIANT_CALLING module
########################################################

if (( VARIANT_CALLING == 1 )); then

    DEP_OPT=""
    export GLOB_PATTERN="${GLOB_PATTERN_BAM[*]}"

    if (( ALIGNMENT == 1 )); then
        export ALIGNMENT_OUT="${ALL_OUT}/${JOB_NAME}_alignment"
        DEP_OPT="--dependency=afterok:${ALIGNMENT_SUBMIT_JOBID}"
    fi

    export OUT="${ALL_OUT}/${JOB_NAME}_variant_calling" \
        GLOB_PATTERN JOB_DIR JOB_NAME MAIL_USER CODE_DIR VCF_NAME \
        REFERENCE_GENOME VAR_MAP_QUAL MAX_READ_DEPTH BASE_QUAL ANNOT_TAGS \
        VARIANT_QUAL MIN_BASE_DEPTH MAX_BASE_DEPTH MIN_ALT_DEPTH \
        VARIANT_CALLING_CORES PRE_FILTERING_CORES VARIANT_CALLING_MEM VARIANT_CALLING_TIME_LIMIT \
        REGION MULTIQC_MEM MULTIQC_TIME_LIMIT MULTIQC_CONFIG \
        STATS_MEM STATS_TIME_LIMIT MERGE_VCF_CORES MERGE_VCF_MEM MERGE_VCF_TIME_LIMIT

    VARIANT_CALLING_SUBMIT_JOBID=$(sbatch \
       --job-name="${JOB_NAME}_submit_variant_calling" \
       $DEP_OPT \
       --output="${JOB_DIR}/${JOB_NAME}_submit_variant_calling_%j.out" \
       --error="${JOB_DIR}/${JOB_NAME}_submit_variant_calling_%j.err" \
       --mail-user="$MAIL_USER" \
       --export=ALL \
       "$CODE_DIR"/variant_calling_module.sbatch | awk '{print $NF}')

fi

########################################################

########################################################
# This Block calls the VARIANT_FILTERING module
########################################################

if (( VARIANT_FILTERING == 1 )); then

    export OUT="${ALL_OUT}/${JOB_NAME}_vcftools_stats" \
        JOB_DIR JOB_NAME CODE_DIR MAIL_USER \
        MIN_DP MAX_DP MAF MAX_MISSING MIN_Q RM_INDELS RM_MULTIALLELIC \
        TAJIMAD_WINDOW PI_WINDOW VCF_IN\
        VARIANT_FILTERING_MEM VARIANT_FILTERING_TIME_LIMIT \
        VARIANT_QC_MEM VARIANT_QC_TIME_LIMIT QC_STATS_FLAG

    VARIANT_FILTERING_SUBMIT_JOBID=$(sbatch \
       --job-name="${JOB_NAME}_submit_variant_filtering" \
       --output="${JOB_DIR}/${JOB_NAME}_submit_variant_filtering_%j.out" \
       --error="${JOB_DIR}/${JOB_NAME}_submit_variant_filtering_%j.err" \
       --mail-user="$MAIL_USER" \
       --export=ALL \
       "$CODE_DIR"/variant_filtering_module.sbatch | awk '{print $NF}')

fi

########################################################

