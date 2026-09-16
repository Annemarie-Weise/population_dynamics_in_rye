#!/bin/bash

# Worker script for batch-wise BAM merging.
# This script is NOT meant to be run directly. It must be launched by the main
# SLURM submission script, which provides all required variables and resources.
#
# Required variables:
#   INPUT              – Directory containing sorted BAM files
#   NAME_MERGED_BAM    – Name of the final merged BAM
#   MERGING_CORES      – Number of CPU threads for samtools
#   MERGING_BATCH_SIZE – Number of BAMs merged per batch
#   MERGED             - Directory to save the merged BAM file
#
# Function:
#   Merges all BAM files in $INPUT in batches to reduce memory usage.
#   Temporary batch files are generated and merged iteratively until only one
#   final BAM remains. A CSI index is created at the end.
#
# Output:
#   ${INPUT}/merged_bam/${NAME_MERGED_BAM}.bam (+ .csi index)
#
# Note:
#   Requires SLURM resources (cpus-per-task, mem, time) set by the main script.

set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load samtools

TMP_DIR="${INPUT}/tmp_merge"
mkdir -p "$TMP_DIR"

BAM_FILES=( "$INPUT"/*.bam )

CURRENT_BATCH=""
BATCH_NUM=1
TOTAL=${#BAM_FILES[@]}
i=0

# batch logic:
while [ "$i" -lt "$TOTAL" ]; do

    CHUNK=( "${BAM_FILES[@]:$i:$MERGING_BATCH_SIZE}" )
    echo "Batch ${BATCH_NUM}: i=$i, nehme ${#CHUNK[@]} Dateien"

    NEW_BATCH="${TMP_DIR}/batch_${BATCH_NUM}.bam"

    if [ -z "$CURRENT_BATCH" ]; then
        samtools merge -@ "$MERGING_CORES" -f -o "$NEW_BATCH" "${CHUNK[@]}"
    else
        samtools merge -@ "$MERGING_CORES" -f -o "$NEW_BATCH" "$CURRENT_BATCH" "${CHUNK[@]}"
        rm -f "$CURRENT_BATCH"
    fi

    CURRENT_BATCH="$NEW_BATCH"
    i=$(( i + MERGING_BATCH_SIZE ))
    BATCH_NUM=$(( BATCH_NUM + 1 ))

done

FINAL_OUT="${MERGED}/${NAME_MERGED_BAM}.bam"
mv "$CURRENT_BATCH" "$FINAL_OUT"
rm -rf "$TMP_DIR"


# indexing the final bam
samtools index -@ "$MERGING_CORES" -c "$FINAL_OUT"
