#!/bin/bash

# Worker script for per-file read alignment.
# Not meant to be run directly — it must be launched by the main SLURM
# submission script, which provides all required variables and resources.
#
# Required variables:
#   GLOB_PATTERN     – Glob expression selecting FASTQ files
#   OUT              – Output directory
#   INDEX            – Path to the minimap2 index (.mmi)
#   ALIGNMENT_CORES  – Number of threads for minimap2
#   MAP_QUAL         – Minimum MAPQ for kept alignments
#   SAMTOOLS_CORES   - Number of threads for the samtools operations in the pipeline
#
# Function:
#   Processes one FASTQ file per array task:
#       - aligns reads with minimap2
#       - computes raw alignment statistics on the unfiltered BAM stream
#       - filters out unmapped, secondary, supplementary and low-quality
#         alignments (MAPQ < MAP_QUAL)
#       - sorts the remaining reads into a final BAM file and creates an index
#       - calculates coverage statistics with samtools
#
# Output:
#   $OUT/<sample>.sorted.bam
#   $OUT/<sample>.sorted.bam.csi
#   $OUT/<sample>.raw.stats
#   $OUT/<sample>.samtools.coverage.txt
#
# Note:
#   SLURM array settings (cpus-per-task, mem, time, array range)
#   must be configured in the main submission script.


set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load minimap2
module load samtools
#module load mosdepth
#module load bbmap


FILES=( $GLOB_PATTERN )

idx=${SLURM_ARRAY_TASK_ID:?}
(( idx < ${#FILES[@]} )) || { echo "Index $idx out of ${#FILES[@]}"; exit 1; }

FILE="${FILES[$idx]}"
BASENAME=$(basename "$FILE")
BASENAME="${BASENAME%.fastq.gz}"
BASENAME="${BASENAME%.fq.gz}"
BASENAME="${BASENAME%.fastq}"
BASENAME="${BASENAME%.fq}"

minimap2 -a -x sr --secondary=yes -t "$ALIGNMENT_CORES" -2 -R "@RG\tID:${BASENAME}\tSM:${BASENAME}" "$INDEX" "$FILE" \
  | samtools view -@ "$SAMTOOLS_CORES" -b \
  | tee >(samtools stats -@ "$SAMTOOLS_CORES" - > "${OUT}/${BASENAME}.stats") \
  | samtools view -@ "$SAMTOOLS_CORES" -b -F 0x904 -q "$MAP_QUAL" \
  | samtools sort -@ "$SAMTOOLS_CORES" -o "${OUT}/${BASENAME}.sorted.bam" -

samtools index -c -@ $(( SAMTOOLS_CORES * 2 )) "${OUT}/${BASENAME}.sorted.bam"
samtools coverage --output "${OUT}/${BASENAME}.samtools.coverage.txt" "${OUT}/${BASENAME}.sorted.bam"
#mosdepth -t $(( SAMTOOLS_CORES * 4 )) "${OUT}/${BASENAME}" "${OUT}/${BASENAME}.sorted.bam"
#pileup.sh in="${OUT}/${BASENAME}.sorted.bam" overwrite=t covstats="${OUT}/${BASENAME}.covstats" threads="$(( SAMTOOLS_CORES * 4 ))"
