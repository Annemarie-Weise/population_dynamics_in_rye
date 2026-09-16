#!/bin/bash

# Worker script for generating a MultiQC report.
# Must be launched by the main SLURM submission script, which provides
# all required variables and resource settings.
#
# Required variables:
#   INPUT        – Directory containing files or subfolders to scan
#   OUT          – Output directory for the MultiQC report
#   GLOB_PATTERN – Glob pattern for selecting subfolders/files under INPUT or ""
#   NAME         – Base name for the MultiQC report
#   CONFIG       - Custom config file for multiqc
#
# Function:
#   Selects matching INPUT subfolders/files (or INPUT itself if the pattern
#   is empty or yields no matches), then runs MultiQC on the result.
#
# Output:
#   MultiQC report written to $OUT
#
# Note:
#   Requires SLURM resources (cpus-per-task, mem, time) set by the main script.

set -euo pipefail
shopt -s nullglob
# init module environment
. /etc/profile.d/modules.sh

module load multiqc

if [[ -n "$GLOB_PATTERN" ]]; then
    FILES=( "$INPUT"/$GLOB_PATTERN )
    if (( ${#FILES[@]} == 0 )); then
        FILES=( "$INPUT" )
    fi
else
    FILES=( "$INPUT" )
fi

multiqc -v -f -o "$OUT" --filename "$NAME" "${FILES[@]}" -c "$CONFIG"
