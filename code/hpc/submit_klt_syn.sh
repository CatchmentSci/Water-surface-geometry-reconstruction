#!/bin/bash

set -euo pipefail

LOGDIR="/mnt/nfs/home/nmp65/klt_slurm_logs"
RUNTIMEDIR="/mnt/nfs/home/nmp65/klt_runtime_jobs"
WORKER_SCRIPT="run_one_klt_case.sh"

mkdir -p "$LOGDIR" "$RUNTIMEDIR"

if [[ ! -f "$WORKER_SCRIPT" ]]; then
    echo "ERROR: Cannot find $WORKER_SCRIPT in the current directory."
    exit 1
fi

# Excel-visible CSV rows to run. Row 1 is the header; the worker converts to MATLAB data-row index.
S_IDS=(4 7 11 14 {16..24})

for S_ID in "${S_IDS[@]}"
do
    echo "Submitting job for Excel row ${S_ID}"

    sbatch \
        --job-name="klt_row${S_ID}" \
        --export=ALL,KLT_S_ID="${S_ID}" \
        "$WORKER_SCRIPT"
done

LAST_INDEX=$(( ${#S_IDS[@]} - 1 ))
echo "Submitted jobs for Excel rows: ${S_IDS[*]}."
echo "Slurm logs will be written to: $LOGDIR"
