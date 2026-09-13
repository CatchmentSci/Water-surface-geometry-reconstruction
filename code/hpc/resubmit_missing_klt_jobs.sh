#!/bin/bash

set -euo pipefail

WORKER_SCRIPT="run_one_klt_case.sh"
USER_NAME="${USER}"
LOGDIR="/mnt/nfs/home/nmp65/klt_slurm_logs"
RUNTIMEDIR="/mnt/nfs/home/nmp65/klt_runtime_jobs"

mkdir -p "$LOGDIR" "$RUNTIMEDIR"

if [[ ! -f "$WORKER_SCRIPT" ]]; then
    echo "ERROR: Cannot find $WORKER_SCRIPT in the current directory."
    exit 1
fi

echo "Checking KLT jobs for user: $USER_NAME"
echo "Worker script: $WORKER_SCRIPT"
echo "Slurm logs: $LOGDIR"
echo

# Excel-visible CSV rows to check/resubmit. Row 1 is the header.
S_IDS=(4 7 11 14 {21..29})

for S_ID in "${S_IDS[@]}"
do
    JOB_NAME="klt_row${S_ID}"

    ACTIVE_JOB_IDS=$(squeue \
        -h \
        -u "$USER_NAME" \
        --name="$JOB_NAME" \
        --states=RUNNING,PENDING,CONFIGURING,COMPLETING \
        -o "%i")

    if [[ -n "$ACTIVE_JOB_IDS" ]]; then
        echo "Excel row ${S_ID}: already active as job(s): $ACTIVE_JOB_IDS"
    else
        echo "Excel row ${S_ID}: not active; submitting new job..."

        sbatch \
            --job-name="$JOB_NAME" \
            --export=ALL,KLT_S_ID="$S_ID" \
            "$WORKER_SCRIPT"
    fi
done

echo
echo "Done. Current matching jobs:"
squeue -u "$USER_NAME" -o "%.18i %.12j %.10T %.10M %.20R"
