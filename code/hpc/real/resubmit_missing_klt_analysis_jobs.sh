#!/bin/bash

set -euo pipefail

USER_NAME="${USER}"
LOGDIR="/mnt/nfs/home/nmp65/klt_slurm_logs"
RUNTIMEDIR="/mnt/nfs/home/nmp65/klt_runtime_jobs"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_one_klt_analysis_case.sh"
SUBMIT_SCRIPT="${SCRIPT_DIR}/submit_klt_analysis_batch.sh"
LOOKUP_FILE="${SCRIPT_DIR}/klt_analysis_case_lookup.tsv"

mkdir -p "$LOGDIR" "$RUNTIMEDIR"

if [[ ! -f "$WORKER_SCRIPT" ]]; then
    echo "ERROR: Cannot find worker script: $WORKER_SCRIPT"
    exit 1
fi

if [[ ! -f "$SUBMIT_SCRIPT" ]]; then
    echo "ERROR: Cannot find submit script: $SUBMIT_SCRIPT"
    exit 1
fi

if [[ ! -f "$LOOKUP_FILE" ]]; then
    echo "ERROR: Cannot find lookup file: $LOOKUP_FILE"
    exit 1
fi

# Reuse CASES, resolve_case_ref, and make_klt_job_name from submit script.
# shellcheck source=/dev/null
source "$SUBMIT_SCRIPT"

# The source command above sets WORKER_SCRIPT and LOOKUP_FILE from the submit
# script directory as absolute paths. Keep using those values here.

echo "Checking KLT analysis jobs for user: $USER_NAME"
echo "Worker script: $WORKER_SCRIPT"
echo "Submit script: $SUBMIT_SCRIPT"
echo "Lookup file: $LOOKUP_FILE"
echo "Slurm logs: $LOGDIR"
echo "Runtime dir: $RUNTIMEDIR"
echo "Job naming: klt_v<videoNumber>"
echo

for CASE_REF in "${CASES[@]}"
do
    CASE_ROW=$(resolve_case_ref "$CASE_REF" "$LOOKUP_FILE")
    IFS=$'\t' read -r FILENAME_IN VIDEO_NUMBER <<< "$CASE_ROW"
    JOB_NAME=$(make_klt_job_name "$VIDEO_NUMBER")

    ACTIVE_JOB_IDS=$(squeue \
        -h \
        -u "$USER_NAME" \
        --name="$JOB_NAME" \
        --states=RUNNING,PENDING,CONFIGURING,COMPLETING \
        -o "%i")

    if [[ -n "$ACTIVE_JOB_IDS" ]]; then
        echo "$CASE_REF / videoNumber $VIDEO_NUMBER: already active as job(s): $ACTIVE_JOB_IDS"
    else
        echo "$CASE_REF / videoNumber $VIDEO_NUMBER: not active; submitting new job as $JOB_NAME..."

        sbatch \
            --job-name="$JOB_NAME" \
            --chdir=/mnt/nfs/home/nmp65 \
            --output="${LOGDIR}/%x-%j.out" \
            --error="${LOGDIR}/%x-%j.err" \
            --export=ALL,KLT_CASE_REF="$CASE_REF",KLT_LOOKUP_FILE="$LOOKUP_FILE" \
            "$WORKER_SCRIPT"
    fi
done

echo
echo "Done. Current matching jobs:"
squeue -u "$USER_NAME" -o "%.18i %.24j %.10T %.10M %.20R"
