#!/bin/bash

set -euo pipefail

LOGDIR="/mnt/nfs/home/nmp65/klt_slurm_logs"
RUNTIMEDIR="/mnt/nfs/home/nmp65/klt_runtime_jobs"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKER_SCRIPT="${SCRIPT_DIR}/run_one_klt_analysis_case.sh"
LOOKUP_FILE="${SCRIPT_DIR}/klt_analysis_case_lookup.tsv"

# -----------------------------------------------------------------------------
# Edit this list to choose which real/video analysis cases to run.
# Each item may be either:
#   1) filenameIn, e.g. devon_dart20181129_13490
#   2) videoNumber / solver value, e.g. 109
# The filenameIn <-> videoNumber mapping is in klt_analysis_case_lookup.tsv.
#
# The Slurm job name is always made from the resolved videoNumber, not the
# filename. For example:
#   devon_dart20181129_13490 -> klt_v109
#   109                     -> klt_v109
# -----------------------------------------------------------------------------
CASES=(
    #15
    #30
    #51
    #61
    #123
    #127
    #141
    #152
    #175
    #182
    #    109
    #    70
    #    81
    #    90
    #    96
    #	20
    #	40
    #	117
    #   123
    #    141
    152
    #    175
    # another_filename_here
)

mkdir -p "$LOGDIR" "$RUNTIMEDIR"

resolve_case_ref() {
    local case_ref="$1"
    local lookup_file="${2:-$LOOKUP_FILE}"
    local matches

    if [[ ! -f "$lookup_file" ]]; then
        echo "ERROR: Cannot find lookup file: $lookup_file" >&2
        return 1
    fi

    matches=$(awk -v ref="$case_ref" '
        BEGIN { FS = "[\t, ]+" }
        /^[[:space:]]*($|#)/ { next }
        {
            filename_in = $1
            video_number = $2

            filename_lower = tolower(filename_in)
            video_lower = tolower(video_number)
            gsub(/_/, "", video_lower)

            if ((filename_lower == "filename" || filename_lower == "filename_in" || filename_lower == "filenamein") && \
                (video_lower == "video" || video_lower == "videonumber" || video_lower == "sweep" || video_lower == "sweepid" || video_lower == "sweepvalue" || video_lower == "value")) {
                next
            }

            if (filename_in == ref || video_number == ref) {
                print filename_in "\t" video_number
            }
        }
    ' "$lookup_file")

    local n_matches
    n_matches=$(printf '%s\n' "$matches" | sed '/^$/d' | wc -l | tr -d ' ')

    if [[ "$n_matches" -eq 0 ]]; then
        echo "ERROR: Case reference '$case_ref' was not found in: $lookup_file" >&2
        echo "       Use either a filenameIn value or a videoNumber from the lookup table." >&2
        return 1
    fi

    if [[ "$n_matches" -gt 1 ]]; then
        echo "ERROR: Case reference '$case_ref' matched more than one row in: $lookup_file" >&2
        printf '%s\n' "$matches" >&2
        echo "       Make filenameIn and videoNumber values unique." >&2
        return 1
    fi

    printf '%s\n' "$matches"
}

make_klt_job_name() {
    local video_number="$1"
    local safe_video
    safe_video=$(printf '%s' "$video_number" | tr -c 'A-Za-z0-9_' '_' | cut -c1-70)
    printf 'klt_v%s' "$safe_video"
}

submit_klt_cases() {
    local dry_run="0"

    if [[ "${1:-}" == "--dry-run" ]]; then
        dry_run="1"
    fi

    if [[ ! -f "$WORKER_SCRIPT" ]]; then
        echo "ERROR: Cannot find worker script: $WORKER_SCRIPT"
        exit 1
    fi

    if [[ ! -f "$LOOKUP_FILE" ]]; then
        echo "ERROR: Cannot find lookup file: $LOOKUP_FILE"
        exit 1
    fi

    echo "Submit script dir: $SCRIPT_DIR"
    echo "Worker script:     $WORKER_SCRIPT"
    echo "Lookup file:       $LOOKUP_FILE"
    echo "Slurm logs:        $LOGDIR"
    echo "Runtime dir:       $RUNTIMEDIR"
    echo

    for CASE_REF in "${CASES[@]}"
    do
        CASE_ROW=$(resolve_case_ref "$CASE_REF" "$LOOKUP_FILE")
        IFS=$'\t' read -r FILENAME_IN VIDEO_NUMBER <<< "$CASE_ROW"

        if [[ -z "$FILENAME_IN" || -z "$VIDEO_NUMBER" ]]; then
            echo "ERROR: Lookup row must contain filenameIn and videoNumber for case: $CASE_REF"
            exit 1
        fi

        JOB_NAME=$(make_klt_job_name "$VIDEO_NUMBER")

        echo "Submitting job for case: $CASE_REF"
        echo "filenameIn:  $FILENAME_IN"
        echo "videoNumber: $VIDEO_NUMBER"
        echo "Job name:    $JOB_NAME"

        if [[ "$dry_run" == "1" ]]; then
            echo "DRY RUN: sbatch --job-name=$JOB_NAME --chdir=/mnt/nfs/home/nmp65 --output=${LOGDIR}/%x-%j.out --error=${LOGDIR}/%x-%j.err --export=ALL,KLT_CASE_REF=$CASE_REF,KLT_LOOKUP_FILE=$LOOKUP_FILE $WORKER_SCRIPT"
        else
            sbatch \
                --job-name="$JOB_NAME" \
                --chdir=/mnt/nfs/home/nmp65 \
                --output="${LOGDIR}/%x-%j.out" \
                --error="${LOGDIR}/%x-%j.err" \
                --export=ALL,KLT_CASE_REF="$CASE_REF",KLT_LOOKUP_FILE="$LOOKUP_FILE" \
                "$WORKER_SCRIPT"
        fi

        echo
    done

    if [[ "$dry_run" == "1" ]]; then
        echo "Dry run complete. No jobs were submitted."
    else
        echo "Submitted ${#CASES[@]} KLT analysis job(s)."
        echo "Slurm logs will be written to: $LOGDIR"
        echo "Runtime files will be written to: $RUNTIMEDIR"
    fi
}

# This allows resubmit_missing_klt_analysis_jobs.sh to source this file and reuse
# the CASES list and helper functions without submitting jobs immediately.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    submit_klt_cases "${1:-}"
fi
