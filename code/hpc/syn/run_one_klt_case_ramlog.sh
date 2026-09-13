#!/bin/bash

#SBATCH --account=comet_wave
#SBATCH --partition=default_paid
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matthew.perks@newcastle.ac.uk
#SBATCH --job-name=klt_case
#SBATCH --chdir=/mnt/nfs/home/nmp65
#SBATCH --output=/mnt/nfs/home/nmp65/klt_slurm_logs/%x-%j.out
#SBATCH --error=/mnt/nfs/home/nmp65/klt_slurm_logs/%x-%j.err

# Run one synthetic KLT case under Slurm using MATLAB/2024a.
#
# Key choices based on testing on Comet:
#   - Explicitly load MATLAB/2024a.
#   - Redirect MATLAB startup/temp/cache/prefs paths away from /tmp and /tmp/slurmd.
#   - Launch MATLAB through stdin rather than -batch/-r inline commands.
#   - Keep Slurm logs in /mnt/nfs/home/nmp65/klt_slurm_logs.
#   - Do not set parcluster.JobStorageLocation in bash.
#   - Do not override HPC MATLAB module/site configuration beyond startup-safe paths.
#
# This script expects KLT_S_ID to be supplied by sbatch --export, e.g.
#   sbatch --job-name=klt_s1 --export=ALL,KLT_S_ID=1 run_one_klt_case_ramlog.sh
#
# RAM logging:
#   - Creates a CSV time series of MATLAB + descendant-process RSS memory.
#   - Prints the peak observed RSS at the end.
#   - Prints Slurm accounting memory information where available.

set -euo pipefail

LOGDIR="/mnt/nfs/home/nmp65/klt_slurm_logs"
PROJECT_ROOT="/mnt/nfs/home/nmp65/Downloads/v1_current"
SYN_DIR="${PROJECT_ROOT}/Syn"
INPUTS_DIR="${PROJECT_ROOT}/Videos/Inputs"

JOBTAG="${SLURM_JOB_ID:-manual}"
JOBNAME="${SLURM_JOB_NAME:-klt_case}"
SAFE_BASE="/mnt/nfs/home/nmp65/klt_runtime_jobs/${JOBNAME}_${JOBTAG}"

WORKDIR="${SAFE_BASE}/work"
TMPROOT="${SAFE_BASE}/tmp"
PREFDIR="${SAFE_BASE}/prefs"
MCRROOT="${SAFE_BASE}/mcr_cache"
XDGCACHE="${SAFE_BASE}/xdg_cache"
XDGCONFIG="${SAFE_BASE}/xdg_config"
XDGRUNTIME="${SAFE_BASE}/xdg_runtime"

mkdir -p "$LOGDIR" "$WORKDIR" "$TMPROOT" "$PREFDIR" "$MCRROOT" "$XDGCACHE" "$XDGCONFIG" "$XDGRUNTIME"
chmod 700 "$SAFE_BASE" "$WORKDIR" "$TMPROOT" "$PREFDIR" "$MCRROOT" "$XDGCACHE" "$XDGCONFIG" "$XDGRUNTIME"

# MATLAB/2024a failed under Slurm when it used the default /tmp paths.
# These are startup/cache/temp redirects only; they do not configure parpool
# or MATLAB parallel cluster profiles.
export TMPDIR="$TMPROOT"
export TMP="$TMPROOT"
export TEMP="$TMPROOT"
export MATLAB_PREFDIR="$PREFDIR"
export MCR_CACHE_ROOT="$MCRROOT"
export XDG_CACHE_HOME="$XDGCACHE"
export XDG_CONFIG_HOME="$XDGCONFIG"
export XDG_RUNTIME_DIR="$XDGRUNTIME"

cd "$WORKDIR"

# -------------------------------------------------------------------------
# RAM logging helper functions
# -------------------------------------------------------------------------

collect_descendant_pids() {
    # Print the root PID and all descendant PIDs recursively.
    # This is used to capture MATLAB plus worker/child processes.
    local root_pid="$1"

    if [[ -z "$root_pid" ]]; then
        return 0
    fi

    local queue=("$root_pid")
    local all_pids=("$root_pid")
    local parent_pid
    local child_pid

    while [[ "${#queue[@]}" -gt 0 ]]; do
        parent_pid="${queue[0]}"
        queue=("${queue[@]:1}")

        while read -r child_pid; do
            if [[ -n "$child_pid" ]]; then
                all_pids+=("$child_pid")
                queue+=("$child_pid")
            fi
        done < <(pgrep -P "$parent_pid" 2>/dev/null || true)
    done

    printf "%s\n" "${all_pids[@]}" | awk 'NF' | sort -n -u
}

start_memory_logger() {
    # Logs RSS memory for MATLAB and all descendant processes.
    #
    # Arguments:
    #   1: root MATLAB launcher PID
    #   2: output CSV file
    #   3: sampling interval in seconds
    local root_pid="$1"
    local memlog="$2"
    local interval_sec="${3:-60}"

    echo "timestamp,total_rss_gb,process_count,pids" > "$memlog"

    while true; do
        # Stop once the root process has gone and no descendant PIDs remain.
        if ! kill -0 "$root_pid" 2>/dev/null; then
            break
        fi

        local timestamp
        timestamp="$(date '+%Y-%m-%d %H:%M:%S')"

        local pids
        pids="$(collect_descendant_pids "$root_pid" | paste -sd, -)"

        local total_rss_kb="0"
        local process_count="0"

        if [[ -n "$pids" ]]; then
            total_rss_kb="$(
                ps -o rss= -p "$pids" 2>/dev/null \
                | awk '{sum += $1} END {print sum + 0}'
            )"

            process_count="$(
                ps -o pid= -p "$pids" 2>/dev/null \
                | wc -l \
                | tr -d ' '
            )"
        fi

        local total_rss_gb
        total_rss_gb="$(
            awk -v kb="$total_rss_kb" 'BEGIN {printf "%.3f", kb / 1024 / 1024}'
        )"

        echo "${timestamp},${total_rss_gb},${process_count},${pids}" >> "$memlog"

        sleep "$interval_sec"
    done

    # Take one final sample if possible.
    local timestamp
    timestamp="$(date '+%Y-%m-%d %H:%M:%S')"

    local pids
    pids="$(collect_descendant_pids "$root_pid" | paste -sd, - || true)"

    local total_rss_kb="0"
    local process_count="0"

    if [[ -n "$pids" ]]; then
        total_rss_kb="$(
            ps -o rss= -p "$pids" 2>/dev/null \
            | awk '{sum += $1} END {print sum + 0}'
        )"

        process_count="$(
            ps -o pid= -p "$pids" 2>/dev/null \
            | wc -l \
            | tr -d ' '
        )"
    fi

    local total_rss_gb
    total_rss_gb="$(
        awk -v kb="$total_rss_kb" 'BEGIN {printf "%.3f", kb / 1024 / 1024}'
    )"

    echo "${timestamp},${total_rss_gb},${process_count},${pids}" >> "$memlog"
}

print_memory_log_summary() {
    local memlog="$1"

    echo
    echo "============================================================"
    echo "Process-based RAM log summary"
    echo "============================================================"

    if [[ -f "$memlog" ]]; then
        echo "Memory CSV log:"
        echo "  $memlog"
        echo

        echo "Peak observed MATLAB + child-process RSS:"
        awk -F, '
            NR > 1 && $2 + 0 > max_rss {
                max_rss = $2 + 0
                max_time = $1
                max_proc = $3
            }
            END {
                if (NR > 1) {
                    printf "  %.3f GB at %s with %s processes\n", max_rss, max_time, max_proc
                } else {
                    print "  No samples recorded."
                }
            }
        ' "$memlog"

        echo
        echo "Last 10 RAM samples:"
        tail -n 10 "$memlog"
    else
        echo "No memory log found:"
        echo "  $memlog"
    fi

    echo "============================================================"
}

print_slurm_memory_summary() {
    echo
    echo "============================================================"
    echo "Slurm memory accounting summary"
    echo "============================================================"

    if command -v sacct >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "sacct memory report for job ${SLURM_JOB_ID}:"
        echo

        sacct -j "$SLURM_JOB_ID" \
            --units=G \
            --format=JobID,JobName%25,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,AveRSS,MaxVMSize \
            || echo "WARNING: sacct command failed or accounting is not available yet."

        echo
        echo "Note: sacct MaxRSS may be most reliable after the whole Slurm job has fully ended."
        echo "You can re-run this manually after completion:"
        echo "  sacct -j ${SLURM_JOB_ID} --units=G --format=JobID,JobName%25,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,AveRSS,MaxVMSize"
    else
        echo "sacct unavailable or not running inside Slurm."
    fi

    echo "============================================================"
}

cleanup_memory_logger() {
    if [[ -n "${MEMLOGGER_PID:-}" ]]; then
        kill "$MEMLOGGER_PID" 2>/dev/null || true
        wait "$MEMLOGGER_PID" 2>/dev/null || true
    fi
}

trap cleanup_memory_logger EXIT

echo "============================================================"
echo "KLT MATLAB/2024a Slurm job"
echo "============================================================"
echo "Started:          $(date)"
echo "Host:             $(hostname)"
echo "User:             $(whoami)"
echo "SLURM job ID:     ${SLURM_JOB_ID:-not_in_slurm}"
echo "SLURM job name:   ${SLURM_JOB_NAME:-not_in_slurm}"
echo "Submit dir:       ${SLURM_SUBMIT_DIR:-not_in_slurm}"
echo "Work dir:         $WORKDIR"
echo "Log dir:          $LOGDIR"
echo "TMPDIR:           $TMPDIR"
echo "MATLAB_PREFDIR:   $MATLAB_PREFDIR"
echo "MCR_CACHE_ROOT:   $MCR_CACHE_ROOT"
echo "XDG_RUNTIME_DIR:  $XDG_RUNTIME_DIR"
echo "CPUs allocated:   ${SLURM_CPUS_PER_TASK:-not_set}"
echo "RAM requested:    ${SLURM_MEM_PER_NODE:-not_set_MB_per_node} MB per node"
echo "============================================================"
echo

# Initialise the module command if needed.
if ! command -v module >/dev/null 2>&1; then
    [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
    [[ -f /usr/share/Modules/init/bash ]] && source /usr/share/Modules/init/bash
fi

if ! command -v module >/dev/null 2>&1; then
    echo "ERROR: module command is unavailable."
    exit 2
fi

echo "Loading MATLAB/2024a..."
module unload MATLAB 2>/dev/null || true
module unload MATLAB/2026a 2>/dev/null || true
module unload MATLAB/2024a 2>/dev/null || true
module load MATLAB/2024a

echo
echo "Loaded modules:"
module list 2>&1
echo

MATLAB_EXE="$(which matlab)"
echo "MATLAB executable: $MATLAB_EXE"

if [[ "$MATLAB_EXE" != *"2024a"* ]]; then
    echo "ERROR: matlab executable does not appear to be MATLAB/2024a."
    echo "Resolved executable was:"
    echo "  $MATLAB_EXE"
    exit 3
fi

# Get Excel-visible CSV row number. Row 1 is assumed to be the header.
# MATLAB readtable uses data-row indexing, so testConfig.caseCsvRow = Excel row - 1.
S_ID="${KLT_S_ID:?KLT_S_ID is not set}"
if ! [[ "$S_ID" =~ ^[0-9]+$ ]]; then
    echo "ERROR: KLT_S_ID must be a positive integer Excel row number. Got: $S_ID"
    exit 1
fi

CASE_CSV_ROW=$((S_ID - 1))
if [[ "$CASE_CSV_ROW" -lt 1 ]]; then
    echo "ERROR: KLT_S_ID=$S_ID is not a valid Excel data row. Row 1 is the header."
    exit 1
fi

echo "Running synthetic Excel row ${S_ID}"
echo "Using MATLAB testConfig.caseCsvRow=${CASE_CSV_ROW} (Excel row ${S_ID} minus header row)"

# Select the correct input MAT file based on Excel-visible CSV row number.
if [[ "$S_ID" -ge 2 && "$S_ID" -le 10 ]]; then
    INPUT_MAT_FILE="${SYN_DIR}/solver_inputs_0pt88.mat"
elif [[ "$S_ID" -ge 11 && "$S_ID" -le 19 ]]; then
    INPUT_MAT_FILE="${SYN_DIR}/solver_inputs_1pt50.mat"
elif [[ "$S_ID" -ge 20 && "$S_ID" -le 31 ]]; then
    INPUT_MAT_FILE="${SYN_DIR}/solver_inputs_1pt70.mat"
elif [[ "$S_ID" -ge 32 && "$S_ID" -le 40 ]]; then
    INPUT_MAT_FILE="${SYN_DIR}/solver_inputs_1pt10.mat"
else
    echo "ERROR: Invalid Excel row S_ID: $S_ID"
    exit 1
fi

KLT_FUNCTION_NAME="KLT_synthetic_truth_test_v3"

echo "MATLAB function: $KLT_FUNCTION_NAME"
echo "Input MAT file:  $INPUT_MAT_FILE"

if [[ ! -f "$INPUT_MAT_FILE" ]]; then
    echo "ERROR: Input MAT file not found: $INPUT_MAT_FILE"
    exit 1
fi

# Pass values into MATLAB.
export KLT_FUNCTION_NAME
export KLT_INPUT_MAT_FILE="$INPUT_MAT_FILE"
export KLT_PROJECT_ROOT="$PROJECT_ROOT"
export KLT_INPUTS_DIR="$INPUTS_DIR"
export KLT_SAVE_DIR="$SYN_DIR"
export KLT_RESULTS_TAG="syn"
export KLT_CASE_CSV_FILE="${SYN_DIR}/synthetic_cases.csv"
export KLT_EXCEL_ROW="$S_ID"
export KLT_CASE_CSV_ROW="$CASE_CSV_ROW"

MATLAB_SCRIPT="${WORKDIR}/run_klt_row${S_ID}_${JOBTAG}.m"
MATLAB_LOG="${LOGDIR}/matlab-${JOBNAME}-${JOBTAG}.log"
MEMLOG="${LOGDIR}/mem-${JOBNAME}-${JOBTAG}.csv"

cat > "$MATLAB_SCRIPT" <<'MATLAB_EOF'
try
    fprintf('\n============================================================\n');
    fprintf('MATLAB KLT driver started.\n');
    fprintf('version: %s\n', version);
    fprintf('pwd: %s\n', pwd);
    fprintf('tempdir: %s\n', tempdir);
    fprintf('prefdir: %s\n', prefdir);
    fprintf('computer: %s\n', computer);
    fprintf('maxNumCompThreads: %d\n', maxNumCompThreads);
    fprintf('============================================================\n');

    projectRoot = getenv('KLT_PROJECT_ROOT');
    inputsDir = getenv('KLT_INPUTS_DIR');
    functionName = getenv('KLT_FUNCTION_NAME');
    inputMatFile = getenv('KLT_INPUT_MAT_FILE');
    caseCsvFile = getenv('KLT_CASE_CSV_FILE');
    excelRowText = getenv('KLT_EXCEL_ROW');
    caseCsvRowText = getenv('KLT_CASE_CSV_ROW');
    caseCsvRow = str2double(caseCsvRowText);

    cd(inputsDir);
    addpath(genpath(projectRoot));

    fprintf('\n============================================================\n');
    fprintf('Running function: %s\n', functionName);
    fprintf('Input MAT file:  %s\n', inputMatFile);
    fprintf('Save dir:        %s\n', getenv('KLT_SAVE_DIR'));
    fprintf('Results tag:     %s\n', getenv('KLT_RESULTS_TAG'));
    fprintf('CSV case file:   %s\n', caseCsvFile);
    fprintf('Excel row:       %s\n', excelRowText);
    fprintf('MATLAB data row: %s\n', caseCsvRowText);
    fprintf('============================================================\n');

    if ~isfile(inputMatFile)
        error('Input MAT file not found: %s', inputMatFile);
    end

    if exist(functionName, 'file') ~= 2
        error('MATLAB function not found on path: %s', functionName);
    end

    testConfig = struct();
    testConfig.saveDir = getenv('KLT_SAVE_DIR');
    testConfig.resultsTag = getenv('KLT_RESULTS_TAG');
    testConfig.verbose = true;
    testConfig.caseCsvFile = caseCsvFile;

    if ~isfinite(caseCsvRow) || caseCsvRow < 1 || caseCsvRow ~= round(caseCsvRow)
        error('Invalid KLT_CASE_CSV_ROW: %s', caseCsvRowText);
    end

    testConfig.caseCsvRow = round(caseCsvRow);

    % Resume from the deterministic per-case checkpoint/case-state files if available.
    % Leave these blank so the MATLAB harness derives them from:
    %   resultsTag + inputMatFile basename + CSV amplitude/wavelength.
    testConfig.resumeCheckpointFile = '';
    testConfig.resumeCaseStateFile = '';
    testConfig.resumeFromCheckpoint = true;
    testConfig.forceFresh = false;

    % Let MATLAB/HPC defaults decide parallel configuration.
    % Do not set parcluster.JobStorageLocation here.
    % Do not force c.NumWorkers here.
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('No parallel pool currently open before KLT function call.\n');
    else
        fprintf('Existing parallel pool before KLT function call: %d workers.\n', pool.NumWorkers);
    end

    results = feval(functionName, inputMatFile, testConfig);

    % Print resume diagnostics.
    if exist('results', 'var') && ~isempty(results)
        for ii = 1:numel(results)
            if isfield(results, 'checkpointFile')
                fprintf('Case %d checkpointFile: %s\n', ii, results(ii).checkpointFile);
            end
            if isfield(results, 'caseStateFile')
                fprintf('Case %d caseStateFile:  %s\n', ii, results(ii).caseStateFile);
            end
            if isfield(results, 'resumedFromCheckpoint')
                fprintf('Case %d resumedFromCheckpoint: %d\n', ii, results(ii).resumedFromCheckpoint);
            end
        end
    end

    pool = gcp('nocreate');
    if ~isempty(pool)
        fprintf('Leaving MATLAB: closing parallel pool with %d workers.\n', pool.NumWorkers);
        delete(pool);
    end

    fprintf('\nMATLAB KLT driver finished successfully.\n');
    exit(0);

catch ME
    fprintf(2, '\n============================================================\n');
    fprintf(2, 'MATLAB KLT driver failed.\n');
    fprintf(2, '============================================================\n');
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    exit(1);
end
MATLAB_EOF

echo "Created MATLAB script:"
echo "$MATLAB_SCRIPT"
echo "MATLAB logfile:"
echo "$MATLAB_LOG"
echo "Memory logfile:"
echo "$MEMLOG"
echo

echo "Launching MATLAB/2024a using stdin..."

# Run MATLAB in the background so the script can sample its memory use.
set +e

matlab -nodisplay -nosplash -logfile "$MATLAB_LOG" < "$MATLAB_SCRIPT" &
MATLAB_PID=$!

echo "MATLAB launcher PID: $MATLAB_PID"
echo "Starting RAM logger at 60-second interval..."

start_memory_logger "$MATLAB_PID" "$MEMLOG" 60 &
MEMLOGGER_PID=$!

wait "$MATLAB_PID"
MATLAB_STATUS=$?

cleanup_memory_logger

set -e

echo
echo "MATLAB exit status: $MATLAB_STATUS"
echo "Finished: $(date)"

print_memory_log_summary "$MEMLOG"
print_slurm_memory_summary

exit "$MATLAB_STATUS"
