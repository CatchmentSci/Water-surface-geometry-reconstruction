#!/bin/bash

#SBATCH --account=comet_wave
#SBATCH --partition=default_paid
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=matthew.perks@newcastle.ac.uk
#SBATCH --job-name=klt_analysis_case
#SBATCH --chdir=/mnt/nfs/home/nmp65
#SBATCH --output=/mnt/nfs/home/nmp65/klt_slurm_logs/%x-%j.out
#SBATCH --error=/mnt/nfs/home/nmp65/klt_slurm_logs/%x-%j.err

# Run one real/video KLT analysis case under Slurm using MATLAB/2024a.
#
# This version follows the working synthetic-case launch pattern:
#   - Explicitly load MATLAB/2024a.
#   - Redirect MATLAB startup/temp/cache/prefs paths away from /tmp and /tmp/slurmd.
#   - Launch MATLAB through stdin rather than -batch/-r inline commands.
#   - Keep Slurm logs in /mnt/nfs/home/nmp65/klt_slurm_logs.
#   - Keep generated runtime files in /mnt/nfs/home/nmp65/klt_runtime_jobs.
#   - Do not set parcluster.JobStorageLocation in bash.
#   - Do not create /tmp/matlab_parallel_* directories.
#
# Required sbatch export:
#   KLT_CASE_REF: either filenameIn, e.g. devon_dart20181129_13490, or videoNumber / solver value, e.g. 109
#
# Optional sbatch exports:
#   KLT_LOOKUP_FILE=/path/to/klt_analysis_case_lookup.tsv
#   KLT_PROJECT_ROOT=/mnt/nfs/home/nmp65/Downloads/v1_current
#   KLT_INPUTS_DIR=/mnt/nfs/home/nmp65/Downloads/v1_current/Videos/Inputs
#   KLT_OUTPUT_DIR=/mnt/nfs/home/nmp65/Downloads/v1_current/Videos/Inputs
#   KLT_SWEEP_LIMITS_FILE=sweep_limits.csv
#   KLT_RESUME_FROM_CHECKPOINT=1   # 1 means resume only if checkpoint exists; otherwise start fresh
#   KLT_SAVE_OUTPUTS=1

set -euo pipefail

LOGDIR="/mnt/nfs/home/nmp65/klt_slurm_logs"
RUNTIMEDIR="/mnt/nfs/home/nmp65/klt_runtime_jobs"
PROJECT_ROOT="${KLT_PROJECT_ROOT:-/mnt/nfs/home/nmp65/Downloads/v1_current}"
INPUTS_DIR="${KLT_INPUTS_DIR:-${PROJECT_ROOT}/Videos/Inputs}"
OUTPUT_DIR="${KLT_OUTPUT_DIR:-$INPUTS_DIR}"
SWEEP_LIMITS_FILE="${KLT_SWEEP_LIMITS_FILE:-sweep_limits.csv}"
RESUME_FROM_CHECKPOINT="${KLT_RESUME_FROM_CHECKPOINT:-1}"
SAVE_OUTPUTS="${KLT_SAVE_OUTPUTS:-1}"

JOBTAG="${SLURM_JOB_ID:-manual}"
JOBNAME="${SLURM_JOB_NAME:-klt_analysis_case}"
# Slurm keeps the original submission directory in SLURM_SUBMIT_DIR.
# Use it to find files kept beside the submit script, even after this
# worker moves into its runtime work directory.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
SAFE_BASE="${RUNTIMEDIR}/${JOBNAME}_${JOBTAG}"

WORKDIR="${SAFE_BASE}/work"
TMPROOT="${SAFE_BASE}/tmp"
PREFDIR="${SAFE_BASE}/prefs"
MCRROOT="${SAFE_BASE}/mcr_cache"
XDGCACHE="${SAFE_BASE}/xdg_cache"
XDGCONFIG="${SAFE_BASE}/xdg_config"
XDGRUNTIME="${SAFE_BASE}/xdg_runtime"

mkdir -p "$LOGDIR" "$RUNTIMEDIR" "$WORKDIR" "$TMPROOT" "$PREFDIR" "$MCRROOT" "$XDGCACHE" "$XDGCONFIG" "$XDGRUNTIME" "$OUTPUT_DIR"
chmod 700 "$SAFE_BASE" "$WORKDIR" "$TMPROOT" "$PREFDIR" "$MCRROOT" "$XDGCACHE" "$XDGCONFIG" "$XDGRUNTIME"

# MATLAB/2024a failed under Slurm when it used default /tmp paths.
# These redirects are startup/cache/temp paths only; they do not configure
# parpool or MATLAB parallel cluster profiles.
export TMPDIR="$TMPROOT"
export TMP="$TMPROOT"
export TEMP="$TMPROOT"
export MATLAB_PREFDIR="$PREFDIR"
export MCR_CACHE_ROOT="$MCRROOT"
export XDG_CACHE_HOME="$XDGCACHE"
export XDG_CONFIG_HOME="$XDGCONFIG"
export XDG_RUNTIME_DIR="$XDGRUNTIME"

cd "$WORKDIR"

echo "============================================================"
echo "KLT real/video analysis MATLAB/2024a Slurm job"
echo "============================================================"
echo "Started:          $(date)"
echo "Host:             $(hostname)"
echo "User:             $(whoami)"
echo "SLURM job ID:     ${SLURM_JOB_ID:-not_in_slurm}"
echo "SLURM job name:   ${SLURM_JOB_NAME:-not_in_slurm}"
echo "Submit dir:       ${SLURM_SUBMIT_DIR:-not_in_slurm}"
echo "Resolved submit:  $SUBMIT_DIR"
echo "Work dir:         $WORKDIR"
echo "Runtime dir:      $RUNTIMEDIR"
echo "Log dir:          $LOGDIR"
echo "TMPDIR:           $TMPDIR"
echo "MATLAB_PREFDIR:   $MATLAB_PREFDIR"
echo "MCR_CACHE_ROOT:   $MCR_CACHE_ROOT"
echo "XDG_RUNTIME_DIR:  $XDG_RUNTIME_DIR"
echo "CPUs allocated:   ${SLURM_CPUS_PER_TASK:-not_set}"
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

# -----------------------------------------------------------------------------
# Case selection
# -----------------------------------------------------------------------------
# KLT_CASE_REF may be either:
#   1) filenameIn, e.g. devon_dart20181129_13490
#   2) videoNumber / solver value, e.g. 109
# The mapping is read from klt_analysis_case_lookup.tsv by default.
# -----------------------------------------------------------------------------

CASE_REF="${KLT_CASE_REF:?KLT_CASE_REF is not set}"
LOOKUP_FILE_RAW="${KLT_LOOKUP_FILE:-klt_analysis_case_lookup.tsv}"

# If KLT_LOOKUP_FILE is absolute, use it as-is. If it is relative, interpret
# it relative to the Slurm submission directory, not the runtime work dir.
if [[ "$LOOKUP_FILE_RAW" = /* ]]; then
    LOOKUP_FILE="$LOOKUP_FILE_RAW"
else
    LOOKUP_FILE="${SUBMIT_DIR}/${LOOKUP_FILE_RAW}"
fi

if [[ ! -f "$LOOKUP_FILE" ]]; then
    echo "ERROR: Cannot find lookup file: $LOOKUP_FILE"
    echo "       KLT_LOOKUP_FILE was: ${KLT_LOOKUP_FILE:-not_set}"
    echo "       SLURM_SUBMIT_DIR was: ${SLURM_SUBMIT_DIR:-not_set}"
    echo "       Runtime pwd is: $(pwd)"
    exit 1
fi

mapfile -t MATCHES < <(
    awk -v ref="$CASE_REF" '
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
    ' "$LOOKUP_FILE"
)

if [[ "${#MATCHES[@]}" -eq 0 ]]; then
    echo "ERROR: Case reference '$CASE_REF' was not found in: $LOOKUP_FILE"
    echo "Use either a filenameIn value or a videoNumber from the lookup table."
    exit 1
fi

if [[ "${#MATCHES[@]}" -gt 1 ]]; then
    echo "ERROR: Case reference '$CASE_REF' matched more than one row in: $LOOKUP_FILE"
    printf '  %s\n' "${MATCHES[@]}"
    echo "Make filenameIn and videoNumber values unique."
    exit 1
fi

IFS=$'\t' read -r FILENAME_IN VIDEO_NUMBER <<< "${MATCHES[0]}"
SWEEP_VALUE="$VIDEO_NUMBER"

if [[ -z "$FILENAME_IN" || -z "$VIDEO_NUMBER" ]]; then
    echo "ERROR: Lookup row must contain both filenameIn and videoNumber."
    exit 1
fi

if [[ ! "$SWEEP_VALUE" =~ ^[0-9]+$ ]]; then
    echo "ERROR: videoNumber must be numeric. Got '$VIDEO_NUMBER' for '$FILENAME_IN'."
    exit 1
fi

INPUT_MAT_FILE="${INPUTS_DIR}/${FILENAME_IN}_solver_inputs.mat"
CHECKPOINT_FILE="${INPUTS_DIR}/${FILENAME_IN}_checkpoint.mat"
SWEEP_LIMITS_PATH="${INPUTS_DIR}/${SWEEP_LIMITS_FILE}"

if [[ ! -f "$INPUT_MAT_FILE" ]]; then
    echo "ERROR: Input MAT file not found: $INPUT_MAT_FILE"
    exit 1
fi

if [[ ! -f "$SWEEP_LIMITS_PATH" ]]; then
    echo "ERROR: Sweep limits file not found: $SWEEP_LIMITS_PATH"
    exit 1
fi

echo "Running analysis case"
echo "Case reference:          $CASE_REF"
echo "filenameIn:              $FILENAME_IN"
echo "videoNumber:             $VIDEO_NUMBER"
echo "Solver value argument:    $SWEEP_VALUE"
echo "Project root:            $PROJECT_ROOT"
echo "Inputs dir:              $INPUTS_DIR"
echo "Output dir:              $OUTPUT_DIR"
echo "Lookup file:             $LOOKUP_FILE"
echo "Input MAT file:          $INPUT_MAT_FILE"
echo "Checkpoint file:         $CHECKPOINT_FILE"
echo "Sweep limits file:       $SWEEP_LIMITS_PATH"
echo "Resume from checkpoint:  $RESUME_FROM_CHECKPOINT"
echo "Save outputs:            $SAVE_OUTPUTS"
echo

# Avoid oversubscription in numerical libraries. This mirrors your original
# synthetic script but does not create a MATLAB parallel JobStorageLocation.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Pass values into MATLAB.
export KLT_PROJECT_ROOT="$PROJECT_ROOT"
export KLT_INPUTS_DIR="$INPUTS_DIR"
export KLT_OUTPUT_DIR="$OUTPUT_DIR"
export KLT_FILENAME_IN="$FILENAME_IN"
export KLT_VIDEO_NUMBER="$VIDEO_NUMBER"
export KLT_SWEEP_VALUE="$SWEEP_VALUE"
export KLT_SWEEP_LIMITS_FILE="$SWEEP_LIMITS_FILE"
export KLT_RESUME_FROM_CHECKPOINT="$RESUME_FROM_CHECKPOINT"
export KLT_SAVE_OUTPUTS="$SAVE_OUTPUTS"

SAFE_VIDEO_NUMBER=$(printf '%s' "$VIDEO_NUMBER" | tr -c 'A-Za-z0-9_' '_')
MATLAB_SCRIPT="${WORKDIR}/run_klt_v${SAFE_VIDEO_NUMBER}_${JOBTAG}.m"
MATLAB_LOG="${LOGDIR}/matlab-${JOBNAME}-${JOBTAG}.log"

cat > "$MATLAB_SCRIPT" <<'EOF_MATLAB'
try
    fprintf('\n============================================================\n');
    fprintf('MATLAB KLT real/video analysis driver started.\n');
    fprintf('version: %s\n', version);
    fprintf('pwd: %s\n', pwd);
    fprintf('tempdir: %s\n', tempdir);
    fprintf('prefdir: %s\n', prefdir);
    fprintf('computer: %s\n', computer);
    fprintf('maxNumCompThreads: %d\n', maxNumCompThreads);
    fprintf('============================================================\n');

    projectRoot = getenv('KLT_PROJECT_ROOT');
    inputsDir = getenv('KLT_INPUTS_DIR');
    outputDir = getenv('KLT_OUTPUT_DIR');
    filenameIn = getenv('KLT_FILENAME_IN');

    % Use absolute paths so the solver can still find files even if it changes pwd internally.
    inputMatFile = fullfile(inputsDir, [filenameIn '_solver_inputs.mat']);
    checkpointFile = fullfile(inputsDir, [filenameIn '_checkpoint.mat']);

    requestedResumeFromCheckpoint = str2double(getenv('KLT_RESUME_FROM_CHECKPOINT'));
    if isnan(requestedResumeFromCheckpoint)
        requestedResumeFromCheckpoint = 1;
    end

    sweepLimitsFileName = getenv('KLT_SWEEP_LIMITS_FILE');
    sweepLimitsFile = fullfile(inputsDir, sweepLimitsFileName);
    videoNumber = str2double(getenv('KLT_VIDEO_NUMBER'));
    sweepValue = str2double(getenv('KLT_SWEEP_VALUE'));
    saveOutputs = str2double(getenv('KLT_SAVE_OUTPUTS'));
    outputMatFile = fullfile(outputDir, [filenameIn '_hpc_outputs.mat']);

    cd(inputsDir);
    addpath(genpath(projectRoot));

    fprintf('\n============================================================\n');
    fprintf('Running KLT analysis case\n');
    fprintf('filenameIn:             %s\n', filenameIn);
    fprintf('inputMatFile:           %s\n', inputMatFile);
    fprintf('checkpointFile:         %s\n', checkpointFile);
    fprintf('requested resume:       %d\n', requestedResumeFromCheckpoint);
    fprintf('sweepLimitsFile:        %s\n', sweepLimitsFile);
    fprintf('videoNumber:            %d\n', videoNumber);
    fprintf('solver value argument:   %d\n', sweepValue);
    fprintf('outputMatFile:          %s\n', outputMatFile);
    fprintf('============================================================\n');

    if ~isfile(inputMatFile)
        error('Input MAT file not found: %s', inputMatFile);
    end

    if ~isfile(sweepLimitsFile)
        error('Sweep limits file not found: %s', sweepLimitsFile);
    end

    % Resume only when a checkpoint actually exists. If the user requested
    % resume but no checkpoint is present, force a fresh run before calling
    % the solver. This prevents KLT_wse_solver_paths_Apoint_block_jacobi_v5
    % from stopping at its missing-checkpoint guard.
    checkpointExists = isfile(checkpointFile);
    resumeFromCheckpoint = (requestedResumeFromCheckpoint ~= 0) && checkpointExists;

    fprintf('checkpoint exists:      %d\n', checkpointExists);
    fprintf('effective resume:       %d\n', resumeFromCheckpoint);

    if requestedResumeFromCheckpoint ~= 0 && checkpointExists
        fprintf('Checkpoint found; resuming from: %s\n', checkpointFile);
    elseif requestedResumeFromCheckpoint ~= 0 && ~checkpointExists
        fprintf('Requested resume, but checkpoint file was not found: %s\n', checkpointFile);
        fprintf('Starting this case from fresh.\n');
    elseif requestedResumeFromCheckpoint == 0 && checkpointExists
        fprintf('Checkpoint exists but resume was not requested; starting fresh.\n');
    else
        fprintf('No checkpoint found; starting fresh.\n');
    end

    if exist('KLT_wse_solver_paths_Apoint_block_jacobi_v5', 'file') ~= 2
        error('MATLAB function not found on path: KLT_wse_solver_paths_Apoint_block_jacobi_v5');
    end

    % Follow the working synthetic-job approach: let MATLAB/HPC defaults or
    % the KLT solver decide parallel configuration. Do not set
    % parcluster.JobStorageLocation here.
    pool = gcp('nocreate');
    if isempty(pool)
        fprintf('No parallel pool currently open before KLT solver call.\n');
    else
        fprintf('Existing parallel pool before KLT solver call: %d workers.\n', pool.NumWorkers);
    end

    [wse_map, metrics, phi] = KLT_wse_solver_paths_Apoint_block_jacobi_v5( ...
        [], [], [], [], [], [], ...
        checkpointFile, resumeFromCheckpoint, inputMatFile, ...
        sweepLimitsFile, sweepValue, [], ...
        []);

    if saveOutputs
        save(outputMatFile, 'wse_map', 'metrics', 'phi', 'filenameIn', ...
            'inputMatFile', 'checkpointFile', 'requestedResumeFromCheckpoint', ...
            'resumeFromCheckpoint', 'sweepLimitsFile', 'videoNumber', ...
            'sweepValue', '-v7.3');
        fprintf('Saved outputs to: %s\n', outputMatFile);
    end

    pool = gcp('nocreate');
    if ~isempty(pool)
        fprintf('Leaving MATLAB: closing parallel pool with %d workers.\n', pool.NumWorkers);
        delete(pool);
    end

    fprintf('\nMATLAB KLT real/video analysis driver finished successfully.\n');
    exit(0);

catch ME
    fprintf(2, '\n============================================================\n');
    fprintf(2, 'MATLAB KLT real/video analysis driver failed.\n');
    fprintf(2, '============================================================\n');
    fprintf(2, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    exit(1);
end
EOF_MATLAB

echo "Created MATLAB script:"
echo "$MATLAB_SCRIPT"
echo "MATLAB logfile:"
echo "$MATLAB_LOG"
echo

echo "Launching MATLAB/2024a using stdin..."
matlab -nodisplay -nosplash -logfile "$MATLAB_LOG" < "$MATLAB_SCRIPT"
MATLAB_STATUS=$?

echo
echo "MATLAB exit status: $MATLAB_STATUS"
echo "Finished: $(date)"

exit "$MATLAB_STATUS"
